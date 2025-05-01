#!/usr/bin/env Rscript
## Script name: search_db.R
## Purpose: Source of functions to search REDCap fields for matching records
## Date Created: April 11, 2024
## Version: 1.0.0
## Author: Jonathan Serrano
## Copyright (c) NYULH Jonathan Serrano, 2025

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) stop("Usage: Rscript search_db.R <input_csv> <api_token>")

# Main Inputs -----------------------------------------------------------------
input_csv <- args[1]
api_token <- args[2]
api_url <- "https://redcap.nyumc.org/apps/redcap/api/"

# Fields in REDCap to check against -------------------------------------------
fields2check <- c(
    "record_id",
    "barcode_and_row_column",
    "tm_number",
    "accession_number",
    "nyu_mrn"
)

# Fields from REDCap to export for matched records ----------------------------
outputFields <- c(
    "record_id","last_name", "first_name", "tm_number", "accession_number",
    "nyu_mrn", "barcode_and_row_column"
) # View all fields available using `redcapAPI::exportFieldNames(rcon)`


# Load required libraries (install if necessary) ------------------------------
loadPacks <- function() {
    options(repos = c(getOption("repos"), CRAN = "https://cran.r-project.org"))
    pkgs <- c("utils", "data.table", "tidyr", "dplyr", "redcapAPI", "stringr")
    for (pkgName in pkgs) {
        if (!pkgName %in% rownames(installed.packages()))
            install.packages(pkgName, dependencies = T, ask = F)
        if (!pkgName %in% loadedNamespaces())
            suppressPackageStartupMessages(library(
                pkgName, character.only = TRUE, quietly = TRUE, logical.return = TRUE))
    }
}

# Function to export all REDCap records for given fields ----------------------
grabAllRecords <- function(fields2check, rcon) {
    params <- list(rcon, fields = fields2check, survey = FALSE, dag = FALSE,
                   factors = FALSE, form_complete_auto = FALSE)
    df <- do.call(redcapAPI::exportRecordsTyped, params)  %>% as.data.frame()
    # Filter out any records that do not have RD-number and no sentrix ID
    redcap_df <- df %>%
        dplyr::filter(stringr::str_detect(record_id, "^RD-"),
                      stringr::str_detect(barcode_and_row_column, "_R0"))
    return(redcap_df)
}

# Function to search for string matches across REDCap fields ------------------
get_matches <- function(input_df, redcap_df) {
    search_values <- input_df %>%
        dplyr::mutate(dplyr::across(dplyr::everything(), as.character)) %>%
        tidyr::pivot_longer(cols = dplyr::everything(), values_to = "value") %>%
        dplyr::filter(!is.na(value), value != "") %>%
        dplyr::distinct(value) %>%
        dplyr::pull(value)

    matched_df <- redcap_df %>%
        dplyr::filter(dplyr::if_any(dplyr::everything(), ~ .x %in% search_values)) %>%
        dplyr::distinct()
    return(matched_df)
}


# Saves output dataframe to a tsv file on the Desktop -------------------------
save_output <- function(matched_df, input_csv, output_suffix) {
    output_file <- sub("\\.csv$", output_suffix, basename(input_csv))
    output_path <- file.path(fs::path_home(), "Desktop", output_file)
    write.table(matched_df, file = output_path, quote = FALSE,
                row.names = FALSE, sep = "\t")
    message("File saved to: ", output_path)
}

# Main function loads packages, connects to REDCap, and finds matches ---------
find_matches <- function(input_csv, fields2check, outputFields) {

    stopifnot("record_id" %in% fields2check,
              "barcode_and_row_column" %in% fields2check)
    loadPacks()
    rcon <- redcapAPI::redcapConnection(url = api_url, token = api_token)

    redcap_df <- grabAllRecords(fields2check, rcon)
    input_df <- as.data.frame(read.csv(input_csv))

    matched_df <- get_matches(input_df, redcap_df)
    save_output(matched_df, input_csv, "_matched.tsv")

    all_matches <- redcapAPI::exportRecordsTyped(
      rcon, records = paste(matched_df$record_id), fields = outputFields)
    save_output(all_matches, input_csv, "_all_matches.tsv")
}

# MAIN
find_matches(input_csv, fields2check, outputFields)
