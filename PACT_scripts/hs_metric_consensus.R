#!/usr/bin/env Rscript
## Script name: hs_metric_consensus.R
## Purpose: Saves HS Metric column values for input
## Date Created: June 2, 2025
## Version: 1.0.0
## Author: Jonathan Serrano
## Copyright (c) NYULH Jonathan Serrano, 2025

args <- commandArgs(TRUE)
args[1] -> run_id       # PACT RUN ID
args[1] -> pact_id      # PACT ID (i.e. "PACT-YY-##")

runYear <- paste0("20", substring(run_id, 1, 2))

MOL_DIR <- "/Volumes/molecular/Molecular/NGS607"
OUTPUT_DIR <- file.path(MOL_DIR, runYear, run_id, "output")

METRICS_DIR <- file.path(OUTPUT_DIR, "CollectHsMetrics")
OUT_FOLDER <- file.path(OUTPUT_DIR, "HS_Metrics_CSV")

# Format: <HS_metric column> = <Column Label>
metric_cols <- list(
    AT_DROPOUT = "AT_DROPOUT",
    GC_DROPOUT = "GC_DROPOUT",
    PCT_250_CVG = "PCT_TARGET_BASES_250X",
    FOLD_80_BASE_PENALTY = "FOLD_80_BASE_PENALTY",
    MEAN_TARG_CVG = "MEAN_TARGET_COVERAGE"
)

options(repos = c(CRAN = "https://cran.rstudio.com/"))
user_lib <- Sys.getenv("R_LIBS_USER", unset = "~/R/library")
dir.create(user_lib, showWarnings = FALSE, recursive = TRUE)
.libPaths(c(user_lib, .libPaths()))

snapshot_date <- "2025-05-01"
bioc_version  <- "3.22"

options(
    download.file.method = "curl",
    repos = c(
        CRAN = sprintf("https://packagemanager.posit.co/cran/%s", snapshot_date),
        BiocManager::repositories(version = bioc_version)
    )
)

if (!"devtools" %in% rownames(installed.packages())) {
    install.packages("devtools", ask = FALSE, dependencies = TRUE)
}

library("devtools")

if (!"Rcpp" %in% rownames(installed.packages())) {
    install.packages("Rcpp", type = "binary", repos = "https://packagemanager.rstudio.com/all/latest")
}

# Installs the pak package wtih a tryCatch
install_pak <- function() {
    tryCatch(
        install.packages(
            "pak", repos = sprintf("https://r-lib.github.io/p/pak/stable/%s/%s/%s",
                        .Platform$pkgType, R.Version()$os, R.Version()$arch)),
        error = function(e) {
            install.packages("pak", ask = FALSE, dependencies = TRUE,
                repos = "https://packagemanager.rstudio.com/all/latest"
            )
        })
}


#' Ensures Required Packages are Loaded the `pak` package installs any missing
#' @param pkgs A character vector of package names
ensure_packages <- function(pkgs) {
    installed_pkgs <- rownames(installed.packages())
    missing_pkgs <- setdiff(pkgs, installed_pkgs)
    if (length(missing_pkgs) > 0) {
        message(missing_pkgs)
        if (!"pak" %in% installed_pkgs) install_pak()
        library("pak")
        for (pkg in missing_pkgs) {
            tryCatch(
                pak::pkg_install(pkg, ask = FALSE),
                error = function(e) {
                    install.packages(pkg, ask = FALSE, dependencies = TRUE)
                }
            )
        }
    }
    stopifnot(all(sapply(pkgs, function(pkg) {
        suppressPackageStartupMessages(library(
            pkg, mask.ok = TRUE, character.only = TRUE, logical.return = TRUE
        ))
    })))
}

# Add all the packages used for subsequent scripts here below
ensure_packages(c("tidyr", "readr", "dplyr"))

#' Retrieve the NGS607 Directory Path for a Given Run ID ----------------------
#' @return A character vector containing the full path to the NGS607 directory
check_pact_path <- function(run_id) {
    runYear <- paste0("20", substring(run_id, 1, 2))
    NGS_607_DIR <- file.path("/Volumes/molecular/Molecular/NGS607", runYear)
    runPath <- dir(path = NGS_607_DIR, pattern = run_id,
                   include.dirs = TRUE, full.names = TRUE)
    if (length(runPath) > 1) {
        message("Multiple paths:\n", paste(runPath, collapse = "\n"))
        pact_path <- grepl("_PACT", basename(runPath))
        runPath <- if (any(pact_path)) runPath[pact_path] else runPath[1]
        message("Path chosen:\n", basename(runPath))
    }
    return(runPath)
}

# Appends the mean and standard deviations to the CSV data frames -------------
get_mean_and_sd <- function(combined_df) {
    mean2 <- round(mean(combined_df[[2]], na.rm = TRUE), 2)
    sd2   <- round(sd(combined_df[[2]], na.rm = TRUE), 2)

    mean_row <- data.frame(
        setNames(list("Mean:"), names(combined_df)[1]),
        setNames(list(mean2)  , names(combined_df)[2]),
        stringsAsFactors = FALSE,
        check.names = FALSE
    )
    sd_row <- data.frame(
        setNames(list("SD:"), names(combined_df)[1]),
        setNames(list(sd2)   , names(combined_df)[2]),
        stringsAsFactors = FALSE,
        check.names = FALSE
    )

    combined_df <- rbind(combined_df, mean_row, sd_row)
    return(combined_df)
}


# Saves data frame list of sample names to a csv file in output folder
metrics_to_csv <- function(df_list, metric_type) {
    padded_list <- lapply(df_list, function(df) {
        colname <- names(df)[1]
        new_df <- data.frame(sample_names = rownames(df), values = df[[1]])
        sams <- stringr::str_split_fixed(new_df$sample_names, "_", 6)[, 6]
        new_df$sample_names <- sams
        names(new_df) <- c(paste(colname, "Sample Names"), metric_type)
        new_df
    })
    combined_df <- do.call(cbind, padded_list)
    combined_df[, 2] <- round(as.numeric(combined_df[, 2]), 2)
    combined_df <- get_mean_and_sd(combined_df)
    output_csv <- paste0(metric_type,"_hs_metric_values.csv")
    output_csv_path <- file.path(OUT_FOLDER, output_csv)
    if (!dir.exists(OUT_FOLDER)) dir.create(OUT_FOLDER, recursive = TRUE)
    message("Saving file:\n", output_csv_path)
    write.csv(combined_df, file = output_csv_path, row.names = F, quote = F)
}


#' Reads the input text file and extracts value from AT_DROPOUT column
#' @param txt_file A file path to the collect HS metric text file
#' @return A single AT dropout value for the input sample file
read_hs_metric <- function(txt_file, col_name = "AT_DROPOUT") {
    sam_name <-
        stringr::str_split_fixed(basename(txt_file), "_hs_metrics", 2)[1, 1]
    if (is.null(sam_name)) stop(txt_file)
    dropout_value <- read.table(
        txt_file, sep = "\t", skip = 6, header = TRUE, check.names = FALSE,
        blank.lines.skip = FALSE, comment.char = "")[1, col_name]
    if (is.null(dropout_value))
        stop(paste("Column not found in file:", col_name, txt_file))
    names(dropout_value) <- sam_name
    return(dropout_value)
}

# Per-run helper: locate tumor & normal hs_metric file paths
get_file_lists <- function(run_id) {
    runPath <- check_pact_path(run_id)
    if (length(runPath) == 0) return(NULL)
    tn_df <- read.csv(file.path(runPath, "samples.tumor.normal.csv"))
    all_files  <- dir(METRICS_DIR,
                      pattern = "_hs_metrics\\.txt$", full.names = TRUE)
    tumor_files <- all_files[grepl(
        paste(tn_df$Tumor, collapse = "|"), all_files)]
    normal_files <- all_files[grepl(
        paste(c(tn_df$Normal, "SC_SERACARE"), collapse = "|"), all_files)]
    return(list(tumor = tumor_files, normal = normal_files))
}

# Saves the hs_metric values for all runs to a csv file and then plots values
save_tumor_normals <- function(metric_res) {
    label   <- metric_res$label
    normals <- metric_res$normals
    tumors  <- metric_res$tumors

    if (length(normals) > 0) {
        metrics_to_csv(df_list = normals, metric_type = paste0(label, "_normals"))
    }
    if (length(tumors) > 0) {
        metrics_to_csv(df_list = tumors, metric_type = paste0(label, "_tumors"))
    }
}


# Extracts one metric column into a single-column data.frame
process_dataframe <- function(hs_files, pact_id, col_name) {
    setNames(data.frame(unlist(
        lapply(hs_files, read_hs_metric, col_name = col_name)
    )), pact_id)
}


# Per-metric helper: for one (column, label) pair, extract the two data.frames
process_metric <- function(hs_files, pact_id, col_name, label) {
    list(label   = label,
         normals = process_dataframe(hs_files$normal, pact_id, col_name),
         tumors  = process_dataframe(hs_files$tumor,  pact_id, col_name)
    )
}

# Helper to append a metric’s data.frames into the results object
append_to_results <- function(results, metric_res) {
    results$normals <- append(results$normals, list(metric_res$normals))
    results$tumors  <- append(results$tumors,  list(metric_res$tumors))
    return(results)
}


#############
# MAIN: loop over runs, reading the hs_metrics files and extracting the metrics
#############
loop_metrics <- function(run_id, pact_id, metric_cols) {
    # iterate through each metric column
    for (metric in names(metric_cols)) {
        metric_label <- metric_cols[[metric]]
        metric_res <- list(label = metric, normals = list(), tumors = list())
        hs_files <- get_file_lists(run_id)
        if (is.null(hs_files)) stop(paste("No HS metrics for:", run_id))
        # Read hs_metric files per metric column and save tumor and normal csv
        metric_data <- process_metric(hs_files, pact_id, metric_label, metric)
        metric_res <- append_to_results(metric_res, metric_data)
        save_tumor_normals(metric_res)
    }
}

loop_metrics(run_id, pact_id, metric_cols)
