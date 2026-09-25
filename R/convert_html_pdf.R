#!/usr/bin/env Rscript
## Script name: convert_html_pdf.R
## Purpose: Functions to save methylation output files as pdfs to upload for REDCap
## Date Created: May 22, 2026
## Version: 1.0.0
## Author: Jonathan Serrano
## Copyright (c) NYULH Jonathan Serrano, 2026

options(stringsAsFactors = FALSE)
API_URL = "https://redcap.nyumc.org/apps/redcap/api/"
CHROME_BIN <- "/Applications/Google Chrome.app/Contents/MacOS/Google Chrome"

TOKEN_PATH <- "/Volumes/CBioinformatics/scripts/METH_DB_API.txt"

if (!file.exists(CHROME_BIN)) {
    stop("Chrome executable does not exist: ", CHROME_BIN, call. = FALSE)
}

stopifnot(file.exists(TOKEN_PATH))
REDCAP_API_TOKEN <- trimws(readLines(TOKEN_PATH, n = 1, warn = FALSE))

legend_css <- paste0(paste(
    readLines("/Volumes/CBioinformatics/scripts/legend_css.txt", warn = FALSE),
    collapse = "\n", sep = ""),"\n")

legend_js <- paste(trimws(readLines(
    "/Volumes/CBioinformatics/scripts/legend_js.txt", warn = FALSE),
    which = "right"), collapse = "\n")

legend_js <- paste0("<script>\n", legend_js, "\n</script>\n")


# Checks if a package is installed
not_installed <- function(pkgName) return(!pkgName %in% rownames(installed.packages()))

if (not_installed("pagedown")) {
    install.packages("pagedown", ask = FALSE, dependencies = c("Depends","Imports"))
}

suppressPackageStartupMessages({requireNamespace("pagedown", quietly = TRUE)})
options(pagedown.remote.maxattempts = 40)
options(pagedown.remote.sleeptime = 2)

sanitize_string <- function(input_file, ending = NULL) {
    input_file <- paste0(tools::file_path_sans_ext(basename(input_file)))
    input_file <- iconv(input_file, from = "", to = "ASCII//TRANSLIT", sub = "")
    input_file <- gsub("[^A-Za-z0-9._-]+", "-", trimws(input_file))
    input_file <- gsub("-+", "-", input_file)
    input_file <- gsub("^[._-]+|[._-]+$", "", input_file)
    if (!is.null(ending))
        return(paste0(input_file, ending))
    return(input_file)
}


inject_plotly_legend <- function(html_file, output_html = NULL) {
    html_text <- paste(readLines(html_file, warn = FALSE), collapse = "\n")

    patch <- paste0(legend_css, "\n", legend_js)

    if (grepl("</head>", html_text, ignore.case = TRUE)) {
        html_text <- sub("</head>", paste0(patch, "</head>"), html_text, ignore.case = TRUE)
    } else {
        html_text <- paste0(patch, "\n", html_text)
    }

    if (is.null(output_html)) {
        output_html <- tempfile(fileext = ".html")
    }

    writeLines(html_text, output_html, useBytes = TRUE)
    return(output_html)
}


reports_to_pdf <- function(input_dir, sam_name = NULL) {

    html_files <- list.files(path = input_dir, pattern = ".html?$", full.names = TRUE,
        recursive = FALSE, ignore.case = TRUE)

    html_files <- html_files[!grepl("_QC\\.html$", html_files)]

    if (!length(html_files) > 0) {
        return(message("No HTML files found in:\n", input_dir))
    }

    if (!is.null(sam_name)) {
        html_files <- html_files[stringr::str_detect(basename(html_files), pattern = sam_name)]
        if (length(html_files) != 1) {
            message("Multiple matches for sample name: ", sam_name)
            message("Files:\n", paste(html_files, collapse = "\n"))
            base_htmls <- stringr::str_split_fixed(basename(html_files), "_", 2)[,
                1]
            base_htmls <- stringr::str_split_fixed(basename(base_htmls), ".html",
                2)[, 1]
            idx_match <- which(sam_name %in% base_htmls)
            html_files <- html_files[idx_match]
            message("Sample: ", sam_name, "\n", "Match found: ", html_files)
        }
    }

    for (html_file in html_files) {
        pdf_file <- file.path(input_dir, sanitize_string(html_file, ".pdf"))
        if (file.exists(pdf_file))
            next
        message("Converting file to PDF: ", basename(html_file))

        chrome_args <- c("--headless=new", "--no-sandbox", "--enable-webgl", "--ignore-gpu-blocklist",
            "--use-angle=swiftshader", "--use-gl=angle", "--window-size=1600,1200",
            "--no-first-run", "--no-default-browser-check", "--disable-extensions",
            "--disable-background-networking", "--disable-sync", "--disable-translate",
            "--disable-component-update", "--disable-client-side-phishing-detection",
            "--disable-popup-blocking", "--metrics-recording-only", "--mute-audio")

        work_dir <- file.path(tempdir(), paste0("pagedown_", sanitize_string(html_file)))
        patched_html <- inject_plotly_legend(html_file)
        tryCatch(expr = suppressWarnings(pagedown::chrome_print(
            input = patched_html,
            output = pdf_file, browser = CHROME_BIN, wait = 4, timeout = 120, extra_args = chrome_args,
            work_dir = work_dir, verbose = FALSE, outline = FALSE)),
            error = function(e) {
            suppressWarnings(pagedown::chrome_print(
                input = patched_html, output = pdf_file,
                browser = CHROME_BIN, wait = 4, timeout = 120, extra_args = chrome_args,
                work_dir = work_dir, verbose = TRUE, outline = FALSE
            ))
        })
    }
}


upload_pdf <- function(recordName, input_dir, fld = "classifier_pdf") {
    message("\nChecking PDF Record Report upload for: ", recordName)

    pdf_path <- dir(path = input_dir, pattern = sprintf("^%s([^0-9].*)?\\.pdf$",
        recordName), full.names = TRUE)

    if (length(pdf_path) == 0) {
        pdf_path <- pdf_path[grepl(pattern = "V13_1\\.pdf$",
                                   pdf_path)]
    }

    if (length(pdf_path) > 1) {
        pdf_path <- pdf_path[grepl(pattern = "V13_1\\.pdf$",
            pdf_path)]
    }

    if (length(pdf_path) != 1) {
        message("REDCap Import error for record: ", recordName)
        return(message("length(pdf_path) != 1"))
    }

    check_payload <- list(
        token = REDCAP_API_TOKEN, content = "record",
        action = "export", format = "json", type = "flat", records = recordName,
        fields = fld, rawOrLabel = "raw", rawOrLabelHeaders = "raw",
        exportCheckboxLabel = "false", exportSurveyFields = "false",
        exportDataAccessGroups = "false"
        )

    res <- httr::POST(url = API_URL, body = check_payload, encode = "form")

    httr::stop_for_status(res)

    record <- jsonlite::fromJSON(httr::content(res, as = "text", encoding = "UTF-8"))

    if (length(record) != 0) {
        return(message("File already exists for: ", recordName))
    }

    message("Uploading file:\n", pdf_path, "\nTo REDCap Record: ",
        recordName)

    api_payload <- list(token = REDCAP_API_TOKEN, content = "file",
        action = "import", record = recordName, field = fld,
        file = httr::upload_file(pdf_path, type = "application/pdf"))

    tryCatch({
        res <- httr::POST(url = API_URL, body = api_payload, encode = "multipart")

        if (httr::status_code(res) == 200) {
            message("Successfully uploaded PDF to REDCap record: ", recordName)
        } else {
            message("REDCap Import error for record: ", recordName,
                "\nError Code: ", httr::status_code(res), "\nDetails: ",
                httr::content(res, as = "text", encoding = "UTF-8"))
        }
    }, error = function(e) {
        message("REDCap Import error for record: ", recordName)
        message("File failed to upload: ", pdf_path)
    })
}


import_files_redcap <- function(input_dir) {
    csv_file <- dir(input_dir, pattern = "samplesheet.csv", full.names = TRUE)
    if (length(csv_file) == 1) {
        sam_sheet <- as.data.frame(read.csv(csv_file))
    } else {
        stop("samplesheet.csv is missng")
    }

    sam_list <- sam_sheet$Sample_Name
    for (sam in sam_list) {
        upload_pdf(sam, input_dir)
    }

}
