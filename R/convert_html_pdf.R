#!/usr/bin/env Rscript
## Script name: convert_html_pdf.R
## Purpose: Functions to save methylation output files as pdfs to upload for REDCap
## Date Created: May 22, 2026
## Version: 1.0.0
## Author: Jonathan Serrano
## Copyright (c) NYULH Jonathan Serrano, 2026

options(stringsAsFactors = FALSE)
gb <- globalenv(); assign("gb", gb)
apiLink = "https://redcap.nyumc.org/apps/redcap/api/"

# Checks if a package is installed
not_installed <- function(pkgName) return(!pkgName %in% rownames(installed.packages()))

if (not_installed("pagedown")) {
    install.packages("pagedown", ask = FALSE, dependencies = c("Depends", "Imports"))
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
    if (!is.null(ending)) return(paste0(input_file, ending))
    return(input_file)
}


inject_plotly_legend <- function(html_file, output_html = NULL) {
    html_text <- paste(readLines(html_file, warn = FALSE), collapse = "\n")

    legend_css <- paste0(
        "<style>\n",
        ".legend-export-page {\n",
        "  break-before: page;\n",
        "  page-break-before: always;\n",
        "  padding: 24px;\n",
        "  font-family: 'Open Sans', Verdana, Arial, sans-serif;\n",
        "}\n",
        ".legend-export-page h2 { margin: 0 0 16px 0; }\n",
        ".legend-export-grid { column-count: 3; column-gap: 32px; }\n",
        ".legend-export-row {\n",
        "  break-inside: avoid;\n",
        "  display: flex;\n",
        "  align-items: center;\n",
        "  margin-bottom: 6px;\n",
        "  font-size: 11px;\n",
        "}\n",
        ".legend-export-symbol {\n",
        "  width: 28px;\n",
        "  height: 16px;\n",
        "  margin-right: 8px;\n",
        "  flex: 0 0 auto;\n",
        "  overflow: visible;\n",
        "}\n",
        ".legend-export-label {\n",
        "  word-break: break-word;\n",
        "}\n",
        "</style>\n"
    )

    legend_js <- paste0(
        "<script>\n",
        "window.addEventListener('load', function() {\n",
        "  setTimeout(function() {\n",
        "    const traces = Array.from(document.querySelectorAll('.legend .traces'));\n",
        "    if (!traces.length) return;\n",
        "\n",
        "    const seen = new Set();\n",
        "    const page = document.createElement('div');\n",
        "    page.className = 'legend-export-page';\n",
        "    page.innerHTML = '<h2>Legend</h2><div class=\"legend-export-grid\"></div>';\n",
        "    const grid = page.querySelector('.legend-export-grid');\n",
        "\n",
        "    function cssValue(el, prop, fallback) {\n",
        "      if (!el) return fallback;\n",
        "      const computed = window.getComputedStyle(el);\n",
        "      const value = computed.getPropertyValue(prop);\n",
        "      return value && value !== 'none' ? value : fallback;\n",
        "    }\n",
        "\n",
        "    function attrOrStyle(el, attr, cssProp, fallback) {\n",
        "      if (!el) return fallback;\n",
        "      const attrValue = el.getAttribute(attr);\n",
        "      if (attrValue && attrValue !== 'none') return attrValue;\n",
        "      return cssValue(el, cssProp, fallback);\n",
        "    }\n",
        "\n",
        "    traces.forEach(function(trace) {\n",
        "      const textNode = trace.querySelector('.legendtext');\n",
        "      if (!textNode) return;\n",
        "\n",
        "      const label = textNode.textContent.trim();\n",
        "      if (!label || seen.has(label)) return;\n",
        "      seen.add(label);\n",
        "\n",
        "      const point = trace.querySelector('.legendsymbols path, .legendpoints path, path.scatterpts, circle, rect');\n",
        "      const line = trace.querySelector('.legendlines path, path.js-line, line');\n",
        "\n",
        "      const fill = attrOrStyle(point, 'fill', 'fill', '#777777');\n",
        "      const stroke = attrOrStyle(point, 'stroke', 'stroke', 'none');\n",
        "      const opacity = attrOrStyle(point, 'opacity', 'opacity', '1');\n",
        "      const fillOpacity = attrOrStyle(point, 'fill-opacity', 'fill-opacity', opacity || '1');\n",
        "      const strokeOpacity = attrOrStyle(point, 'stroke-opacity', 'stroke-opacity', opacity || '1');\n",
        "      const strokeWidth = attrOrStyle(point, 'stroke-width', 'stroke-width', '0');\n",
        "\n",
        "      const lineStroke = attrOrStyle(line, 'stroke', 'stroke', fill);\n",
        "      const lineWidth = attrOrStyle(line, 'stroke-width', 'stroke-width', '0');\n",
        "\n",
        "      let pathD = null;\n",
        "      if (point && point.tagName.toLowerCase() === 'path') {\n",
        "        pathD = point.getAttribute('d');\n",
        "      }\n",
        "\n",
        "      const row = document.createElement('div');\n",
        "      row.className = 'legend-export-row';\n",
        "\n",
        "      const symbolSvg = document.createElementNS('http://www.w3.org/2000/svg', 'svg');\n",
        "      symbolSvg.setAttribute('class', 'legend-export-symbol');\n",
        "      symbolSvg.setAttribute('viewBox', '0 0 28 16');\n",
        "\n",
        "      if (line && parseFloat(lineWidth) > 0) {\n",
        "        const staticLine = document.createElementNS('http://www.w3.org/2000/svg', 'line');\n",
        "        staticLine.setAttribute('x1', '2');\n",
        "        staticLine.setAttribute('x2', '26');\n",
        "        staticLine.setAttribute('y1', '8');\n",
        "        staticLine.setAttribute('y2', '8');\n",
        "        staticLine.setAttribute('stroke', lineStroke);\n",
        "        staticLine.setAttribute('stroke-width', lineWidth);\n",
        "        symbolSvg.appendChild(staticLine);\n",
        "      }\n",
        "\n",
        "      if (pathD) {\n",
        "        const staticPath = document.createElementNS('http://www.w3.org/2000/svg', 'path');\n",
        "        staticPath.setAttribute('d', pathD);\n",
        "        staticPath.setAttribute('transform', 'translate(14,8)');\n",
        "        staticPath.setAttribute('fill', fill);\n",
        "        staticPath.setAttribute('fill-opacity', fillOpacity);\n",
        "        staticPath.setAttribute('stroke', stroke);\n",
        "        staticPath.setAttribute('stroke-opacity', strokeOpacity);\n",
        "        staticPath.setAttribute('stroke-width', strokeWidth);\n",
        "        symbolSvg.appendChild(staticPath);\n",
        "      } else {\n",
        "        const staticCircle = document.createElementNS('http://www.w3.org/2000/svg', 'circle');\n",
        "        staticCircle.setAttribute('cx', '14');\n",
        "        staticCircle.setAttribute('cy', '8');\n",
        "        staticCircle.setAttribute('r', '5');\n",
        "        staticCircle.setAttribute('fill', fill);\n",
        "        staticCircle.setAttribute('fill-opacity', fillOpacity);\n",
        "        staticCircle.setAttribute('stroke', stroke);\n",
        "        staticCircle.setAttribute('stroke-opacity', strokeOpacity);\n",
        "        staticCircle.setAttribute('stroke-width', strokeWidth);\n",
        "        symbolSvg.appendChild(staticCircle);\n",
        "      }\n",
        "\n",
        "      const labelSpan = document.createElement('span');\n",
        "      labelSpan.className = 'legend-export-label';\n",
        "      labelSpan.textContent = label;\n",
        "\n",
        "      row.appendChild(symbolSvg);\n",
        "      row.appendChild(labelSpan);\n",
        "      grid.appendChild(row);\n",
        "    });\n",
        "\n",
        "    document.body.appendChild(page);\n",
        "  }, 2500);\n",
        "});\n",
        "</script>\n"
    )

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

reports_to_pdf <- function(input_dir){

    html_files <- list.files(
        path = input_dir,
        pattern = "\\.html?$",
        full.names = TRUE,
        recursive = FALSE,
        ignore.case = TRUE
    )

    html_files <- html_files[!grepl("_QC", html_files)]

    if (!length(html_files) > 0) {
        return(message("No HTML files found in:\n", input_dir))
    }
    chrome_bin <- "/Applications/Google Chrome.app/Contents/MacOS/Google Chrome"
    if (!file.exists(chrome_bin)) {
        stop("Chrome executable does not exist: ", chrome_bin, call. = FALSE)
    }

    for (html_file in html_files) {
        pdf_file <- file.path(input_dir, sanitize_string(html_file, ".pdf"))
        message("Converting file to PDF: ", basename(html_file))

        chrome_args <- c(
            "--headless=new",
            "--no-sandbox",
            "--enable-webgl",
            "--ignore-gpu-blocklist",
            "--use-angle=swiftshader",
            "--use-gl=angle",
            "--window-size=1600,1200",
            "--no-first-run",
            "--no-default-browser-check",
            "--disable-extensions",
            "--disable-background-networking",
            "--disable-sync",
            "--disable-translate",
            "--disable-component-update",
            "--disable-client-side-phishing-detection",
            "--disable-popup-blocking",
            "--metrics-recording-only",
            "--mute-audio"
        )

        work_dir <- file.path(tempdir(), paste0("pagedown_", sanitize_string(html_file)))
        patched_html <- inject_plotly_legend(html_file)
        tryCatch(
            expr = suppressWarnings(pagedown::chrome_print(
                input = patched_html, output = pdf_file, browser = chrome_bin,
                wait = 4, timeout = 120, extra_args = chrome_args,
                work_dir = work_dir, verbose = TRUE, outline = FALSE
            )),
        error = function(e){
            suppressWarnings(pagedown::chrome_print(
                input = patched_html, output = pdf_file, browser = chrome_bin,
                wait = 4, timeout = 120, extra_args = chrome_args,
                work_dir = work_dir, verbose = TRUE, outline = FALSE
            ))
        })
    }
}


# Finds html reports in current directory and uploads to matching record_id
upload_pdf <- function(recordName, input_dir, fld = "classifier_pdf") {
    rcon <- redcapAPI::redcapConnection(url = gb$apiLink, gb$ApiToken)

    message("\n", "Importing Record Report for:", " ", recordName)
    pdf_path <- dir(path = input_dir,
                  pattern = sprintf("^%s.*\\.pdf$", recordName), full.names = TRUE)

    if (length(pdf_path) == 1) {
        message("Uploading file:\n", pdf_path, "\nTo REDCap Record: ", recordName)

        tryCatch(
            expr = redcapAPI::importFiles(rcon, file = pdf_path, record = recordName, field = fld),
            error = function(e) {
                message("REDCap Import error for record: ", recordName)
                message("File failed to upload: ", pdf_path)
            }
        )
    } else{
        message("REDCap Import error for record: ", recordName)
    }
}


import_files_redcap <- function(input_dir) {
    library("redcapAPI")
    csv_file <- dir(input_dir, pattern = "samplesheet.csv", full.names = TRUE)
    if (length(csv_file) == 1) {
        sam_sheet <- as.data.frame(read.csv(csv_file))
    }

    sam_list <- sam_sheet$Sample_Name
    for (sam in sam_list) {
        upload_pdf(sam, input_dir)
    }

}
