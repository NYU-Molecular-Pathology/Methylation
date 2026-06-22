## RMD rendering / knitting (pipelineHelper.R make_knit_report + generateQCreport).
##
## The default report/QC templates are the paths hard-coded in pipelineHelper.R
## (gb$reportMd, gb$QC_file). Override either via an environment variable to test
## a locally modified template:
##   REPORT_RMD=/path/My_new_report.Rmd
##   QC_RMD=/path/My_new_QC.Rmd
##
## Rendering is heavy (needs mnp.v12epicv2, the template, and EPIC v2 demo data)
## and is OPT-IN: set RUN_RENDER_TESTS=1 to run the two render tests below. All
## sample values are derived from the generated samplesheet.csv.

# --- Configuration checks (always run, no heavy deps) -------------------------

options(warn = -1)

test_that("default report/QC RMD paths match the pipelineHelper hard-coded values", {
  expect_identical(gb$reportMd,
                   "/Volumes/CBioinformatics/Methylation/EPIC_V2_report_2.Rmd")
  expect_identical(basename(gb$QC_file), "Methyl_QC.Rmd")
})

test_that("REPORT_RMD / QC_RMD environment variables override the defaults", {
  withr::local_envvar(REPORT_RMD = "/tmp/custom_report.Rmd",
                      QC_RMD     = "/tmp/custom_QC.Rmd")
  expect_identical(get_report_rmd(), "/tmp/custom_report.Rmd")
  expect_identical(get_qc_rmd(),     "/tmp/custom_QC.Rmd")
})


# --- Report rendering (opt-in) -----------------------------------------------
# Writes straight to stderr() so output is NOT muffled by testthat's reporter
# (message()/print() are captured during tests; a direct cat() to stderr is not).
emit <- function(...) cat(..., "\n", sep = "", file = stderr())

brick_msg <- function(main_msg) {
  border <- base::strrep("~", nchar(main_msg))
  emit("\n", border, "\n", main_msg, "\n", border)
}

log_list_values <- function(list_item, label = "None") {
  brick_msg(paste(label, "values"))
  for (item_name in names(list_item)) {
    value <- list_item[[item_name]]

    emit("\n[x] ", item_name, " <", paste(class(value), collapse = ", "), ">")

    if (is.data.frame(value)) {
      emit("    dimensions: ", nrow(value), " rows ", ncol(value), " columns")
      emit(paste(utils::capture.output(utils::head(value, 6L)), collapse = "\n"))
    } else if (length(value) == 1L) {
      emit("    value: ", as.character(value))
    } else {
      emit("    length: ", length(value))
      emit(paste(utils::capture.output(utils::head(value, 6L)), collapse = "\n"))
    }
  }
}

test_that("report RMD knits one EPIC v2 demo sample to HTML", {
  skip_if(!render_enabled(), "Set RUN_RENDER_TESTS=1 to run rendering tests")
  skip_if_not_installed("rmarkdown")
  skip_if_not_installed("minfi")
  skip_if_not_installed("Biobase")
  skip_if_not_installed("RnBeads")
  skip_if_not_installed("MethylAid")
  skip_if_not_installed("mnp.v12epicv2")

  report_rmd <- get_report_rmd()
  skip_if(!file.exists(report_rmd), paste("Report RMD not accessible:", report_rmd))

  ss <- setup_demo_run(copy_idats = TRUE)
  skip_if(is.null(ss), "EPIC v2 demo data / required packages unavailable")
  source_report_scripts()  # defines gb$SuppInfoTable, gb$GetRedcapDF, etc.
  # Don't auto-open the rendered HTML; just save it for inspection
  withr::local_options(browser = function(url) invisible(url), viewer = NULL)

  # Pick the first real (non-control) sample straight from the samplesheet.
  real <- ss[ss$Sample_Name != "control", , drop = FALSE]
  skip_if(nrow(real) == 0, "No RD samples in the generated samplesheet")
  single <- real[1, , drop = FALSE]

  dat <- suppressMessages(gb$getRunData(single))
  idat_dir <- demo_idat_dir_for(single$SentrixID_Pos)
  rg  <- suppressWarnings(suppressMessages(gb$getRGset(idat_dir, single$SentrixID_Pos)))
  params_init <- list(token = gb$ApiToken, rundata = dat, RGsetEpic = rg, knitDir = getwd())

  log_list_values(params_init, "params_init")
  suppressWarnings(suppressMessages(gb$make_knit_report(dat, report_rmd, params_init)))
  html_file <- file.path(getwd(), dat$outFi)

  expect_true(file.exists(html_file))

  if (identical(Sys.getenv("OPEN_TEST_HTML"), "true")) {
    utils::browseURL(normalizePath(html_file))
  }
  message("\n~~~~~Expected file output:\n", dat$outFi)
})

# --- QC rendering (opt-in) ---------------------------------------------------

test_that("QC RMD knits to <RUN_ID>_QC.html", {
  skip_if(!render_enabled(), "Set RUN_RENDER_TESTS=1 to run rendering tests")
  skip_if_not_installed("rmarkdown")
  skip_if_not_installed("mnp.v12epicv2")

  qc_rmd <- get_qc_rmd()
  skip_if(!file.exists(qc_rmd), paste("QC RMD not accessible:", qc_rmd))

  ss <- setup_demo_run(copy_idats = TRUE)
  skip_if(is.null(ss), "EPIC v2 demo data / required packages unavailable")
  source_report_scripts()  # report helper functions used by the QC RMD
  # Don't auto-open the rendered HTML; just save it for inspection
  withr::local_options(browser = function(url) invisible(url), viewer = NULL)

  gb$QC_file <- qc_rmd
  suppressWarnings(suppressMessages(gb$generateQCreport(TEST_RUN_ID)))
  expect_true(length(dir(getwd(), pattern = "_QC\\.html$")) >= 1)
})
