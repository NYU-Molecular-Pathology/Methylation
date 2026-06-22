#!/usr/bin/env Rscript
## ---------------------------
## Script name: run_tests.R
## Purpose: Run the Methylation pipeline unit/integration test suite (testthat)
## Usage:    Rscript Testing/run_tests.R        # from the repo root
##           # or in RStudio: setwd(<repo root>); source("Testing/run_tests.R")
## Notes:    Tests are offline — no /Volumes mounts, REDCap, or patient data.
##           Tests needing optional packages (minfi, minfiDataEPIC, mockery, ...)
##           SKIP rather than fail when those packages are absent.
## ---------------------------

suppressWarnings(suppressPackageStartupMessages(library("testthat")))

muffle_test_warning <- function(expr) {
  withCallingHandlers(expr, warning = function(w) {
    warn_msg <- conditionMessage(w)
    if (
      grepl("was built under R version", warn_msg, fixed = TRUE) ||
      grepl("replacing previous import", warn_msg, fixed = TRUE)
    ) invokeRestart("muffleWarning")
  })
}

muffle_test_warning({
  testthat::test_dir("testthat", reporter = "progress", stop_on_failure = TRUE, stop_on_warning = FALSE)
})


skip_if_not_installed <- function(pkg, minimum_version = NULL) {
  tryCatch(find.package(pkg), error = function(e)
    testthat::skip(paste0("{", pkg, "} is not installed")))
  if (!suppressPackageStartupMessages(suppressWarnings(requireNamespace(pkg, quietly = TRUE))))
    testthat::skip(paste0("{", pkg, "} cannot be loaded"))
  if (!is.null(minimum_version)) {
    installed_version <- package_version(pkg)
    if (installed_version < minimum_version) {
      testthat::skip(paste0("Installed ", pkg, " is version ", installed_version, "; ", "but ",
                            minimum_version, " is required"))
    }
  }
  invisible()
}

unlockBinding("skip_if_not_installed", asNamespace("testthat"))

assign("skip_if_not_installed", skip_if_not_installed, envir = asNamespace("testthat"))

lockBinding("skip_if_not_installed", asNamespace("testthat"))

# Resolve the path to Testing/testthat whether run via Rscript or source() ------
args_all <- commandArgs(FALSE)
file_arg <- grep("^--file=", args_all, value = TRUE)
this_dir <- if (length(file_arg) == 1) {
  dirname(normalizePath(sub("^--file=", "", file_arg)))
} else if (!is.null(sys.frames()[[1]]$ofile)) {
  dirname(normalizePath(sys.frames()[[1]]$ofile))
} else {
  "Testing"
}

test_dir_path <- file.path(this_dir, "testthat")
if (!dir.exists(test_dir_path)) test_dir_path <- file.path("Testing", "testthat")
stopifnot("Cannot locate Testing/testthat" = dir.exists(test_dir_path))

cat("Running tests in:", normalizePath(test_dir_path), "\n\n")

res <- test_dir(test_dir_path, reporter = "summary", stop_on_failure = FALSE)

df <- as.data.frame(res)
n_failed  <- sum(df$failed)
n_error   <- sum(as.integer(df$error))
n_skipped <- sum(df$skipped)
n_warn    <- sum(df$warning)

cat(sprintf(
  "\n==== RESULT: %d failed, %d errors, %d skipped, %d warnings (%d test blocks) ====\n",
  n_failed, n_error, n_skipped, n_warn, nrow(df)
))

if (!interactive()) quit(status = if ((n_failed + n_error) > 0) 1L else 0L)
