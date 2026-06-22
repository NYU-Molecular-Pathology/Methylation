#!/usr/bin/env Rscript
## prepare_test_data.R — used by `make test-data`.
## 1. Locate (or download) the Illumina EPIC v2 demo idats.
## 2. Generate a demo <RUN_ID>.xlsm worksheet in the Test_Data directory, built
##    from the real sentrix ids of the demo BeadChips (DEMO_BARCODES).

suppressWarnings({
  source(file.path("testthat", "helper-fixtures.R"))       # TEST_RUN_ID
  source(file.path("testthat", "helper-make_dummy_xlsm.R")) # make_dummy_xlsm()
  source(file.path("testthat", "helper-idats.R"))          # demo data helpers
})

data_dir <- test_data_dir()
cat("TEST_DATA_DIR: ", data_dir, "\n", sep = "")

dir <- ensure_epicv2_demo(data_dir)
if (is.null(dir)) {
  stop(
    "Could not locate or download the EPIC v2 demo idats.\n",
    "Looked in: ", data_dir, "\n",
    "Point TEST_DATA_DIR at the folder containing DemoDataEPICv2/, e.g.:\n",
    "  export TEST_DATA_DIR=\"$HOME/Documents/Development/Test_Data\"\n",
    "or run `unset TEST_DATA_DIR` to use the default next to the repo.",
    call. = FALSE
  )
}

samples <- demo_samples(dir)
if (nrow(samples) == 0) stop("No idat files detected under: ", dir, call. = FALSE)

cat(sprintf("\nDetected %d sample(s) across barcode(s): %s\n",
            nrow(samples), paste(unique(samples$barcode), collapse = ", ")))
print(samples[, c("barcode", "sentrix")], row.names = FALSE)

xlsm_path <- file.path(data_dir, paste0(TEST_RUN_ID, ".xlsm"))
make_dummy_xlsm(xlsm_path, run_id = TEST_RUN_ID, sentrix = samples$sentrix)

cat("\nGenerated demo worksheet:\n  ", xlsm_path, "\n", sep = "")
cat("Test data ready.\n")
