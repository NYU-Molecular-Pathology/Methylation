## helper-idats.R — EPIC v2 demo idat data (Illumina DemoDataEPIC_v2), offline.
##
## The pipeline's classifier (mnp.v12epicv2) processes both EPIC v1 and v2; the
## tests use only the public Illumina EPIC v2 demo dataset. Idats are expected at
##   <TEST_DATA_DIR>/.../<barcode>/<barcode>_RxxCxx_{Grn,Red}.idat
## TEST_DATA_DIR defaults to "<repo parent>/Test_Data" and can be overridden with
## the TEST_DATA_DIR environment variable. If no idats are present, the dataset is
## downloaded and unzipped on first use; tests SKIP if that is not possible.
##
## By default only two BeadChips are used to keep runs fast (the demo ships four
## arrays; these two together provide eight samples).

DEMO_EPICV2_URL <-
  "https://support.illumina.com/content/dam/illumina-support/documents/downloads/productfiles/methylationepic/DemoDataEPIC_v2.zip"

## BeadChips (barcodes) used by the tests.
DEMO_BARCODES <- c("206891110001", "206891110002")

test_data_dir <- function() {
  v <- Sys.getenv("TEST_DATA_DIR", unset = "")
  if (nzchar(v)) return(v)
  # Walk up to the repo root (folder with R/pipelineHelper.R); demo data lives in
  # a sibling "Test_Data" directory. Works whether called from Testing/ or testthat/.
  d <- normalizePath(getwd(), mustWork = FALSE)
  for (i in 1:8) {
    if (file.exists(file.path(d, "R", "pipelineHelper.R"))) {
      return(file.path(dirname(d), "Test_Data"))
    }
    parent <- dirname(d)
    if (identical(parent, d)) break
    d <- parent
  }
  if (exists("gb") && exists("REPO_ROOT", envir = gb)) {
    return(file.path(dirname(gb$REPO_ROOT), "Test_Data"))
  }
  file.path(dirname(getwd()), "Test_Data")
}

## All Grn idats found under `dir` (recursive), optionally limited to `barcodes`.
list_demo_grn <- function(dir = test_data_dir(), barcodes = NULL) {
  if (!dir.exists(dir)) return(character(0))
  grn <- list.files(dir, pattern = "_Grn\\.idat$", recursive = TRUE, full.names = TRUE)
  if (!is.null(barcodes) && length(grn) > 0) {
    keep <- vapply(basename(grn), function(b) any(startsWith(b, barcodes)), logical(1))
    grn <- grn[keep]
  }
  grn
}

## Ensure the EPIC v2 demo data exists locally; download + unzip on first use.
## Returns the data directory, or NULL if it could not be made available.
ensure_epicv2_demo <- function(dir = test_data_dir()) {
  if (length(list_demo_grn(dir)) > 0) return(dir)
  ok <- tryCatch({
    dir.create(dir, recursive = TRUE, showWarnings = FALSE)
    zip_path <- file.path(tempdir(), "DemoDataEPIC_v2.zip")
    utils::download.file(DEMO_EPICV2_URL, zip_path, mode = "wb", quiet = TRUE)
    utils::unzip(zip_path, exdir = dir)
    TRUE
  }, error = function(e) FALSE)
  if (isTRUE(ok) && length(list_demo_grn(dir)) > 0) dir else NULL
}

## Data frame of complete (Grn + Red) samples for the chosen barcodes:
##   sentrix, barcode, grn, red, idat_root (folder holding the <barcode>/ dirs).
## Falls back to all available barcodes if the requested ones are not present.
demo_samples <- function(dir = test_data_dir(), barcodes = DEMO_BARCODES) {
  grn <- list_demo_grn(dir, barcodes)
  if (length(grn) == 0) grn <- list_demo_grn(dir)
  if (length(grn) == 0) return(data.frame())
  red <- sub("_Grn\\.idat$", "_Red.idat", grn)
  keep <- file.exists(red)
  grn <- grn[keep]; red <- red[keep]
  sentrix <- sub("_Grn\\.idat$", "", basename(grn))
  ord <- order(sentrix)
  data.frame(
    sentrix   = sentrix[ord],
    barcode   = sub("_.*$", "", sentrix[ord]),
    grn       = grn[ord],
    red       = red[ord],
    idat_root = dirname(dirname(grn[ord])),
    stringsAsFactors = FALSE
  )
}

## One complete sample (the first under the chosen barcodes), as a list.
pick_one_epicv2 <- function(dir = test_data_dir()) {
  s <- demo_samples(dir)
  if (nrow(s) == 0) return(NULL)
  list(grn = s$grn[1], red = s$red[1], sentrix = s$sentrix[1],
       barcode = s$barcode[1], idat_root = s$idat_root[1])
}

## Copy one sample's Grn + Red idats (flat) into `dest`; return its sentrix.
copy_one_epicv2_sample <- function(dest, dir = test_data_dir()) {
  s <- pick_one_epicv2(dir)
  if (is.null(s)) return(NULL)
  file.copy(c(s$grn, s$red), dest, overwrite = TRUE)
  s$sentrix
}

## Copy all idats for the chosen barcodes (flat) into `dest`; return the sentrix
## vector (one per sample).
copy_demo_samples <- function(dest, dir = test_data_dir(), barcodes = DEMO_BARCODES) {
  s <- demo_samples(dir, barcodes)
  if (nrow(s) == 0) return(character(0))
  file.copy(c(s$grn, s$red), dest, overwrite = TRUE)
  s$sentrix
}

## Build a demo <RUN_ID>.xlsm at `path`, populated with the real sentrix ids of
## the chosen demo barcodes (delegates to make_dummy_xlsm()).
make_demo_xlsm <- function(path, dir = test_data_dir(), run_id = TEST_RUN_ID,
                           barcodes = DEMO_BARCODES) {
  s <- demo_samples(dir, barcodes)
  if (nrow(s) == 0) stop("No EPIC v2 demo samples found under: ", dir)
  make_dummy_xlsm(path, run_id = run_id, sentrix = s$sentrix)
  invisible(path)
}

## Return the directory holding a given sentrix's idats (the <barcode>/ folder),
## so getRGset() can read them in place without copying. NULL if not found.
demo_idat_dir_for <- function(sentrix, dir = test_data_dir()) {
  s <- demo_samples(dir)
  i <- match(sentrix, s$sentrix)
  if (length(i) == 0 || is.na(i)) return(NULL)
  dirname(s$grn[i])
}

## End-to-end demo setup: generate <RUN_ID>.xlsm from the real demo sentrix ids,
## run the pipeline's readSheetWrite() to produce samplesheet.csv, and return the
## parsed samplesheet (data.frame). Tests then derive ALL sample names / sentrix
## from this output. The working directory is left as the run dir for the test.
## Returns NULL when the demo data or required packages are unavailable.
## `copy_idats = TRUE` also copies the demo idats (flat) into the run directory,
## so functions that read idats from the run dir (e.g. the QC RMD via MethylAid /
## minfi::read.metharray) find them. Off by default to keep light tests fast.
setup_demo_run <- function(dir = test_data_dir(), env = parent.frame(),
                           copy_idats = FALSE) {
  if (!requireNamespace("openxlsx2", quietly = TRUE)) return(NULL)
  if (!requireNamespace("readxl", quietly = TRUE)) return(NULL)
  data_dir <- ensure_epicv2_demo(dir)
  if (is.null(data_dir)) return(NULL)
  s <- demo_samples(data_dir)
  if (nrow(s) < 8) return(NULL)  # readSheetWrite requires >= 8 samples
  rd <- local_run_dir(env = env, persist = keep_outputs())
  make_dummy_xlsm(file.path(rd, paste0(TEST_RUN_ID, ".xlsm")),
                  run_id = TEST_RUN_ID, sentrix = s$sentrix)
  gb$runID <- TEST_RUN_ID
  suppressWarnings(suppressMessages(
    gb$readSheetWrite(sampleNumb = nrow(s), runID = TEST_RUN_ID)
  ))
  if (!file.exists("samplesheet.csv")) return(NULL)
  if (isTRUE(copy_idats)) {
    file.copy(c(s$grn, s$red), rd, overwrite = TRUE)
  }
  utils::read.csv("samplesheet.csv", stringsAsFactors = FALSE)
}

## Minimal samplesheet.csv describing a single demo sentrix, for getAllFiles()
## / get_idats() path resolution against the demo idat root.
write_demo_samplesheet <- function(path, sentrix, sample_name = "RD-26-1") {
  barcode <- sub("_.*$", "", sentrix)
  pos <- sub("^[^_]*_", "", sentrix)
  df <- data.frame(
    Sample_Name = sample_name, DNA_Number = "B1",
    Sentrix_ID = barcode, Sentrix_Position = pos, SentrixID_Pos = sentrix,
    Basename = "bn", RunID = TEST_RUN_ID, MP_num = "MP-1",
    tech = "T1", tech2 = "T2", Date = "2026-01-01",
    stringsAsFactors = FALSE
  )
  utils::write.csv(df, path, row.names = FALSE, quote = FALSE)
  invisible(df)
}
