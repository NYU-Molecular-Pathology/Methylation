## helper-fixtures.R — shared constants, temp run-dir builders, RMD config.

## RUN_ID used across the suite (also the worksheet / run-directory name).
TEST_RUN_ID <- "26-MGDM-Test"

## Single source of truth for sample names, matching the worksheet builder
## (make_dummy_xlsm): first sample is the control, the rest are RD-26-<i>-test.
test_sample_names <- function(n = 8L) {
  c("control", paste0("RD-26-", seq_len(n - 1L), "-test"))
}

## Where persisted test outputs go for inspection: <repo>/Testing/test_outputs
## (add `Testing/test_outputs/` to .gitignore).
test_outputs_dir <- function() {
  root <- if (exists("REPO_ROOT", envir = gb)) gb$REPO_ROOT else getwd()
  file.path(root, "Testing", "test_outputs")
}

## Run outputs (samplesheet.csv, rendered HTML, ...) are kept on disk under
## Testing/test_outputs for inspection by default. Set KEEP_OUTPUTS=0 (or false/no)
## to use a throwaway temp dir instead.
keep_outputs <- function() {
  !tolower(Sys.getenv("KEEP_OUTPUTS", unset = "1")) %in% c("0", "false", "no")
}

## Create a <base>/<RUN_ID> run directory, setwd() into it for the duration of the
## calling test, and set the gb globals the pipeline relies on. With persist =
## TRUE the run dir is the single, predictable Testing/test_outputs/<RUN_ID> and is
## kept after the test (cleared once per R session, so tests within a session - e.g.
## the report + QC renders - accumulate into the same folder); otherwise it is a
## throwaway temp dir.
local_run_dir <- function(run_id = TEST_RUN_ID, env = parent.frame(),
                          persist = FALSE) {
  if (isTRUE(persist)) {
    base <- test_outputs_dir()
    rd <- file.path(base, run_id)
    # Clear once per session so stale files from prior runs don't linger
    cleared <- get0(".test_outputs_cleared", envir = gb, ifnotfound = character(0))
    if (!run_id %in% cleared) {
      if (dir.exists(rd)) unlink(rd, recursive = TRUE)
      assign(".test_outputs_cleared", union(cleared, run_id), envir = gb)
    }
    dir.create(rd, recursive = TRUE, showWarnings = FALSE)
    withr::local_dir(rd, .local_envir = env)  # cwd restored at test end; dir kept
    message("Test outputs kept at: ", rd)
  } else {
    base <- withr::local_tempdir(.local_envir = env)
    rd <- file.path(base, run_id)
    dir.create(rd, recursive = TRUE, showWarnings = FALSE)
    withr::local_dir(rd, .local_envir = env)
  }
  gb$runID      <- run_id
  gb$methDir    <- base
  gb$baseDir    <- base
  gb$baseFolder <- base
  rd
}

## Minimal minfi-style samplesheet.csv (the central hand-off artifact).
mk_samplesheet <- function(path = "samplesheet.csv",
                           sample_names = test_sample_names(3L),
                           mp = NULL) {
  n <- length(sample_names)
  if (is.null(mp)) mp <- paste0("MP-", seq_len(n))
  barcode <- sprintf("204%09d", seq_len(n))
  pos <- sprintf("R0%dC01", ((seq_len(n) - 1) %% 8) + 1)
  df <- data.frame(
    Sample_Name      = sample_names,
    DNA_Number       = paste0("B", seq_len(n)),
    Sentrix_ID       = barcode,
    Sentrix_Position = pos,
    SentrixID_Pos    = paste0(barcode, "_", pos),
    Basename         = "basename",
    RunID            = TEST_RUN_ID,
    MP_num           = mp,
    tech             = "T1",
    tech2            = "T2",
    Date             = "2026-01-01",
    stringsAsFactors = FALSE
  )
  utils::write.csv(df, path, row.names = FALSE, quote = FALSE)
  invisible(df)
}

## <RUN_ID>_qc_data.csv (the per-sample metrics the QC RMD emits). Rows are
## engineered so the number of in-range metrics (cutoffs: SI<=9.0, BP<=11.0,
## BS<=10.0, HC<=12.75, NC[Pvalue]<=0.95) is exactly 5, 2, 1, 0 respectively.
## Per CheckSampleQCmetrics: qc_passed == "no" when >= 2 metrics are in range,
## "yes" when < 2 are in range.
mk_qc_csv <- function(path, run_id = TEST_RUN_ID) {
  nm <- test_sample_names(5L)  # control + RD-26-1-test .. RD-26-4-test
  df <- data.frame(
    RD.number        = c(paste0(run_id, "_control"), nm[-1]),
    #                    control  rd1   rd2(2 in)  rd3(1 in)  rd4(0 in)
    Log2sqrt.M.U.    = c(8.0,     8.0,  8.0,       8.0,       12.0),  # <= 9.0
    log2sqrt.R.G.    = c(10.0,    10.0, 10.0,      13.0,      13.0),  # <= 11.0
    BS_log2sqrt.R.G. = c(9.0,     9.0,  13.0,      13.0,      13.0),  # <= 10.0
    log2sqrt.H.L.    = c(11.0,    11.0, 14.0,      14.0,      14.0),  # <= 12.75
    Pvalue           = c(0.50,    0.50, 0.99,      0.99,      0.99),  # <= 0.95
    check.names = FALSE, stringsAsFactors = FALSE
  )
  utils::write.csv(df, path, row.names = FALSE, quote = FALSE)
  invisible(df)
}

## --- RMD render configuration ------------------------------------------------
## Rendering tests default to the paths hard-coded in pipelineHelper.R
## (gb$reportMd, gb$QC_file). Override either with an environment variable to
## test a locally modified RMD, e.g.:
##   REPORT_RMD=/path/to/My_new_report.Rmd
##   QC_RMD=/path/to/My_new_QC.Rmd
get_report_rmd <- function() {
  v <- Sys.getenv("REPORT_RMD", unset = "")
  if (nzchar(v)) return(v)
  if (!is.null(gb$reportMd) && nzchar(gb$reportMd)) return(gb$reportMd)
  "/Volumes/CBioinformatics/Methylation/EPIC_V2_report_2.Rmd"
}

get_qc_rmd <- function() {
  v <- Sys.getenv("QC_RMD", unset = "")
  if (nzchar(v)) return(v)
  if (!is.null(gb$QC_file) && nzchar(gb$QC_file)) return(path.expand(gb$QC_file))
  path.expand("~/Methyl_QC.Rmd")
}

## Heavy rendering tests are opt-in: set RUN_RENDER_TESTS=1 to enable them.
render_enabled <- function() {
  tolower(Sys.getenv("RUN_RENDER_TESTS", unset = "")) %in% c("1", "true", "yes")
}

## Sources the report helper scripts into gb, as the production pipeline does, so
## the report/QC RMDs can call gb$SuppInfoTable(), gb$GetClassProbTables(),
## gb$GetRedcapDF(), etc. Used only by the opt-in render tests (ClassTables.R
## loads heavy report packages at source time).
source_report_scripts <- function() {
  rs_dir <- file.path(gb$R_DIR, "Report-Scripts")
  for (s in c("ClassTables.R", "MLH1_Functions.R", "RedcapOutput.R")) {
    path <- file.path(rs_dir, s)
    if (file.exists(path)) {
      try(suppressWarnings(sys.source(path, envir = gb)), silent = TRUE)
    }
  }
}
