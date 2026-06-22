## setup.R — runs once before the test files.
## Sources ONLY the pipeline scripts that are definition-only (no install /
## checkMounts / report-package loading at source time), in the same order
## methylExpress.R uses, so shadowed functions resolve to the production copy.

gb <- globalenv()
assign("gb", gb)
options(stringsAsFactors = FALSE)

# --- Locate the repo root (folder containing R/pipelineHelper.R) --------------
find_repo_root <- function(start = getwd()) {
  d <- normalizePath(start, mustWork = FALSE)
  for (i in 1:8) {
    if (file.exists(file.path(d, "R", "pipelineHelper.R"))) return(d)
    parent <- dirname(d)
    if (identical(parent, d)) break
    d <- parent
  }
  tp <- tryCatch(testthat::test_path(), error = function(e) NA_character_)
  if (!is.na(tp)) {
    cand <- normalizePath(file.path(tp, "..", ".."), mustWork = FALSE)
    if (file.exists(file.path(cand, "R", "pipelineHelper.R"))) return(cand)
  }
  stop("setup.R: could not locate repo root containing R/pipelineHelper.R")
}

REPO_ROOT <- find_repo_root()
R_DIR <- file.path(REPO_ROOT, "R")
assign("REPO_ROOT", REPO_ROOT, envir = gb)
assign("R_DIR", R_DIR, envir = gb)

# --- Source the safe (definition-only) scripts, in production order -----------
safe_scripts <- c(
  "SetRunParams.R",
  "MakeSampleSheet.R",
  "CopyInputs.R",
  "CopyOutput.R",
  "pipelineHelper.R",
  "CustomRuns.R",
  file.path("Report-Scripts", "MLH1_Functions.R")
)
for (s in safe_scripts) {
  sys.source(file.path(R_DIR, s), envir = gb)
}

# convert_html_pdf.R installs pagedown at source time; only source it (for
# sanitize_string) when pagedown is already present, so tests never trigger an
# install. Tests that use sanitize_string skip if it is unavailable.
if (requireNamespace("pagedown", quietly = TRUE)) {
  try(sys.source(file.path(R_DIR, "convert_html_pdf.R"), envir = gb), silent = TRUE)
}

# --- Minimal global state the functions expect -------------------------------
gb$ApiToken <- gb$token <- "FAKE_TOKEN_0000000000000000000000000000"
gb$apiLink  <- "https://example.invalid/redcap/api/"  # neutralize any REDCap call
gb$redcapUp <- FALSE

# Pin the production report/QC template defaults. methylExpress.R sets reportMd
# last (after CustomRuns.R reassigns gb$reportMd at source time), so mirror that
# here. gb$QC_file is set by pipelineHelper.R (~/Methyl_QC.Rmd); keep it.
gb$reportMd <- "/Volumes/CBioinformatics/Methylation/EPIC_V2_report_2.Rmd"
if (is.null(gb$QC_file) || !nzchar(gb$QC_file)) {
  gb$QC_file <- file.path(path.expand("~"), "Methyl_QC.Rmd")
}
