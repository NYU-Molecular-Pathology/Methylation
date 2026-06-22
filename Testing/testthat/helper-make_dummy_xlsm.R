## helper-make_dummy_xlsm.R — build a fake <RUN_ID>.xlsm worksheet for parsing.
##
## Written with openxlsx2, which produces valid macro-enabled (.xlsm) workbooks
## directly (Excel opens them cleanly; the pipeline reads them via readxl).
##
## Mirrors the real lab template (blank_template.xlsm):
##   - Sheet "worksheet" (1st sheet): cell B4 = total samples, F4 = run date
##   - Sheet "raw_labels" (matches grep "raw"): header columns in the exact order
##       1 Sentrix_ID | 2 Slide | 3 Project | 4 Array | 5 Sample_Plate |
##       6 Sample_Well | 7 Basename | 8 MP_number | 9 Sample_Name | 10 b_number |
##       11 Tech | 12 Tech2 | 13 Notes | 14 Batch | 15 FileName
##     (checkHeaders() keys off positions 1, 8, 9, 10)
##   - Sheet "REDCap_Import" (matches grep "REDCap"): record cross-check sheet
##
## `sentrix` (optional): a vector of real "<barcode>_RxxCxx" ids. When supplied,
## the worksheet is built from those (n is derived); otherwise synthetic ids are
## generated. Use the real demo sentrix to drive an end-to-end pipeline test.

## `mp` (optional): a vector of MP numbers (recycled to length n). Use it to
## inject "BN0" values for testing the BN run-order prioritization.
make_dummy_xlsm <- function(path, run_id = TEST_RUN_ID, n = 8L,
                            sentrix = NULL, mp = NULL) {
  stopifnot(requireNamespace("openxlsx2", quietly = TRUE))

  if (!is.null(sentrix)) {
    n <- length(sentrix)
    barcode <- sub("_.*$", "", sentrix)
    pos     <- sub("^[^_]*_", "", sentrix)
  } else {
    n <- as.integer(n)
    barcode <- sprintf("20689111%04d", seq_len(n))
    pos     <- sprintf("R0%dC01", ((seq_len(n) - 1L) %% 8L) + 1L)
    sentrix <- paste0(barcode, "_", pos)
  }
  stopifnot(n >= 1L)
  sample_names <- c("control", paste0("RD-26-", seq_len(n - 1L), "-test"))
  mp_nums <- if (is.null(mp)) paste0("MP-26-", seq_len(n), "-test") else rep_len(mp, n)
  b_nums  <- paste0("B26-", seq_len(n), "-test")

  raw <- data.frame(
    Sentrix_ID   = sentrix,     # 1
    Slide        = barcode,     # 2
    Project      = run_id,      # 3
    Array        = pos,         # 4
    Sample_Plate = "",          # 5
    Sample_Well  = "",          # 6
    Basename     = "",          # 7
    MP_number    = mp_nums,     # 8
    Sample_Name  = sample_names,# 9
    b_number     = b_nums,      # 10
    Tech         = "T1",        # 11
    Tech2        = "T2",        # 12
    Notes        = "",          # 13
    Batch        = run_id,      # 14
    FileName     = "",          # 15
    check.names = FALSE, stringsAsFactors = FALSE
  )

  redcap <- data.frame(
    record_id              = sample_names,
    tm_number              = mp_nums,
    b_number               = b_nums,
    barcode_and_row_column = sentrix,
    run_number             = run_id,
    check.names = FALSE, stringsAsFactors = FALSE
  )

  wb <- openxlsx2::wb_workbook()

  wb <- openxlsx2::wb_add_worksheet(wb, "worksheet")
  wb <- openxlsx2::wb_add_data(wb, "worksheet", "Total # samples:", dims = "A4", col_names = FALSE)
  wb <- openxlsx2::wb_add_data(wb, "worksheet", n,                  dims = "B4", col_names = FALSE) # total
  wb <- openxlsx2::wb_add_data(wb, "worksheet", "Start Date:",      dims = "E4", col_names = FALSE)
  wb <- openxlsx2::wb_add_data(wb, "worksheet", "2026-01-01",       dims = "F4", col_names = FALSE) # date

  wb <- openxlsx2::wb_add_worksheet(wb, "raw_labels")
  wb <- openxlsx2::wb_add_data(wb, "raw_labels", raw)

  wb <- openxlsx2::wb_add_worksheet(wb, "REDCap_Import")
  wb <- openxlsx2::wb_add_data(wb, "REDCap_Import", redcap)

  wb <- openxlsx2::wb_add_worksheet(wb, "EPIC")
  wb <- openxlsx2::wb_add_worksheet(wb, "DNA_Input")

  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  openxlsx2::wb_save(wb, file = path, overwrite = TRUE)
  invisible(path)
}
