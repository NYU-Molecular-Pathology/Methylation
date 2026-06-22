## QC metric computation and failed-report renaming (pipelineHelper.R).
##
## CheckSampleQCmetrics counts how many of the five metrics are within their
## cutoffs (SI<=9.0, BP<=11.0, BS<=10.0, HC<=12.75, NC[Pvalue]<=0.95) and sets
##   qc_passed = "no"  when >= 2 metrics are in range
##   qc_passed = "yes" when <  2 metrics are in range
## The fixture rows have 5, 2, 1, 0 metrics in range (see mk_qc_csv).

test_that("CheckSampleQCmetrics applies the >=2-in-range criteria", {
  rd <- local_run_dir()
  mk_qc_csv(file.path(rd, paste0(TEST_RUN_ID, "_qc_data.csv")))

  res <- suppressWarnings(suppressMessages(gb$CheckSampleQCmetrics(TEST_RUN_ID)))

  expect_true(all(c("record_id", "qc_passed") %in% colnames(res)))
  expect_identical(res$qc_passed[res$record_id == "RD-26-1-test"], "no")   # 5 in range
  expect_identical(res$qc_passed[res$record_id == "RD-26-2-test"], "no")   # 2 in range (boundary)
  expect_identical(res$qc_passed[res$record_id == "RD-26-3-test"], "yes")  # 1 in range
  expect_identical(res$qc_passed[res$record_id == "RD-26-4-test"], "yes")  # 0 in range
})

test_that("RenameFailed renames the HTML of records flagged 'no'", {
  local_run_dir()
  flagged <- "RD-26-2-test"   # qc_passed == "no" -> renamed to _QC_FAILED.html
  other   <- "RD-26-1-test"   # qc_passed == "yes" -> left untouched
  file.create(paste0(flagged, ".html"))
  file.create(paste0(other, ".html"))

  qc <- data.frame(record_id = c(other, flagged),
                   qc_passed = c("yes", "no"),
                   stringsAsFactors = FALSE)

  suppressWarnings(suppressMessages(gb$RenameFailed(qc)))

  expect_true(file.exists(paste0(flagged, "_QC_FAILED.html")))
  expect_false(file.exists(paste0(flagged, ".html")))
  expect_true(file.exists(paste0(other, ".html")))
})
