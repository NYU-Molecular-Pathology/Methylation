## Pure data-reshaping helpers (CopyOutput.R, pipelineHelper.R) — no network.

test_that("check_validation appends _VAL only for validation runs", {
  x <- data.frame(run_number = "26-MGDM-VAL", record_id = "RD-26-1", stringsAsFactors = FALSE)
  expect_identical(gb$check_validation(x)$record_id, "RD-26-1_VAL")

  y <- data.frame(run_number = "26-MGDM-Test", record_id = "RD-26-1", stringsAsFactors = FALSE)
  expect_identical(gb$check_validation(y)$record_id, "RD-26-1")
})

test_that("getRunData reshapes a single samplesheet row", {
  local_run_dir()
  mk_samplesheet("samplesheet.csv")                       # control + RD-26-N-test rows
  ss <- read.csv("samplesheet.csv", stringsAsFactors = FALSE)
  single <- ss[ss$Sample_Name != "control", , drop = FALSE][1, , drop = FALSE]

  dat <- suppressMessages(gb$getRunData(single))
  expect_identical(dat$sampleID, single$Sample_Name)
  expect_identical(dat$senLi,    single$SentrixID_Pos)
  expect_identical(dat$outFi,    paste0(single$Sample_Name, ".html"))
  expect_identical(dat$tech,     single$tech)
})

test_that("rename_reorder_output renames and reorders QC/classifier columns", {
  skip_if_not_installed("dplyr")
  input <- data.frame(
    RunID = TEST_RUN_ID, RD.number = "RD-26-1", B.number = "B1", TM.number = "MP-1",
    Log2sqrt.M.U. = 8.0, Log2.M.U. = 0.1, log2sqrt.R.G. = 10.0, log2.R.G. = 0.2,
    BS_log2sqrt.R.G. = 9.0, BS_log2.R.G. = 0.3, log2sqrt.H.L. = 11.0, log2.H.L. = 0.4,
    Pvalue = 0.5, record_id = "RD-26-1", subgroup_score = "0.9",
    classifier_score = "0.95", classifier_value = "GBM", subgroup = "x",
    b_number = "B1", primary_tech = "T1", run_number = TEST_RUN_ID,
    accession_number = "A1", block = "BL1", tm_number = "MP-1", diagnosis = "dx",
    check.names = FALSE, stringsAsFactors = FALSE
  )
  out <- suppressWarnings(suppressMessages(gb$rename_reorder_output(input)))
  expect_identical(colnames(out)[1], "RunID")
  expect_true("RD-number" %in% colnames(out))
  expect_equal(nrow(out), 1L)
})
