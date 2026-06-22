## Priority ordering and run-list filtering (pipelineHelper.R).

test_that("reOrderRun puts selected RDs first and returns all indices", {
  local_run_dir()
  nm <- test_sample_names(4L)            # control, RD-26-1-test, RD-26-2-test, RD-26-3-test
  mk_samplesheet("samplesheet.csv", nm)
  sel <- nm[3]                           # RD-26-2-test (row 3)
  ord <- suppressMessages(gb$reOrderRun(sel))
  expect_equal(ord[1], 3L)
  expect_setequal(ord, seq_along(nm))
})

test_that("reOrderRun returns NULL when no selection is given", {
  local_run_dir()
  mk_samplesheet("samplesheet.csv")
  expect_null(suppressMessages(gb$reOrderRun(NULL)))
})

test_that("GetPriorityCases promotes samples whose MP_num matches BN0", {
  local_run_dir()
  nm <- test_sample_names(3L)            # control, RD-26-1-test, RD-26-2-test
  mk_samplesheet("samplesheet.csv", nm, mp = c("MP-1", "BN0-99", "MP-3"))
  res <- suppressMessages(gb$GetPriorityCases(NULL))
  expect_true(nm[2] %in% res)            # RD-26-1-test (its MP_num holds BN0)
})

test_that("GetPriorityCases returns input unchanged with no BN0 match", {
  local_run_dir()
  mk_samplesheet("samplesheet.csv", test_sample_names(3L))
  expect_null(suppressMessages(gb$GetPriorityCases(NULL)))
})

test_that("BN0 cases are auto-prioritized to the front of the run order", {
  skip_if_not_installed("openxlsx2")
  skip_if_not_installed("readxl")
  rd <- local_run_dir()
  nm <- test_sample_names(8L)                  # control + RD-26-1-test .. RD-26-7-test
  mp <- paste0("MP-26-", seq_len(8L), "-test")
  bn_idx <- 4L                                 # flag one sample as a BN0 priority case
  mp[bn_idx] <- "BN0-123"

  make_dummy_xlsm(file.path(rd, paste0(TEST_RUN_ID, ".xlsm")),
                  run_id = TEST_RUN_ID, n = 8L, mp = mp)
  gb$runID <- TEST_RUN_ID
  suppressWarnings(suppressMessages(
    gb$readSheetWrite(sampleNumb = 8L, runID = TEST_RUN_ID)
  ))

  # GetPriorityCases detects the BN0 sample from the generated samplesheet
  priority <- suppressMessages(gb$GetPriorityCases(NULL))
  expect_true(nm[bn_idx] %in% priority)

  # reOrderRun moves the BN0 sample to the front of the run order
  run_order <- suppressMessages(gb$reOrderRun(priority))
  ss <- read.csv("samplesheet.csv", stringsAsFactors = FALSE)
  expect_identical(ss$Sample_Name[run_order[1]], nm[bn_idx])
})

test_that("getRunList skips samples that already have an HTML report", {
  local_run_dir()
  nm <- test_sample_names(4L)
  data <- data.frame(Sample_Name = nm, stringsAsFactors = FALSE)

  # No HTML present -> nothing skipped
  expect_equal(suppressMessages(gb$getRunList(data, seq_along(nm))), seq_along(nm))

  # Add a report for one sample -> that index is dropped
  done <- nm[3]                          # RD-26-2-test
  file.create(paste0(done, ".html"))
  res <- suppressMessages(gb$getRunList(data, seq_along(nm)))
  expect_false(3L %in% res)
  expect_true(all(c(1L, 2L, 4L) %in% res))
})
