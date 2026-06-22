## Core integration test: dummy <RUN_ID>.xlsm  ->  samplesheet.csv (for minfi).

test_that("getTotalSamples reads cell B4 of the worksheet", {
  skip_if_not_installed("openxlsx2")
  skip_if_not_installed("readxl")
  rd <- local_run_dir()
  make_dummy_xlsm(file.path(rd, paste0(TEST_RUN_ID, ".xlsm")), TEST_RUN_ID, n = 8L)
  expect_equal(suppressMessages(gb$getTotalSamples()), 8L)
})

test_that("readSheetWrite generates a samplesheet.csv from the dummy xlsm", {
  skip_if_not_installed("openxlsx2")
  skip_if_not_installed("readxl")
  rd <- local_run_dir()
  make_dummy_xlsm(file.path(rd, paste0(TEST_RUN_ID, ".xlsm")), TEST_RUN_ID, n = 8L)
  gb$runID <- TEST_RUN_ID

  suppressWarnings(suppressMessages(
    gb$readSheetWrite(sampleNumb = 8L, runID = TEST_RUN_ID)
  ))

  expect_true(file.exists("samplesheet.csv"))
  ss <- read.csv("samplesheet.csv", stringsAsFactors = FALSE)

  expect_setequal(
    colnames(ss),
    c("Sample_Name", "DNA_Number", "Sentrix_ID", "Sentrix_Position",
      "SentrixID_Pos", "Basename", "RunID", "MP_num", "tech", "tech2", "Date")
  )
  expect_equal(nrow(ss), 8L)
  expect_true(all(ss$RunID == TEST_RUN_ID))
  expect_identical(ss$Sample_Name[1], "control")
  # Sentrix split is consistent: <barcode>_<position>
  expect_true(all(paste0(ss$Sentrix_ID, "_", ss$Sentrix_Position) == ss$SentrixID_Pos))
})

test_that("readSheetWrite does not overwrite an existing samplesheet.csv", {
  skip_if_not_installed("openxlsx2")
  rd <- local_run_dir()
  make_dummy_xlsm(file.path(rd, paste0(TEST_RUN_ID, ".xlsm")), TEST_RUN_ID, n = 8L)
  writeLines("PRE-EXISTING", "samplesheet.csv")
  suppressWarnings(suppressMessages(gb$readSheetWrite(sampleNumb = 8L, runID = TEST_RUN_ID)))
  expect_identical(readLines("samplesheet.csv"), "PRE-EXISTING")
})
