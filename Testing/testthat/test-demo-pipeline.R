## End-to-end on the real EPIC v2 demo data (two BeadChips):
##   worksheet (.xlsm) -> samplesheet.csv -> idat path resolution.
## All sample values are read from the generated samplesheet.csv. SKIPS when the
## demo data or required packages are unavailable.

test_that("demo_samples finds complete EPIC v2 demo samples on disk", {
  data_dir <- ensure_epicv2_demo()
  skip_if(is.null(data_dir), "EPIC v2 demo data unavailable")
  s <- demo_samples(data_dir)
  skip_if(nrow(s) == 0, "No demo samples found")
  expect_gte(nrow(s), 8L)
  expect_true(all(file.exists(s$grn)))
  expect_true(all(file.exists(s$red)))
})

test_that("readSheetWrite builds a samplesheet matching the demo sentrix ids", {
  skip_if_not_installed("openxlsx2")
  skip_if_not_installed("readxl")

  ss <- setup_demo_run()
  skip_if(is.null(ss), "EPIC v2 demo data / required packages unavailable")

  samples <- demo_samples()
  expect_equal(nrow(ss), nrow(samples))
  expect_setequal(ss$SentrixID_Pos, samples$sentrix)
  expect_identical(ss$Sample_Name[1], "control")
  expect_true(all(ss$RunID == TEST_RUN_ID))
})

test_that("getAllFiles resolves the generated samplesheet to idats that exist", {
  skip_if_not_installed("openxlsx2")
  skip_if_not_installed("readxl")

  ss <- setup_demo_run()
  skip_if(is.null(ss), "EPIC v2 demo data / required packages unavailable")

  idat_root <- demo_samples()$idat_root[1]
  allFi <- suppressMessages(
    gb$getAllFiles(idatDir = idat_root, csvNam = "samplesheet.csv")
  )
  expect_equal(length(allFi), 2L * nrow(ss))
  expect_true(all(file.exists(allFi)))
})
