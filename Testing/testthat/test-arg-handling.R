## Command-line argument validation / normalization (SetRunParams.R).
## These cover the input-simplification refactor (RUN_ID-only entry contract).

test_that("CheckInputArg maps NULL-like values to the default", {
  token <- "NULL"
  suppressMessages(gb$CheckInputArg(token, gb))
  expect_null(gb$token)

  selectRDs <- "NULL"
  suppressMessages(gb$CheckInputArg(selectRDs, gb))
  expect_null(gb$selectRDs)
})

test_that("CheckInputArg preserves a real value and assigns into gb", {
  runID <- "26-MGDM-Test"
  suppressMessages(gb$CheckInputArg(runID, gb))
  expect_identical(gb$runID, "26-MGDM-Test")
})

test_that("CheckInputArg applies the supplied default", {
  redcapUp <- "NULL"
  suppressMessages(gb$CheckInputArg(redcapUp, gb, FALSE))
  expect_false(gb$redcapUp)
})

test_that("CheckInputArg currently keeps 'T'/'F' as STRINGS (pre-refactor baseline)", {
  ## NOTE: current code does NOT coerce "T"->TRUE. After the planned T/F
  ## refactor, change this to expect_true(isTRUE(gb$redcapUp)).
  redcapUp <- "T"
  suppressMessages(gb$CheckInputArg(redcapUp, gb, FALSE))
  expect_identical(gb$redcapUp, "T")
})

test_that("CheckBaseFolderInput accepts existing dirs and nullifies the rest", {
  d <- withr::local_tempdir()
  expect_identical(gb$CheckBaseFolderInput(d), d)
  expect_null(gb$CheckBaseFolderInput("NULL"))
  expect_null(gb$CheckBaseFolderInput(NULL))
  expect_error(gb$CheckBaseFolderInput(file.path(tempdir(), "no_such_dir_xyz_123")))
})
