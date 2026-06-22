## REDCap-touching logic exercised WITHOUT a live server.
## The early-return path needs no network; the export path is mocked.

test_that("checkRedcapRecord short-circuits for validation-control records", {
  expect_false(suppressMessages(gb$checkRedcapRecord("26_control_sample")))
})

test_that("checkRedcapRecord returns is.na(field) via mocked redcapAPI", {
  skip_if_not_installed("mockery")
  skip_if_not_installed("redcapAPI")

  fn <- gb$checkRedcapRecord  # local copy to stub

  mockery::stub(fn, "redcapAPI::redcapConnection", function(...) NULL)
  mockery::stub(fn, "redcapAPI::exportRecordsTyped",
                function(...) data.frame(classifier_pdf = NA, stringsAsFactors = FALSE))
  expect_true(suppressMessages(fn("RD-26-1")))

  mockery::stub(fn, "redcapAPI::exportRecordsTyped",
                function(...) data.frame(classifier_pdf = "report.pdf", stringsAsFactors = FALSE))
  expect_false(suppressMessages(fn("RD-26-1")))
})
