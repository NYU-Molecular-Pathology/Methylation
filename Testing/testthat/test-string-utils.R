## Pure string / logic utilities — fast, no fixtures, no network.

test_that("grabYear derives the 4-digit year from a RUN_ID", {
  expect_identical(suppressMessages(gb$grabYear("26-MGDM-Test")), "2026")
  expect_identical(suppressMessages(gb$grabYear("21-MGDM17")), "2021")
  expect_identical(suppressMessages(gb$grabYear(" 22-MR-5 ")), "2022")  # spaces stripped
})

test_that("get_next_rd increments the highest suffix", {
  expect_identical(gb$get_next_rd(c("RD-26-1", "RD-26-2", "RD-26-3")), "RD-26-4")
  expect_identical(gb$get_next_rd(c("RD-26-9", "RD-26-2")), "RD-26-10")
})

test_that("get_next_rd returns NA when no numeric suffix exists", {
  expect_true(is.na(gb$get_next_rd(c("abc", "def"))))
})

test_that("get_newest_rd selects the newest RD across years", {
  expect_identical(
    gb$get_newest_rd(c("RD-25-100", "RD-26-1", "RD-26-50")),
    "RD-26-50"
  )
})

test_that("is_present treats NULL/NA/empty/zero as absent", {
  expect_false(gb$is_present(NULL))
  expect_false(gb$is_present(NA))
  expect_false(gb$is_present(""))
  expect_false(gb$is_present("0"))
  expect_true(gb$is_present("RD-26-1"))
})

test_that("CheckSarcRDnumber strips the _sarc suffix only when present", {
  expect_identical(suppressMessages(gb$CheckSarcRDnumber("RD-26-1")), "RD-26-1")
  expect_identical(suppressMessages(gb$CheckSarcRDnumber("RD-26-1_sarc")), "RD-26-1")
})

test_that("sanitize_string normalizes a filename", {
  skip_if_not(exists("sanitize_string", envir = gb, inherits = FALSE),
              "convert_html_pdf.R not sourced (pagedown absent)")
  expect_identical(gb$sanitize_string("My Report (v2).html"), "My-Report-v2")
  expect_identical(gb$sanitize_string("RD-26-1", ".pdf"), "RD-26-1.pdf")
})
