## Samplesheet validation helpers (MakeSampleSheet.R, pipelineHelper.R).

mk_val_df <- function(sample_names) {
  n <- length(sample_names)
  data.frame(
    Sentrix_ID   = sprintf("204%09d_R0%dC01", seq_len(n), ((seq_len(n) - 1) %% 8) + 1),
    c2 = "", c3 = "", c4 = "", c5 = "", c6 = "",
    MP_number    = paste0("MP-", seq_len(n)),
    Sample_Name  = sample_names,
    b_number     = paste0("B", seq_len(n)),
    Project      = TEST_RUN_ID, Batch = TEST_RUN_ID,
    Tech = "T1", Tech2 = "T2", Notes = "",
    check.names = FALSE, stringsAsFactors = FALSE
  )
}

valid_names <- test_sample_names(8L)  # control + RD-26-1-test .. RD-26-7-test

test_that("checkHeaders errors when required headers are missing", {
  bad <- data.frame(A = 1, B = 2, C = 3)
  expect_error(suppressMessages(gb$checkHeaders(bad)))
})

test_that("checkSampleSheet accepts a well-formed sheet", {
  df <- suppressWarnings(suppressMessages(gb$checkSampleSheet(mk_val_df(valid_names))))
  expect_equal(nrow(df), 8L)
  expect_identical(df$Sample_Name[1], "control")
})

test_that("checkSampleSheet rejects duplicated sample names", {
  dupe <- c("control", "RD-26-1", "RD-26-1", paste0("RD-26-", 3:7))
  expect_error(suppressWarnings(suppressMessages(gb$checkSampleSheet(mk_val_df(dupe)))))
})

test_that("checkSampleSheet rejects a capitalized 'Control' first row", {
  caps <- c("Control", paste0("RD-26-", 1:7))
  expect_error(suppressWarnings(suppressMessages(gb$checkSampleSheet(mk_val_df(caps)))))
})

test_that("NameControl prefixes unnamed control rows with the RunID", {
  d <- data.frame(V1 = c("control", "RD-26-1"), stringsAsFactors = FALSE)
  out <- suppressWarnings(suppressMessages(gb$NameControl(d, TEST_RUN_ID)))
  expect_identical(out$V1[1], paste0(TEST_RUN_ID, "_control"))
  expect_identical(out$V1[2], "RD-26-1")
})

test_that("NameControl leaves already-named controls unchanged", {
  d <- data.frame(V1 = c(paste0(TEST_RUN_ID, "_control"), "RD-26-1"), stringsAsFactors = FALSE)
  out <- suppressWarnings(suppressMessages(gb$NameControl(d, TEST_RUN_ID)))
  expect_identical(out$V1[1], paste0(TEST_RUN_ID, "_control"))
})
