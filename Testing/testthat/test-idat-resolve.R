## idat path resolution, copying, and completeness checks (CopyInputs.R).

test_that("getAllFiles builds Grn/Red idat paths for each sample", {
  local_run_dir()
  mk_samplesheet("samplesheet.csv", c("RD-26-1", "RD-26-2"))
  allFi <- suppressMessages(gb$getAllFiles(idatDir = "/fake/idat", csvNam = "samplesheet.csv"))
  expect_equal(length(allFi), 4L)                       # 2 samples x (Grn + Red)
  expect_true(all(grepl("\\.idat$", allFi)))
  expect_equal(sum(grepl("_Grn\\.idat$", allFi)), 2L)
  expect_equal(sum(grepl("_Red\\.idat$", allFi)), 2L)
})

test_that("getAllFiles resolves real EPIC v2 demo idat paths that exist on disk", {
  data_dir <- ensure_epicv2_demo()
  skip_if(is.null(data_dir), "EPIC v2 demo data unavailable")
  s <- pick_one_epicv2(data_dir)
  skip_if(is.null(s), "No EPIC v2 demo sample available")

  local_run_dir()
  write_demo_samplesheet("samplesheet.csv", s$sentrix)
  allFi <- suppressMessages(gb$getAllFiles(idatDir = s$idat_root, csvNam = "samplesheet.csv"))
  expect_equal(length(allFi), 2L)
  expect_true(all(file.exists(allFi)))
})

test_that("copyBaseIdats copies real demo idats into the run directory", {
  skip_if_not_installed("fs")
  skip_if_not_installed("cli")
  data_dir <- ensure_epicv2_demo()
  skip_if(is.null(data_dir), "EPIC v2 demo data unavailable")
  s <- pick_one_epicv2(data_dir)
  skip_if(is.null(s), "No EPIC v2 demo sample available")

  rd <- local_run_dir()
  suppressWarnings(suppressMessages(gb$copyBaseIdats(c(s$grn, s$red))))
  expect_true(file.exists(file.path(rd, basename(s$grn))))
  expect_true(file.exists(file.path(rd, basename(s$red))))
})

test_that("check_idat_sizes warns when idat files differ in size", {
  rd <- local_run_dir()
  writeBin(raw(200000), "206891110001_R01C01_Grn.idat")  # ~0.2 MB
  writeBin(raw(500000), "206891110002_R01C01_Grn.idat")  # ~0.5 MB
  expect_message(
    suppressWarnings(gb$check_idat_sizes(rd)),
    "Not all idat files are the same file size"
  )
})

test_that("CheckIdatsReal logs missing idats to a CSV", {
  local_run_dir()
  ssheet <- data.frame(SentrixID_Pos = c("206891110001_R01C01", "206891110001_R02C01"),
                       stringsAsFactors = FALSE)
  allFi <- c("206891110001_R01C01_Grn.idat", "206891110001_R01C01_Red.idat")  # 1 of 2
  suppressWarnings(suppressMessages(gb$CheckIdatsReal(ssheet, allFi)))
  expect_true(file.exists("missing_idats_log.csv"))
})
