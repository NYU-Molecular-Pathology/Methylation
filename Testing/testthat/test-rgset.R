## getRGset on the EPIC v2 demo idats (pipelineHelper.R).
## The sentrix is taken from the generated samplesheet.csv; idats are read in
## place from the demo directory. Requires minfi + the EPIC v2 manifest (installed
## alongside mnp.v12epicv2), openxlsx2/readxl, and the demo data; SKIPS otherwise.

test_that("getRGset reads EPIC v2 demo idats and assigns EPICv2 annotation", {
  skip_if_not_installed("minfi")
  skip_if_not_installed("IlluminaHumanMethylationEPICv2manifest")

  ss <- setup_demo_run()
  skip_if(is.null(ss), "EPIC v2 demo data / required packages unavailable")

  sentrix  <- ss$SentrixID_Pos[1]
  idat_dir <- demo_idat_dir_for(sentrix)
  skip_if(is.null(idat_dir), "Could not locate demo idats for samplesheet sentrix")

  rg <- suppressWarnings(suppressMessages(gb$getRGset(idat_dir, sentrix)))
  expect_s4_class(rg, "RGChannelSet")
  expect_match(minfi::annotation(rg)[["array"]], "EPICv2")
})
