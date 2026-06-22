#!/usr/bin/env Rscript
## install_test_deps.R — ensure test-only packages are available.
## Pipeline runtime packages are installed/loaded by R/LoadInstallPackages.R when
## the pipeline runs; this only covers the harness packages.

required_packages <- c("testthat", "withr", "openxlsx2")

optional_packages <- c(
  "rmarkdown",
  "dplyr",
  "mockery",
  "fs",
  "cli",
  "crayon",
  "stringr",
  "stringi",
  "sjmisc",
  "readxl",
  "BiocManager",
  "devtools"
)

local_packages <- c("minfi", "IlluminaHumanMethylationEPICv2manifest", "mnp.v12epicv2")

print_package_status <- function(packages, installed_packages) {
  for (package in packages) {
    if (package %in% installed_packages) {
      cat(sprintf(
        "  [✓] %s (%s)\n", package, as.character(utils::packageVersion(package))
      ))
    } else {
      cat(sprintf("  [ ] %s\n", package))
    }
  }
}

installed_packages <- row.names(utils::installed.packages())

missing_required <- setdiff(required_packages, installed_packages)

if (length(missing_required) > 0) {
  message("Installing required packages: ", paste(missing_required, collapse = ", "))

  utils::install.packages(
    missing_required,
    repos = "https://cloud.r-project.org",
    ask = FALSE,
    dependencies = TRUE
  )
}

installed_packages <- row.names(utils::installed.packages())
missing_required <- setdiff(required_packages, installed_packages)

if (length(missing_required) > 0) {
  stop("Failed to install required packages: ", paste(missing_required, collapse = ", "))
}

if (utils::compareVersion(as.character(utils::packageVersion("testthat")), "3.0.0") < 0) {
  stop("testthat >= 3.0.0 is required")
}

missing_optional <- setdiff(optional_packages, installed_packages)

if (length(missing_optional) > 0) {
  message("Installing optional CRAN packages: ", paste(missing_optional, collapse = ", "))

  utils::install.packages(
    missing_optional,
    repos = "https://cloud.r-project.org",
    ask = FALSE,
    dependencies = TRUE
  )
}

installed_packages <- row.names(utils::installed.packages())

missing_optional <- setdiff(optional_packages, installed_packages)
missing_local <- setdiff(local_packages, installed_packages)

cat("\nRequired packages:\n")
print_package_status(required_packages, installed_packages)

cat("\nOptional packages:\n")
print_package_status(optional_packages, installed_packages)

cat("\nBioconductor/private packages:\n")
print_package_status(local_packages, installed_packages)

if (length(missing_optional) > 0) {
  warning("Optional packages missing: ", paste(missing_optional, collapse = ", "))
}

if (length(missing_local) > 0) {
  message("Bioconductor/private packages missing: ", paste(missing_local, collapse = ", "))
}

message("\nPackage check complete!\n")
