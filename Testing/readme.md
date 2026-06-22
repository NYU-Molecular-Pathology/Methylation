# Methylation Pipeline — Testing

`testthat` unit/integration tests for the R workflow. They source the **local**
pipeline scripts in `../R` and exercise the same functions the production
pipeline runs, so you can confirm that a refactor (or a change merged to `main`)
still behaves the same.

The suite runs **offline**: no `/Volumes/` mounts, no REDCap network calls, no
patient data. idat-based tests use the public **Illumina EPIC v2 demo dataset**.
The classifier package `mnp.v12epicv2` (which processes both EPIC v1 and v2 idats)
is assumed to be installed locally — it is not public and is never downloaded by
these tests.

---

## 1. Prerequisites

Always required:

- R (the version the pipeline targets, ≥ 4.3.2)
- `testthat` (>= 3), `withr`, `openxlsx2`
- Pipeline runtime deps already used by the scripts: `crayon`, `fs`, `stringr`,
  `stringi`, `sjmisc`, `readxl`

Required only for the idat / RGset / rendering tests (these tests **skip**
automatically if missing):

- `minfi` and `IlluminaHumanMethylationEPICv2manifest` (installed with `mnp.v12epicv2`)
- `mnp.v12epicv2` (local, non-public) — rendering tests only
- `rmarkdown` — rendering tests only

Optional (related tests skip if absent): `dplyr`, `mockery`, `pagedown`.

Install the test-only helpers if needed:

```r
install.packages(c("testthat", "withr", "openxlsx2"))
```

---

## 2. Test data (EPIC v2 demo)

idat tests look for the Illumina demo idats under a `Test_Data` directory with
the standard layout:

```
Test_Data/.../<barcode>/<barcode>_RxxCxx_Grn.idat
Test_Data/.../<barcode>/<barcode>_RxxCxx_Red.idat
```

Resolution order:

1. `TEST_DATA_DIR` environment variable, if set.
2. Default: `Test_Data` next to the repo (i.e. `<repo>/../Test_Data`).

`make test-data` (runs `prepare_test_data.R`):

1. Finds the demo idats — and if absent, downloads + unzips `DemoDataEPIC_v2.zip`
   from Illumina into `Test_Data/`.
2. Generates a demo worksheet `Test_Data/<RUN_ID>.xlsm` (RUN_ID `26-MGDM-Test`)
   built from the **real** sentrix ids of two BeadChips: `206891110001` and
   `206891110002` (the demo ships four arrays; these two give eight samples).

By default tests use those two BeadChips; single-sample tests (`rgset`) use just
the first sample to stay fast. If the data cannot be made available, the
idat tests skip rather than fail.

Point the tests at an existing copy with your **real** path (the `you` below is a
placeholder — substitute your username, or just `unset TEST_DATA_DIR` to use the
default next to the repo):

```bash
export TEST_DATA_DIR="$HOME/Documents/Development/Test_Data"
```

> Troubleshooting `make test-data` → "Could not obtain EPIC v2 demo data": this
> almost always means `TEST_DATA_DIR` points at a non-existent folder (e.g. the
> literal `/Users/you/...` from a copied example). The error prints the exact
> path it checked — fix the export or `unset TEST_DATA_DIR`.

---

## 3. Running the tests

### With `make` (recommended — from the `Testing/` directory)

```bash
cd Testing
make setup        # install/verify test-only packages (testthat, withr, openxlsx2)
make test-data    # locate or download the EPIC v2 demo idats into Test_Data/
make test-all     # run the whole suite
```

Run a single group with `make test-<type>`:

```bash
make test-string-utils
make test-samplesheet-parsing
make test-rgset
make test-render          # opt-in; see section 4
```

| Target | What it does |
|---|---|
| `make help` | List all targets (default goal) |
| `make setup` | Install/verify `testthat (>=3)`, `withr`, `openxlsx2` (pipeline pkgs are handled by `LoadInstallPackages.R`) |
| `make test-data` | Ensure the EPIC v2 demo idats exist (download + unzip if missing) **and** build a demo `<RUN_ID>.xlsm` in `Test_Data/` |
| `make test-all` (or `make test`) | Run the entire suite; non-zero exit on any failure |
| `make test-<type>` | Run one group: `string-utils`, `arg-handling`, `samplesheet-parsing`, `samplesheet-validation`, `priority-ordering`, `qc-metrics`, `idat-resolve`, `rgset`, `pure-helpers`, `redcap-mocked`, `api`, `demo-pipeline`, `render` |
| `make clean` | Remove transient artifacts (HTML/PDF/CSV/`.Rout`/`.RData`); keeps `Test_Data` |
| `make clean-data` | Remove the downloaded demo data (`Test_Data`) |

Environment variables are forwarded to the tests, e.g.:

```bash
TEST_DATA_DIR=/path/to/Test_Data make test-rgset
RUN_RENDER_TESTS=1 REPORT_RMD=/path/My_new_report.Rmd make test-render
API_TOKEN=<32-char REDCap token> make test-api
```

### From a terminal without make

```bash
Rscript Testing/run_tests.R     # from the repo root
```

Exit status is `0` when everything passes (skips are OK) and `1` on any
failure/error — convenient for CI or a pre-merge check.

### From RStudio / an R console

```r
setwd("/path/to/Methylation")   # repo root
source("Testing/run_tests.R")
# or, for the live testthat UI:
testthat::test_dir("Testing/testthat")
# a single group:
testthat::test_dir("Testing/testthat", filter = "samplesheet-parsing")
```

---

## 4. Rendering / knitting tests (opt-in)

The report and QC render tests are heavy and disabled by default. They knit a
real `.Rmd` template and need `mnp.v12epicv2`, `rmarkdown`, and the demo data.

Enable them and choose which templates to render:

```bash
export RUN_RENDER_TESTS=1
# Test a locally modified report template (defaults to the pipeline's hard-coded
# /Volumes/.../EPIC_V2_report_2.Rmd):
export REPORT_RMD="/Users/you/Documents/Development/My_new_test_report.Rmd"
# Test a locally modified QC template (defaults to ~/Methyl_QC.Rmd):
export QC_RMD="/Users/you/Documents/Development/My_new_QC_test.Rmd"

Rscript Testing/run_tests.R
```

Notes:

- With `RUN_RENDER_TESTS` unset, these two tests skip; the rest of the suite
  still runs.
- The defaults are exactly the paths hard-coded in `pipelineHelper.R`
  (`gb$reportMd`, `gb$QC_file`). A lightweight test (always run) asserts those
  defaults so a change to the hard-coded paths is flagged.
- The QC template is discovered by the pipeline via a `*QC.Rmd` pattern, so a
  custom QC template filename should contain `QC.Rmd` (e.g. `My_new_QC.Rmd`).
- The suite uses a placeholder REDCap token; if your template performs live
  REDCap calls during knitting, set `gb$ApiToken` accordingly before running.

---

## 5. Interpreting results

`testthat` reports three outcomes:

- **PASS** — the function behaved as expected.
- **FAIL / ERROR** — a behavior changed (or a dependency is genuinely broken).
- **SKIP** — a prerequisite was absent (e.g. `mnp.v12epicv2`, the demo data, or
  `RUN_RENDER_TESTS` not set). Skips are expected and not failures.

A clean baseline looks like `0 failed, 0 errors` with some skips depending on
what is installed.

---

## 6. Using the suite when refactoring

1. On the current/unchanged code, run the suite and confirm a green baseline.
2. Make the refactor (e.g. flatten `get_idats`, simplify inputs, `T`→`TRUE`).
3. Run the suite again. Any new FAIL points to a behavior change.
4. If the change was intentional, update the affected expectation deliberately
   (and note why). If not, fix the code.

The tests assert **outcomes** (files produced, returned values, samplesheet
contents), not internal call structure, so they remain valid as functions are
reorganized.

### QC pass/fail criteria

`test-qc-metrics.R` asserts the `CheckSampleQCmetrics` criteria: it counts how
many of the five metrics are within their cutoffs (SI<=9.0, BP<=11.0, BS<=10.0,
HC<=12.75, NC[Pvalue]<=0.95) and sets `qc_passed = "no"` when **>= 2** metrics
are in range, `"yes"` when **< 2** are in range. The fixture rows exercise the
boundary (5, 2, 1, and 0 metrics in range).

### REDCap API tests

`test-api.R` reads the token from the `API_TOKEN` environment variable so it is
never committed. If `API_TOKEN` is not exported, the API tests **skip** with a
message (a plain `make test-all` stays green). To run them:

```bash
export API_TOKEN="<your 32-char REDCap token>"
make test-api
```

---

## 7. Coverage map

| Test file | Area | Notable deps |
|---|---|---|
| `test-string-utils.R` | `grabYear`, `get_next_rd`, `get_newest_rd`, `is_present`, `CheckSarcRDnumber`, `sanitize_string` | none (pagedown for sanitize_string) |
| `test-arg-handling.R` | `CheckInputArg`, `CheckBaseFolderInput` | none |
| `test-samplesheet-parsing.R` | `getTotalSamples`, `readSheetWrite` (xlsm → samplesheet.csv) | `openxlsx2`, `readxl` |
| `test-samplesheet-validation.R` | `checkHeaders`, `checkSampleSheet`, `NameControl` | none |
| `test-priority-ordering.R` | `reOrderRun`, `GetPriorityCases`, `getRunList` | none |
| `test-qc-metrics.R` | `CheckSampleQCmetrics`, `RenameFailed` | none |
| `test-idat-resolve.R` | `getAllFiles`, `copyBaseIdats`, `check_idat_sizes`, `CheckIdatsReal` | EPIC v2 demo, `fs`, `cli` |
| `test-rgset.R` | `getRGset` (EPIC v2) | `minfi`, EPIC v2 manifest, demo data |
| `test-api.R` | REDCap token present + connection | `API_TOKEN` env var; `redcapAPI` |
| `test-demo-pipeline.R` | real demo worksheet → samplesheet → idat resolution (2 BeadChips) | `openxlsx2`, `readxl`, demo data |
| `test-pure-helpers.R` | `check_validation`, `getRunData`, `rename_reorder_output` | `dplyr` (last one) |
| `test-redcap-mocked.R` | `checkRedcapRecord` (stubbed REDCap) | `mockery`, `redcapAPI` |
| `test-render.R` | RMD config + `make_knit_report` / `generateQCreport` | opt-in; `mnp.v12epicv2`, `rmarkdown`, demo data |

---

## 8. Files

```
Testing/
├── Makefile                        # make setup / test-data / test-all / test-<type> / clean
├── install_test_deps.R             # used by `make setup`
├── prepare_test_data.R             # used by `make test-data` (locate idats + build demo .xlsm)
├── run_tests.R                     # runner (Rscript or source())
├── README.md                       # this file
├── TESTING_PLAN.md                 # design notes / rationale
└── testthat/
    ├── setup.R                     # sources local R/ scripts into gb (no side effects)
    ├── helper-fixtures.R           # RUN_ID, temp run dirs, samplesheet/QC builders, RMD config
    ├── helper-make_dummy_xlsm.R    # builds a fake <RUN_ID>.xlsm matching the lab template
    ├── helper-idats.R              # locate/download EPIC v2 demo idats; pick one sample
    └── test-*.R                    # the tests above
```
