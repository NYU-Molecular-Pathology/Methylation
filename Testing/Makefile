# Makefile — Methylation pipeline test suite
# Run from the Testing/ directory:  make <target>
#
# Common usage:
#   make setup        # install/verify test-only packages (testthat, withr, openxlsx)
#   make test-data    # download/locate the EPIC v2 demo idats (Test_Data/)
#   make test-all     # run the full suite
#   make test-rgset   # run a single test group (see TEST_TYPES below)
#   make clean        # remove transient test artifacts
#
# Env vars honored (exported to R): TEST_DATA_DIR, RUN_RENDER_TESTS,
#   REPORT_RMD, QC_RMD  (e.g. `RUN_RENDER_TESTS=1 make test-render`)

RSCRIPT  ?= Rscript
TESTTHAT := testthat

# Demo data lives outside the repo and is git-ignored. Override with TEST_DATA_DIR.
TEST_DATA_DIR ?= $(abspath $(CURDIR)/../../Test_Data)
export TEST_DATA_DIR

# Test groups that `make test-<type>` can run individually.
TEST_TYPES := string-utils arg-handling samplesheet-parsing samplesheet-validation \
              priority-ordering qc-metrics idat-resolve rgset pure-helpers \
              redcap-mocked api demo-pipeline render

.DEFAULT_GOAL := help
# Note: the per-group `test-<type>` targets are intentionally NOT in .PHONY —
# GNU make skips pattern-rule (test-%) matching for phony targets.
.PHONY: help setup test-data test test-all clean clean-data

help:
	@echo "Methylation test suite — available targets:"
	@echo ""
	@echo "  make setup         Install/verify test packages (testthat, withr, openxlsx)"
	@echo "  make test-data     Ensure EPIC v2 demo idats exist + build demo .xlsm in Test_Data/"
	@echo "  make test-all      Run the entire suite"
	@echo "  make test-<type>   Run one group; <type> is one of:"
	@echo "                     $(TEST_TYPES)"
	@echo "  make clean         Remove transient test artifacts (keeps Test_Data)"
	@echo "  make clean-data    Remove the downloaded demo data (Test_Data)"
	@echo ""
	@echo "  Env: TEST_DATA_DIR=$(TEST_DATA_DIR)"
	@echo "       RUN_RENDER_TESTS=1 enables the opt-in render tests"
	@echo "       REPORT_RMD=/path.Rmd  QC_RMD=/path.Rmd  override default templates"
	@echo "       API_TOKEN=<token> enables the REDCap API tests (else they skip)"

setup:
	$(RSCRIPT) install_test_deps.R

test-data:
	$(RSCRIPT) prepare_test_data.R

# Full suite (exit status non-zero on any failure).
test test-all:
	$(RSCRIPT) run_tests.R

# Single group: `make test-rgset` -> filter "rgset" against testthat/test-rgset.R
test-%:
	$(RSCRIPT) -e 'testthat::test_dir("$(TESTTHAT)", filter="$*", stop_on_failure=TRUE)'

clean:
	@echo "Removing transient test artifacts (Test_Data preserved)..."
	@rm -f $(TESTTHAT)/*.Rout *.Rout .RData .Rhistory
	@rm -f *.html *.pdf *_QC.html samplesheet.csv missing_idats_log.csv
	@rm -rf test_outputs
	@find . -name '.DS_Store' -delete 2>/dev/null || true
	@echo "Done. Use 'make clean-data' to also remove $(TEST_DATA_DIR)."

clean-data:
	@echo "Removing demo data at: $(TEST_DATA_DIR)"
	@rm -rf "$(TEST_DATA_DIR)"
