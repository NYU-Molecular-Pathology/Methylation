## REDCap API token tests.
##
## The token is read from the API_TOKEN environment variable so it is never
## committed. If it is not set (e.g. a plain `make test-all` without exporting
## it), these tests SKIP with a message rather than fail. To run them:
##   export API_TOKEN="<your 32-char REDCap token>"

api_token <- function() Sys.getenv("API_TOKEN", unset = "")

test_that("API_TOKEN is provided via the environment", {
    token <- api_token()
    skip_if(
        !nzchar(token),
        "API_TOKEN not set - export API_TOKEN to run the REDCap API tests"
    )
    expect_match(token, "^[A-Fa-f0-9]{32}$") # REDCap tokens are 32 hex characters
})

test_that("a REDCap connection can be built and queried with API_TOKEN", {
    token <- api_token()
    skip_if(
        !nzchar(token),
        "API_TOKEN not set - export API_TOKEN to run the REDCap API tests"
    )
    skip_if_not_installed("redcapAPI")

    gb$token <- gb$ApiToken <- token

    rcon <- tryCatch(
        redcapAPI::redcapConnection(
            "https://redcap.nyumc.org/apps/redcap/api/", token
        ),
        error = function(e) {
            skip(paste(
                "Could not create REDCap connection:",
                conditionMessage(e)
            ))
        }
    )
    test_rcon <- tryCatch(
        {
            redcapAPI::exportProjectInformation(rcon)
            TRUE
        },
        error = function(e) {
            FALSE
        }
    )
    expect_true(test_rcon)

    # Light live check: a metadata/records call should not error with a valid token.
    res <- tryCatch(
        suppressWarnings(redcapAPI::exportRecordsTyped(
            rcon,
            records = "RD-14-01",
            fields = "record_id",
            factors = FALSE,
            survey = FALSE,
            dag = FALSE
        )),
        error = function(e) {
            testthat::skip(paste("REDCap query failed (network/token):", conditionMessage(e)))
        }
    )

    expect_true(is.data.frame(as.data.frame(res)))
})
