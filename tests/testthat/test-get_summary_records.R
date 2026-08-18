# get_summary_records ----
## Test 1: deprecation error if function is called ----
test_that("get_summary_records Test 1: deprecation error if function is called", {
  expect_error(
    get_summary_records(
      tibble::tribble(~USUBJID, ~PARAM, ~AVAL),
      by_vars = exprs(USUBJID, PARAM),
      set_values_to = exprs(AVAL = mean(AVAL, na.rm = TRUE))
    ),
    class = "lifecycle_error_deprecated"
  )
})
