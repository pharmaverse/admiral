# derive_param_extreme_record ----
## Test 1: deprecation error if function is called ----
test_that("derive_param_extreme_record Test 1: deprecation error if function is called", {
  expect_error(
    derive_param_extreme_record(
      dataset = tibble::tribble(
        ~STUDYID, ~USUBJID,
        "1001", "1"
      ),
      sources = list(),
      source_datasets = list(),
      order = exprs(ADT),
      mode = "first",
      set_values_to = exprs(PARAMCD = "TEST")
    ),
    class = "lifecycle_error_deprecated"
  )
})
