# date_source ----
## Test 1: deprecation error if date_source is called ----
test_that("date_source Test 1: deprecation error if function is called", {
  expect_error(
    date_source(
      dataset_name = "adsl",
      date = TRTEDT
    ),
    class = "lifecycle_error_deprecated"
  )
})

# derive_var_extreme_dt ----
## Test 2: deprecation error if derive_var_extreme_dt is called ----
test_that("derive_var_extreme_dt Test 2: deprecation error if function is called", {
  expect_error(
    derive_var_extreme_dt(
      tibble::tribble(~STUDYID, ~USUBJID),
      new_var = LSTALVDT,
      source_datasets = list(),
      mode = "last"
    ),
    class = "lifecycle_error_deprecated"
  )
})

# derive_var_extreme_dtm ----
## Test 3: deprecation error if derive_var_extreme_dtm is called ----
test_that("derive_var_extreme_dtm Test 3: deprecation error if function is called", {
  expect_error(
    derive_var_extreme_dtm(
      tibble::tribble(~STUDYID, ~USUBJID),
      new_var = LSTALVDTM,
      source_datasets = list(),
      mode = "last"
    ),
    class = "lifecycle_error_deprecated"
  )
})
