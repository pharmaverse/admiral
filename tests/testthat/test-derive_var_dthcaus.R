# dthcaus_source ----
## Test 1: deprecation error if dthcaus_source is called ----
test_that("dthcaus_source Test 1: deprecation error if function is called", {
  expect_error(
    dthcaus_source(
      dataset_name = "ae",
      filter = AEOUT == "FATAL",
      date = AEDTHDT,
      mode = "first",
      dthcaus = AEDECOD
    ),
    class = "lifecycle_error_deprecated"
  )
})

# derive_var_dthcaus ----
## Test 2: deprecation error if derive_var_dthcaus is called ----
test_that("derive_var_dthcaus Test 2: deprecation error if function is called", {
  expect_error(
    derive_var_dthcaus(
      tibble::tribble(~STUDYID, ~USUBJID),
      source_datasets = list()
    ),
    class = "lifecycle_error_deprecated"
  )
})
