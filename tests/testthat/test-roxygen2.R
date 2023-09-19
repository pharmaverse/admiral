# roxygen_param_dataset ----
## Test 1: Input dataset standardized text works ----
test_that("roxygen_param_dataset Test 1: Input dataset standardized text works", {
  expect_equal("Input dataset", roxygen_param_dataset(expected_vars = NULL))
})

## Test 2: Input dataset standardized text works ----
test_that("roxygen_param_dataset Test 2: Input dataset standardized text works", {
  expect_equal(
    paste0(
      "Input dataset \n \n",
      "The variables specified by the `by_vars` argument(s) to be expected."
    ),
    roxygen_param_dataset(expected_vars = c("by_vars")))
})
