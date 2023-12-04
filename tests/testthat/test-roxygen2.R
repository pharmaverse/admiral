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
      "The variables specified by the `by_vars` argument are expected to be in the dataset."
    ),
    roxygen_param_dataset(expected_vars = c("by_vars"))
  )
})


test_that("roxygen_param_by_vars Test 1: Text variations", {
  expect_snapshot(
    roxygen_param_by_vars()
  )
  expect_snapshot(
    roxygen_param_by_vars(rename = TRUE)
  )
})

# roxygen_order_na_handling ----
test_that("Standardized text works", {
  expect_equal(
    roxygen_order_na_handling(),
    paste(
      "For handling of `NA`s in sorting variables see",
      "[Sort Order](../articles/generic.html#sort_order)."
    )
  )
})
