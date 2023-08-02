test_that("roxygen2 documentation helper functions work", {
  expect_true(inherits(roxygen_param_by_var(), "character"))
})
