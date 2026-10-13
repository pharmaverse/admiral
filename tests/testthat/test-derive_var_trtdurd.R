test_that("TRTDURD is added", {
  input <- tibble::tribble(
    ~TRTSDT, ~TRTEDT,
    ymd("2020-01-01"), ymd("2020-02-24")
  )
  expected_output <- mutate(input, TRTDURD = 55)
  actual_output <- derive_var_trtdurd(input)

  expect_equal(actual_output, as_admiral_df(expected_output))
  expect_s3_class(actual_output, "admiral_df")
})
