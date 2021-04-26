context("test-derive_var_trtdurd")


test_that("TRTDURD is added", {
  input <- tibble::tribble(
    ~TRTSDT, ~TRTEDT,
    ymd("2020-01-01"), ymd("2020-02-24")
  )
  expected_output <- mutate(input, TRTDURD = 55)

  expect_equal(derive_var_trtdurd(input), expected_output)
})
