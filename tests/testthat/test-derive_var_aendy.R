context("test-derive_var_aendy")


test_that("AENDY is added", {
  input <- tibble::tribble(
    ~TRTSDT, ~AENDT,
    ymd("2020-01-01"), ymd("2020-02-24"),
    ymd("2020-01-01"), ymd("2020-01-01"),
    ymd("2020-02-24"), ymd("2020-01-01")
  )
  expected_output <- mutate(input, AENDY = c(55, 1, -54))

  expect_equal(derive_var_aendy(input), expected_output)
})
