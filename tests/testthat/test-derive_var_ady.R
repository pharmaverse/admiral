test_that("ADY is added", {
  input <- tibble::tribble(
    ~TRTSDT, ~ADT,
    ymd("2020-01-01"), ymd("2020-02-24"),
    ymd("2020-01-01"), ymd("2020-01-01"),
    ymd("2020-02-24"), ymd("2020-01-01")
  )

  expected_output <- mutate(input, ADY = c(55, 1, -54))

  expect_equal(derive_var_ady(input), expected_output)
})
