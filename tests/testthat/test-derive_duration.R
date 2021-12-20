test_that("duration and unit variable are added", {
  input <- tibble::tribble(
    ~BRTHDT, ~RANDDT,
    ymd("1999-09-09"), ymd("2020-02-20")
  )
  expected_output <- mutate(input, AGE = 20, AGEU = "YEARS")
  actual_output <- derive_vars_duration(
    input,
    new_var = AGE,
    start_date = BRTHDT,
    end_date = RANDDT,
    new_var_unit = AGEU,
    out_unit = "years",
    trunc_out = TRUE
  )

  expect_equal(actual_output, expected_output)
})
