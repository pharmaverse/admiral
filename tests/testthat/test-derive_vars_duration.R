test_that("duration and unit variable are added", {
  input <- tibble::tribble(
    ~BRTHDT, ~RANDDT,
    lubridate::ymd("1984-09-06"), lubridate::ymd("2020-02-24"),
    lubridate::ymd("1985-01-01"), NA,
    NA, lubridate::ymd("2021-03-10"),
    NA, NA
  )
  expected_output <- mutate(input, AGE = c(35, NA, NA, NA), AGEU = c("YEARS", NA, NA, NA))
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
