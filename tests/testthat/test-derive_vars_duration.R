library(tibble)
library(lubridate)

test_that("duration and unit variable are added", {
  input <- tribble(
    ~BRTHDT, ~RANDDT,
    ymd("1984-09-06"), ymd("2020-02-24"),
    ymd("1985-01-01"), NA,
    NA, ymd("2021-03-10"),
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

test_that("duration and unit variable are added", {
  input <- tribble(
    ~ASTDT, ~AENDT,
    ymd("2021-03-05"), ymd("2021-03-02"),
    ymd("2019-09-18"), ymd("2019-09-18"),
    ymd("1985-01-01"), NA,
    NA, NA
  )
  expected_output <- mutate(input, ADURN = c(-3, 1, NA, NA), ADURU = c("DAYS", "DAYS", NA, NA))
  actual_output <- derive_vars_duration(
    input,
    new_var = ADURN,
    new_var_unit = ADURU,
    start_date = ASTDT,
    end_date = AENDT,
    out_unit = "days"
  )

  expect_equal(actual_output, expected_output)
})

test_that("duration and unit variable are added", {
  input <- tribble(
    ~ADTM, ~TRTSDTM,
    ymd_hms("2019-08-09T04:30:56"), ymd_hms("2019-08-09T05:00:00"),
    ymd_hms("2019-11-11T10:30:00"), ymd_hms("2019-11-11T11:30:00"),
    ymd("2019-11-11"), ymd_hms("2019-11-11T04:00:00"),
    NA, ymd_hms("2019-11-11T12:34:56"),
  )
  expected_output <- mutate(input, ADURN = c(30, 60, 240, NA), ADURU = c("MINUTES", "MINUTES", "MINUTES", NA))
  actual_output <- derive_vars_duration(
    input,
    new_var = ADURN,
    new_var_unit = ADURU,
    start_date = ADTM,
    end_date = TRTSDTM,
    in_unit = "minutes",
    out_unit = "minutes",
    add_one = FALSE
  )

  expect_equal(actual_output, expected_output)
})
