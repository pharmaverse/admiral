library(tibble)
library(lubridate)

test_that("derive_vars_duration Test 1: Duration and unit variable are added", {
  input <- tribble(
    ~USUBJID, ~BRTHDT, ~RANDDT,
    "P01", ymd("1984-09-06"), ymd("2020-02-24"),
    "P02", ymd("1985-01-01"), NA,
    "P03", NA, ymd("2021-03-10"),
    "P04", NA, NA
  )
  expected_output <- mutate(input,
    AGE = c(35, NA, NA, NA),
    AGEU = c("YEARS", NA_character_, NA_character_, NA_character_)
  )
  actual_output <- derive_vars_duration(
    input,
    new_var = AGE,
    start_date = BRTHDT,
    end_date = RANDDT,
    new_var_unit = AGEU,
    out_unit = "years",
    trunc_out = TRUE
  )

  expect_dfs_equal(actual_output, expected_output, keys = "USUBJID")
})

test_that("derive_vars_duration Test 2: Duration and unit variable are added", {
  input <- tribble(
    ~USUBJID, ~ASTDT, ~AENDT,
    "P01", ymd("2021-03-05"), ymd("2021-03-02"),
    "P02", ymd("2019-09-18"), ymd("2019-09-18"),
    "P03", ymd("1985-01-01"), NA,
    "P04", NA, NA
  )
  expected_output <- mutate(
    input,
    ADURN = c(-3, 1, NA, NA),
    ADURU = c("DAYS", "DAYS", NA_character_, NA_character_)
  )
  actual_output <- derive_vars_duration(
    input,
    new_var = ADURN,
    new_var_unit = ADURU,
    start_date = ASTDT,
    end_date = AENDT,
    out_unit = "days"
  )

  expect_dfs_equal(actual_output, expected_output, keys = "USUBJID")
})

test_that("derive_vars_duration Test 3: Duration and unit variable are added", {
  input <- tribble(
    ~USUBJID, ~ADTM, ~TRTSDTM,
    "P01", ymd_hms("2019-08-09T04:30:56"), ymd_hms("2019-08-09T05:00:00"),
    "P02", ymd_hms("2019-11-11T10:30:00"), ymd_hms("2019-11-11T11:30:00"),
    "P03", ymd("2019-11-11"), ymd_hms("2019-11-11T04:00:00"),
    "P04", NA, ymd_hms("2019-11-11T12:34:56"),
  )
  expected_output <- mutate(
    input,
    ADURN = c(30, 60, 240, NA),
    ADURU = c("MINUTES", "MINUTES", "MINUTES", NA_character_)
  )
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

  expect_dfs_equal(actual_output, expected_output, keys = "USUBJID")
})
