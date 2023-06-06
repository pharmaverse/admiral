## Test 1: compute_age_years() works when `age_unit` is a string ----
test_that("compute_age_years Test 1: compute_age_years() works when `age_unit` is a string", {
  age_input <- c(240, 360, 480, NA)
  age_unit_input <- "MONTHS"

  expected_output <- c(20, 30, 40, NA)

  expect_equal(
    compute_age_years(
      age_input,
      age_unit_input
    ),
    expected_output
  )
})

## Test 2: compute_age_years() works when `age_unit` is a vector ----
test_that("compute_age_years Test 2: compute_age_years() works when `age_unit` is a vector", {
  age_input <- c(28, 1461, 10227)
  age_unit_input <- c("YEARS", "WEEKS", "DAYS")

  expected_output <- rep(28, 3)

  expect_equal(
    compute_age_years(
      age_input,
      age_unit_input
    ),
    expected_output
  )
})

## Test 3: Error is issued when `age_unit` has invalid length ----
test_that("compute_age_years Test 3: Error is issued when `age_unit` has invalid length", {
  age_input <- c(28, 1461, 10227)
  age_unit_input <- c("YEARS", "WEEKS")

  expected_output <- rep(28, 3)

  expect_error(
    compute_age_years(age_input, age_unit_input),
    paste0(
      "`age_unit` must be a single string or a vector of the same length as",
      "`age`, but there are 3 values in `age` and 2 values in `age_unit`."
    )
  )
})

## Test 4: `age_unit` processes values in a case insensitive manner ----
test_that("compute_age_years Test 4: `age_unit` processes values in a case insensitive manner", {
  age_input <- c(240, 360, 480)
  age_unit_input <- c("MONTHS", "Months", "months")

  expected_output <- c(20, 30, 40)

  expect_equal(
    compute_age_years(
      age_input,
      age_unit_input
    ),
    expected_output
  )
})
