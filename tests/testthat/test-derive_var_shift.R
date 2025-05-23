# derive_var_shift ----

## Test 1: Shift based on character variables ----
test_that("derive_var_shift Test 1: Shift based on character variables", {
  input <- tibble::tribble(
    ~USUBJID, ~PARAMCD, ~AVAL, ~ABLFL, ~BNRIND, ~ANRIND,
    "P01", "ALB", 33, "Y", "LOW", "LOW",
    "P01", "ALB", 38, NA_character_, "LOW", "NORMAL",
    "P02", "ALB", 37, "Y", "NORMAL", "NORMAL",
    "P02", "ALB", 49, NA_character_, "NORMAL", "HIGH",
    "P02", "SODIUM", 147, "Y", "HIGH", "HIGH"
  )
  expected_output <- tibble::tribble(
    ~USUBJID, ~PARAMCD, ~AVAL, ~ABLFL, ~BNRIND, ~ANRIND, ~SHIFT1,
    "P01", "ALB", 33, "Y", "LOW", "LOW", "LOW to LOW",
    "P01", "ALB", 38, NA_character_, "LOW", "NORMAL", "LOW to NORMAL",
    "P02", "ALB", 37, "Y", "NORMAL", "NORMAL", "NORMAL to NORMAL",
    "P02", "ALB", 49, NA_character_, "NORMAL", "HIGH", "NORMAL to HIGH",
    "P02", "SODIUM", 147, "Y", "HIGH", "HIGH", "HIGH to HIGH"
  )

  expect_equal(
    derive_var_shift(
      input,
      new_var = SHIFT1,
      from_var = BNRIND,
      to_var = ANRIND
    ),
    expected_output
  )
})


## Test 2: Shift based on character variables with missing values ----
test_that("derive_var_shift Test 2: Shift based on character variables with missing values", {
  input <- tibble::tribble(
    ~USUBJID, ~PARAMCD, ~AVAL, ~ABLFL, ~BNRIND, ~ANRIND,
    "P01", "ALB", 33, "Y", "LOW", "LOW",
    "P01", "ALB", 38, NA_character_, "LOW", "NORMAL",
    "P01", "ALB", NA, NA_character_, "LOW", NA,
    "P02", "ALB", NA, "Y", NA, NA,
    "P02", "ALB", 49, NA_character_, NA, "HIGH",
    "P02", "SODIUM", 147, "Y", "HIGH", "HIGH"
  )
  expected_output <- tibble::tribble(
    ~USUBJID, ~PARAMCD, ~AVAL, ~ABLFL, ~BNRIND, ~ANRIND, ~SHIFT1,
    "P01", "ALB", 33, "Y", "LOW", "LOW", "LOW to LOW",
    "P01", "ALB", 38, NA_character_, "LOW", "NORMAL", "LOW to NORMAL",
    "P01", "ALB", NA, NA_character_, "LOW", NA, "LOW to NULL",
    "P02", "ALB", NA, "Y", NA, NA, "NULL to NULL",
    "P02", "ALB", 49, NA_character_, NA, "HIGH", "NULL to HIGH",
    "P02", "SODIUM", 147, "Y", "HIGH", "HIGH", "HIGH to HIGH"
  )

  expect_equal(
    derive_var_shift(
      input,
      new_var = SHIFT1,
      from_var = BNRIND,
      to_var = ANRIND
    ),
    expected_output
  )
})


## Test 3: Shift based on numeric variables with missing values ----
test_that("derive_var_shift Test 3: Shift based on numeric variables with missing values", {
  input <- tibble::tribble(
    ~USUBJID, ~PARAMCD, ~AVAL, ~ABLFL, ~BASE,
    "P01", "ALB", 33.1, "Y", 33.1,
    "P01", "ALB", 38.5, NA_character_, 33.1,
    "P01", "ALB", NA, NA_character_, 33.1,
    "P02", "ALB", NA, "Y", NA,
    "P02", "ALB", 49.0, NA_character_, NA,
    "P02", "SODIUM", 147.5, "Y", 147.5
  ) %>% convert_blanks_to_na()
  expected_output <- tibble::tribble(
    ~USUBJID, ~PARAMCD, ~AVAL, ~ABLFL, ~BASE, ~SHIFT2,
    "P01", "ALB", 33.1, "Y", 33.1, "33.1 to 33.1",
    "P01", "ALB", 38.5, NA_character_, 33.1, "33.1 to 38.5",
    "P01", "ALB", NA, NA_character_, 33.1, "33.1 to NULL",
    "P02", "ALB", NA, "Y", NA, "NULL to NULL",
    "P02", "ALB", 49, NA_character_, NA, "NULL to 49",
    "P02", "SODIUM", 147.5, "Y", 147.5, "147.5 to 147.5"
  )

  expect_equal(
    derive_var_shift(
      input,
      new_var = SHIFT2,
      from_var = BASE,
      to_var = AVAL
    ),
    expected_output
  )
})

## Test 4: Shift with user-specified missing_value and sep_val ----
test_that("derive_var_shift Test 4: Shift with user-specified missing_value and sep_val", {
  input <- tibble::tribble(
    ~USUBJID, ~PARAMCD, ~AVAL, ~ABLFL, ~BNRIND, ~ANRIND,
    "P01", "ALB", 33, "Y", "LOW", "LOW",
    "P01", "ALB", 38, NA_character_, "LOW", "NORMAL",
    "P01", "ALB", NA, NA_character_, "LOW", NA,
    "P02", "ALB", NA, "Y", NA, NA,
    "P02", "ALB", 49, NA_character_, NA, "HIGH",
    "P02", "SODIUM", 147, "Y", "HIGH", "HIGH"
  )
  expected_output <- tibble::tribble(
    ~USUBJID, ~PARAMCD, ~AVAL, ~ABLFL, ~BNRIND, ~ANRIND, ~SHIFT1,
    "P01", "ALB", 33, "Y", "LOW", "LOW", "LOW - LOW",
    "P01", "ALB", 38, NA_character_, "LOW", "NORMAL", "LOW - NORMAL",
    "P01", "ALB", NA, NA_character_, "LOW", NA, "LOW - MISSING",
    "P02", "ALB", NA, "Y", NA, NA, "MISSING - MISSING",
    "P02", "ALB", 49, NA_character_, NA, "HIGH", "MISSING - HIGH",
    "P02", "SODIUM", 147, "Y", "HIGH", "HIGH", "HIGH - HIGH"
  )

  expect_equal(
    derive_var_shift(
      input,
      new_var = SHIFT1,
      from_var = BNRIND,
      to_var = ANRIND,
      missing_value = "MISSING",
      sep_val = " - "
    ),
    expected_output
  )
})
