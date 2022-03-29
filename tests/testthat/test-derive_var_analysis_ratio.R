test_that("Test 1: All Ratio Variables are Created", {
  expected_data <- tibble::tribble(
    ~USUBJID, ~PARAMCD, ~SEQ, ~AVAL, ~BASE, ~ANRLO, ~ANRHI, ~R2BASE, ~R2ANRLO, ~R2ANRHI,
    "P01", "ALT", 1, 27, 27, 6, 34, 1, 4.5, 0.794,
    "P01", "ALT", 2, 41, 27, 6, 34, 1.52, 6.83, 1.21,
    "P01", "ALT", 3, 17, 27, 6, 34, 0.63, 2.83, 0.5,
    "P02", "ALB", 1, 38, 38, 33, 49, 1, 1.15, 0.776,
    "P02", "ALB", 2, 39, 38, 33, 49, 1.03, 1.18, 0.796,
    "P02", "ALB", 3, 37, 38, 33, 49, 0.974, 1.12, 0.755
  )

  data <- tibble::tribble(
    ~USUBJID, ~PARAMCD, ~SEQ, ~AVAL, ~BASE, ~ANRLO, ~ANRHI,
    "P01", "ALT", 1, 27, 27, 6, 34,
    "P01", "ALT", 2, 41, 27, 6, 34,
    "P01", "ALT", 3, 17, 27, 6, 34,
    "P02", "ALB", 1, 38, 38, 33, 49,
    "P02", "ALB", 2, 39, 38, 33, 49,
    "P02", "ALB", 3, 37, 38, 33, 49
  )

  actual_data <- data %>%
    derive_var_analysis_ratio(
      numer_var = AVAL,
      denom_var = BASE
    ) %>%
    derive_var_analysis_ratio(
      numer_var = AVAL,
      denom_var = ANRLO
    ) %>%
    derive_var_analysis_ratio(
      numer_var = AVAL,
      denom_var = ANRHI
    )

  expect_dfs_equal(expected_data,
    actual_data,
    keys = c("USUBJID", "PARAMCD", "SEQ"),
    tolerance = 0.2
  )
})

test_that("Test 2: All Ratio Variables are Created while NAs present", {
  expected_data <- tibble::tribble(
    ~USUBJID, ~PARAMCD, ~SEQ, ~AVAL, ~BASE, ~ANRLO, ~ANRHI, ~R2BASE, ~R2ANRLO, ~R2ANRHI,
    "P01", "ALT", 1, 27, 27, 6, 34, 1, 4.5, 0.794,
    "P01", "ALT", 2, 41, 27, 6, 34, 1.52, 6.83, 1.21,
    "P01", "ALT", 3, NA_real_, 27, 6, 34, NA_real_, NA_real_, NA_real_,
    "P02", "ALB", 1, 38, 38, 33, 49, 1, 1.15, 0.776,
    "P02", "ALB", 2, 39, 38, 33, 49, 1.03, 1.18, 0.796,
    "P02", "ALB", 3, 37, 38, 33, 49, 0.974, 1.12, 0.755
  )

  data <- tibble::tribble(
    ~USUBJID, ~PARAMCD, ~SEQ, ~AVAL, ~BASE, ~ANRLO, ~ANRHI,
    "P01", "ALT", 1, 27, 27, 6, 34,
    "P01", "ALT", 2, 41, 27, 6, 34,
    "P01", "ALT", 3, NA_real_, 27, 6, 34,
    "P02", "ALB", 1, 38, 38, 33, 49,
    "P02", "ALB", 2, 39, 38, 33, 49,
    "P02", "ALB", 3, 37, 38, 33, 49
  )

  actual_data <- data %>%
    derive_var_analysis_ratio(
      numer_var = AVAL,
      denom_var = BASE
    ) %>%
    derive_var_analysis_ratio(
      numer_var = AVAL,
      denom_var = ANRLO
    ) %>%
    derive_var_analysis_ratio(
      numer_var = AVAL,
      denom_var = ANRHI
    )

  expect_dfs_equal(expected_data,
    actual_data,
    keys = c("USUBJID", "PARAMCD", "SEQ"),
    tolerance = 0.2
  )
})

test_that("Test 3: User can supply custom variable by invoking override = TRUE
          and supplying a new_var name", {
  expected_data <- tibble::tribble(
    ~USUBJID, ~PARAMCD, ~SEQ, ~AVAL, ~BASE, ~ANRLO, ~ANRHI, ~custom1, ~custom2, ~custom3,
    "P01", "ALT", 1, 27, 27, 6, 34, 1, 4.5, 0.794,
    "P01", "ALT", 2, 41, 27, 6, 34, 1.52, 6.83, 1.21,
    "P01", "ALT", 3, NA_real_, 27, 6, 34, NA_real_, NA_real_, NA_real_,
    "P02", "ALB", 1, 38, 38, 33, 49, 1, 1.15, 0.776,
    "P02", "ALB", 2, 39, 38, 33, 49, 1.03, 1.18, 0.796,
    "P02", "ALB", 3, 37, 38, 33, 49, 0.974, 1.12, 0.755
  )

  data <- tibble::tribble(
    ~USUBJID, ~PARAMCD, ~SEQ, ~AVAL, ~BASE, ~ANRLO, ~ANRHI,
    "P01", "ALT", 1, 27, 27, 6, 34,
    "P01", "ALT", 2, 41, 27, 6, 34,
    "P01", "ALT", 3, NA_real_, 27, 6, 34,
    "P02", "ALB", 1, 38, 38, 33, 49,
    "P02", "ALB", 2, 39, 38, 33, 49,
    "P02", "ALB", 3, 37, 38, 33, 49
  )

  actual_data <- data %>%
    derive_var_analysis_ratio(
      numer_var = AVAL,
      denom_var = BASE,
      new_var = custom1
    ) %>%
    derive_var_analysis_ratio(
      numer_var = AVAL,
      denom_var = ANRLO,
      new_var = custom2
    ) %>%
    derive_var_analysis_ratio(
      numer_var = AVAL,
      denom_var = ANRHI,
      new_var = custom3
    )

  expect_dfs_equal(expected_data,
    actual_data,
    keys = c("USUBJID", "PARAMCD", "SEQ"),
    tolerance = 0.2
  )
})
