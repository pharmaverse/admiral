test_that("duration and unit variable are added", {
  input <- tibble::tribble(
    ~BRTHDT, ~RANDDT,
    ymd("1999-09-09"), ymd("2020-02-20")
  )
  expected_output <- mutate(input, AAGE = 20, AAGEU = "YEARS")

  expect_dfs_equal(derive_vars_aage(input), expected_output, keys = c("BRTHDT", "RANDDT"))
})


test_that("derive_var_age_years works as expected", {

  input <- tibble::tibble(AGE = c(12, 24, 36, 48, 60),
                          AGEU = c("months", "months", "months", "months", "months"))

  expected_output <- mutate(
    input,
    AAGE = c(1, 2, 3, 4, 5)
    )

  expect_dfs_equal(derive_var_age_years(input, AGE, new_var = AAGE), expected_output, keys = "AGE")

})

test_that("derive_var_age_years works as expected", {

  input <- tibble::tibble(AGE = c(12, 24, 36, 48, 60))

  expected_output <- mutate(
    input,
    AAGE = c(1, 2, 3, 4, 5)
  )

  expect_dfs_equal(derive_var_age_years(input, AGE, new_var = AAGE, age_unit = "months"),
                   expected_output, keys = "AGE")

})

test_that("derive_var_agegr_fda works as expected", {

  input <- tibble::tibble(AGE = c(10, 17, 18, 50, 64, 65, 80))

  expected_output <- mutate(
    input,
    AGEGR_EXP = factor(
      c("<18", "<18", "18-64", "18-64", "18-64", ">=65", ">=65"),
      levels = c("<18", "18-64", ">=65"),
      exclude = NULL
    )
  )

  expect_dfs_equal(derive_var_agegr_fda(input, AGE, age_unit = "years", AGEGR_EXP), expected_output,
                   keys = "AGE")

})

test_that("derive_var_agegr_fda works with age_unit missing and multiple units in AGEU", {

  input <- tibble::tibble(AGE = c(10, 17, 18, 50, 64, 65, 80, 85),
                          AGEU = c("years", "years", "years", "years", "years", "years", "months",
                                   "months"))

  expected_output <- mutate(
    input,
    AGEGR_EXP = factor(
      c("<18", "<18", "18-64", "18-64", "18-64", ">=65", "<18", "<18"),
      levels = c("<18", "18-64", ">=65"),
      exclude = NULL
    )
  )

  expect_dfs_equal(derive_var_agegr_fda(input, AGE, age_unit = NULL, AGEGR_EXP), expected_output,
                   keys = "AGE")

})

test_that("derive_var_agegr_ema works as expected", {

  input <- tibble::tibble(AGE = c(10, 18, 19, 50, 64, 65, 80, 85))

  expected_output <- mutate(
    input,
    AGEGR_EXP = factor(
      c("2-11 (Children)", "18-64", "18-64", "18-64", "18-64", "65-84", "65-84", ">=85"),
      levels = c("0-27 days (Newborns)", "28 days to 23 months (Infants and Toddlers)",
                 "2-11 (Children)", "12-17 (Adolescents)", "18-64", "65-84", ">=85"),
      exclude = NULL
    )
  )

  expect_dfs_equal(derive_var_agegr_ema(input, AGE, age_unit = "years", AGEGR_EXP), expected_output,
                   keys = "AGE")

})

test_that("derive_var_agegr_ema - works as expected", {

  input <- tibble::tibble(AGE = c(1, 2, 11, 12, 17, 18))

  expected_output <- mutate(
    input,
    AGEGR_EXP = factor(
      c("28 days to 23 months (Infants and Toddlers)", "2-11 (Children)",
        "2-11 (Children)", "12-17 (Adolescents)", "12-17 (Adolescents)",
        "18-64"),
      levels = c("0-27 days (Newborns)", "28 days to 23 months (Infants and Toddlers)",
                 "2-11 (Children)", "12-17 (Adolescents)", "18-64", "65-84", ">=85"),
      exclude = NULL
    )
  )

  expect_dfs_equal(
    derive_var_agegr_ema(input, AGE, age_unit = "years", AGEGR_EXP),
    expected_output,
    keys = "AGE"
  )
})


test_that("derive_var_agegr_ema works with age_unit missing and multiple units in AGEU (adults)", {

  input <- tibble::tibble(AGE = c(10, 18, 19, 50, 64, 65, 80, 85),
                          AGEU = c("years", "years", "years", "years", "years", "years",
                                   "months", "years"))

  expected_output <- mutate(
    input,
    AGEGR_EXP = factor(
      c("2-11 (Children)", "18-64", "18-64", "18-64", "18-64", "65-84", "2-11 (Children)", ">=85"),
      levels = c("0-27 days (Newborns)", "28 days to 23 months (Infants and Toddlers)",
                 "2-11 (Children)", "12-17 (Adolescents)", "18-64", "65-84", ">=85"),
      exclude = NULL
    )
  )

  expect_dfs_equal(derive_var_agegr_ema(input, AGE, new_var = AGEGR_EXP), expected_output,
                   keys = "AGE")
})

test_that("derive_var_agegr_ema - works with age_unit missing and multiple units in AGEU (all)", {

  input <- tibble::tibble(AGE = c(1, 2, 11, 12, 17, 18, 36, 72, 3),
                          AGEU = c("years", "years", "years", "years", "years", "years", "months",
                                   "months", "weeks"))

  expected_output <- mutate(
    input,
    AGEGR_EXP = factor(
      c("28 days to 23 months (Infants and Toddlers)", "2-11 (Children)",
        "2-11 (Children)", "12-17 (Adolescents)", "12-17 (Adolescents)",
        "18-64", "2-11 (Children)", "2-11 (Children)", "0-27 days (Newborns)"),
      levels = c("0-27 days (Newborns)", "28 days to 23 months (Infants and Toddlers)",
                 "2-11 (Children)", "12-17 (Adolescents)", "18-64", "65-84", ">=85"),
      exclude = NULL
    )
  )

  expect_dfs_equal(
    derive_var_agegr_ema(input, AGE, new_var = AGEGR_EXP),
    expected_output,
    keys = "AGE"
  )
})

test_that("derive_var_age_years - Error is thrown when age_unit is not proper unit ", {
  input <- data.frame(AGE = c(12, 24, 36, 48))
  expect_error(
    derive_var_age_years(input, AGE, age_unit = "month", new_var = AAGE),
  "`age_unit` must be one of 'years', 'months', 'weeks', 'days', 'hours', 'minutes' or 'seconds' but is 'month'" # nolint
  )
})

test_that("An Error is issued age_unit is missing", {
  input <- data.frame(AGE = c(12, 24, 36, 48))
  expect_error(
    derive_var_age_years(input, AGE, new_var = AAGE)
  )
  })

















