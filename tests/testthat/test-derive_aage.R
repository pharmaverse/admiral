context("test-derive_vars_aage")


test_that("duration and unit variable are added", {
  input <- tibble::tribble(
    ~BRTHDT, ~RANDDT,
    ymd("1999-09-09"), ymd("2020-02-20")
  )
  expected_output <- mutate(input, AAGE = 20, AAGEU = "YEARS")

  expect_dfs_equal(derive_vars_aage(input), expected_output, keys = c("BRTHDT", "RANDDT"))
})


test_that("derive_agegr_fda works as expected", {

  input <- tibble::tibble(AGE = c(10, 18, 19, 50, 64, 65, 80))

  expected_output <- mutate(
    input,
    AGEGR_EXP = factor(
      c("<=18", "<=18", "19-64", "19-64", "19-64", ">=65", ">=65"),
      levels = c("<=18", "19-64", ">=65"),
      exclude = NULL
    )
  )

  expect_dfs_equal(derive_agegr_fda(input, AGE, AGEGR_EXP), expected_output, keys = "AGE")

})

test_that("derive_agegr_ema works as expected", {

  input <- tibble::tibble(AGE = c(10, 18, 19, 50, 64, 65, 80, 85))

  expected_output <- mutate(
    input,
    AGEGR_EXP = factor(
      c(NA, "18-64", "18-64", "18-64", "18-64", "65-84", "65-84", ">=85"),
      levels = c("18-64", "65-84", ">=85", NA_character_),
      exclude = NULL
    )
  )

  expect_dfs_equal(derive_agegr_ema(input, AGE, AGEGR_EXP), expected_output, keys = "AGE")

})

test_that("derive_agegr_ema - pediatric version - works as expected", {

  input <- tibble::tibble(AGE = c(1, 2, 11, 12, 17, 18))

  expected_output <- mutate(
    input,
    AGEGR_EXP = factor(
      c("0-1 (Newborns / Infants / Toddlers)", "2-11 (Children)",
        "2-11 (Children)", "12-17 (Adolescents)", "12-17 (Adolescents)",
        NA),
      levels = c("0-1 (Newborns / Infants / Toddlers)", "2-11 (Children)",
                 "12-17 (Adolescents)", NA_character_),
      exclude = NULL
    )
  )

  expect_dfs_equal(
    derive_agegr_ema(input, AGE, AGEGR_EXP, adults = FALSE),
    expected_output,
    keys = "AGE"
  )

})
