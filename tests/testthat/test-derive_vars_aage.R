test_that("duration and unit variable are added", {
  input <- tibble::tribble(
    ~BRTHDT, ~RANDDT,
    ymd("1999-09-09"), ymd("2020-02-20")
  )
  expected_output <- mutate(input, AAGE = 20, AAGEU = "YEARS")

  expect_dfs_equal(derive_vars_aage(input), expected_output, keys = c("BRTHDT", "RANDDT"))
})


test_that("derive_var_age_years works as expected", {
  input <- tibble::tibble(
    AGE = c(12, 24, 36, 48, 60),
    AGEU = c("months", "months", "months", "months", "months")
  )

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

test_that("derive_var_age_years - Error is issued if age_unit is missing", {
  input <- data.frame(AGE = c(12, 24, 36, 48))
  expect_error(
    derive_var_age_years(input, AGE, new_var = AAGE)
  )
})

test_that(paste0("derive_var_age_years - Warning is issued if age_unit is not null, but the ",
                 "'unit' variable corresponding to age_var stores more than one unique ",
                 "value."), {
  input <- tibble::tribble(
    ~AGE,   ~AGEU,
    #-------/---------
    25,     "years",
    312,    "months",
    51,     "years",
    402,    "months",
    432,    "months"
  )

  expect_warning(
    derive_var_age_years(input, AGE, age_unit = "months", new_var = AAGE)
  )
})

test_that(paste0("derive_var_age_years - Error is issued if age_unit consists of more than ",
                 "one unique value."), {
  input <- tibble::tribble(
    ~AGE,   ~AGEU,
    #-------/---------
    459,    "months",
    312,    "months",
    510,    "months",
    402,    "months",
    432,    "months"
  )

  expect_error(
    derive_var_age_years(input, AGE, age_unit = c("months", "years"), new_var = AAGE)
  )
})

test_that(paste0("derive_var_age_years - The 'unit' variable corresponding to age_var will ",
                 "be considered as storing one unique unit, if values only differ by case, ",
                 "i.e. 'months', 'Months', 'MONTHS' considered same unit, etc."), {

  # The tibbles "input" and "input2" differ only in the third row: "Months" versus "months".

  input <- tibble::tribble(
    ~AGE,   ~AGEU,
    #-------/---------
    459,    "months",
    312,    "months",
    510,    "Months",
    402,    "months",
    432,    "months"
  )

  input2 <- tibble::tribble(
    ~AGE,   ~AGEU,
    #-------/---------
    459,    "months",
    312,    "months",
    510,    "months",
    402,    "months",
    432,    "months"
  )

  expect_equal(
    derive_var_age_years(input, AGE, age_unit = "months", new_var = AAGE)$AAGE,
    derive_var_age_years(input2, AGE, age_unit = "months", new_var = AAGE)$AAGE
  )
})

test_that(paste0("derive_var_age_years - Warning is issued if age_unit is not null, but the ",
                 "'unit' variable corresponding to age_var stores one unique unit that is ",
                 "not equivalent to age_unit."), {
  input <- tibble::tribble(
    ~AGE,   ~AGEU,
    #-------/---------
    459,    "months",
    312,    "months",
    510,    "months",
    402,    "months",
    432,    "months"
  )

  expect_warning(
    derive_var_age_years(input, AGE, age_unit = "years", new_var = AAGE)
  )
})

