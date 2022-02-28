test_that("Shift based on character variables", {
  input <- tibble::tribble(
    ~USUBJID, ~PARAMCD, ~AVAL, ~ABLFL, ~BNRIND, ~ANRIND,
    "P01",    "ALB",     33,    "Y",    "LOW",   "LOW",
    "P01",    "ALB",     38,    "",     "LOW",   "NORMAL",
    "P02",    "ALB",     37,    "Y",    "NORMAL", "NORMAL",
    "P02",    "ALB",     49,    "",     "NORMAL", "HIGH",
    "P02",    "SODIUM",  147,   "Y",    "HIGH",  "HIGH"
  )
  expected_output <- tibble::tribble(
    ~USUBJID, ~PARAMCD, ~AVAL, ~ABLFL, ~BNRIND, ~ANRIND, ~SHIFT1,
    "P01",    "ALB",     33,    "Y",    "LOW",   "LOW",   "LOW to LOW",
    "P01",    "ALB",     38,    "",     "LOW",   "NORMAL", "LOW to NORMAL",
    "P02",    "ALB",     37,    "Y",    "NORMAL", "NORMAL", "NORMAL to NORMAL",
    "P02",    "ALB",     49,    "",     "NORMAL", "HIGH",  "NORMAL to HIGH",
    "P02",    "SODIUM",  147,   "Y",    "HIGH",  "HIGH",  "HIGH to HIGH"
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


test_that("Shift based on character variables with missing values", {
  input <- tibble::tribble(
    ~USUBJID, ~PARAMCD, ~AVAL, ~ABLFL, ~BNRIND, ~ANRIND,
    "P01",    "ALB",     33,    "Y",    "LOW",   "LOW",
    "P01",    "ALB",     38,    "",     "LOW",   "NORMAL",
    "P01",    "ALB",     NA,    "",     "LOW",   "",
    "P02",    "ALB",     NA,    "Y",    "",      "",
    "P02",    "ALB",     49,    "",     "",      "HIGH",
    "P02",    "SODIUM",  147,   "Y",    "HIGH",  "HIGH"
  ) %>% convert_blanks_to_na
  expected_output <- tibble::tribble(
    ~USUBJID, ~PARAMCD, ~AVAL, ~ABLFL, ~BNRIND, ~ANRIND, ~SHIFT1,
    "P01",    "ALB",     33,    "Y",    "LOW",   "LOW",   "LOW to LOW",
    "P01",    "ALB",     38,    "",     "LOW",   "NORMAL", "LOW to NORMAL",
    "P01",    "ALB",     NA,    "",     "LOW",   "",      "LOW to NULL",
    "P02",    "ALB",     NA,    "Y",    "",      "",      "NULL to NULL",
    "P02",    "ALB",     49,    "",     "",      "HIGH",  "NULL to HIGH",
    "P02",    "SODIUM",  147,   "Y",    "HIGH",  "HIGH",  "HIGH to HIGH"
  ) %>% convert_blanks_to_na

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


test_that("Shift based on numeric variables with missing values", {
  input <- tibble::tribble(
    ~USUBJID, ~PARAMCD, ~AVAL, ~ABLFL, ~BASE,
    "P01",    "ALB",     33.1,  "Y",  33.1,
    "P01",    "ALB",     38.5,  "",   33.1,
    "P01",    "ALB",     NA,    "",   33.1,
    "P02",    "ALB",     NA,    "Y",  NA,
    "P02",    "ALB",     49.0,  "",   NA,
    "P02",    "SODIUM",  147.5, "Y",  147.5
  ) %>% convert_blanks_to_na
  expected_output <- tibble::tribble(
    ~USUBJID, ~PARAMCD, ~AVAL, ~ABLFL, ~BASE, ~SHIFT2,
    "P01",    "ALB",     33.1,  "Y",  33.1,    "33.1 to 33.1",
    "P01",    "ALB",     38.5,  "",   33.1,    "33.1 to 38.5",
    "P01",    "ALB",     NA,    "",   33.1,    "33.1 to NULL",
    "P02",    "ALB",     NA,    "Y",  NA,      "NULL to NULL",
    "P02",    "ALB",     49,    "",   NA,      "NULL to 49",
    "P02",    "SODIUM",  147.5, "Y",  147.5,   "147.5 to 147.5"
  ) %>% convert_blanks_to_na

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


test_that("Only populate post-baseline", {
  input <- tibble::tribble(
    ~USUBJID, ~PARAMCD, ~AVAL, ~ABLFL, ~BASE,
    "P01",    "ALB",     33.1,  "Y",  33.1,
    "P01",    "ALB",     38.5,  "",   33.1,
    "P01",    "ALB",     NA,    "",   33.1,
    "P02",    "ALB",     NA,    "Y",  NA,
    "P02",    "ALB",     49.0,  "",   NA,
    "P02",    "SODIUM",  147.5, "Y",  147.5
  ) %>% convert_blanks_to_na
  expected_output <- tibble::tribble(
    ~USUBJID, ~PARAMCD, ~AVAL, ~ABLFL, ~BASE, ~SHIFT2,
    "P01",    "ALB",     33.1,  "Y",  33.1,    NA_character_,
    "P01",    "ALB",     38.5,  "",   33.1,    "33.1 to 38.5",
    "P01",    "ALB",     NA,    "",   33.1,    "33.1 to NULL",
    "P02",    "ALB",     NA,    "Y",  NA,      NA_character_,
    "P02",    "ALB",     49,    "",   NA,      "NULL to 49",
    "P02",    "SODIUM",  147.5, "Y",  147.5,   NA_character_
  ) %>% convert_blanks_to_na

  expect_equal(
    derive_var_shift(
      input,
      new_var = SHIFT2,
      from_var = BASE,
      to_var = AVAL,
      post_baseline_condition = (ABLFL != "Y")
    ),
    expected_output
  )
})
