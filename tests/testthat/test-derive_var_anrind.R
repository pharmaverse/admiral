# derive_var_anrind ---

## Test 1: two-sided reference ranges work ----
test_that("derive_var_anrind Test 1: two-sided reference ranges work", {
  expected_output <- tibble::tribble(
    ~USUBJID, ~PARAMCD, ~ASEQ, ~AVAL, ~ANRLO, ~ANRHI, ~A1LO, ~A1HI,     ~ANRIND,
    "P01",    "PUL",        1,    70,     60,    100,    40,   110,    "NORMAL",
    "P01",    "PUL",        2,    57,     60,    100,    40,   110,       "LOW",
    "P01",    "PUL",        3,    60,     60,    100,    40,   110,    "NORMAL",
    "P01",    "DIABP",      1,   102,     60,     80,    40,    90, "HIGH HIGH",
    "P02",    "PUL",        1,   109,     60,    100,    40,   110,      "HIGH",
    "P02",    "PUL",        2,   100,     60,    100,    40,   110,    "NORMAL",
    "P02",    "DIABP",      1,    80,     60,     80,    40,    90,    "NORMAL",
    "P03",    "PUL",        1,    39,     60,    100,    40,   110,   "LOW LOW",
    "P03",    "PUL",        2,    40,     60,    100,    40,   110,       "LOW"
  )
  input <- select(expected_output, USUBJID:A1HI)

  expect_dfs_equal(
    derive_var_anrind(input, use_a1hia1lo = TRUE),
    expected_output,
    keys = c("USUBJID", "PARAMCD", "ASEQ")
  )
})

## Test 2: explicitly requesting to use `A1LO` and `A1HI` works ----
test_that("derive_var_anrind Test 2: explicitly requesting to use `A1LO` and `A1HI` works", {
  expected_output <- tibble::tribble(
    ~USUBJID, ~PARAMCD, ~ASEQ, ~AVAL, ~ANRLO, ~ANRHI, ~A1LO, ~A1HI,     ~ANRIND,
    "P01",    "PUL",        1,    69,     60,    100,    40,   110,    "NORMAL",
    "P01",    "PUL",        2,    55,     60,    100,    40,   110,       "LOW",
    "P01",    "PUL",        3,    60,     60,    100,    40,   110,    "NORMAL",
    "P01",    "DIABP",      1,   102,     60,     80,    40,    90, "HIGH HIGH",
    "P02",    "PUL",        1,   107,     60,    100,    40,   110,      "HIGH",
    "P02",    "PUL",        2,   100,     60,    100,    40,   110,    "NORMAL",
    "P02",    "DIABP",      1,    51,     60,     80,    40,    90,       "LOW",
    "P03",    "PUL",        1,    32,     60,    100,    40,   110,   "LOW LOW",
    "P03",    "PUL",        2,   107,     60,    100,    40,   110,      "HIGH"
  )
  input <- select(expected_output, USUBJID:A1HI)

  expect_dfs_equal(
    derive_var_anrind(input, use_a1hia1lo = TRUE),
    expected_output,
    keys = c("USUBJID", "PARAMCD", "ASEQ")
  )
})

## Test 3: implicitly missing extreme ranges are supported ----
test_that("derive_var_anrind Test 3: implicitly missing extreme ranges are supported", {
  expected_output <- tibble::tribble(
    ~USUBJID, ~PARAMCD, ~ASEQ, ~AVAL, ~ANRLO, ~ANRHI,  ~ANRIND,
    "P01",    "PUL",        1,    70,     60,    100, "NORMAL",
    "P01",    "PUL",        2,    57,     60,    100,    "LOW",
    "P01",    "PUL",        3,    60,     60,    100, "NORMAL",
    "P01",    "DIABP",      1,   102,     60,     80,   "HIGH",
    "P02",    "PUL",        1,   109,     60,    100,   "HIGH",
    "P02",    "PUL",        2,   100,     60,    100, "NORMAL",
    "P02",    "DIABP",      1,    80,     60,     80, "NORMAL",
    "P03",    "PUL",        1,    39,     60,    100,    "LOW",
    "P03",    "PUL",        2,    40,     60,    100,    "LOW"
  )
  input <- select(expected_output, USUBJID:ANRHI)

  expect_dfs_equal(
    derive_var_anrind(input, use_a1hia1lo = FALSE),
    expected_output,
    keys = c("USUBJID", "PARAMCD", "ASEQ")
  )
})

## Test 4: explicitly missing extreme ranges are supported ----
test_that("derive_var_anrind Test 4: explicitly missing extreme ranges are supported", {
  expected_output <- tibble::tribble(
    ~USUBJID, ~PARAMCD, ~ASEQ, ~AVAL, ~ANRLO, ~ANRHI, ~A1LO, ~A1HI,     ~ANRIND,
    "P01",    "PUL",        1,    70,     60,    100,    NA,    NA,    "NORMAL",
    "P01",    "PUL",        2,    57,     60,    100,    NA,    NA,       "LOW",
    "P01",    "PUL",        3,    60,     60,    100,    NA,    NA,    "NORMAL",
    "P01",    "DIABP",      1,   102,     60,     80,    40,    90, "HIGH HIGH",
    "P02",    "PUL",        1,   109,     60,    100,    NA,    NA,      "HIGH",
    "P02",    "PUL",        2,   100,     60,    100,    NA,    NA,    "NORMAL",
    "P02",    "DIABP",      1,    80,     60,     80,    40,    90,    "NORMAL",
    "P03",    "PUL",        1,    39,     60,    100,    NA,    NA,       "LOW",
    "P03",    "PUL",        2,    40,     60,    100,    NA,    NA,       "LOW"
  )
  input <- select(expected_output, USUBJID:A1HI)

  expect_dfs_equal(
    derive_var_anrind(input, use_a1hia1lo = TRUE),
    expected_output,
    keys = c("USUBJID", "PARAMCD", "ASEQ")
  )
})

## Test 5: one-sided reference ranges work ----
test_that("derive_var_anrind Test 5: one-sided reference ranges work", {
  expected_output <- tibble::tribble(
    ~USUBJID, ~PARAMCD, ~ASEQ, ~AVAL, ~ANRLO, ~ANRHI, ~A1LO, ~A1HI,     ~ANRIND,
    "P01",    "PUL",        1,   101,     NA,    100,    NA,   120,      "HIGH",
    "P01",    "PUL",        2,    99,     NA,    100,    NA,   120,    "NORMAL",
    "P01",    "PUL",        3,   123,     NA,    100,    NA,   120, "HIGH HIGH",
    "P01",    "DIABP",      1,   102,     60,     NA,    40,    NA,    "NORMAL",
    "P02",    "PUL",        1,   109,     NA,    100,    NA,   120,      "HIGH",
    "P02",    "PUL",        2,   100,     NA,    100,    NA,   120,    "NORMAL",
    "P02",    "DIABP",      1,    58,     60,     NA,    40,    NA,       "LOW",
    "P03",    "PUL",        1,    39,     NA,    100,    NA,   120,    "NORMAL",
    "P03",    "PUL",        2,    40,     NA,    100,    NA,   120,    "NORMAL"
  )
  input <- select(expected_output, USUBJID:A1HI)

  expect_dfs_equal(
    derive_var_anrind(input, use_a1hia1lo = TRUE),
    expected_output,
    keys = c("USUBJID", "PARAMCD", "ASEQ")
  )
})

## Test 6: missing `AVAL` is handled properly ----
test_that("derive_var_anrind Test 6: missing `AVAL` is handled properly", {
  expected_output <- tibble::tribble(
    ~USUBJID, ~PARAMCD, ~ASEQ,    ~AVAL, ~ANRLO, ~ANRHI,          ~ANRIND,
    "P01",    "PUL",    1,     NA_real_,     60,    100,    NA_character_
  )
  input <- select(expected_output, USUBJID:ANRHI)

  expect_dfs_equal(
    derive_var_anrind(input),
    expected_output,
    keys = c("USUBJID", "PARAMCD", "ASEQ")
  )
})

# Show floating point issue ----
expected_output_fp <- tibble::tribble(
  ~USUBJID, ~PARAMCD, ~ASEQ, ~AVAL, ~ANRLO, ~ANRHI,  ~ANRIND,
  "P01",       "PUL",     1,   100,     60,    110, "NORMAL",
) %>%
  mutate(AVAL = 1.1 * AVAL)

input_fp <- select(expected_output_fp, USUBJID:ANRHI)


## Test 7: Show floating points handled properly ----
test_that("derive_var_anrind Test 7: missing `AVAL` is handled properly", {
  expect_dfs_equal(
    derive_var_anrind(input_fp),
    expected_output_fp,
    keys = c("USUBJID", "PARAMCD", "ASEQ")
  )
})

## Test 8: Show floating points not handled correctly ----
test_that("derive_var_anrind Test 8: Show floating points not handled correctly", {
  # when SIGNIFICANT DIGITS = 17 then AVAL > ANRHI
  # even though 1.1 * 100 should equal 110
  expected_output_fpx <- expected_output_fp %>%
    mutate(ANRIND = "HIGH")

  expect_dfs_equal(
    derive_var_anrind(input_fp, signif_dig = 17),
    expected_output_fpx,
    keys = c("USUBJID", "PARAMCD", "ASEQ")
  )
})
