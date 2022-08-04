test_that("two-sided reference ranges work", {
  library(tibble)

  expected_output <- tribble(
    ~USUBJID, ~PARAMCD, ~ASEQ, ~AVAL, ~ANRLO, ~ANRHI, ~A1LO, ~A1HI, ~ANRIND,
    "P01", "PUL", 1, 70, 60, 100, 40, 110, "NORMAL",
    "P01", "PUL", 2, 57, 60, 100, 40, 110, "LOW",
    "P01", "PUL", 3, 60, 60, 100, 40, 110, "NORMAL",
    "P01", "DIABP", 1, 102, 60, 80, 40, 90, "HIGH HIGH",
    "P02", "PUL", 1, 109, 60, 100, 40, 110, "HIGH",
    "P02", "PUL", 2, 100, 60, 100, 40, 110, "NORMAL",
    "P02", "DIABP", 1, 80, 60, 80, 40, 90, "NORMAL",
    "P03", "PUL", 1, 39, 60, 100, 40, 110, "LOW LOW",
    "P03", "PUL", 2, 40, 60, 100, 40, 110, "LOW"
  )
  input <- select(expected_output, USUBJID:A1HI)

  expect_dfs_equal(
    derive_var_anrind(input),
    expected_output,
    keys = c("USUBJID", "PARAMCD", "ASEQ")
  )
})

test_that("implicitly missing extreme ranges are supported", {
  library(tibble)

  expected_output <- tribble(
    ~USUBJID, ~PARAMCD, ~ASEQ, ~AVAL, ~ANRLO, ~ANRHI, ~ANRIND,
    "P01", "PUL", 1, 70, 60, 100, "NORMAL",
    "P01", "PUL", 2, 57, 60, 100, "LOW",
    "P01", "PUL", 3, 60, 60, 100, "NORMAL",
    "P01", "DIABP", 1, 102, 60, 80, "HIGH",
    "P02", "PUL", 1, 109, 60, 100, "HIGH",
    "P02", "PUL", 2, 100, 60, 100, "NORMAL",
    "P02", "DIABP", 1, 80, 60, 80, "NORMAL",
    "P03", "PUL", 1, 39, 60, 100, "LOW",
    "P03", "PUL", 2, 40, 60, 100, "LOW"
  )
  input <- select(expected_output, USUBJID:ANRHI)

  expect_dfs_equal(
    derive_var_anrind(input),
    expected_output,
    keys = c("USUBJID", "PARAMCD", "ASEQ")
  )
})

test_that("explicitly missing extreme ranges are supported", {
  library(tibble)

  expected_output <- tribble(
    ~USUBJID, ~PARAMCD, ~ASEQ, ~AVAL, ~ANRLO, ~ANRHI, ~A1LO, ~A1HI, ~ANRIND,
    "P01", "PUL", 1, 70, 60, 100, NA, NA, "NORMAL",
    "P01", "PUL", 2, 57, 60, 100, NA, NA, "LOW",
    "P01", "PUL", 3, 60, 60, 100, NA, NA, "NORMAL",
    "P01", "DIABP", 1, 102, 60, 80, 40, 90, "HIGH HIGH",
    "P02", "PUL", 1, 109, 60, 100, NA, NA, "HIGH",
    "P02", "PUL", 2, 100, 60, 100, NA, NA, "NORMAL",
    "P02", "DIABP", 1, 80, 60, 80, 40, 90, "NORMAL",
    "P03", "PUL", 1, 39, 60, 100, NA, NA, "LOW",
    "P03", "PUL", 2, 40, 60, 100, NA, NA, "LOW"
  )
  input <- select(expected_output, USUBJID:A1HI)

  expect_dfs_equal(
    derive_var_anrind(input),
    expected_output,
    keys = c("USUBJID", "PARAMCD", "ASEQ")
  )
})

test_that("one-sided reference ranges work", {
  library(tibble)

  expected_output <- tribble(
    ~USUBJID, ~PARAMCD, ~ASEQ, ~AVAL, ~ANRLO, ~ANRHI, ~A1LO, ~A1HI, ~ANRIND,
    "P01", "PUL", 1, 101, NA, 100, NA, 120, "HIGH",
    "P01", "PUL", 2, 99, NA, 100, NA, 120, "NORMAL",
    "P01", "PUL", 3, 123, NA, 100, NA, 120, "HIGH HIGH",
    "P01", "DIABP", 1, 102, 60, NA, 40, NA, "NORMAL",
    "P02", "PUL", 1, 109, NA, 100, NA, 120, "HIGH",
    "P02", "PUL", 2, 100, NA, 100, NA, 120, "NORMAL",
    "P02", "DIABP", 1, 58, 60, NA, 40, NA, "LOW",
    "P03", "PUL", 1, 39, NA, 100, NA, 120, "NORMAL",
    "P03", "PUL", 2, 40, NA, 100, NA, 120, "NORMAL"
  )
  input <- select(expected_output, USUBJID:A1HI)

  expect_dfs_equal(
    derive_var_anrind(input),
    expected_output,
    keys = c("USUBJID", "PARAMCD", "ASEQ")
  )
})

test_that("missing `AVAL` is handled properly", {
  library(tibble)

  expected_output <- tribble(
    ~USUBJID, ~PARAMCD, ~ASEQ, ~AVAL,    ~ANRLO, ~ANRHI, ~ANRIND,
    "P01",    "PUL",    1,     NA_real_, 60,     100,    NA_character_
  )
  input <- select(expected_output, USUBJID:ANRHI)

  expect_dfs_equal(
    derive_var_anrind(input),
    expected_output,
    keys = c("USUBJID", "PARAMCD", "ASEQ")
  )
})
