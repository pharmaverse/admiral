context("test-derive_reference_ranges")

test_that("reference ranges are mapped correctly", {
  ref_ranges <- tibble::tribble(
    ~PARAMCD, ~ANRLO, ~ANRHI, ~A1LO, ~A1HI,
    "DIABP",  60,      80,    40,     90,
    "PUL",    60,     100,    40,    110
  )
  expected_output <- tibble::tribble(
    ~USUBJID, ~PARAMCD, ~ASEQ, ~AVAL, ~ANRLO, ~ANRHI, ~A1LO, ~A1HI,
    "P01",    "PUL",    1,      70,   60,     100,    40,    110,
    "P01",    "PUL",    2,      67,   60,     100,    40,    110,
    "P01",    "DIABP",  1,     102,   60,      80,    40,     90,
    "P02",    "DIABP",  1,      80,   60,      80,    40,     90
  )
  input <- select(expected_output, USUBJID, PARAMCD, ASEQ, AVAL)

  expect_dfs_equal(
    derive_reference_ranges(input, ref_ranges),
    expected_output,
    keys = c("USUBJID", "PARAMCD", "ASEQ")
  )
})

test_that("a warning is issued if reference ranges are missing for a parameter", {
  ref_ranges <- tibble::tribble(
    ~PARAMCD, ~ANRLO, ~ANRHI, ~A1LO, ~A1HI,
    "DIABP",  60,      80,    40,     90
  )
  input <- tibble::tribble(
    ~USUBJID, ~PARAMCD, ~ASEQ, ~AVAL,
    "P01",    "PUL",    1,      70,
    "P01",    "TEMP",   1,      36.8,
    "P02",    "DIABP",  1,      80
  )

  expect_warning(
    derive_reference_ranges(input, ref_ranges),
    "Reference ranges are missing for the following `PARAMCD`: 'PUL' and 'TEMP'"
  )
})
