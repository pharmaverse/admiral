context("test-derive_reference_ranges")

test_that("reference ranges are mapped correctly", {
  ref_ranges <- tibble::tribble(
    ~PARAMCD, ~ANRLO, ~ANRHI, ~A1LO, ~A1HI,
    "DIABP",  60,      80,    40,     90,
    "PUL",    60,     100,    40,    110
  )
  expected_output <- tibble::tribble(
    ~USUBJID, ~PARAMCD, ~ASEQ, ~AVAL, ~ANRLO, ~ANRHI, ~A1LO, ~A1HI, ~ANRIND,
    "P01",    "PUL",    1,      70,   60,     100,    40,    110,   "NORMAL",
    "P01",    "PUL",    2,      57,   60,     100,    40,    110,   "LOW",
    "P01",    "PUL",    3,      60,   60,     100,    40,    110,   "NORMAL",
    "P01",    "DIABP",  1,     102,   60,      80,    40,     90,   "HIGH HIGH",
    "P02",    "PUL",    1,     109,   60,     100,    40,    110,   "HIGH",
    "P02",    "PUL",    2,     100,   60,     100,    40,    110,   "NORMAL",
    "P02",    "DIABP",  1,      80,   60,      80,    40,     90,   "NORMAL",
    "P03",    "PUL",    1,      39,   60,     100,    40,    110,   "LOW LOW",
    "P03",    "PUL",    2,      40,   60,     100,    40,    110,   "LOW"
  )
  input <- select(expected_output, USUBJID:AVAL)

  expect_dfs_equal(
    derive_reference_ranges(input, ref_ranges),
    expected_output,
    keys = c("USUBJID", "PARAMCD", "ASEQ")
  )
})

test_that("reference ranges are mapped correctly 2", {
  ref_ranges <- tibble::tribble(
    ~PARAMCD, ~ANRLO, ~ANRHI,
    "DIABP",  60,      80,
    "PUL",    60,     100
  )
  expected_output <- tibble::tribble(
    ~USUBJID, ~PARAMCD, ~ASEQ, ~AVAL, ~ANRLO, ~ANRHI, ~ANRIND,
    "P01",    "PUL",    1,      70,   60,     100,    "NORMAL",
    "P01",    "PUL",    2,      57,   60,     100,    "LOW",
    "P01",    "PUL",    3,      60,   60,     100,    "NORMAL",
    "P01",    "DIABP",  1,     102,   60,      80,    "HIGH",
    "P02",    "PUL",    1,     109,   60,     100,    "HIGH",
    "P02",    "PUL",    2,     100,   60,     100,    "NORMAL",
    "P02",    "DIABP",  1,      80,   60,      80,    "NORMAL",
    "P03",    "PUL",    1,      39,   60,     100,    "LOW",
    "P03",    "PUL",    2,      40,   60,     100,    "LOW"
  )
  input <- select(expected_output, USUBJID:AVAL)

  expect_dfs_equal(
    derive_reference_ranges(input, ref_ranges),
    expected_output,
    keys = c("USUBJID", "PARAMCD", "ASEQ")
  )
})

test_that("reference ranges are mapped correctly 3", {
  ref_ranges <- tibble::tribble(
    ~PARAMCD, ~ANRLO, ~ANRHI, ~A1LO, ~A1HI,
    "DIABP",  60,      80,    40,    90,
    "PUL",    60,     100,    NA,    NA
  )
  expected_output <- tibble::tribble(
    ~USUBJID, ~PARAMCD, ~ASEQ, ~AVAL, ~ANRLO, ~ANRHI, ~A1LO, ~A1HI, ~ANRIND,
    "P01",    "PUL",    1,      70,   60,     100,    NA,    NA,    "NORMAL",
    "P01",    "PUL",    2,      57,   60,     100,    NA,    NA,    "LOW",
    "P01",    "PUL",    3,      60,   60,     100,    NA,    NA,    "NORMAL",
    "P01",    "DIABP",  1,     102,   60,      80,    40,    90,    "HIGH HIGH",
    "P02",    "PUL",    1,     109,   60,     100,    NA,    NA,    "HIGH",
    "P02",    "PUL",    2,     100,   60,     100,    NA,    NA,    "NORMAL",
    "P02",    "DIABP",  1,      80,   60,      80,    40,    90,    "NORMAL",
    "P03",    "PUL",    1,      39,   60,     100,    NA,    NA,    "LOW",
    "P03",    "PUL",    2,      40,   60,     100,    NA,    NA,    "LOW"
  )
  input <- select(expected_output, USUBJID:AVAL)

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
