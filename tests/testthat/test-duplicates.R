# extract_duplicate_records ----
## Test 1: duplicate records are extracted ----
test_that("extract_duplicate_records Test 1: duplicate records are extracted", {
  input <- tibble::tribble(
    ~USUBJID, ~COUNTRY, ~AAGE,
    "P01",    "GER",    22,
    "P01",    "JPN",    34,
    "P02",    "CZE",    41,
    "P03",    "AUS",    39,
    "P04",    "BRA",    21,
    "P04",    "BRA",    21
  )
  expected_ouput <- input[c(1:2, 5:6), ]

  expect_equal(
    expected_ouput,
    extract_duplicate_records(input, exprs(USUBJID))
  )
})

## Test 2: duplicate records for all variables ----
test_that("extract_duplicate_records Test 2: duplicate records for all variables", {
  input <- tibble::tribble(
    ~USUBJID, ~COUNTRY, ~AAGE,
    "P01",    "GER",    22,
    "P01",    "JPN",    34,
    "P02",    "CZE",    41,
    "P03",    "AUS",    39,
    "P04",    "BRA",    21,
    "P04",    "BRA",    21
  )
  expected_ouput <- input[c(5:6), ]

  expect_equal(
    expected_ouput,
    extract_duplicate_records(input)
  )
})


# signal_duplicate_records ----
## Test 3: dataset of duplicate records can be accessed using `get_duplicates_dataset()` ----
test_that("signal_duplicate_records Test 3: dataset of duplicate records can be accessed using `get_duplicates_dataset()`", { # nolint
  input <- tibble::tribble(
    ~USUBJID, ~COUNTRY, ~AAGE,
    "P01",    "GER",    22,
    "P01",    "JPN",    34,
    "P02",    "CZE",    41,
    "P03",    "AUS",    39,
    "P04",    "BRA",    21,
    "P04",    "BRA",    21
  )
  expected_ouput <- input[c(1:2, 5:6), ]

  expect_error(
    signal_duplicate_records(input, exprs(USUBJID)),
    "Dataset contains duplicate records with respect to `USUBJID`",
    fixed = TRUE
  )

  expect_true(all(expected_ouput == get_duplicates_dataset()))

  expect_snapshot(get_duplicates_dataset())
})
