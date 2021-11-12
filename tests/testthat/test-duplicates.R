test_that("duplicate records are extracted", {
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
    extract_duplicate_records(input, vars(USUBJID))
  )
})

test_that("dataset of duplicate records can be accessed using `get_duplicates_dataset()`", {
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
    signal_duplicate_records(input, vars(USUBJID)),
    "Dataset contains duplicate records with respect to `USUBJID`",
    fixed = TRUE
  )

  expect_true(expected_ouput == get_duplicates_dataset())
})
