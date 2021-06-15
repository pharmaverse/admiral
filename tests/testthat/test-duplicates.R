context("test-duplicates")

test_that("multiplication works", {
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
    get_duplicate_records(input, vars(USUBJID))
  )
})
