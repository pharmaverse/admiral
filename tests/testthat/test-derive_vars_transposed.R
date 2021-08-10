dataset <- tibble::tribble(
  ~USUBJID, ~VAR1,
  "P01",    3,
  "P02",    31,
  "P03",    42
)
dataset_merge <- tibble::tribble(
  ~USUBJID, ~TESTCD, ~VALUE,
  "P01",    "T01",   31,
  "P01",    "T02",   5,
  "P02",    "T01",   3,
  "P03",    "T02",   9
)

test_that("multiplication works", {
  expected_output <- tibble::tribble(
    ~USUBJID, ~VAR1, ~T01, ~T02,
    "P01",     3,    31,   5,
    "P02",    31,     3,   NA,
    "P03",    42,    NA,   9
  )
  actual_output <- derive_vars_transposed(
    dataset,
    dataset_merge,
    by_vars = vars(USUBJID),
    key_var = TESTCD,
    value_var = VALUE
  )

  expect_dfs_equal(expected_output, actual_output, keys = "USUBJID")
})


test_that("multiplication works", {
  expected_output <- tibble::tribble(
    ~USUBJID, ~VAR1, ~T01,
    "P01",     3,    31,
    "P02",    31,     3,
    "P03",    42,    NA
  )
  actual_output <- derive_vars_transposed(
    dataset,
    dataset_merge,
    by_vars = vars(USUBJID),
    key_var = TESTCD,
    value_var = VALUE,
    filter = TESTCD == "T01"
  )

  expect_dfs_equal(expected_output, actual_output, keys = "USUBJID")
})
