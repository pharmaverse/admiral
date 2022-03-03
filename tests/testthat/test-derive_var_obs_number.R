test_that("first observation for each group are selected", {
  input <- tibble::tribble(
    ~USUBJID, ~AVISITN, ~AVAL,
    1, 1, 12,
    1, 3, 9,
    2, 2, 42,
    3, 3, 14,
    3, 3, 10
  )

  expected_output <- tibble::tribble(
    ~USUBJID, ~AVISITN, ~AVAL,
    1, 1, 12,
    2, 2, 42,
    3, 3, 10
  )

  actual_output <- filter_extreme(
    input,
    by_vars = vars(USUBJID),
    order = vars(AVISITN, AVAL),
    mode = "first"
  )

  expect_dfs_equal(
    base = expected_output,
    compare = actual_output,
    keys = "USUBJID"
  )
})

test_that("first observation is selected without grouping", {
  input <- tibble::tribble(
    ~USUBJID, ~AVISITN, ~AVAL,
    1, 1, 12,
    1, 3, 9,
    2, 2, 42,
    3, 3, 14,
    3, 3, 10
  )

  expected_output <- tibble::tribble(
    ~USUBJID, ~AVISITN, ~AVAL,
    1, 1, 12
  )

  actual_output <- filter_extreme(
    input,
    order = vars(AVISITN, AVAL),
    mode = "first"
  )

  expect_dfs_equal(
    base = expected_output,
    compare = actual_output,
    keys = "USUBJID"
  )
})
