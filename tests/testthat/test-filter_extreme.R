context("test-filter_extreme")


test_that("first observation for each group are selected", {
  input <- tibble::tribble(
    ~USUBJID, ~AVISITN, ~AVAL,
    1, 1, 12,
    1, 3, 9,
    2, 2, 42,
    3, 3, 14,
    3, 3, 10)

  expected_output <- tibble::tribble(
    ~USUBJID, ~AVISITN, ~AVAL,
    1, 1, 12,
    2, 2, 42,
    3, 3, 10)

  expect_equal(filter_extreme(input,
                              order = exprs(AVISITN, AVAL),
                              by_vars = exprs(USUBJID),
                              mode = 'first'),
               expected_output)
})

test_that("first observation is selected without grouping", {
  input <- tibble::tribble(
    ~USUBJID, ~AVISITN, ~AVAL,
    1, 1, 12,
    1, 3, 9,
    2, 2, 42,
    3, 3, 14,
    3, 3, 10)

  expected_output <- tibble::tribble(
    ~USUBJID, ~AVISITN, ~AVAL,
    1, 1, 12)

  expect_equal(filter_extreme(input,
                              order = exprs(AVISITN, AVAL),
                              mode = 'first'),
               expected_output)
})
