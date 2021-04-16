context("test-derive_extreme_flag")


test_that("first observation for each group is flagged", {
  input <- tibble::tribble(
    ~USUBJID, ~AVISITN, ~AVAL,
    1, 1, 12,
    1, 3, 9,
    2, 2, 42,
    3, 3, 14,
    3, 3, 10)

  expected_output <- input %>% mutate(firstfl := c("Y", "", "Y", "Y", ""))

  expect_equal(derive_extreme_flag(input,
                                   new_var = firstfl,
                                   order = exprs(AVISITN, desc(AVAL)),
                                   by_vars = exprs(USUBJID),
                                   mode = "first"),
               expected_output)
})

test_that("last observation for each group is flagged", {
  input <- tibble::tribble(
    ~USUBJID, ~AVISITN, ~AVAL,
    1, 1, 12,
    1, 3, 9,
    2, 2, 42,
    3, 3, 14,
    3, 3, 10)

  expected_output <- input %>% mutate(lastfl := c("", "Y", "Y", "", "Y"))

  expect_equal(derive_extreme_flag(input,
                                   new_var = lastfl,
                                   order = exprs(AVISITN, desc(AVAL)),
                                   by_vars = exprs(USUBJID)),
               expected_output)
})
