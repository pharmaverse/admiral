context("test-derive_extreme_flag")


test_that("first observation for each group is flagged", {
  input <- tibble::tribble(
    ~USUBJID, ~AVISITN, ~AVAL,
    1, 1, 12,
    1, 3, 9,
    2, 2, 42,
    3, 3, 14,
    3, 3, 10)

  expected_output <- input %>% mutate(firstfl := c("Y", NA, "Y", "Y", NA))

  actual_output <- derive_extreme_flag(input,
                                       new_var = firstfl,
                                       order = exprs(AVISITN, desc(AVAL)),
                                       by_vars = exprs(USUBJID),
                                       mode = "first")

  expect_dfs_equal(base = expected_output,
                   compare = actual_output,
                   keys = c("USUBJID", "AVISITN", "AVAL"))
})

test_that("last observation for each group is flagged, flag_filter works", {
  input <- tibble::tribble(
    ~USUBJID, ~AVISITN, ~AVAL,
    1, 1, 12,
    1, 3, 9,
    2, 2, 42,
    3, 3, 14,
    3, 3, 10)

  expected_output <- input %>% mutate(lastfl := c(NA, "Y", NA, NA, "Y"))

  actual_output <- derive_extreme_flag(input,
                                       new_var = lastfl,
                                       order = exprs(AVISITN, desc(AVAL)),
                                       by_vars = exprs(USUBJID),
                                       mode = "last",
                                       flag_filter = expr(USUBJID != 2))
  expect_dfs_equal(base = expected_output,
                   compare = actual_output,
                   keys = c("USUBJID", "AVISITN", "AVAL"))
})
