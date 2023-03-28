input_worst_flag <- tibble::tribble(
  ~STUDYID, ~USUBJID, ~PARAMCD, ~AVISIT, ~ADT, ~AVAL,
  "TEST01", "PAT01", "PARAM01", "BASELINE", as.Date("2021-04-27"), 15.0,
  "TEST01", "PAT01", "PARAM01", "BASELINE", as.Date("2021-04-25"), 14.0,
  "TEST01", "PAT01", "PARAM01", "BASELINE", as.Date("2021-04-23"), 15.0,
  "TEST01", "PAT01", "PARAM01", "WEEK 1", as.Date("2021-04-27"), 10.0,
  "TEST01", "PAT01", "PARAM01", "WEEK 2", as.Date("2021-04-30"), 12.0,
  "TEST01", "PAT02", "PARAM01", "SCREENING", as.Date("2021-04-27"), 15.0,
  "TEST01", "PAT02", "PARAM01", "BASELINE", as.Date("2021-04-25"), 14.0,
  "TEST01", "PAT02", "PARAM01", "BASELINE", as.Date("2021-04-23"), 15.0,
  "TEST01", "PAT02", "PARAM01", "WEEK 1", as.Date("2021-04-27"), 10.0,
  "TEST01", "PAT02", "PARAM01", "WEEK 2", as.Date("2021-04-30"), 12.0,
  "TEST01", "PAT01", "PARAM02", "SCREENING", as.Date("2021-04-27"), 15.0,
  "TEST01", "PAT01", "PARAM02", "SCREENING", as.Date("2021-04-25"), 14.0,
  "TEST01", "PAT01", "PARAM02", "SCREENING", as.Date("2021-04-23"), 15.0,
  "TEST01", "PAT01", "PARAM02", "BASELINE", as.Date("2021-04-27"), 10.0,
  "TEST01", "PAT01", "PARAM02", "WEEK 2", as.Date("2021-04-30"), 12.0,
  "TEST01", "PAT02", "PARAM02", "SCREENING", as.Date("2021-04-27"), 15.0,
  "TEST01", "PAT02", "PARAM02", "BASELINE", as.Date("2021-04-25"), 14.0,
  "TEST01", "PAT02", "PARAM02", "WEEK 1", as.Date("2021-04-23"), 15.0,
  "TEST01", "PAT02", "PARAM02", "WEEK 1", as.Date("2021-04-27"), 10.0,
  "TEST01", "PAT02", "PARAM02", "BASELINE", as.Date("2021-04-30"), 12.0,
  "TEST01", "PAT02", "PARAM03", "SCREENING", as.Date("2021-04-27"), 15.0,
  "TEST01", "PAT02", "PARAM03", "BASELINE", as.Date("2021-04-25"), 14.0,
  "TEST01", "PAT02", "PARAM03", "WEEK 1", as.Date("2021-04-23"), 15.0,
  "TEST01", "PAT02", "PARAM03", "WEEK 1", as.Date("2021-04-27"), 10.0,
  "TEST01", "PAT02", "PARAM03", "BASELINE", as.Date("2021-04-30"), 12.0
)

test_that("first observation for each group is flagged", {
  input <- tibble::tribble(
    ~USUBJID, ~AVISITN, ~AVAL,
    1, 1, 12,
    1, 3, 9,
    2, 2, 42,
    3, 3, 14,
    3, 3, 10
  )

  expected_output <- input %>% mutate(firstfl = c("Y", NA, "Y", "Y", NA))

  actual_output <- derive_var_extreme_flag(
    input,
    by_vars = exprs(USUBJID),
    order = exprs(AVISITN, desc(AVAL)),
    new_var = firstfl,
    mode = "first"
  )

  expect_dfs_equal(
    base = expected_output,
    compare = actual_output,
    keys = c("USUBJID", "AVISITN", "AVAL")
  )
})

test_that("last observation for each group is flagged", {
  input <- tibble::tribble(
    ~USUBJID, ~AVISITN, ~AVAL,
    1, 1, 12,
    1, 3, 9,
    2, 2, 42,
    3, 3, 14,
    3, 3, 10
  )

  expected_output <- input %>% mutate(lastfl = c(NA, "Y", "Y", NA, "Y"))

  actual_output <- derive_var_extreme_flag(
    input,
    by_vars = exprs(USUBJID),
    order = exprs(AVISITN, desc(AVAL)),
    new_var = lastfl,
    mode = "last"
  )

  expect_dfs_equal(
    base = expected_output,
    compare = actual_output,
    keys = c("USUBJID", "AVISITN", "AVAL")
  )
})

## Test 7: An error is issued if `derive_var_worst_flag()` is called ----
test_that("deprecation Test 7: An error is issued if Derive worst flag is called", {
  expect_error(
    derive_var_worst_flag(
      input_worst_flag,
      by_vars = exprs(USUBJID, PARAMCD, AVISIT),
      order = exprs(desc(ADT)),
      new_var = WORSTFL,
      param_var = PARAMCD,
      analysis_var = AVAL,
      worst_high = c("PARAM01", "PARAM03"),
      worst_low = "PARAM02"
    ),
    class = "lifecycle_error_deprecated"
  )
})
