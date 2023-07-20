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

## Test 1: first observation for each group is flagged ----
test_that("derive_var_extreme_flag Test 1: first observation for each group is flagged", {
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

## Test 2: last observation for each group is flagged ----
test_that("derive_var_extreme_flag Test 2: last observation for each group is flagged", {
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


test_flag_all <- tibble::tribble(
  ~STUDYID, ~USUBJID,               ~ADTM, ~AVISITN, ~BASETYPE,
  "TEST01",  "PAT01",  "2020-02-01T12:00",        1,   "ontrt",
  "TEST01",  "PAT01",  "2020-02-01T12:00",        1,   "ontrt",
  "TEST01",  "PAT01",  "2020-02-01T12:01",        1,   "ontrt",
  "TEST01",  "PAT01",  "2020-02-01T13:00",        1,   "ontrt",
  "TEST01",  "PAT01",  "2020-02-01T13:00",        1,   "ontrt"
)

## Test 3: flag_all = FALSE when mode is first ----
test_that("derive_var_extreme_flag Test 3: flag_all = FALSE when mode is first", {
  expected_output <- test_flag_all %>%
    mutate(FIRSTFL = c("Y", NA, NA, NA, NA))

  actual_output <- derive_var_extreme_flag(
    test_flag_all,
    by_vars = exprs(STUDYID, USUBJID, BASETYPE, AVISITN),
    order = exprs(ADTM),
    new_var = FIRSTFL,
    mode = "first",
    check_type = "none",
    flag_all = FALSE
  )

  expect_true(all.equal(expected_output, actual_output))
})

## Test 4: flag_all = TRUE when mode is first ----
test_that("derive_var_extreme_flag Test 4: flag_all = TRUE when mode is first", {
  expected_output <- test_flag_all %>%
    mutate(FIRSTFL = c("Y", "Y", NA, NA, NA))

  actual_output <- derive_var_extreme_flag(
    test_flag_all,
    by_vars = exprs(STUDYID, USUBJID, BASETYPE, AVISITN),
    order = exprs(ADTM),
    new_var = FIRSTFL,
    mode = "first",
    flag_all = TRUE
  )

  expect_true(all.equal(expected_output, actual_output))
})

## Test 5: flag_all = FALSE when mode is last ----
test_that("derive_var_extreme_flag Test 5: flag_all = FALSE when mode is last", {
  expected_output <- test_flag_all %>%
    mutate(LASTFL = c(NA, NA, NA, NA, "Y"))

  actual_output <- derive_var_extreme_flag(
    test_flag_all,
    by_vars = exprs(STUDYID, USUBJID, BASETYPE, AVISITN),
    order = exprs(ADTM),
    new_var = LASTFL,
    mode = "last",
    check_type = "none",
    flag_all = FALSE
  )

  expect_true(all.equal(expected_output, actual_output))
})

## Test 6: flag_all = TRUE when mode is last ----
test_that("derive_var_extreme_flag Test 6: flag_all = TRUE when mode is last", {
  expected_output <- test_flag_all %>%
    mutate(LASTFL = c(NA, NA, NA, "Y", "Y"))

  actual_output <- derive_var_extreme_flag(
    test_flag_all,
    by_vars = exprs(STUDYID, USUBJID, BASETYPE, AVISITN),
    order = exprs(ADTM),
    new_var = LASTFL,
    mode = "last",
    flag_all = TRUE
  )

  expect_true(all.equal(expected_output, actual_output))
})
