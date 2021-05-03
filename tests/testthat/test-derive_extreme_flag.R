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

test_that("ABLFL = Y using last observation within a subset", {
  input <- tibble::tribble(
    ~STUDYID, ~USUBJID, ~PARAMCD,  ~AVISIT,    ~ADT,                 ~AVAL,
    "TEST01", "PAT01",  "PARAM01", "BASELINE", as.Date("2021-04-27"), 15.0,
    "TEST01", "PAT01",  "PARAM01", "BASELINE", as.Date("2021-04-25"), 14.0,
    "TEST01", "PAT01",  "PARAM01", "BASELINE", as.Date("2021-04-23"), 15.0,
    "TEST01", "PAT01",  "PARAM01", "WEEK 1",   as.Date("2021-04-27"), 10.0,
    "TEST01", "PAT01",  "PARAM01", "WEEK 2",   as.Date("2021-04-30"), 12.0,

    "TEST01", "PAT02",  "PARAM01", "SCREEN",   as.Date("2021-04-27"), 15.0,
    "TEST01", "PAT02",  "PARAM01", "BASELINE", as.Date("2021-04-25"), 14.0,
    "TEST01", "PAT02",  "PARAM01", "BASELINE", as.Date("2021-04-23"), 15.0,
    "TEST01", "PAT02",  "PARAM01", "WEEK 1",   as.Date("2021-04-27"), 10.0,
    "TEST01", "PAT02",  "PARAM01", "WEEK 2",   as.Date("2021-04-30"), 12.0,

    "TEST01", "PAT01",  "PARAM02", "SCREEN",   as.Date("2021-04-27"), 15.0,
    "TEST01", "PAT01",  "PARAM02", "SCREEN",   as.Date("2021-04-25"), 14.0,
    "TEST01", "PAT01",  "PARAM02", "SCREEN",   as.Date("2021-04-23"), 15.0,
    "TEST01", "PAT01",  "PARAM02", "BASELINE", as.Date("2021-04-27"), 10.0,
    "TEST01", "PAT01",  "PARAM02", "WEEK 2",   as.Date("2021-04-30"), 12.0,

    "TEST01", "PAT02",  "PARAM02", "SCREEN",   as.Date("2021-04-27"), 15.0,
    "TEST01", "PAT02",  "PARAM02", "BASELINE", as.Date("2021-04-25"), 14.0,
    "TEST01", "PAT02",  "PARAM02", "WEEK 1",   as.Date("2021-04-23"), 15.0,
    "TEST01", "PAT02",  "PARAM02", "WEEK 1",   as.Date("2021-04-27"), 10.0,
    "TEST01", "PAT02",  "PARAM02", "BASELINE", as.Date("2021-04-30"), 12.0

  )
  expected_output <- tibble::tribble(
    ~STUDYID, ~USUBJID, ~PARAMCD,  ~AVISIT,    ~ADT,                 ~AVAL, ~ABLFL,
    "TEST01", "PAT01",  "PARAM01", "BASELINE", as.Date("2021-04-27"), 15.0, "Y",
    "TEST01", "PAT01",  "PARAM01", "BASELINE", as.Date("2021-04-25"), 14.0, NA,
    "TEST01", "PAT01",  "PARAM01", "BASELINE", as.Date("2021-04-23"), 15.0, NA,
    "TEST01", "PAT01",  "PARAM01", "WEEK 1",   as.Date("2021-04-27"), 10.0, NA,
    "TEST01", "PAT01",  "PARAM01", "WEEK 2",   as.Date("2021-04-30"), 12.0, NA,

    "TEST01", "PAT02",  "PARAM01", "SCREEN",   as.Date("2021-04-27"), 15.0, NA,
    "TEST01", "PAT02",  "PARAM01", "BASELINE", as.Date("2021-04-25"), 14.0, "Y",
    "TEST01", "PAT02",  "PARAM01", "BASELINE", as.Date("2021-04-23"), 15.0, NA,
    "TEST01", "PAT02",  "PARAM01", "WEEK 1",   as.Date("2021-04-27"), 10.0, NA,
    "TEST01", "PAT02",  "PARAM01", "WEEK 2",   as.Date("2021-04-30"), 12.0, NA,

    "TEST01", "PAT01",  "PARAM02", "SCREEN",   as.Date("2021-04-27"), 15.0, NA,
    "TEST01", "PAT01",  "PARAM02", "SCREEN",   as.Date("2021-04-25"), 14.0, NA,
    "TEST01", "PAT01",  "PARAM02", "SCREEN",   as.Date("2021-04-23"), 15.0, NA,
    "TEST01", "PAT01",  "PARAM02", "BASELINE", as.Date("2021-04-27"), 10.0, "Y",
    "TEST01", "PAT01",  "PARAM02", "WEEK 2",   as.Date("2021-04-30"), 12.0, NA,

    "TEST01", "PAT02",  "PARAM02", "SCREEN",   as.Date("2021-04-27"), 15.0, NA,
    "TEST01", "PAT02",  "PARAM02", "BASELINE", as.Date("2021-04-25"), 14.0, NA,
    "TEST01", "PAT02",  "PARAM02", "WEEK 1",   as.Date("2021-04-23"), 15.0, NA,
    "TEST01", "PAT02",  "PARAM02", "WEEK 1",   as.Date("2021-04-27"), 10.0, NA,
    "TEST01", "PAT02",  "PARAM02", "BASELINE", as.Date("2021-04-30"), 12.0, "Y"
  )

  actual_output <- derive_extreme_flag(
    input,
    new_var = ABLFL,
    by_vars = exprs(USUBJID, PARAMCD),
    order = exprs(ADT),
    flag_filter = expr(AVISIT == "BASELINE")
  )

  expect_dfs_equal(
    expected_output,
    actual_output,
    keys = c("STUDYID", "USUBJID", "PARAMCD", "AVISIT", "ADT")
  )
})

test_that("ABLFL = Y worst observation = HI within a subset", {
  input <- tibble::tribble(
    ~STUDYID, ~USUBJID, ~PARAMCD,  ~AVISIT,    ~ADT,                 ~AVAL,
    "TEST01", "PAT01",  "PARAM01", "BASELINE", as.Date("2021-04-27"), 15.0,
    "TEST01", "PAT01",  "PARAM01", "BASELINE", as.Date("2021-04-25"), 14.0,
    "TEST01", "PAT01",  "PARAM01", "BASELINE", as.Date("2021-04-23"), 15.0,
    "TEST01", "PAT01",  "PARAM01", "WEEK 1",   as.Date("2021-04-27"), 10.0,
    "TEST01", "PAT01",  "PARAM01", "WEEK 2",   as.Date("2021-04-30"), 12.0,

    "TEST01", "PAT02",  "PARAM01", "SCREEN",   as.Date("2021-04-27"), 15.0,
    "TEST01", "PAT02",  "PARAM01", "BASELINE", as.Date("2021-04-25"), 14.0,
    "TEST01", "PAT02",  "PARAM01", "BASELINE", as.Date("2021-04-23"), 15.0,
    "TEST01", "PAT02",  "PARAM01", "WEEK 1",   as.Date("2021-04-27"), 10.0,
    "TEST01", "PAT02",  "PARAM01", "WEEK 2",   as.Date("2021-04-30"), 12.0,

    "TEST01", "PAT01",  "PARAM02", "SCREEN",   as.Date("2021-04-27"), 15.0,
    "TEST01", "PAT01",  "PARAM02", "SCREEN",   as.Date("2021-04-25"), 14.0,
    "TEST01", "PAT01",  "PARAM02", "SCREEN",   as.Date("2021-04-23"), 15.0,
    "TEST01", "PAT01",  "PARAM02", "BASELINE", as.Date("2021-04-27"), 10.0,
    "TEST01", "PAT01",  "PARAM02", "WEEK 2",   as.Date("2021-04-30"), 12.0,

    "TEST01", "PAT02",  "PARAM02", "SCREEN",   as.Date("2021-04-27"), 15.0,
    "TEST01", "PAT02",  "PARAM02", "BASELINE", as.Date("2021-04-25"), 14.0,
    "TEST01", "PAT02",  "PARAM02", "WEEK 1",   as.Date("2021-04-23"), 15.0,
    "TEST01", "PAT02",  "PARAM02", "WEEK 1",   as.Date("2021-04-27"), 10.0,
    "TEST01", "PAT02",  "PARAM02", "BASELINE", as.Date("2021-04-30"), 12.0

  )
  expected_output <- tibble::tribble(
    ~STUDYID, ~USUBJID, ~PARAMCD,  ~AVISIT,    ~ADT,                 ~AVAL, ~ABLFL,
    "TEST01", "PAT01",  "PARAM01", "BASELINE", as.Date("2021-04-27"), 15.0, "Y",
    "TEST01", "PAT01",  "PARAM01", "BASELINE", as.Date("2021-04-25"), 14.0, NA,
    "TEST01", "PAT01",  "PARAM01", "BASELINE", as.Date("2021-04-23"), 15.0, NA,
    "TEST01", "PAT01",  "PARAM01", "WEEK 1",   as.Date("2021-04-27"), 10.0, NA,
    "TEST01", "PAT01",  "PARAM01", "WEEK 2",   as.Date("2021-04-30"), 12.0, NA,

    "TEST01", "PAT02",  "PARAM01", "SCREEN",   as.Date("2021-04-27"), 15.0, NA,
    "TEST01", "PAT02",  "PARAM01", "BASELINE", as.Date("2021-04-25"), 14.0, NA,
    "TEST01", "PAT02",  "PARAM01", "BASELINE", as.Date("2021-04-23"), 15.0, "Y",
    "TEST01", "PAT02",  "PARAM01", "WEEK 1",   as.Date("2021-04-27"), 10.0, NA,
    "TEST01", "PAT02",  "PARAM01", "WEEK 2",   as.Date("2021-04-30"), 12.0, NA,

    "TEST01", "PAT01",  "PARAM02", "SCREEN",   as.Date("2021-04-27"), 15.0, NA,
    "TEST01", "PAT01",  "PARAM02", "SCREEN",   as.Date("2021-04-25"), 14.0, NA,
    "TEST01", "PAT01",  "PARAM02", "SCREEN",   as.Date("2021-04-23"), 15.0, NA,
    "TEST01", "PAT01",  "PARAM02", "BASELINE", as.Date("2021-04-27"), 10.0, "Y",
    "TEST01", "PAT01",  "PARAM02", "WEEK 2",   as.Date("2021-04-30"), 12.0, NA,

    "TEST01", "PAT02",  "PARAM02", "SCREEN",   as.Date("2021-04-27"), 15.0, NA,
    "TEST01", "PAT02",  "PARAM02", "BASELINE", as.Date("2021-04-25"), 14.0, "Y",
    "TEST01", "PAT02",  "PARAM02", "WEEK 1",   as.Date("2021-04-23"), 15.0, NA,
    "TEST01", "PAT02",  "PARAM02", "WEEK 1",   as.Date("2021-04-27"), 10.0, NA,
    "TEST01", "PAT02",  "PARAM02", "BASELINE", as.Date("2021-04-30"), 12.0, NA
  )

  actual_output <- derive_extreme_flag(
    input,
    new_var = ABLFL,
    by_vars = exprs(USUBJID, PARAMCD),
    order = exprs(AVAL, ADT),
    flag_filter = expr(AVISIT == "BASELINE")
  )

  expect_dfs_equal(
    expected_output,
    actual_output,
    keys = c("STUDYID", "USUBJID", "PARAMCD", "AVISIT", "ADT")
  )
})

test_that("ABLFL = Y worst observation = LO within a subset", {
  input <- tibble::tribble(
    ~STUDYID, ~USUBJID, ~PARAMCD,  ~AVISIT,    ~ADT,                 ~AVAL,
    "TEST01", "PAT01",  "PARAM01", "BASELINE", as.Date("2021-04-27"), 15.0,
    "TEST01", "PAT01",  "PARAM01", "BASELINE", as.Date("2021-04-25"), 14.0,
    "TEST01", "PAT01",  "PARAM01", "BASELINE", as.Date("2021-04-23"), 15.0,
    "TEST01", "PAT01",  "PARAM01", "WEEK 1",   as.Date("2021-04-27"), 10.0,
    "TEST01", "PAT01",  "PARAM01", "WEEK 2",   as.Date("2021-04-30"), 12.0,

    "TEST01", "PAT02",  "PARAM01", "SCREEN",   as.Date("2021-04-27"), 15.0,
    "TEST01", "PAT02",  "PARAM01", "BASELINE", as.Date("2021-04-25"), 14.0,
    "TEST01", "PAT02",  "PARAM01", "BASELINE", as.Date("2021-04-23"), 15.0,
    "TEST01", "PAT02",  "PARAM01", "WEEK 1",   as.Date("2021-04-27"), 10.0,
    "TEST01", "PAT02",  "PARAM01", "WEEK 2",   as.Date("2021-04-30"), 12.0,

    "TEST01", "PAT01",  "PARAM02", "SCREEN",   as.Date("2021-04-27"), 15.0,
    "TEST01", "PAT01",  "PARAM02", "SCREEN",   as.Date("2021-04-25"), 14.0,
    "TEST01", "PAT01",  "PARAM02", "SCREEN",   as.Date("2021-04-23"), 15.0,
    "TEST01", "PAT01",  "PARAM02", "BASELINE", as.Date("2021-04-27"), 10.0,
    "TEST01", "PAT01",  "PARAM02", "WEEK 2",   as.Date("2021-04-30"), 12.0,

    "TEST01", "PAT02",  "PARAM02", "SCREEN",   as.Date("2021-04-27"), 15.0,
    "TEST01", "PAT02",  "PARAM02", "BASELINE", as.Date("2021-04-25"), 14.0,
    "TEST01", "PAT02",  "PARAM02", "WEEK 1",   as.Date("2021-04-23"), 15.0,
    "TEST01", "PAT02",  "PARAM02", "WEEK 1",   as.Date("2021-04-27"), 10.0,
    "TEST01", "PAT02",  "PARAM02", "BASELINE", as.Date("2021-04-30"), 12.0

  )
  expected_output <- tibble::tribble(
    ~STUDYID, ~USUBJID, ~PARAMCD,  ~AVISIT,    ~ADT,                 ~AVAL, ~ABLFL,
    "TEST01", "PAT01",  "PARAM01", "BASELINE", as.Date("2021-04-27"), 15.0, NA,
    "TEST01", "PAT01",  "PARAM01", "BASELINE", as.Date("2021-04-25"), 14.0, "Y",
    "TEST01", "PAT01",  "PARAM01", "BASELINE", as.Date("2021-04-23"), 15.0, NA,
    "TEST01", "PAT01",  "PARAM01", "WEEK 1",   as.Date("2021-04-27"), 10.0, NA,
    "TEST01", "PAT01",  "PARAM01", "WEEK 2",   as.Date("2021-04-30"), 12.0, NA,

    "TEST01", "PAT02",  "PARAM01", "SCREEN",   as.Date("2021-04-27"), 15.0, NA,
    "TEST01", "PAT02",  "PARAM01", "BASELINE", as.Date("2021-04-25"), 14.0, "Y",
    "TEST01", "PAT02",  "PARAM01", "BASELINE", as.Date("2021-04-23"), 15.0, NA,
    "TEST01", "PAT02",  "PARAM01", "WEEK 1",   as.Date("2021-04-27"), 10.0, NA,
    "TEST01", "PAT02",  "PARAM01", "WEEK 2",   as.Date("2021-04-30"), 12.0, NA,

    "TEST01", "PAT01",  "PARAM02", "SCREEN",   as.Date("2021-04-27"), 15.0, NA,
    "TEST01", "PAT01",  "PARAM02", "SCREEN",   as.Date("2021-04-25"), 14.0, NA,
    "TEST01", "PAT01",  "PARAM02", "SCREEN",   as.Date("2021-04-23"), 15.0, NA,
    "TEST01", "PAT01",  "PARAM02", "BASELINE", as.Date("2021-04-27"), 10.0, "Y",
    "TEST01", "PAT01",  "PARAM02", "WEEK 2",   as.Date("2021-04-30"), 12.0, NA,

    "TEST01", "PAT02",  "PARAM02", "SCREEN",   as.Date("2021-04-27"), 15.0, NA,
    "TEST01", "PAT02",  "PARAM02", "BASELINE", as.Date("2021-04-25"), 14.0, NA,
    "TEST01", "PAT02",  "PARAM02", "WEEK 1",   as.Date("2021-04-23"), 15.0, NA,
    "TEST01", "PAT02",  "PARAM02", "WEEK 1",   as.Date("2021-04-27"), 10.0, NA,
    "TEST01", "PAT02",  "PARAM02", "BASELINE", as.Date("2021-04-30"), 12.0, "Y"
  )

  actual_output <- derive_extreme_flag(
    input,
    new_var = ABLFL,
    by_vars = exprs(USUBJID, PARAMCD),
    order = exprs(desc(AVAL), ADT),
    flag_filter = expr(AVISIT == "BASELINE")
  )

  expect_dfs_equal(
    expected_output,
    actual_output,
    keys = c("STUDYID", "USUBJID", "PARAMCD", "AVISIT", "ADT")
  )
})

test_that("ABLFL = Y average records within a subset", {
  input <- tibble::tribble(
    ~STUDYID, ~USUBJID, ~PARAMCD,  ~AVISIT,    ~ADT,                 ~AVAL, ~DTYPE,
    "TEST01", "PAT01",  "PARAM01", "BASELINE", as.Date("2021-04-27"), 15.0, NA,
    "TEST01", "PAT01",  "PARAM01", "BASELINE", as.Date("2021-04-25"), 14.0, NA,
    "TEST01", "PAT01",  "PARAM01", "BASELINE", as.Date("2021-04-23"), 15.0, "AVERAGE",
    "TEST01", "PAT01",  "PARAM01", "WEEK 1",   as.Date("2021-04-27"), 10.0, "AVERAGE",
    "TEST01", "PAT01",  "PARAM01", "WEEK 2",   as.Date("2021-04-30"), 12.0, NA,

    "TEST01", "PAT02",  "PARAM01", "SCREEN",   as.Date("2021-04-27"), 15.0, "AVERAGE",
    "TEST01", "PAT02",  "PARAM01", "BASELINE", as.Date("2021-04-25"), 14.0, "AVERAGE",
    "TEST01", "PAT02",  "PARAM01", "BASELINE", as.Date("2021-04-23"), 15.0, "AVERAGE",
    "TEST01", "PAT02",  "PARAM01", "WEEK 1",   as.Date("2021-04-27"), 10.0, "AVERAGE",
    "TEST01", "PAT02",  "PARAM01", "WEEK 2",   as.Date("2021-04-30"), 12.0, "AVERAGE",

    "TEST01", "PAT01",  "PARAM02", "SCREEN",   as.Date("2021-04-27"), 15.0, "AVERAGE",
    "TEST01", "PAT01",  "PARAM02", "SCREEN",   as.Date("2021-04-25"), 14.0, "AVERAGE",
    "TEST01", "PAT01",  "PARAM02", "SCREEN",   as.Date("2021-04-23"), 15.0, NA,
    "TEST01", "PAT01",  "PARAM02", "BASELINE", as.Date("2021-04-27"), 10.0, "AVERAGE",
    "TEST01", "PAT01",  "PARAM02", "WEEK 2",   as.Date("2021-04-30"), 12.0, NA,

    "TEST01", "PAT02",  "PARAM02", "SCREEN",   as.Date("2021-04-27"), 15.0, NA,
    "TEST01", "PAT02",  "PARAM02", "BASELINE", as.Date("2021-04-25"), 14.0, NA,
    "TEST01", "PAT02",  "PARAM02", "WEEK 1",   as.Date("2021-04-23"), 15.0, NA,
    "TEST01", "PAT02",  "PARAM02", "WEEK 1",   as.Date("2021-04-27"), 10.0, NA,
    "TEST01", "PAT02",  "PARAM02", "BASELINE", as.Date("2021-04-30"), 12.0, NA

  )
  expected_output <- tibble::tribble(
    ~STUDYID, ~USUBJID, ~PARAMCD,  ~AVISIT,    ~ADT,                 ~AVAL, ~DTYPE,    ~ABLFL,
    "TEST01", "PAT01",  "PARAM01", "BASELINE", as.Date("2021-04-27"), 15.0, NA,        NA,
    "TEST01", "PAT01",  "PARAM01", "BASELINE", as.Date("2021-04-25"), 14.0, NA,        NA,
    "TEST01", "PAT01",  "PARAM01", "BASELINE", as.Date("2021-04-23"), 15.0, "AVERAGE", "Y",
    "TEST01", "PAT01",  "PARAM01", "WEEK 1",   as.Date("2021-04-27"), 10.0, "AVERAGE", NA,
    "TEST01", "PAT01",  "PARAM01", "WEEK 2",   as.Date("2021-04-30"), 12.0, NA,        NA,

    "TEST01", "PAT02",  "PARAM01", "SCREEN",   as.Date("2021-04-27"), 15.0, "AVERAGE", NA,
    "TEST01", "PAT02",  "PARAM01", "BASELINE", as.Date("2021-04-25"), 14.0, "AVERAGE", "Y",
    "TEST01", "PAT02",  "PARAM01", "BASELINE", as.Date("2021-04-23"), 15.0, "AVERAGE", NA,
    "TEST01", "PAT02",  "PARAM01", "WEEK 1",   as.Date("2021-04-27"), 10.0, "AVERAGE", NA,
    "TEST01", "PAT02",  "PARAM01", "WEEK 2",   as.Date("2021-04-30"), 12.0, "AVERAGE", NA,

    "TEST01", "PAT01",  "PARAM02", "SCREEN",   as.Date("2021-04-27"), 15.0, "AVERAGE", NA,
    "TEST01", "PAT01",  "PARAM02", "SCREEN",   as.Date("2021-04-25"), 14.0, "AVERAGE", NA,
    "TEST01", "PAT01",  "PARAM02", "SCREEN",   as.Date("2021-04-23"), 15.0, NA,        NA,
    "TEST01", "PAT01",  "PARAM02", "BASELINE", as.Date("2021-04-27"), 10.0, "AVERAGE", "Y",
    "TEST01", "PAT01",  "PARAM02", "WEEK 2",   as.Date("2021-04-30"), 12.0, NA,        NA,

    "TEST01", "PAT02",  "PARAM02", "SCREEN",   as.Date("2021-04-27"), 15.0, NA,        NA,
    "TEST01", "PAT02",  "PARAM02", "BASELINE", as.Date("2021-04-25"), 14.0, NA,        NA,
    "TEST01", "PAT02",  "PARAM02", "WEEK 1",   as.Date("2021-04-23"), 15.0, NA,        NA,
    "TEST01", "PAT02",  "PARAM02", "WEEK 1",   as.Date("2021-04-27"), 10.0, NA,        NA,
    "TEST01", "PAT02",  "PARAM02", "BASELINE", as.Date("2021-04-30"), 12.0, NA,        NA
  )

  actual_output <- derive_extreme_flag(
    input,
    new_var = ABLFL,
    by_vars = exprs(USUBJID, PARAMCD),
    order = exprs(ADT, desc(AVAL)),
    flag_filter = expr(AVISIT == "BASELINE" & DTYPE == "AVERAGE")
  )

  expect_dfs_equal(
    expected_output,
    actual_output,
    keys = c("STUDYID", "USUBJID", "PARAMCD", "AVISIT", "ADT", "AVAL")
  )
})

test_that("ABLFL = Y using last observation within a subset and multiple baselines
          possible", {
  input <- tibble::tribble(
    ~STUDYID, ~USUBJID, ~PARAMCD,  ~AVISIT,    ~ADT,                 ~AVAL,
    "TEST01", "PAT01",  "PARAM01", "BASELINE", as.Date("2021-04-27"), 15.0,
    "TEST01", "PAT01",  "PARAM01", "BASELINE", as.Date("2021-04-25"), 14.0,
    "TEST01", "PAT01",  "PARAM01", "BASELINE", as.Date("2021-04-23"), 15.0,
    "TEST01", "PAT01",  "PARAM01", "WEEK 1",   as.Date("2021-04-27"), 10.0,
    "TEST01", "PAT01",  "PARAM01", "WEEK 2",   as.Date("2021-04-30"), 12.0,

    "TEST01", "PAT02",  "PARAM01", "SCREEN",   as.Date("2021-04-27"), 15.0,
    "TEST01", "PAT02",  "PARAM01", "BASELINE", as.Date("2021-04-25"), 14.0,
    "TEST01", "PAT02",  "PARAM01", "BASELINE", as.Date("2021-04-23"), 15.0,
    "TEST01", "PAT02",  "PARAM01", "WEEK 1",   as.Date("2021-04-27"), 10.0,
    "TEST01", "PAT02",  "PARAM01", "WEEK 2",   as.Date("2021-04-30"), 12.0,

    "TEST01", "PAT01",  "PARAM02", "SCREEN",   as.Date("2021-04-27"), 15.0,
    "TEST01", "PAT01",  "PARAM02", "SCREEN",   as.Date("2021-04-25"), 14.0,
    "TEST01", "PAT01",  "PARAM02", "SCREEN",   as.Date("2021-04-23"), 15.0,
    "TEST01", "PAT01",  "PARAM02", "BASELINE", as.Date("2021-04-27"), 10.0,
    "TEST01", "PAT01",  "PARAM02", "WEEK 2",   as.Date("2021-04-30"), 12.0,

    "TEST01", "PAT02",  "PARAM02", "SCREEN",   as.Date("2021-04-27"), 15.0,
    "TEST01", "PAT02",  "PARAM02", "BASELINE", as.Date("2021-04-25"), 14.0,
    "TEST01", "PAT02",  "PARAM02", "WEEK 1",   as.Date("2021-04-23"), 15.0,
    "TEST01", "PAT02",  "PARAM02", "WEEK 1",   as.Date("2021-04-27"), 10.0,
    "TEST01", "PAT02",  "PARAM02", "BASELINE", as.Date("2021-04-30"), 12.0
  )

  expected_output <- tibble::tribble(
    ~STUDYID, ~USUBJID, ~PARAMCD,  ~AVISIT,    ~ADT,                 ~AVAL, ~ABLFL,
    "TEST01", "PAT01",  "PARAM01", "BASELINE", as.Date("2021-04-27"), 15.0, "Y",
    "TEST01", "PAT01",  "PARAM01", "BASELINE", as.Date("2021-04-25"), 14.0, NA,
    "TEST01", "PAT01",  "PARAM01", "BASELINE", as.Date("2021-04-23"), 15.0, NA,
    "TEST01", "PAT01",  "PARAM01", "WEEK 1",   as.Date("2021-04-27"), 10.0, "Y",
    "TEST01", "PAT01",  "PARAM01", "WEEK 2",   as.Date("2021-04-30"), 12.0, NA,

    "TEST01", "PAT02",  "PARAM01", "SCREEN",   as.Date("2021-04-27"), 15.0, NA,
    "TEST01", "PAT02",  "PARAM01", "BASELINE", as.Date("2021-04-25"), 14.0, "Y",
    "TEST01", "PAT02",  "PARAM01", "BASELINE", as.Date("2021-04-23"), 15.0, NA,
    "TEST01", "PAT02",  "PARAM01", "WEEK 1",   as.Date("2021-04-27"), 10.0, "Y",
    "TEST01", "PAT02",  "PARAM01", "WEEK 2",   as.Date("2021-04-30"), 12.0, NA,

    "TEST01", "PAT01",  "PARAM02", "SCREEN",   as.Date("2021-04-27"), 15.0, NA,
    "TEST01", "PAT01",  "PARAM02", "SCREEN",   as.Date("2021-04-25"), 14.0, NA,
    "TEST01", "PAT01",  "PARAM02", "SCREEN",   as.Date("2021-04-23"), 15.0, NA,
    "TEST01", "PAT01",  "PARAM02", "BASELINE", as.Date("2021-04-27"), 10.0, "Y",
    "TEST01", "PAT01",  "PARAM02", "WEEK 2",   as.Date("2021-04-30"), 12.0, NA,

    "TEST01", "PAT02",  "PARAM02", "SCREEN",   as.Date("2021-04-27"), 15.0, NA,
    "TEST01", "PAT02",  "PARAM02", "BASELINE", as.Date("2021-04-25"), 14.0, NA,
    "TEST01", "PAT02",  "PARAM02", "WEEK 1",   as.Date("2021-04-23"), 15.0, NA,
    "TEST01", "PAT02",  "PARAM02", "WEEK 1",   as.Date("2021-04-27"), 10.0, "Y",
    "TEST01", "PAT02",  "PARAM02", "BASELINE", as.Date("2021-04-30"), 12.0, "Y"
  )

  actual_output <- derive_extreme_flag(
    input,
    new_var = ABLFL,
    by_vars = exprs(USUBJID, PARAMCD, AVISIT),
    order = exprs(ADT),
    flag_filter = expr(AVISIT %in% c("BASELINE","WEEK 1"))
  )

  expect_dfs_equal(
    expected_output,
    actual_output,
    keys = c("STUDYID", "USUBJID", "PARAMCD", "AVISIT", "ADT")
  )
})
