context("test-derive_var_ablfl")

test_that("ABLFL = Y using last observation within a subset", {
  input <- tibble::tribble(
    ~STUDYID, ~USUBJID, ~PARAMCD,  ~AVISIT,    ~ADT,                 ~AVAL,
    "TEST01", "PAT01",  "PARAM01", "BASELINE", as.Date("2021-04-27"), 15.0,
    "TEST01", "PAT01",  "PARAM01", "BASELINE", as.Date("2021-04-25"), 14.0,
    "TEST01", "PAT01",  "PARAM01", "BASELINE", as.Date("2021-04-23"), 15.0,
    "TEST01", "PAT01",  "PARAM01", "WEEK 1",   as.Date("2021-04-27"), 10.0,
    "TEST01", "PAT01",  "PARAM01", "WEEK 2",   as.Date("2021-04-30"), 12.0,

    "TEST01", "PAT02",  "PARAM01", "SCREENING", as.Date("2021-04-27"), 15.0,
    "TEST01", "PAT02",  "PARAM01", "BASELINE", as.Date("2021-04-25"), 14.0,
    "TEST01", "PAT02",  "PARAM01", "BASELINE", as.Date("2021-04-23"), 15.0,
    "TEST01", "PAT02",  "PARAM01", "WEEK 1",   as.Date("2021-04-27"), 10.0,
    "TEST01", "PAT02",  "PARAM01", "WEEK 2",   as.Date("2021-04-30"), 12.0,

    "TEST01", "PAT01",  "PARAM02", "SCREENING", as.Date("2021-04-27"), 15.0,
    "TEST01", "PAT01",  "PARAM02", "SCREENING", as.Date("2021-04-25"), 14.0,
    "TEST01", "PAT01",  "PARAM02", "SCREENING", as.Date("2021-04-23"), 15.0,
    "TEST01", "PAT01",  "PARAM02", "BASELINE",  as.Date("2021-04-27"), 10.0,
    "TEST01", "PAT01",  "PARAM02", "WEEK 2",    as.Date("2021-04-30"), 12.0,

    "TEST01", "PAT02",  "PARAM02", "SCREENING",as.Date("2021-04-27"), 15.0,
    "TEST01", "PAT02",  "PARAM02", "BASELINE", as.Date("2021-04-25"), 14.0,
    "TEST01", "PAT02",  "PARAM02", "WEEK 1",   as.Date("2021-04-23"), 15.0,
    "TEST01", "PAT02",  "PARAM02", "WEEK 1",   as.Date("2021-04-27"), 10.0,
    "TEST01", "PAT02",  "PARAM02", "BASELINE", as.Date("2021-04-30"), 12.0

  )
  expected_output <- tibble::tribble(
    ~STUDYID, ~USUBJID, ~PARAMCD,  ~AVISIT,    ~ADT,                 ~AVAL, ~ABLFL,
    "TEST01", "PAT01",  "PARAM01", "BASELINE", as.Date("2021-04-27"), 15.0, "Y",
    "TEST01", "PAT01",  "PARAM01", "BASELINE", as.Date("2021-04-25"), 14.0, "",
    "TEST01", "PAT01",  "PARAM01", "BASELINE", as.Date("2021-04-23"), 15.0, "",
    "TEST01", "PAT01",  "PARAM01", "WEEK 1",   as.Date("2021-04-27"), 10.0, "",
    "TEST01", "PAT01",  "PARAM01", "WEEK 2",   as.Date("2021-04-30"), 12.0, "",

    "TEST01", "PAT02",  "PARAM01", "SCREENING", as.Date("2021-04-27"), 15.0, "",
    "TEST01", "PAT02",  "PARAM01", "BASELINE",  as.Date("2021-04-25"), 14.0, "Y",
    "TEST01", "PAT02",  "PARAM01", "BASELINE", as.Date("2021-04-23"), 15.0, "",
    "TEST01", "PAT02",  "PARAM01", "WEEK 1",    as.Date("2021-04-27"), 10.0, "",
    "TEST01", "PAT02",  "PARAM01", "WEEK 2",    as.Date("2021-04-30"), 12.0, "",

    "TEST01", "PAT01",  "PARAM02", "SCREENING", as.Date("2021-04-27"), 15.0, "",
    "TEST01", "PAT01",  "PARAM02", "SCREENING", as.Date("2021-04-25"), 14.0, "",
    "TEST01", "PAT01",  "PARAM02", "SCREENING", as.Date("2021-04-23"), 15.0, "",
    "TEST01", "PAT01",  "PARAM02", "BASELINE",  as.Date("2021-04-27"), 10.0, "Y",
    "TEST01", "PAT01",  "PARAM02", "WEEK 2",    as.Date("2021-04-30"), 12.0, "",

    "TEST01", "PAT02",  "PARAM02", "SCREENING", as.Date("2021-04-27"), 15.0, "",
    "TEST01", "PAT02",  "PARAM02", "BASELINE",  as.Date("2021-04-25"), 14.0, "",
    "TEST01", "PAT02",  "PARAM02", "WEEK 1",    as.Date("2021-04-23"), 15.0, "",
    "TEST01", "PAT02",  "PARAM02", "WEEK 1",    as.Date("2021-04-27"), 10.0, "",
    "TEST01", "PAT02",  "PARAM02", "BASELINE",  as.Date("2021-04-30"), 12.0, "Y"
    )

    actual_output <- derive_var_ablfl(
    input,
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

    "TEST01", "PAT02",  "PARAM01", "SCREENING", as.Date("2021-04-27"), 15.0,
    "TEST01", "PAT02",  "PARAM01", "BASELINE", as.Date("2021-04-25"), 14.0,
    "TEST01", "PAT02",  "PARAM01", "BASELINE", as.Date("2021-04-23"), 15.0,
    "TEST01", "PAT02",  "PARAM01", "WEEK 1",   as.Date("2021-04-27"), 10.0,
    "TEST01", "PAT02",  "PARAM01", "WEEK 2",   as.Date("2021-04-30"), 12.0,

    "TEST01", "PAT01",  "PARAM02", "SCREENING", as.Date("2021-04-27"), 15.0,
    "TEST01", "PAT01",  "PARAM02", "SCREENING", as.Date("2021-04-25"), 14.0,
    "TEST01", "PAT01",  "PARAM02", "SCREENING", as.Date("2021-04-23"), 15.0,
    "TEST01", "PAT01",  "PARAM02", "BASELINE",  as.Date("2021-04-27"), 10.0,
    "TEST01", "PAT01",  "PARAM02", "WEEK 2",    as.Date("2021-04-30"), 12.0,

    "TEST01", "PAT02",  "PARAM02", "SCREENING",as.Date("2021-04-27"), 15.0,
    "TEST01", "PAT02",  "PARAM02", "BASELINE", as.Date("2021-04-25"), 14.0,
    "TEST01", "PAT02",  "PARAM02", "WEEK 1",   as.Date("2021-04-23"), 15.0,
    "TEST01", "PAT02",  "PARAM02", "WEEK 1",   as.Date("2021-04-27"), 10.0,
    "TEST01", "PAT02",  "PARAM02", "BASELINE", as.Date("2021-04-30"), 12.0

  )
  expected_output <- tibble::tribble(
    ~STUDYID, ~USUBJID, ~PARAMCD,  ~AVISIT,    ~ADT,                 ~AVAL, ~ABLFL,
    "TEST01", "PAT01",  "PARAM01", "BASELINE", as.Date("2021-04-27"), 15.0, "Y",
    "TEST01", "PAT01",  "PARAM01", "BASELINE", as.Date("2021-04-25"), 14.0, "",
    "TEST01", "PAT01",  "PARAM01", "BASELINE", as.Date("2021-04-23"), 15.0, "",
    "TEST01", "PAT01",  "PARAM01", "WEEK 1",   as.Date("2021-04-27"), 10.0, "",
    "TEST01", "PAT01",  "PARAM01", "WEEK 2",   as.Date("2021-04-30"), 12.0, "",

    "TEST01", "PAT02",  "PARAM01", "SCREENING", as.Date("2021-04-27"), 15.0, "",
    "TEST01", "PAT02",  "PARAM01", "BASELINE",  as.Date("2021-04-25"), 14.0, "",
    "TEST01", "PAT02",  "PARAM01", "BASELINE", as.Date("2021-04-23"), 15.0, "Y",
    "TEST01", "PAT02",  "PARAM01", "WEEK 1",    as.Date("2021-04-27"), 10.0, "",
    "TEST01", "PAT02",  "PARAM01", "WEEK 2",    as.Date("2021-04-30"), 12.0, "",

    "TEST01", "PAT01",  "PARAM02", "SCREENING", as.Date("2021-04-27"), 15.0, "",
    "TEST01", "PAT01",  "PARAM02", "SCREENING", as.Date("2021-04-25"), 14.0, "",
    "TEST01", "PAT01",  "PARAM02", "SCREENING", as.Date("2021-04-23"), 15.0, "",
    "TEST01", "PAT01",  "PARAM02", "BASELINE",  as.Date("2021-04-27"), 10.0, "Y",
    "TEST01", "PAT01",  "PARAM02", "WEEK 2",    as.Date("2021-04-30"), 12.0, "",

    "TEST01", "PAT02",  "PARAM02", "SCREENING", as.Date("2021-04-27"), 15.0, "",
    "TEST01", "PAT02",  "PARAM02", "BASELINE",  as.Date("2021-04-25"), 14.0, "Y",
    "TEST01", "PAT02",  "PARAM02", "WEEK 1",    as.Date("2021-04-23"), 15.0, "",
    "TEST01", "PAT02",  "PARAM02", "WEEK 1",    as.Date("2021-04-27"), 10.0, "",
    "TEST01", "PAT02",  "PARAM02", "BASELINE",  as.Date("2021-04-30"), 12.0, ""
  )

  actual_output <- derive_var_ablfl(
    input,
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

    "TEST01", "PAT02",  "PARAM01", "SCREENING", as.Date("2021-04-27"), 15.0,
    "TEST01", "PAT02",  "PARAM01", "BASELINE", as.Date("2021-04-25"), 14.0,
    "TEST01", "PAT02",  "PARAM01", "BASELINE", as.Date("2021-04-23"), 15.0,
    "TEST01", "PAT02",  "PARAM01", "WEEK 1",   as.Date("2021-04-27"), 10.0,
    "TEST01", "PAT02",  "PARAM01", "WEEK 2",   as.Date("2021-04-30"), 12.0,

    "TEST01", "PAT01",  "PARAM02", "SCREENING", as.Date("2021-04-27"), 15.0,
    "TEST01", "PAT01",  "PARAM02", "SCREENING", as.Date("2021-04-25"), 14.0,
    "TEST01", "PAT01",  "PARAM02", "SCREENING", as.Date("2021-04-23"), 15.0,
    "TEST01", "PAT01",  "PARAM02", "BASELINE",  as.Date("2021-04-27"), 10.0,
    "TEST01", "PAT01",  "PARAM02", "WEEK 2",    as.Date("2021-04-30"), 12.0,

    "TEST01", "PAT02",  "PARAM02", "SCREENING",as.Date("2021-04-27"), 15.0,
    "TEST01", "PAT02",  "PARAM02", "BASELINE", as.Date("2021-04-25"), 14.0,
    "TEST01", "PAT02",  "PARAM02", "WEEK 1",   as.Date("2021-04-23"), 15.0,
    "TEST01", "PAT02",  "PARAM02", "WEEK 1",   as.Date("2021-04-27"), 10.0,
    "TEST01", "PAT02",  "PARAM02", "BASELINE", as.Date("2021-04-30"), 12.0

  )
  expected_output <- tibble::tribble(
    ~STUDYID, ~USUBJID, ~PARAMCD,  ~AVISIT,    ~ADT,                 ~AVAL, ~ABLFL,
    "TEST01", "PAT01",  "PARAM01", "BASELINE", as.Date("2021-04-27"), 15.0, "",
    "TEST01", "PAT01",  "PARAM01", "BASELINE", as.Date("2021-04-25"), 14.0, "Y",
    "TEST01", "PAT01",  "PARAM01", "BASELINE", as.Date("2021-04-23"), 15.0, "",
    "TEST01", "PAT01",  "PARAM01", "WEEK 1",   as.Date("2021-04-27"), 10.0, "",
    "TEST01", "PAT01",  "PARAM01", "WEEK 2",   as.Date("2021-04-30"), 12.0, "",

    "TEST01", "PAT02",  "PARAM01", "SCREENING", as.Date("2021-04-27"), 15.0, "",
    "TEST01", "PAT02",  "PARAM01", "BASELINE",  as.Date("2021-04-25"), 14.0, "Y",
    "TEST01", "PAT02",  "PARAM01", "BASELINE", as.Date("2021-04-23"), 15.0, "",
    "TEST01", "PAT02",  "PARAM01", "WEEK 1",    as.Date("2021-04-27"), 10.0, "",
    "TEST01", "PAT02",  "PARAM01", "WEEK 2",    as.Date("2021-04-30"), 12.0, "",

    "TEST01", "PAT01",  "PARAM02", "SCREENING", as.Date("2021-04-27"), 15.0, "",
    "TEST01", "PAT01",  "PARAM02", "SCREENING", as.Date("2021-04-25"), 14.0, "",
    "TEST01", "PAT01",  "PARAM02", "SCREENING", as.Date("2021-04-23"), 15.0, "",
    "TEST01", "PAT01",  "PARAM02", "BASELINE",  as.Date("2021-04-27"), 10.0, "Y",
    "TEST01", "PAT01",  "PARAM02", "WEEK 2",    as.Date("2021-04-30"), 12.0, "",

    "TEST01", "PAT02",  "PARAM02", "SCREENING", as.Date("2021-04-27"), 15.0, "",
    "TEST01", "PAT02",  "PARAM02", "BASELINE",  as.Date("2021-04-25"), 14.0, "",
    "TEST01", "PAT02",  "PARAM02", "WEEK 1",    as.Date("2021-04-23"), 15.0, "",
    "TEST01", "PAT02",  "PARAM02", "WEEK 1",    as.Date("2021-04-27"), 10.0, "",
    "TEST01", "PAT02",  "PARAM02", "BASELINE",  as.Date("2021-04-30"), 12.0, "Y"
  )

  actual_output <- derive_var_ablfl(
    input,
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

    "TEST01", "PAT02",  "PARAM01", "SCREENING", as.Date("2021-04-27"), 15.0, "AVERAGE",
    "TEST01", "PAT02",  "PARAM01", "BASELINE", as.Date("2021-04-25"), 14.0, "AVERAGE",
    "TEST01", "PAT02",  "PARAM01", "BASELINE", as.Date("2021-04-23"), 15.0, "AVERAGE",
    "TEST01", "PAT02",  "PARAM01", "WEEK 1",   as.Date("2021-04-27"), 10.0, "AVERAGE",
    "TEST01", "PAT02",  "PARAM01", "WEEK 2",   as.Date("2021-04-30"), 12.0, "AVERAGE",

    "TEST01", "PAT01",  "PARAM02", "SCREENING", as.Date("2021-04-27"), 15.0, "AVERAGE",
    "TEST01", "PAT01",  "PARAM02", "SCREENING", as.Date("2021-04-25"), 14.0, "AVERAGE",
    "TEST01", "PAT01",  "PARAM02", "SCREENING", as.Date("2021-04-23"), 15.0, NA,
    "TEST01", "PAT01",  "PARAM02", "BASELINE",  as.Date("2021-04-27"), 10.0, "AVERAGE",
    "TEST01", "PAT01",  "PARAM02", "WEEK 2",    as.Date("2021-04-30"), 12.0, NA,

    "TEST01", "PAT02",  "PARAM02", "SCREENING",as.Date("2021-04-27"), 15.0, NA,
    "TEST01", "PAT02",  "PARAM02", "BASELINE", as.Date("2021-04-25"), 14.0, NA,
    "TEST01", "PAT02",  "PARAM02", "WEEK 1",   as.Date("2021-04-23"), 15.0, NA,
    "TEST01", "PAT02",  "PARAM02", "WEEK 1",   as.Date("2021-04-27"), 10.0, NA,
    "TEST01", "PAT02",  "PARAM02", "BASELINE", as.Date("2021-04-30"), 12.0, NA

  )
  expected_output <- tibble::tribble(
    ~STUDYID, ~USUBJID, ~PARAMCD,  ~AVISIT,    ~ADT,                 ~AVAL, ~DTYPE,    ~ABLFL,
    "TEST01", "PAT01",  "PARAM01", "BASELINE", as.Date("2021-04-27"), 15.0, NA,        "",
    "TEST01", "PAT01",  "PARAM01", "BASELINE", as.Date("2021-04-25"), 14.0, NA,        "",
    "TEST01", "PAT01",  "PARAM01", "BASELINE", as.Date("2021-04-23"), 15.0, "AVERAGE", "Y",
    "TEST01", "PAT01",  "PARAM01", "WEEK 1",   as.Date("2021-04-27"), 10.0, "AVERAGE", "",
    "TEST01", "PAT01",  "PARAM01", "WEEK 2",   as.Date("2021-04-30"), 12.0, NA,        "",

    "TEST01", "PAT02",  "PARAM01", "SCREENING", as.Date("2021-04-27"), 15.0, "AVERAGE", "",
    "TEST01", "PAT02",  "PARAM01", "BASELINE", as.Date("2021-04-25"), 14.0, "AVERAGE", "Y",
    "TEST01", "PAT02",  "PARAM01", "BASELINE", as.Date("2021-04-23"), 15.0, "AVERAGE", "",
    "TEST01", "PAT02",  "PARAM01", "WEEK 1",   as.Date("2021-04-27"), 10.0, "AVERAGE", "",
    "TEST01", "PAT02",  "PARAM01", "WEEK 2",   as.Date("2021-04-30"), 12.0, "AVERAGE", "",

    "TEST01", "PAT01",  "PARAM02", "SCREENING", as.Date("2021-04-27"), 15.0, "AVERAGE", "",
    "TEST01", "PAT01",  "PARAM02", "SCREENING", as.Date("2021-04-25"), 14.0, "AVERAGE", "",
    "TEST01", "PAT01",  "PARAM02", "SCREENING", as.Date("2021-04-23"), 15.0, NA,        "",
    "TEST01", "PAT01",  "PARAM02", "BASELINE",  as.Date("2021-04-27"), 10.0, "AVERAGE", "Y",
    "TEST01", "PAT01",  "PARAM02", "WEEK 2",    as.Date("2021-04-30"), 12.0, NA,        "",

    "TEST01", "PAT02",  "PARAM02", "SCREENING",as.Date("2021-04-27"), 15.0, NA,         "",
    "TEST01", "PAT02",  "PARAM02", "BASELINE", as.Date("2021-04-25"), 14.0, NA,         "",
    "TEST01", "PAT02",  "PARAM02", "WEEK 1",   as.Date("2021-04-23"), 15.0, NA,         "",
    "TEST01", "PAT02",  "PARAM02", "WEEK 1",   as.Date("2021-04-27"), 10.0, NA,         "",
    "TEST01", "PAT02",  "PARAM02", "BASELINE", as.Date("2021-04-30"), 12.0, NA,         ""
  )

  actual_output <- derive_var_ablfl(
    input,
    by_vars = exprs(USUBJID, PARAMCD),
    order = exprs(ADT, desc(AVAL)),
    flag_filter = expr(AVISIT == "BASELINE" & DTYPE == "AVERAGE" & is.na(DTYPE) == FALSE)
  )

  expect_dfs_equal(
    expected_output,
    actual_output,
    keys = c("STUDYID", "USUBJID", "PARAMCD", "AVISIT", "ADT", "AVAL")
  )
})
