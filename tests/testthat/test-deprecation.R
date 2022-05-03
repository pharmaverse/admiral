test_that("a warning is issued when using `derive_var_extreme_flag()` with `filter` argument", {
  input <- tibble::tribble(
    ~USUBJID, ~AVISITN, ~AVAL,
    1, 1, 12,
    1, 3, 9,
    2, 2, 42,
    3, 3, 14,
    3, 3, 10
  )

  expect_warning(
    derive_var_extreme_flag(
      input,
      by_vars = vars(USUBJID),
      order = vars(AVISITN, desc(AVAL)),
      new_var = firstfl,
      mode = "first",
      filter = !is.na(AVAL)
    ),
    paste(
      "`filter` is deprecated as of admiral 0.7.0.",
      "Please use `restrict_derivation()` instead (see examples).",
      sep = "\n"
    ),
    fixed = TRUE
  )
}
)

test_that("a warning is issued when using `derive_var_worst_flag()` with `filter` argument", {
  input_worst_flag <- tibble::tribble(
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
    "TEST01", "PAT01",  "PARAM02", "BASELINE", as.Date("2021-04-27"), 10.0,
    "TEST01", "PAT01",  "PARAM02", "WEEK 2",   as.Date("2021-04-30"), 12.0,

    "TEST01", "PAT02",  "PARAM02", "SCREENING", as.Date("2021-04-27"), 15.0,
    "TEST01", "PAT02",  "PARAM02", "BASELINE", as.Date("2021-04-25"), 14.0,
    "TEST01", "PAT02",  "PARAM02", "WEEK 1",   as.Date("2021-04-23"), 15.0,
    "TEST01", "PAT02",  "PARAM02", "WEEK 1",   as.Date("2021-04-27"), 10.0,
    "TEST01", "PAT02",  "PARAM02", "BASELINE", as.Date("2021-04-30"), 12.0,

    "TEST01", "PAT02",  "PARAM03", "SCREENING", as.Date("2021-04-27"), 15.0,
    "TEST01", "PAT02",  "PARAM03", "BASELINE", as.Date("2021-04-25"), 14.0,
    "TEST01", "PAT02",  "PARAM03", "WEEK 1",   as.Date("2021-04-23"), 15.0,
    "TEST01", "PAT02",  "PARAM03", "WEEK 1",   as.Date("2021-04-27"), 10.0,
    "TEST01", "PAT02",  "PARAM03", "BASELINE", as.Date("2021-04-30"), 12.0
  )

  expect_warning(
    derive_var_worst_flag(
      input_worst_flag,
      by_vars = vars(USUBJID, PARAMCD, AVISIT),
      order = vars(desc(ADT)),
      new_var = WORSTFL,
      param_var = PARAMCD,
      analysis_var = AVAL,
      worst_high = c("PARAM01", "PARAM03"),
      worst_low = "PARAM02",
      filter = !is.na(AVAL)
    ),
    paste(
      "`filter` is deprecated as of admiral 0.7.0.",
      "Please use `restrict_derivation()` instead (see examples).",
      sep = "\n"
    ),
    fixed = TRUE
  )
}
)
