adsl <- tibble::tribble(
  ~USUBJID, ~SEX, ~COUNTRY,
  "ST42-1", "F",  "AUT",
  "ST42-2", "M",  "MWI",
  "ST42-3", "M",  "NOR",
  "ST42-4", "F",  "UGA"
) %>% mutate(STUDYID = "ST42")

ex <- tibble::tribble(
  ~USUBJID, ~EXSTDTC,
  "ST42-1", "2020-12-07",
  "ST42-1", "2020-12-14",
  "ST42-2", "2021-01-12T12:00:00",
  "ST42-2", "2021-01-26T13:21",
  "ST42-3", "2021-03-02"
) %>% mutate(STUDYID = "ST42")

input_worst_flag <- tibble::tribble(
  ~STUDYID, ~USUBJID, ~PARAMCD, ~AVISIT, ~ADT, ~AVAL,
  "TEST01", "PAT01", "PARAM01", "BASELINE", as.Date("2021-04-27"), 15.0,
  "TEST01", "PAT01", "PARAM01", "BASELINE", as.Date("2021-04-25"), 14.0,
  "TEST01", "PAT01", "PARAM01", "WEEK 2", as.Date("2021-04-30"), 12.0,
  "TEST01", "PAT02", "PARAM01", "WEEK 2", as.Date("2021-04-30"), 12.0,
  "TEST01", "PAT01", "PARAM02", "SCREENING", as.Date("2021-04-27"), 15.0,
  "TEST01", "PAT01", "PARAM02", "BASELINE", as.Date("2021-04-27"), 10.0,
  "TEST01", "PAT02", "PARAM02", "BASELINE", as.Date("2021-04-30"), 12.0,
  "TEST01", "PAT02", "PARAM03", "SCREENING", as.Date("2021-04-27"), 15.0,
  "TEST01", "PAT02", "PARAM03", "WEEK 1", as.Date("2021-04-27"), 10.0,
  "TEST01", "PAT02", "PARAM03", "BASELINE", as.Date("2021-04-30"), 12.0
)

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
