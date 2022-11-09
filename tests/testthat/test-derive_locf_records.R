test_that("derive_locf_records Test 1: new locf observations are derived correctly when visits are missing", {
  library(tibble)
  library(tidyr)


  input <- tribble(
      ~STUDYID,       ~USUBJID,     ~PARAMCD, ~PARAM,                          ~AVAL, ~AVISITN, ~AVISIT,
      "CDISCPILOT01", "01-701-1015", "DIABP", "Diastolic Blood Pressure (mmHg)", 51,   0,        "BASELINE",
      "CDISCPILOT01", "01-701-1015", "DIABP", "Diastolic Blood Pressure (mmHg)", 50,   2,        "WEEK 2",
      "CDISCPILOT01", "01-701-1015", "SYSBP", "Systolic Blood Pressure (mmHg)", 121,   0,        "BASELINE",
      "CDISCPILOT01", "01-701-1015", "SYSBP", "Systolic Blood Pressure (mmHg)", 121,   2,        "WEEK 2",
      "CDISCPILOT01", "01-701-1028", "DIABP", "Diastolic Blood Pressure (mmHg)", 79,   0,        "BASELINE",
      "CDISCPILOT01", "01-701-1028", "SYSBP", "Systolic Blood Pressure (mmHg)", 130,   0,        "BASELINE"
    )

  advs_expected_obsv <- tibble::tribble(
   ~PARAMCD, ~PARAM,                           ~AVISITN, ~AVISIT,
    "DIABP", "Diastolic Blood Pressure (mmHg)", 0,       "BASELINE",
    "DIABP", "Diastolic Blood Pressure (mmHg)", 2,       "WEEK 2",
    "SYSBP", "Systolic Blood Pressure (mmHg)",  0,       "BASELINE",
    "SYSBP", "Systolic Blood Pressure (mmHg)",  2,       "WEEK 2",
  )


  expected_output <- bind_rows(
    input,
    tibble::tribble(
      ~STUDYID,       ~USUBJID,     ~PARAMCD, ~PARAM,                           ~AVAL, ~AVISITN, ~AVISIT,
      "CDISCPILOT01", "01-701-1028", "DIABP", "Diastolic Blood Pressure (mmHg)", 79,    2,        "WEEK 2",
      "CDISCPILOT01", "01-701-1028", "SYSBP", "Systolic Blood Pressure (mmHg)", 130,    2,        "WEEK 2"
  ) %>%
    mutate(DTYPE = "LOCF"))


  actual_output <- derive_locf_records(
    input,
    dataset_expected_obs=advs_expected_obsv,
    by_vars=vars(STUDYID, USUBJID, PARAM, PARAMCD),
    order=vars(STUDYID, USUBJID, PARAM, PARAMCD, AVISITN, AVISIT))


  expect_dfs_equal(
    base = expected_output,
    compare = actual_output,
    keys = c("STUDYID", "USUBJID", "PARAMCD", "PARAM", "AVISITN", "AVISIT", "DTYPE")
  )

})


test_that("derive_locf_records Test 2: new locf observations are derived correctly when some visits have missing AVAL", {
  library(tibble)

  input <- tribble(
    ~STUDYID,       ~USUBJID,     ~PARAMCD, ~PARAM,                          ~AVAL, ~AVISITN, ~AVISIT,
    "CDISCPILOT01", "01-701-1015", "DIABP", "Diastolic Blood Pressure (mmHg)", 51,   0,        "BASELINE",
    "CDISCPILOT01", "01-701-1015", "DIABP", "Diastolic Blood Pressure (mmHg)", 50,   2,        "WEEK 2",
    "CDISCPILOT01", "01-701-1015", "SYSBP", "Systolic Blood Pressure (mmHg)", 121,   0,        "BASELINE",
    "CDISCPILOT01", "01-701-1015", "SYSBP", "Systolic Blood Pressure (mmHg)", 121,   2,        "WEEK 2",
    "CDISCPILOT01", "01-701-1028", "DIABP", "Diastolic Blood Pressure (mmHg)", 79,   0,        "BASELINE",
    "CDISCPILOT01", "01-701-1028", "DIABP", "Diastolic Blood Pressure (mmHg)", NA,   2,        "WEEK 2",
    "CDISCPILOT01", "01-701-1028", "SYSBP", "Systolic Blood Pressure (mmHg)", 130,   0,        "BASELINE",
    "CDISCPILOT01", "01-701-1028", "SYSBP", "Systolic Blood Pressure (mmHg)", NA,    2,        "WEEK 2"

  )

  advs_expected_obsv <- tibble::tribble(
    ~PARAMCD, ~PARAM,                           ~AVISITN, ~AVISIT,
    "DIABP", "Diastolic Blood Pressure (mmHg)", 0,       "BASELINE",
    "DIABP", "Diastolic Blood Pressure (mmHg)", 2,       "WEEK 2",
    "SYSBP", "Systolic Blood Pressure (mmHg)",  0,       "BASELINE",
    "SYSBP", "Systolic Blood Pressure (mmHg)",  2,       "WEEK 2",
  )


  expected_output <- bind_rows(
    input,
    tibble::tribble(
      ~STUDYID,       ~USUBJID,     ~PARAMCD, ~PARAM,                           ~AVAL, ~AVISITN, ~AVISIT,
      "CDISCPILOT01", "01-701-1028", "DIABP", "Diastolic Blood Pressure (mmHg)", 79,    2,        "WEEK 2",
      "CDISCPILOT01", "01-701-1028", "SYSBP", "Systolic Blood Pressure (mmHg)", 130,    2,        "WEEK 2"
    ) %>%
      mutate(DTYPE = "LOCF"))


  actual_output <- derive_locf_records(
    input,
    dataset_expected_obs=advs_expected_obsv,
    by_vars=vars(STUDYID, USUBJID, PARAM, PARAMCD),
    order=vars(STUDYID, USUBJID, PARAM, PARAMCD, AVISITN, AVISIT))


  expect_dfs_equal(
    base = expected_output,
    compare = actual_output,
    keys = c("STUDYID", "USUBJID", "PARAMCD", "PARAM", "AVISITN", "AVISIT", "DTYPE")
  )

})

