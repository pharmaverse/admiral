## Test 1: visits are missing ----
test_that("derive_locf_records Test 1: visits are missing", {
  input <- tibble::tribble(
    ~STUDYID, ~USUBJID, ~PARAMCD, ~PARAM, ~AVAL, ~AVISITN, ~AVISIT,
    "TEST01", "01-701-1015", "DIABP", "Diastolic Blood Pressure (mmHg)", 51, 0, "BASELINE",
    "TEST01", "01-701-1015", "DIABP", "Diastolic Blood Pressure (mmHg)", 50, 2, "WEEK 2",
    "TEST01", "01-701-1015", "SYSBP", "Systolic Blood Pressure (mmHg)", 121, 0, "BASELINE",
    "TEST01", "01-701-1015", "SYSBP", "Systolic Blood Pressure (mmHg)", 121, 2, "WEEK 2",
    "TEST01", "01-701-1028", "DIABP", "Diastolic Blood Pressure (mmHg)", 79, 0, "BASELINE",
    "TEST01", "01-701-1028", "SYSBP", "Systolic Blood Pressure (mmHg)", 130, 0, "BASELINE"
  )

  advs_expected_obsv <- tibble::tribble(
    ~PARAMCD, ~PARAM, ~AVISITN, ~AVISIT,
    "DIABP", "Diastolic Blood Pressure (mmHg)", 0, "BASELINE",
    "DIABP", "Diastolic Blood Pressure (mmHg)", 2, "WEEK 2",
    "SYSBP", "Systolic Blood Pressure (mmHg)", 0, "BASELINE",
    "SYSBP", "Systolic Blood Pressure (mmHg)", 2, "WEEK 2",
  )


  expected_output <- bind_rows(
    input,
    tibble::tribble(
      ~STUDYID, ~USUBJID, ~PARAMCD, ~PARAM, ~AVAL, ~AVISITN, ~AVISIT,
      "TEST01", "01-701-1028", "DIABP", "Diastolic Blood Pressure (mmHg)", 79, 2, "WEEK 2",
      "TEST01", "01-701-1028", "SYSBP", "Systolic Blood Pressure (mmHg)", 130, 2, "WEEK 2"
    ) %>%
      mutate(DTYPE = "LOCF")
  )


  actual_output <- derive_locf_records(
    input,
    dataset_expected_obs = advs_expected_obsv,
    by_vars = vars(STUDYID, USUBJID, PARAM, PARAMCD),
    order = vars(AVISITN, AVISIT)
  )


  expect_dfs_equal(
    base = expected_output,
    compare = actual_output,
    keys = c("STUDYID", "USUBJID", "PARAMCD", "PARAM", "AVISITN", "AVISIT", "DTYPE")
  )
})


## Test 2: some visits have missing AVAL ----
test_that("derive_locf_records Test 2: some visits have missing AVAL", {
  input <- tibble::tribble(
    ~STUDYID, ~USUBJID, ~PARAMCD, ~PARAM, ~AVAL, ~AVISITN, ~AVISIT,
    "TEST01", "01-701-1015", "DIABP", "Diastolic Blood Pressure (mmHg)", 51, 0, "BASELINE",
    "TEST01", "01-701-1015", "DIABP", "Diastolic Blood Pressure (mmHg)", 50, 2, "WEEK 2",
    "TEST01", "01-701-1015", "SYSBP", "Systolic Blood Pressure (mmHg)", 121, 0, "BASELINE",
    "TEST01", "01-701-1015", "SYSBP", "Systolic Blood Pressure (mmHg)", 121, 2, "WEEK 2",
    "TEST01", "01-701-1028", "DIABP", "Diastolic Blood Pressure (mmHg)", 79, 0, "BASELINE",
    "TEST01", "01-701-1028", "DIABP", "Diastolic Blood Pressure (mmHg)", NA, 2, "WEEK 2",
    "TEST01", "01-701-1028", "SYSBP", "Systolic Blood Pressure (mmHg)", 130, 0, "BASELINE",
    "TEST01", "01-701-1028", "SYSBP", "Systolic Blood Pressure (mmHg)", NA, 2, "WEEK 2"
  )

  advs_expected_obsv <- tibble::tribble(
    ~PARAMCD, ~PARAM, ~AVISITN, ~AVISIT,
    "DIABP", "Diastolic Blood Pressure (mmHg)", 0, "BASELINE",
    "DIABP", "Diastolic Blood Pressure (mmHg)", 2, "WEEK 2",
    "SYSBP", "Systolic Blood Pressure (mmHg)", 0, "BASELINE",
    "SYSBP", "Systolic Blood Pressure (mmHg)", 2, "WEEK 2",
  )


  expected_output <- bind_rows(
    input,
    tibble::tribble(
      ~STUDYID, ~USUBJID, ~PARAMCD, ~PARAM, ~AVAL, ~AVISITN, ~AVISIT,
      "TEST01", "01-701-1028", "DIABP", "Diastolic Blood Pressure (mmHg)", 79, 2, "WEEK 2",
      "TEST01", "01-701-1028", "SYSBP", "Systolic Blood Pressure (mmHg)", 130, 2, "WEEK 2"
    ) %>%
      mutate(DTYPE = "LOCF")
  )


  actual_output <- derive_locf_records(
    input,
    dataset_expected_obs = advs_expected_obsv,
    by_vars = vars(STUDYID, USUBJID, PARAM, PARAMCD),
    order = vars(AVISITN, AVISIT)
  )


  expect_dfs_equal(
    base = expected_output,
    compare = actual_output,
    keys = c("STUDYID", "USUBJID", "PARAMCD", "PARAM", "AVISITN", "AVISIT", "DTYPE")
  )
})


## Test 3: visits are missing - and DTYPE already exits ----
test_that("derive_locf_records Test 3: visits are missing - and DTYPE already exits", {
  input <- tibble::tribble(
    ~STUDYID, ~USUBJID, ~PARAMCD, ~PARAM, ~AVAL, ~AVISITN, ~AVISIT, ~DTYPE,
    "TEST01", "1015", "DIABP", "Diastolic Blood Pressure", 51, 0, "BASELINE", NA,
    "TEST01", "1015", "DIABP", "Diastolic Blood Pressure", 50, 2, "WEEK 2", NA,
    "TEST01", "1015", "SYSBP", "Systolic Blood Pressure", 121, 0, "BASELINE", NA,
    "TEST01", "1015", "SYSBP", "Systolic Blood Pressure", 121, 2, "WEEK 2", NA,
    "TEST01", "1015", "LTDIABP", "Log(Diastolic Blood Pressure)", 1.71, 0, "BASELINE", "LOG",
    "TEST01", "1015", "LTDIABP", "Log(Diastolic Blood Pressure)", 1.69, 2, "WEEK 2", "LOG",
    "TEST01", "1015", "LTSYSBP", "Log(Systolic Blood Pressure)", 2.08, 0, "BASELINE", "LOG",
    "TEST01", "1015", "LTSYSBP", "Log(Systolic Blood Pressure)", 2.08, 2, "WEEK 2", "LOG",
    "TEST01", "1028", "DIABP", "Diastolic Blood Pressure", 79, 0, "BASELINE", NA,
    "TEST01", "1028", "SYSBP", "Systolic Blood Pressure", 130, 0, "BASELINE", NA,
    "TEST01", "1028", "LTDIABP", "Log(Diastolic Blood Pressure)", 1.89, 0, "BASELINE", "LOG",
    "TEST01", "1028", "LTSYSBP", "Log(Systolic Blood Pressure)", 2.11, 0, "BASELINE", "LOG"
  )

  advs_expected_obsv <- tibble::tribble(
    ~PARAMCD, ~PARAM, ~AVISITN, ~AVISIT,
    "DIABP", "Diastolic Blood Pressure", 0, "BASELINE",
    "DIABP", "Diastolic Blood Pressure", 2, "WEEK 2",
    "LTDIABP", "Log(Diastolic Blood Pressure)", 0, "BASELINE",
    "LTDIABP", "Log(Diastolic Blood Pressure)", 2, "WEEK 2",
    "SYSBP", "Systolic Blood Pressure", 0, "BASELINE",
    "SYSBP", "Systolic Blood Pressure", 2, "WEEK 2",
    "LTSYSBP", "Log(Systolic Blood Pressure)", 0, "BASELINE",
    "LTSYSBP", "Log(Systolic Blood Pressure)", 2, "WEEK 2"
  )


  expected_output <- bind_rows(
    input,
    tibble::tribble(
      ~STUDYID, ~USUBJID, ~PARAMCD, ~PARAM, ~AVAL, ~AVISITN, ~AVISIT,
      "TEST01", "1028", "DIABP", "Diastolic Blood Pressure", 79, 2, "WEEK 2",
      "TEST01", "1028", "LTDIABP", "Log(Diastolic Blood Pressure)", 1.89, 2, "WEEK 2",
      "TEST01", "1028", "SYSBP", "Systolic Blood Pressure", 130, 2, "WEEK 2",
      "TEST01", "1028", "LTSYSBP", "Log(Systolic Blood Pressure)", 2.11, 2, "WEEK 2"
    ) %>%
      mutate(DTYPE = "LOCF")
  )


  actual_output <- derive_locf_records(
    input,
    dataset_expected_obs = advs_expected_obsv,
    by_vars = vars(STUDYID, USUBJID, PARAM, PARAMCD),
    order = vars(AVISITN, AVISIT)
  )


  expect_dfs_equal(
    base = expected_output,
    compare = actual_output,
    keys = c("STUDYID", "USUBJID", "PARAMCD", "PARAM", "AVISITN", "AVISIT", "DTYPE")
  )
})


## Test 4: visit variables are parameter independent ----
test_that("derive_locf_records Test 4: visit variables are parameter independent", {
  input <- tibble::tribble(
    ~STUDYID, ~USUBJID, ~PARAMCD, ~PARAM, ~AVAL, ~AVISITN, ~AVISIT,
    "TEST01", "01-701-1015", "DIABP", "Diastolic Blood Pressure (mmHg)", 51, 0, "BASELINE",
    "TEST01", "01-701-1015", "DIABP", "Diastolic Blood Pressure (mmHg)", 50, 2, "WEEK 2",
    "TEST01", "01-701-1015", "SYSBP", "Systolic Blood Pressure (mmHg)", 121, 0, "BASELINE",
    "TEST01", "01-701-1015", "SYSBP", "Systolic Blood Pressure (mmHg)", 121, 2, "WEEK 2",
    "TEST01", "01-701-1028", "DIABP", "Diastolic Blood Pressure (mmHg)", 79, 0, "BASELINE",
    "TEST01", "01-701-1028", "DIABP", "Diastolic Blood Pressure (mmHg)", NA, 2, "WEEK 2",
    "TEST01", "01-701-1028", "SYSBP", "Systolic Blood Pressure (mmHg)", 130, 0, "BASELINE",
    "TEST01", "01-701-1028", "SYSBP", "Systolic Blood Pressure (mmHg)", NA, 2, "WEEK 2"
  )

  advs_expected_obsv <- tibble::tribble(
    ~AVISITN, ~AVISIT,
    0,        "BASELINE",
    2,        "WEEK 2"
  )


  expected_output <- bind_rows(
    input,
    tibble::tribble(
      ~STUDYID, ~USUBJID, ~PARAMCD, ~PARAM, ~AVAL, ~AVISITN, ~AVISIT,
      "TEST01", "01-701-1028", "DIABP", "Diastolic Blood Pressure (mmHg)", 79, 2, "WEEK 2",
      "TEST01", "01-701-1028", "SYSBP", "Systolic Blood Pressure (mmHg)", 130, 2, "WEEK 2"
    ) %>%
      mutate(DTYPE = "LOCF")
  )


  actual_output <- derive_locf_records(
    input,
    dataset_expected_obs = advs_expected_obsv,
    by_vars = vars(STUDYID, USUBJID, PARAM, PARAMCD),
    order = vars(AVISITN, AVISIT)
  )


  expect_dfs_equal(
    base = expected_output,
    compare = actual_output,
    keys = c("STUDYID", "USUBJID", "PARAMCD", "PARAM", "AVISITN", "AVISIT", "DTYPE")
  )
})


## Test 5: visit variables are parameter dependent ----
test_that("derive_locf_records Test 5: visit variables are parameter dependent", {
  input <- tibble::tribble(
    ~STUDYID, ~USUBJID, ~PARAMCD, ~PARAM, ~AVAL, ~AVISITN, ~AVISIT,
    "TEST01", "01-701-1015", "DIABP", "Diastolic Blood Pressure (mmHg)", 51, 0, "BASELINE",
    "TEST01", "01-701-1015", "DIABP", "Diastolic Blood Pressure (mmHg)", 50, 2, "WEEK 2",
    "TEST01", "01-701-1015", "DIABP", "Diastolic Blood Pressure (mmHg)", 52, 4, "WEEK 4",
    "TEST01", "01-701-1015", "DIABP", "Diastolic Blood Pressure (mmHg)", 54, 6, "WEEK 6",
    "TEST01", "01-701-1015", "SYSBP", "Systolic Blood Pressure (mmHg)", 121, 0, "BASELINE",
    "TEST01", "01-701-1015", "SYSBP", "Systolic Blood Pressure (mmHg)", 121, 2, "WEEK 2",
    "TEST01", "01-701-1028", "DIABP", "Diastolic Blood Pressure (mmHg)", 79, 0, "BASELINE",
    "TEST01", "01-701-1028", "DIABP", "Diastolic Blood Pressure (mmHg)", 80, 2, "WEEK 2",
    "TEST01", "01-701-1028", "DIABP", "Diastolic Blood Pressure (mmHg)", NA, 4, "WEEK 4",
    "TEST01", "01-701-1028", "DIABP", "Diastolic Blood Pressure (mmHg)", NA, 6, "WEEK 6",
    "TEST01", "01-701-1028", "SYSBP", "Systolic Blood Pressure (mmHg)", 130, 0, "BASELINE"
  )

  advs_expected_obsv <- tibble::tribble(
    ~PARAMCD, ~PARAM, ~AVISITN, ~AVISIT,
    "DIABP", "Diastolic Blood Pressure (mmHg)", 0, "BASELINE",
    "DIABP", "Diastolic Blood Pressure (mmHg)", 2, "WEEK 2",
    "DIABP", "Diastolic Blood Pressure (mmHg)", 4, "WEEK 4",
    "DIABP", "Diastolic Blood Pressure (mmHg)", 6, "WEEK 6",
    "SYSBP", "Systolic Blood Pressure (mmHg)", 0, "BASELINE",
    "SYSBP", "Systolic Blood Pressure (mmHg)", 2, "WEEK 2"
  )


  expected_output <- bind_rows(
    input,
    tibble::tribble(
      ~STUDYID, ~USUBJID, ~PARAMCD, ~PARAM, ~AVAL, ~AVISITN, ~AVISIT,
      "TEST01", "01-701-1028", "DIABP", "Diastolic Blood Pressure (mmHg)", 80, 4, "WEEK 4",
      "TEST01", "01-701-1028", "DIABP", "Diastolic Blood Pressure (mmHg)", 80, 6, "WEEK 6",
      "TEST01", "01-701-1028", "SYSBP", "Systolic Blood Pressure (mmHg)", 130, 2, "WEEK 2"
    ) %>%
      mutate(DTYPE = "LOCF")
  )


  actual_output <- derive_locf_records(
    input,
    dataset_expected_obs = advs_expected_obsv,
    by_vars = vars(STUDYID, USUBJID, PARAM, PARAMCD),
    order = vars(AVISITN, AVISIT)
  )


  expect_dfs_equal(
    base = expected_output,
    compare = actual_output,
    keys = c("STUDYID", "USUBJID", "PARAMCD", "PARAM", "AVISITN", "AVISIT", "DTYPE")
  )
})
