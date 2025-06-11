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
    dataset_ref = advs_expected_obsv,
    by_vars = exprs(STUDYID, USUBJID, PARAM, PARAMCD),
    order = exprs(AVISITN, AVISIT)
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
    dataset_ref = advs_expected_obsv,
    by_vars = exprs(STUDYID, USUBJID, PARAM, PARAMCD),
    order = exprs(AVISITN, AVISIT)
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
    dataset_ref = advs_expected_obsv,
    by_vars = exprs(STUDYID, USUBJID, PARAM, PARAMCD),
    order = exprs(AVISITN, AVISIT)
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
    dataset_ref = advs_expected_obsv,
    by_vars = exprs(STUDYID, USUBJID, PARAM, PARAMCD),
    order = exprs(AVISITN, AVISIT)
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
    dataset_ref = advs_expected_obsv,
    by_vars = exprs(STUDYID, USUBJID, PARAM, PARAMCD),
    order = exprs(AVISITN, AVISIT)
  )

  expect_dfs_equal(
    base = expected_output,
    compare = actual_output,
    keys = c("STUDYID", "USUBJID", "PARAMCD", "PARAM", "AVISITN", "AVISIT", "DTYPE")
  )
})


## Test 6: analysis visits are assigned based on time windows ----
test_that("derive_locf_records Test 6: analysis visits are assigned based on time windows", {
  input <- tibble::tribble(
    ~USUBJID, ~PARAMCD, ~AVAL, ~AVISITN, ~AVISIT, ~ADY,
    "1", "DIABP", 51, 0, "BASELINE", 0,
    "1", "DIABP", 50, 2, "WEEK 2", 14,
    "1", "DIABP", 52, 4, "WEEK 4", 28,
    "1", "DIABP", 54, 6, "WEEK 6", 42,
    "1", "SYSBP", 21, 0, "BASELINE", 0,
    "1", "SYSBP", 121, 2, "WEEK 2", 14,
    "2", "DIABP", 79, 0, "BASELINE", 0,
    "2", "DIABP", 80, 2, "WEEK 2", 12,
    "2", "DIABP", NA, 4, "WEEK 4", 26,
    "2", "DIABP", NA, 6, "WEEK 6", 44,
    "2", "SYSBP", 130, 0, "BASELINE", 0
  )

  advs_expected_obsv <- tibble::tribble(
    ~PARAMCD, ~AVISITN, ~AVISIT, ~ADY,
    "DIABP", 0, "BASELINE", 0,
    "DIABP", 2, "WEEK 2", 14,
    "DIABP", 4, "WEEK 4", 28,
    "DIABP", 6, "WEEK 6", 42,
    "SYSBP", 0, "BASELINE", 0,
    "SYSBP", 2, "WEEK 2", 14
  )

  expected_output <- bind_rows(
    input,
    tibble::tribble(
      ~USUBJID, ~PARAMCD, ~AVAL, ~AVISITN, ~AVISIT, ~ADY,
      "2", "DIABP", 80, 4, "WEEK 4", 28,
      "2", "DIABP", 80, 6, "WEEK 6", 42,
      "2", "SYSBP", 130, 2, "WEEK 2", 14
    ) %>%
      mutate(DTYPE = "LOCF")
  )

  actual_output <- derive_locf_records(
    input,
    dataset_ref = advs_expected_obsv,
    analysis_var = AVAL,
    by_vars = exprs(USUBJID, PARAMCD),
    id_vars_ref = exprs(USUBJID, PARAMCD, AVISITN, AVISIT),
    order = exprs(AVISITN, AVISIT, ADY)
  )

  expect_dfs_equal(
    base = expected_output,
    compare = actual_output,
    keys = c("USUBJID", "PARAMCD", "AVISITN", "AVISIT", "DTYPE")
  )
})


## Test 7: imputation = update ----
test_that("derive_locf_records Test 7: imputation = update", {
  input <- tibble::tribble(
    ~USUBJID, ~PARAMCD, ~AVAL, ~AVISITN, ~AVISIT, ~ADY,
    "1", "DIABP", 51, 0, "BASELINE", 0,
    "1", "DIABP", 50, 2, "WEEK 2", 14,
    "1", "DIABP", 52, 4, "WEEK 4", 28,
    "1", "DIABP", 54, 6, "WEEK 6", 42,
    "1", "SYSBP", 121, 0, "BASELINE", 0,
    "1", "SYSBP", 121, 2, "WEEK 2", 14,
    "2", "DIABP", 79, 0, "BASELINE", 0,
    "2", "DIABP", 80, 2, "WEEK 2", 12,
    "2", "DIABP", NA, 4, "WEEK 4", 28,
    "2", "DIABP", NA, 6, "WEEK 6", 44,
    "2", "SYSBP", 130, 0, "BASELINE", 0
  )

  advs_expected_obsv <- tibble::tribble(
    ~PARAMCD, ~AVISITN, ~AVISIT, ~ADY,
    "DIABP", 0, "BASELINE", 0,
    "DIABP", 2, "WEEK 2", 14,
    "DIABP", 4, "WEEK 4", 28,
    "DIABP", 6, "WEEK 6", 42,
    "SYSBP", 0, "BASELINE", 0,
    "SYSBP", 2, "WEEK 2", 14
  )

  expected_output <- bind_rows(
    input,
    tibble::tribble(
      ~USUBJID, ~PARAMCD, ~AVAL, ~AVISITN, ~AVISIT, ~ADY,
      "2", "DIABP", 80, 4, "WEEK 4", 28,
      "2", "DIABP", 80, 6, "WEEK 6", 44,
      "2", "SYSBP", 130, 2, "WEEK 2", 14
    ) %>%
      mutate(DTYPE = "LOCF")
  ) %>%
    filter(!(USUBJID == "2" & PARAMCD == "DIABP" & AVISITN %in% c(4, 6) & is.na(DTYPE)))

  actual_output <- derive_locf_records(
    input,
    dataset_ref = advs_expected_obsv,
    analysis_var = AVAL,
    by_vars = exprs(USUBJID, PARAMCD),
    id_vars_ref = exprs(USUBJID, PARAMCD, AVISITN, AVISIT),
    imputation = "update",
    order = exprs(AVISITN, AVISIT, ADY)
  )

  expect_dfs_equal(
    base = expected_output,
    compare = actual_output,
    keys = c("USUBJID", "PARAMCD", "AVISITN", "AVISIT", "DTYPE")
  )
})


## Test 8: imputation = update_add ----
test_that("derive_locf_records Test 8: imputation = update_add", {
  input <- tibble::tribble(
    ~USUBJID, ~PARAMN, ~PARAMCD, ~AVAL, ~AVISITN, ~AVISIT,    ~ADY,
    "1",      1,       "DIABP",  51,    0,        "BASELINE", 0,
    "1",      1,       "DIABP",  50,    2,        "WEEK 2",   14,
    "1",      1,       "DIABP",  52,    4,        "WEEK 4",   28,
    "1",      1,       "DIABP",  54,    6,        "WEEK 6",   42,
    "1",      2,       "SYSBP",  121,   0,        "BASELINE", 0,
    "1",      2,       "SYSBP",  121,   2,        "WEEK 2",   14,
    "2",      1,       "DIABP",  79,    0,        "BASELINE", 0,
    "2",      1,       "DIABP",  80,    2,        "WEEK 2",   12,
    "2",      1,       "DIABP",  NA,    4,        "WEEK 4",   28,
    "2",      1,       "DIABP",  NA,    6,        "WEEK 6",   44,
    "2",      2,       "SYSBP",  130,   0,        "BASELINE", 0
  )

  advs_expected_obsv <- tibble::tribble(
    ~PARAMCD, ~AVISITN, ~AVISIT,    ~ADY,
    "DIABP",  0,        "BASELINE", 0,
    "DIABP",  2,        "WEEK 2",   14,
    "DIABP",  4,        "WEEK 4",   28,
    "DIABP",  6,        "WEEK 6",   42,
    "SYSBP",  0,        "BASELINE", 0,
    "SYSBP",  2,        "WEEK 2",   14
  )

  expected_output <- bind_rows(
    input,
    tibble::tribble(
      ~USUBJID, ~PARAMN, ~PARAMCD, ~AVAL, ~AVISITN, ~AVISIT,  ~ADY,
      "2",      1,       "DIABP",  80,    4,        "WEEK 4", 28,
      "2",      1,       "DIABP",  80,    6,        "WEEK 6", 44,
      "2",      2,       "SYSBP",  130,   2,        "WEEK 2", 14
    ) %>%
      mutate(DTYPE = "LOCF")
  )

  actual_output <- derive_locf_records(
    input,
    dataset_ref = advs_expected_obsv,
    analysis_var = AVAL,
    by_vars = exprs(USUBJID, PARAMCD),
    id_vars_ref = exprs(USUBJID, PARAMCD, AVISITN, AVISIT),
    imputation = "update_add",
    order = exprs(AVISITN, AVISIT, ADY),
    keep_vars = exprs(PARAMN)
  )

  expect_dfs_equal(
    base = expected_output,
    compare = actual_output,
    keys = c("USUBJID", "PARAMCD", "AVISITN", "AVISIT", "DTYPE")
  )
})


## Test 9: fill variables other than analysis_var for imputed records ----
test_that("derive_locf_records Test 9: fill variables other than analysis_var for imputed
          records", {
  input <- tibble::tribble(
    ~USUBJID, ~PARAMCD, ~PARAMN, ~AVAL, ~AVISITN, ~AVISIT,    ~DATEC,      ~DAY,
    "1",      "DIABP",  2,       85,    0,        "BASELINE", "February",  "day a",
    "1",      "DIABP",  2,       50,    4,        "VISIT 4",  "April",     "day c",
    "1",      "DIABP",  2,       20,    6,        "VISIT 6",  "May",       "day d",
    "1",      "DIABP",  2,       35,    8,        "VISIT 8",  "June",      "day e",
    "1",      "DIABP",  2,       NA,    10,       "VISIT 10", "July",      "day f",
    "1",      "DIABP",  2,       20,    12,       "VISIT 12", "August",    "day g",
    "1",      "DIABP",  2,       NA,    14,       "VISIT 14", "September", "day h"
  )

  advs_expected_obsv <- tibble::tribble(
    ~PARAMCD, ~AVISITN, ~AVISIT,
    "DIABP",  0,        "BASELINE",
    "DIABP",  2,        "VISIT 2",
    "DIABP",  4,        "VISIT 4",
    "DIABP",  6,        "VISIT 6",
    "DIABP",  8,        "VISIT 8",
    "DIABP",  10,       "VISIT 10",
    "DIABP",  12,       "VISIT 12",
    "DIABP",  14,       "VISIT 14"
  )

  expected_output <- bind_rows(
    input,
    tibble::tribble(
      ~USUBJID, ~PARAMCD, ~PARAMN, ~AVAL, ~AVISITN, ~AVISIT, ~DATEC, ~DAY,
      "1", "DIABP", 2, 85, 2, "VISIT 2", "February", "day a",
      "1", "DIABP", 2, 35, 10, "VISIT 10", "July", "day f",
      "1", "DIABP", 2, 20, 14, "VISIT 14", "September", "day h"
    ) %>%
      mutate(DTYPE = "LOCF")
  )

  actual_output <- derive_locf_records(
    dataset = input,
    dataset_ref = advs_expected_obsv,
    by_vars = exprs(USUBJID, PARAMCD),
    imputation = "update_add",
    order = exprs(AVISITN, AVISIT),
    keep_vars = exprs(PARAMN, DATEC, DAY),
  ) |>
    arrange(USUBJID, PARAMCD, AVISITN, desc(is.na(DTYPE)), DTYPE)

  expect_dfs_equal(
    base = expected_output,
    compare = actual_output,
    keys = c("USUBJID", "PARAMCD", "PARAMN", "AVISITN", "AVISIT", "DTYPE")
  )
})
