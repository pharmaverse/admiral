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


## Test 1: An error is issued if `derive_derived_param()` is called ----
test_that("deprecation Test 1: An error is issued if `derive_derived_param()`
          is called", {
  input <- tibble::tribble(
    ~USUBJID, ~PARAMCD, ~PARAM, ~AVAL, ~AVALU, ~VISIT,
    "01-701-1015", "DIABP", "Diastolic Blood Pressure (mmHg)", 51, "mmHg", "BASELINE",
    "01-701-1015", "DIABP", "Diastolic Blood Pressure (mmHg)", 50, "mmHg", "WEEK 2",
    "01-701-1015", "SYSBP", "Systolic Blood Pressure (mmHg)", 121, "mmHg", "BASELINE",
    "01-701-1015", "SYSBP", "Systolic Blood Pressure (mmHg)", 121, "mmHg", "WEEK 2",
    "01-701-1028", "DIABP", "Diastolic Blood Pressure (mmHg)", 79, "mmHg", "BASELINE",
    "01-701-1028", "DIABP", "Diastolic Blood Pressure (mmHg)", 80, "mmHg", "WEEK 2",
    "01-701-1028", "SYSBP", "Systolic Blood Pressure (mmHg)", 130, "mmHg", "BASELINE",
    "01-701-1028", "SYSBP", "Systolic Blood Pressure (mmHg)", 132, "mmHg", "WEEK 2"
  )

  expect_error(
    derive_derived_param(
      input,
      parameters = c("SYSBP", "DIABP"),
      by_vars = exprs(USUBJID, VISIT),
      analysis_value = (AVAL.SYSBP + 2 * AVAL.DIABP) / 3,
      set_values_to = exprs(
        PARAMCD = "MAP",
        PARAM = "Mean arterial pressure (mmHg)",
        AVALU = "mmHg"
      )
    ),
    class = "lifecycle_error_deprecated"
  )
})

## Test 2: derive_vars_merged_dt: a deprecation error is issued ----
test_that("deprecation Test 2: derive_vars_merged_dt: a deprecation error
          is issued", {
  expect_error(
    derive_vars_merged_dt(
      adsl,
      dataset_add = ex,
      order = exprs(TRTSDT),
      flag_imputation = "date",
      by_vars = exprs(STUDYID, USUBJID),
      dtc = EXSTDTC,
      new_vars_prefix = "TRTS",
      mode = "first"
    ),
    class = "lifecycle_error_deprecated"
  )
})

## Test 3: derive_vars_merged_dtm: a deprecation error is issued ----
test_that("deprecation Test 3: derive_vars_merged_dtm: a deprecation error
          is issued", {
  expect_error(
    derive_vars_merged_dtm(
      adsl,
      dataset_add = ex,
      order = exprs(TRTSDTM),
      by_vars = exprs(STUDYID, USUBJID),
      dtc = EXSTDTC,
      new_vars_prefix = "TRTS",
      time_imputation = "first",
      mode = "first"
    ),
    class = "lifecycle_error_deprecated"
  )
})


## Test 4: An error is issued if `derive_var_agegr_ema()` is called ----
test_that("deprecation Test 4: An error is issued if `derive_var_agegr_ema()`
          is called", {
  rlang::with_options(lifecycle_verbosity = "error", {
    expect_error(
      derive_var_agegr_ema(admiral.test::admiral_dm, age_var = AGE, new_var = AGEGR1),
      class = "lifecycle_error_deprecated"
    )
  })
})

## Test 5: An error is issued if `derive_var_agegr_fda()` is called ----
test_that("deprecation Test 5: An error is issued if `derive_var_agegr_fda()`
          is called", {
  rlang::with_options(lifecycle_verbosity = "error", {
    expect_error(
      derive_var_agegr_fda(admiral.test::admiral_dm, age_var = AGE, new_var = AGEGR1),
      class = "lifecycle_error_deprecated"
    )
  })
})

## Test 6: An error is issued if `derive_param_first_event()` is called ----
test_that("deprecation Test 6: An error is issued if `derive_param_first_event()`
          is called", {
  rlang::with_options(lifecycle_verbosity = "error", {
    adsl <- tibble::tribble(
      ~STUDYID, ~USUBJID, ~DTHDT,
      "XX1234", "1",      ymd("2022-05-13"),
      "XX1234", "2",      ymd(""),
      "XX1234", "3",      ymd(""),
    )

    adrs <- tibble::tribble(
      ~USUBJID, ~ADTC,        ~AVALC, ~PARAMCD,
      "1",      "2020-01-02", "PR",   "OVR",
      "1",      "2020-02-01", "CR",   "OVR",
      "1",      "2020-03-01", "CR",   "OVR",
      "1",      "2020-04-01", "SD",   "OVR",
      "2",      "2021-06-15", "SD",   "OVR",
      "2",      "2021-07-16", "PD",   "OVR",
      "2",      "2021-09-14", "PD",   "OVR",
    ) %>%
      mutate(
        STUDYID = "XX1234",
        ADT = ymd(ADTC)
      ) %>%
      select(-ADTC)

    expect_error(
      derive_param_first_event(
        adrs,
        dataset_adsl = adsl,
        dataset_source = adrs,
        filter_source = PARAMCD == "OVR" & AVALC == "PD",
        date_var = ADT,
        set_values_to = exprs(
          PARAMCD = "PD",
          ANL01FL = "Y"
        )
      ),
      class = "lifecycle_error_deprecated"
    )
  })
})

## Test 7: An warning is issued if `derive_var_worst_flag()` is called ----
test_that("deprecation Test 7: A warning is issued if Derive worst flag is called", {
  expect_warning(
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
    class = "lifecycle_warning_deprecated"
  )
})

## Test 8: A warning is issued if derive confirmation flag is called ----
test_that("deprecation Test 8: A warning is issued if derive confirmation flag is called", {
  data <- tibble::tribble(
    ~USUBJID, ~AVISITN, ~AVALC,
    "1",      1,        "PR",
    "1",      2,        "CR",
    "1",      3,        "CR",
    "1",      4,        "SD",
    "1",      5,        "NE",
    "2",      1,        "SD",
    "2",      2,        "PR",
    "2",      3,        "PD",
    "3",      1,        "SD",
    "4",      1,        "PR",
    "4",      2,        "PD",
    "4",      3,        "SD",
    "4",      4,        "SD",
    "4",      5,        "PR"
  )
  expect_error(
    derive_var_confirmation_flag(
      data,
      new_var = CONFFL,
      by_vars = exprs(USUBJID),
      join_vars = exprs(AVALC),
      join_type = "after",
      order = exprs(AVISITN),
      filter = AVALC == "PR" & AVALC.join %in% c("CR", "PR")
    ),
    class = "lifecycle_error_deprecated"
  )
})


## Test 9: A warning is issued if filter_joined is called ----
test_that("deprecation Test 9: A warning is issued if filter_confirmation is called", {
  data <- tibble::tribble(
    ~USUBJID, ~AVISITN, ~AVALC,
    "1",      1,        "PR",
    "1",      2,        "CR",
    "1",      3,        "CR",
    "1",      4,        "SD",
    "1",      5,        "NE",
    "2",      1,        "SD",
    "2",      2,        "PR",
    "2",      3,        "PD",
    "3",      1,        "SD",
    "4",      1,        "PR",
    "4",      2,        "PD",
    "4",      3,        "SD",
    "4",      4,        "SD",
    "4",      5,        "PR"
  )
  expect_error(
    filter_confirmation(
      data,
      by_vars = exprs(USUBJID),
      join_vars = exprs(AVISITN, AVALC),
      join_type = "after",
      order = exprs(AVISITN),
      filter = AVALC == "PR" & AVALC.join %in% c("CR", "PR") &
        AVISITN < AVISITN.join
    ),
    class = "lifecycle_error_deprecated"
  )
})
