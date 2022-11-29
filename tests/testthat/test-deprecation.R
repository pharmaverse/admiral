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

## Test 1: An error is thrown if `derive_var_extreme_flag()` with `filter` argument is called ----
test_that("deprecation Test 1: An error is thrown if `derive_var_extreme_flag()`
          with `filter` argument is called", {
  expect_error(
    derive_var_extreme_flag(
      filter = !is.na(AVAL)
    ),
    class = "lifecycle_error_deprecated"
  )
})

## Test 2: An error is thrown if `derive_var_worst_flag()` with `filter` argument is called ----
test_that("deprecation Test 2: An error is thrown if `derive_var_worst_flag()`
          with `filter` argument is called", {
  expect_error(
    derive_var_worst_flag(
      filter = !is.na(AVAL)
    ),
    class = "lifecycle_error_deprecated"
  )
})

## Test 3: derive_var_ady() An error is thrown if `derive_var_ady()` is called ----
test_that("deprecation Test 3: derive_var_ady() An error is thrown if
          `derive_var_ady()` is called", {
  expect_error(
    derive_var_ady(),
    class = "lifecycle_error_deprecated"
  )
})

## Test 4: An error is thrown if `derive_var_aendy()` is called ----
test_that("deprecation Test 4: An error is thrown if `derive_var_aendy()`
          is called", {
  expect_error(
    derive_var_aendy(),
    class = "lifecycle_error_deprecated"
  )
})

## Test 5: An error is thrown if `derive_var_astdy()` is called ----
test_that("deprecation Test 5: An error is thrown if `derive_var_astdy()`
          is called", {
  expect_error(
    derive_var_astdy(),
    class = "lifecycle_error_deprecated"
  )
})

## Test 6: An error is thrown if `derive_var_atirel()` is called ----
test_that("deprecation Test 6: An error is thrown if `derive_var_atirel()`
          is called", {
  expect_error(
    derive_var_atirel(),
    class = "lifecycle_error_deprecated"
  )
})

## Test 7: An error is thrown if `derive_vars_suppqual()` is called ----
test_that("deprecation Test 7: An error is thrown if `derive_vars_suppqual()`
          is called", {
  expect_error(
    derive_vars_suppqual(),
    class = "lifecycle_error_deprecated"
  )
})

## Test 8: A warning is issued if `derive_derived_param()` is called ----
test_that("deprecation Test 8: A warning is issued if `derive_derived_param()`
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

  expect_warning(
    derive_derived_param(
      input,
      parameters = c("SYSBP", "DIABP"),
      by_vars = vars(USUBJID, VISIT),
      analysis_value = (AVAL.SYSBP + 2 * AVAL.DIABP) / 3,
      set_values_to = vars(
        PARAMCD = "MAP",
        PARAM = "Mean arterial pressure (mmHg)",
        AVALU = "mmHg"
      )
    ),
    class = "lifecycle_warning_deprecated"
  )
})

## Test 9: derive_vars_merged_dt: a deprecation warning is issued ----
test_that("deprecation Test 9: derive_vars_merged_dt: a deprecation warning
          is issued", {
  expect_warning(
    derive_vars_merged_dt(
      adsl,
      dataset_add = ex,
      order = vars(TRTSDT),
      flag_imputation = "date",
      by_vars = vars(STUDYID, USUBJID),
      dtc = EXSTDTC,
      new_vars_prefix = "TRTS",
      mode = "first"
    ),
    class = "lifecycle_warning_deprecated"
  )
})

## Test 10: derive_vars_merged_dtm: a deprecation warning is issued ----
test_that("deprecation Test 10: derive_vars_merged_dtm: a deprecation warning
          is issued", {
  expect_warning(
    derive_vars_merged_dtm(
      adsl,
      dataset_add = ex,
      order = vars(TRTSDTM),
      by_vars = vars(STUDYID, USUBJID),
      dtc = EXSTDTC,
      new_vars_prefix = "TRTS",
      time_imputation = "first",
      mode = "first"
    ),
    class = "lifecycle_warning_deprecated"
  )
})

## Test 11: date_source: errors when date_imputation is specified ----
test_that("deprecation Test 11: date_source: errors when date_imputation
          is specified", {
  expect_error(
    date_source(
      dataset_name = "ae",
      date = ASTDTM,
      date_imputation = "first"
    ),
    class = "lifecycle_error_deprecated"
  )
})

## Test 12: date_source: errors when time_imputation is specified ----
test_that("deprecation Test 12: date_source: errors when time_imputation
          is specified", {
  expect_error(
    date_source(
      dataset_name = "ae",
      date = ASTDTM,
      time_imputation = "first"
    ),
    class = "lifecycle_error_deprecated"
  )
})

## Test 13: date_source: errors when preserve is specified ----
test_that("deprecation Test 13: date_source: errors when preserve
          is specified", {
  expect_error(
    date_source(
      dataset_name = "ae",
      date = ASTDTM,
      preserve = TRUE
    ),
    class = "lifecycle_error_deprecated"
  )
})

## Test 14: A warning is issued if `derive_var_agegr_ema()` is called ----
test_that("deprecation Test 14: A warning is issued if `derive_var_agegr_ema()`
          is called", {
  rlang::with_options(lifecycle_verbosity = "warning", {
    expect_warning(
      derive_var_agegr_ema(admiral.test::admiral_dm, age_var = AGE, new_var = AGEGR1),
      class = "lifecycle_warning_deprecated"
    )
  })
})

## Test 15: A warning is issued if `derive_var_agegr_fda()` is called ----
test_that("deprecation Test 15: A warning is issued if `derive_var_agegr_fda()`
          is called", {
  rlang::with_options(lifecycle_verbosity = "warning", {
    expect_warning(
      derive_var_agegr_fda(admiral.test::admiral_dm, age_var = AGE, new_var = AGEGR1),
      class = "lifecycle_warning_deprecated"
    )
  })
})

## Test 16: A warning is issued if `derive_param_first_event()` is called ----
test_that("deprecation Test 16: A warning is issued if `derive_param_first_event()`
          is called", {
  rlang::with_options(lifecycle_verbosity = "warning", {
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

    expect_warning(
      derive_param_first_event(
        adrs,
        dataset_adsl = adsl,
        dataset_source = adrs,
        filter_source = PARAMCD == "OVR" & AVALC == "PD",
        date_var = ADT,
        set_values_to = vars(
          PARAMCD = "PD",
          ANL01FL = "Y"
        )
      ),
      class = "lifecycle_warning_deprecated"
    )
  })
})

## Test 17: An error is thrown if `smq_select()` is called ----
test_that("deprecation Test 17: An error is thrown if `smq_select()`
          is called", {
  expect_error(
    smq_select(),
    class = "lifecycle_error_deprecated"
  )
})

## Test 18: An error is thrown if `sdg_select()` is called ----
test_that("deprecation Test 18: An error is thrown if `sdg_select()`
          is called", {
  expect_error(
    sdg_select(),
    class = "lifecycle_error_deprecated"
  )
})

## Test 19: An error is thrown if `create_query_data()` and `meddra_version` argument is called ----
test_that("deprecation Test 19: An error is thrown if `create_query_data()`
           with `meddra_version` argument is called", {
  expect_error(
    create_query_data(
      meddra_version = "20.1"
    ),
    class = "lifecycle_error_deprecated"
  )
})

## Test 20: An error is thrown if `create_query_data()` and `whodd_version` argument is called ----
test_that("deprecation Test 20: An error is thrown if `create_query_data()`
           with `whodd_version` argument is called", {
  expect_error(
    create_query_data(
      whodd_version = "2019-09"
    ),
    class = "lifecycle_error_deprecated"
  )
})

## Test 21: An error is thrown if `create_query_data()` and `get_smq_fun` argument is called ----
test_that("deprecation Test 21: An error is thrown if `create_query_data()`
           with `get_smq_fun` argument is called", {
  expect_error(
    create_query_data(
      get_smq_fun = admiral.test:::get_smq_terms
    ),
    class = "lifecycle_error_deprecated"
  )
})

## Test 22: An error is thrown if `create_query_data()` and `get_sdg_fun` argument is called ----
test_that("deprecation Test 22: An error is thrown if `create_query_data()`
           with `get_sdg_fun` argument is called", {
  expect_error(
    create_query_data(
      get_sdg_fun = admiral.test:::get_sdg_terms
    ),
    class = "lifecycle_error_deprecated"
  )
})
