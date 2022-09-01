library(admiral.test)
library(rlang)
library(tibble)
data("admiral_dm")

adsl <- tribble(
  ~USUBJID, ~SEX, ~COUNTRY,
  "ST42-1", "F",  "AUT",
  "ST42-2", "M",  "MWI",
  "ST42-3", "M",  "NOR",
  "ST42-4", "F",  "UGA"
) %>% mutate(STUDYID = "ST42")

ex <- tribble(
  ~USUBJID, ~EXSTDTC,
  "ST42-1", "2020-12-07",
  "ST42-1", "2020-12-14",
  "ST42-2", "2021-01-12T12:00:00",
  "ST42-2", "2021-01-26T13:21",
  "ST42-3", "2021-03-02"
) %>% mutate(STUDYID = "ST42")

test_that("An error is thrown if `derive_var_extreme_flag()` with `filter` argument is called", {
  expect_error(
    derive_var_extreme_flag(
      filter = !is.na(AVAL)
    ),
    "deprecated",
    fixed = TRUE,
    class = "lifecycle_error_deprecated"
  )
})

test_that("An error is thrown if `derive_var_worst_flag()` with `filter` argument is called", {
  expect_error(
    derive_var_worst_flag(
      filter = !is.na(AVAL)
    ),
    "deprecated",
    fixed = TRUE,
    class = "lifecycle_error_deprecated"
  )
})

test_that("derive_var_ady() Test 1: An error is thrown if `derive_var_ady()` is called", {
  expect_error(
    derive_var_ady(),
    "deprecated",
    fixed = TRUE,
    class = "lifecycle_error_deprecated"
  )
})

test_that("derive_var_aendy Test 1: An error is thrown if `derive_var_aendy()` is called", {
  expect_error(
    derive_var_aendy(),
    "deprecated",
    fixed = TRUE,
    class = "lifecycle_error_deprecated"
  )
})

test_that("derive_var_astdy Test 1: An error is thrown if `derive_var_astdy()` is called", {
  expect_error(
    derive_var_astdy(),
    "deprecated",
    fixed = TRUE,
    class = "lifecycle_error_deprecated"
  )
})

test_that("derive_var_atirel Test 1: An error is thrown if `derive_var_atirel()` is called", {
  expect_error(
    derive_var_atirel(),
    "deprecated",
    fixed = TRUE,
    class = "lifecycle_error_deprecated"
  )
})

test_that("derive_vars_suppqual Test 1: An error is thrown if `derive_vars_suppqual()` is called", {
  expect_error(
    derive_vars_suppqual(),
    "deprecated",
    fixed = TRUE
  )
})

test_that("derive_derived_param Test 1: A warning is issued if `derive_derived_param()` is called", { # nolint
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
    "deprecated",
    fixed = TRUE
  )
})

test_that("derive_vars_merged_dt: a deprecation warning is issued", {
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
    "deprecated"
  )
})

test_that("derive_vars_merged_dtm: a deprecation warning is issued", {
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
    "deprecated"
  )
})

test_that("date_source: errors when date_imputation is specified", {
  expect_error(
    date_source(
      dataset_name = "ae",
      date = ASTDTM,
      date_imputation = "first"
    ),
    paste(
      "The `date_imputation` argument of `date_source\\(\\)` .* deprecated .* admiral 0.8.0.*",
      "Please use `derive_vars_dtm\\(\\)` to convert DTC variables to datetime variables in the dataset.", # nolint
      sep = "\n"
    )
  )
})

test_that("date_source: errors when time_imputation is specified", {
  expect_error(
    date_source(
      dataset_name = "ae",
      date = ASTDTM,
      time_imputation = "first"
    ),
    paste(
      "The `time_imputation` argument of `date_source\\(\\)` .* deprecated .* admiral 0.8.0.*",
      "Please use `derive_vars_dtm\\(\\)` to convert DTC variables to datetime variables in the dataset.", # nolint
      sep = "\n"
    )
  )
})

test_that("date_source: errors when preserve is specified", {
  expect_error(
    date_source(
      dataset_name = "ae",
      date = ASTDTM,
      preserve = TRUE
    ),
    paste(
      "The `preserve` argument of `date_source\\(\\)` .* deprecated .* admiral 0.8.0.*",
      "Please use `derive_vars_dtm\\(\\)` to convert DTC variables to datetime variables in the dataset.", # nolint
      sep = "\n"
    )
  )
})

test_that("derive_var_agegr_ema Test 1: A warning is issued if `derive_var_agegr_ema()` is called", { # nolint
  with_options(lifecycle_verbosity = "warning", {
    expect_warning(
      derive_var_agegr_ema(admiral_dm, age_var = AGE, new_var = AGEGR1),
      "deprecated",
      fixed = TRUE
    )
  })
})

test_that("derive_var_agegr_fda Test 1: A warning is issued if `derive_var_agegr_fda()` is called", { # nolint
  with_options(lifecycle_verbosity = "warning", {
    expect_warning(
      derive_var_agegr_fda(admiral_dm, age_var = AGE, new_var = AGEGR1),
      "deprecated",
      fixed = TRUE
    )
  })
})
