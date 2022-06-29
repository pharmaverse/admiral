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
})

test_that("a warning is issued when using `derive_var_worst_flag()` with `filter` argument", {
  input_worst_flag <- tibble::tribble(
    ~STUDYID, ~USUBJID, ~PARAMCD, ~AVISIT, ~ADT, ~AVAL,
    "TEST01", "PAT01", "PARAM01", "BASELINE", as.Date("2021-04-27"), 15.0,
    "TEST01", "PAT01", "PARAM01", "BASELINE", as.Date("2021-04-25"), 14.0,
    "TEST01", "PAT01", "PARAM01", "BASELINE", as.Date("2021-04-23"), 15.0,
    "TEST01", "PAT01", "PARAM01", "WEEK 1", as.Date("2021-04-27"), 10.0,
    "TEST01", "PAT01", "PARAM01", "WEEK 2", as.Date("2021-04-30"), 12.0,
    "TEST01", "PAT02", "PARAM01", "SCREENING", as.Date("2021-04-27"), 15.0,
    "TEST01", "PAT02", "PARAM01", "BASELINE", as.Date("2021-04-25"), 14.0,
    "TEST01", "PAT02", "PARAM01", "BASELINE", as.Date("2021-04-23"), 15.0,
    "TEST01", "PAT02", "PARAM01", "WEEK 1", as.Date("2021-04-27"), 10.0,
    "TEST01", "PAT02", "PARAM01", "WEEK 2", as.Date("2021-04-30"), 12.0,
    "TEST01", "PAT01", "PARAM02", "SCREENING", as.Date("2021-04-27"), 15.0,
    "TEST01", "PAT01", "PARAM02", "SCREENING", as.Date("2021-04-25"), 14.0,
    "TEST01", "PAT01", "PARAM02", "SCREENING", as.Date("2021-04-23"), 15.0,
    "TEST01", "PAT01", "PARAM02", "BASELINE", as.Date("2021-04-27"), 10.0,
    "TEST01", "PAT01", "PARAM02", "WEEK 2", as.Date("2021-04-30"), 12.0,
    "TEST01", "PAT02", "PARAM02", "SCREENING", as.Date("2021-04-27"), 15.0,
    "TEST01", "PAT02", "PARAM02", "BASELINE", as.Date("2021-04-25"), 14.0,
    "TEST01", "PAT02", "PARAM02", "WEEK 1", as.Date("2021-04-23"), 15.0,
    "TEST01", "PAT02", "PARAM02", "WEEK 1", as.Date("2021-04-27"), 10.0,
    "TEST01", "PAT02", "PARAM02", "BASELINE", as.Date("2021-04-30"), 12.0,
    "TEST01", "PAT02", "PARAM03", "SCREENING", as.Date("2021-04-27"), 15.0,
    "TEST01", "PAT02", "PARAM03", "BASELINE", as.Date("2021-04-25"), 14.0,
    "TEST01", "PAT02", "PARAM03", "WEEK 1", as.Date("2021-04-23"), 15.0,
    "TEST01", "PAT02", "PARAM03", "WEEK 1", as.Date("2021-04-27"), 10.0,
    "TEST01", "PAT02", "PARAM03", "BASELINE", as.Date("2021-04-30"), 12.0
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
})

test_that("derive_var_ady() Test 1: A warning is issued when using `derive_var_ady()`", {
  input <- tibble::tribble(
    ~TRTSDT, ~ADT,
    ymd("2020-01-01"), ymd("2020-02-24"),
    ymd("2020-01-01"), ymd("2020-01-01"),
    ymd("2020-02-24"), ymd("2020-01-01")
  )

  expect_warning(
    derive_var_ady(input),
    "deprecated",
    fixed = TRUE
  )
})

test_that("derive_var_aendy Test 1: A warning is issued when using `derive_var_aendy()`", {
  input <- tibble::tribble(
    ~TRTSDT, ~AENDT,
    ymd("2020-01-01"), ymd("2020-02-24"),
    ymd("2020-01-01"), ymd("2020-01-01"),
    ymd("2020-02-24"), ymd("2020-01-01")
  )

  expect_warning(
    derive_var_aendy(input),
    "deprecated",
    fixed = TRUE
  )
})

test_that("derive_var_astdy Test 1: A warning is issued when using `derive_var_astdy()`", {
  input <- tibble::tribble(
    ~TRTSDT, ~ASTDT,
    ymd("2020-01-01"), ymd("2020-02-24"),
    ymd("2020-01-01"), ymd("2020-01-01"),
    ymd("2020-02-24"), ymd("2020-01-01")
  )

  expect_warning(
    derive_var_astdy(input),
    "deprecated",
    fixed = TRUE
  )
})

test_that("derive_var_trtedtm Test 1: A warning is issued when using `derive_var_trtedtm()`", {
  adsl <- tibble::tibble(STUDYID = "STUDY", USUBJID = 1:3)
  ex <- tibble::tribble(
    ~USUBJID, ~EXENDTC, ~EXSEQ, ~EXDOSE, ~EXTRT,
    1L, "2020-01-01", 1, 12, "ACTIVE",
    1L, "2020-02-03", 2, 9, "ACTIVE",
    2L, "2020-01-02", 1, 0, "PLACEBO",
    3L, "2020-03-13", 1, 14, "ACTIVE",
    3L, "2020-03-21", 2, 0, "ACTIVE"
  )

  expect_warning(derive_var_trtedtm(
    adsl,
    dataset_ex = ex,
    subject_keys = vars(USUBJID)
  ),
  "deprecated",
  fixed = TRUE
  )
})

test_that("derive_var_trtsdtm Test 1: A warning is issued when using `derive_var_trtsdtm()`", {
  adsl <- tibble::tibble(STUDYID = "STUDY", USUBJID = 1:3)
  ex <- tibble::tribble(
    ~USUBJID, ~EXSTDTC, ~EXSEQ, ~EXDOSE, ~EXTRT,
    1L, "2020-01-01", 1, 12, "ACTIVE",
    1L, "2020-02-03", 2, 9, "ACTIVE",
    2L, "2020-01-02", 1, 0, "PLACEBO",
    3L, "2020-03-13", 1, 14, "ACTIVE",
    3L, "2020-03-21", 2, 0, "ACTIVE"
  )

  expect_warning(derive_var_trtsdtm(
    adsl,
    dataset_ex = ex,
    subject_keys = vars(USUBJID)
  ),
  "deprecated",
  fixed = TRUE
  )
})

test_that("derive_var_disposition Test 1: A warning is issued when using `derive_var_disposition_dt()`", { # nolint
  adsl <- tibble::tribble(
    ~STUDYID, ~USUBJID,
    "TEST01", "PAT01",
    "TEST01", "PAT02"
  )

  ds <- tibble::tribble(
    ~STUDYID, ~USUBJID, ~DSCAT, ~DSDECOD, ~DSSTDTC,
    "TEST01", "PAT01", "PROTOCOL MILESTONE", "INFORMED CONSENT OBTAINED", "2021-04-01",
    "TEST01", "PAT01", "PROTOCOL MILESTONE", "RANDOMIZATION", "2021-04-11",
    "TEST01", "PAT01", "DISPOSITION EVENT", "ADVERSE EVENT", "2021-12-01",
    "TEST01", "PAT01", "OTHER EVENT", "DEATH", "2022-02-01",
    "TEST01", "PAT02", "PROTOCOL MILESTONE", "INFORMED CONSENT OBTAINED", "2021-04-02",
    "TEST01", "PAT02", "PROTOCOL MILESTONE", "RANDOMIZATION", "2021-04-11",
    "TEST01", "PAT02", "DISPOSITION EVENT", "COMPLETED", "2021-12-01",
    "TEST01", "PAT02", "OTHER EVENT", "DEATH", "2022-04"
  )

  expect_warning(
    derive_var_disposition_dt(
      adsl,
      dataset_ds = ds,
      new_var = RFICDT,
      dtc = DSSTDTC,
      filter = DSCAT == "PROTOCOL MILESTONE" &
        DSDECOD == "INFORMED CONSENT OBTAINED",
      date_imputation = NULL
    ),
    "deprecated",
    fixed = TRUE
  )
})

test_that("derive_var_atirel Test 1: A warning is issued when using `derive_var_atirel()`", {
  input <- tibble::tribble(
    ~STUDYID, ~USUBJID, ~TRTSDTM, ~ASTDTM, ~AENDTM, ~ASTTMF, # nolint
    "TEST01", "PAT01", "2012-02-25 23:00:00", "2012-03-28 19:00:00", "2012-05-25 23:00:00", "",
    "TEST01", "PAT01", "", "2012-02-28 19:00:00", "", "",
    "TEST01", "PAT01", "2017-02-25 23:00:00", "2013-02-25 19:00:00", "2014-02-25 19:00:00", "",
    "TEST01", "PAT01", "2017-02-25 23:00:00", "2017-02-25 14:00:00", "2017-03-25 23:00:00", "m",
    "TEST01", "PAT01", "2017-02-25 23:00:00", "2017-01-25 14:00:00", "2018-04-29 14:00:00", "",
  ) %>% mutate(
    TRTSDTM = lubridate::as_datetime(TRTSDTM),
    ASTDTM = lubridate::as_datetime(ASTDTM),
    AENDTM = lubridate::as_datetime(AENDTM)
  )

  expect_warning(
    derive_var_atirel(
      dataset = input,
      flag_var = ASTTMF,
      new_var = ATIREL
    ),
    "deprecated",
    fixed = TRUE
  )
})

test_that("derive_vars_suppqual Test 1: An error is thrown if `derive_vars_suppqual()` is called", {
  expect_error(
    derive_vars_suppqual(),
    "deprecated",
    fixed = TRUE,
    class = "lifecycle_error_deprecated"
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
