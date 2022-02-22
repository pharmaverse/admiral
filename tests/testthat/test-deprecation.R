test_that("a warning is issued when specifying `dthcaus_source(dataset = )", {
  expect_warning(
    dthcaus_source(
      dataset = ae,
      filter = AEOUT == "FATAL",
      date = AEDTHDTC,
      mode = "first",
      dthcaus = AEDECOD
    ),
    "deprecated",
    fixed = TRUE
  )
})

test_that("a warning is issued when specifying `lstalvdt_source(dataset = )", {
  expect_warning(
    lstalvdt_source(
      dataset = lb,
      date = LBDTC,
      filter = nchar(LBDTC) >= 10
    ),
    "deprecated",
    fixed = TRUE
  )
})

test_that("a warning is issued when using `derive_var_basec()", {
  dataset <- tibble::tribble(
    ~STUDYID, ~USUBJID, ~PARAMCD,  ~AVAL, ~AVALC,   ~AVISIT,    ~ABLFL,
    "TEST01", "PAT01",  "PARAM03", NA,    "LOW",    "Baseline", "Y",
    "TEST01", "PAT01",  "PARAM03", NA,    "LOW",    "Day 7",    "N",
    "TEST01", "PAT01",  "PARAM03", NA,    "MEDIUM", "Day 14",   "N",
    "TEST01", "PAT01",  "PARAM04", NA,    "HIGH",   "Baseline", "Y",
    "TEST01", "PAT01",  "PARAM04", NA,    "HIGH",   "Day 7",    "N",
    "TEST01", "PAT01",  "PARAM04", NA,    "MEDIUM", "Day 14",   "N"
  )

  expect_warning(
    derive_var_basec(dataset, by_vars = vars(STUDYID, USUBJID, PARAMCD)),
    "deprecated",
    fixed = TRUE
  )
})

test_that("a warning is issued when using `derive_baseline()", {
  dataset <- tibble::tribble(
    ~STUDYID, ~USUBJID, ~PARAMCD,  ~AVAL,  ~AVALC, ~AVISIT,    ~ABLFL,
    "TEST01", "PAT01",  "PARAM01", 10.12,  NA,     "Baseline", "Y",
    "TEST01", "PAT01",  "PARAM01",  9.7,   NA,     "Day 7",    "N",
    "TEST01", "PAT01",  "PARAM01", 15.01,  NA,     "Day 14",   "N",
    "TEST01", "PAT01",  "PARAM02",  8.35,  NA,     "Baseline", "Y",
    "TEST01", "PAT01",  "PARAM02", NA,     NA,     "Day 7",    "N",
    "TEST01", "PAT01",  "PARAM02",  8.35,  NA,     "Day 14",   "N"
  )

  expect_warning(
    derive_baseline(
      dataset,
      by_vars = vars(STUDYID, USUBJID, PARAMCD),
      source_var = AVAL,
      new_var = BASE
    ),
    "deprecated",
    fixed = TRUE
  )
})


test_that("a warning is issued when using `derive_disposition_dt()`", {
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
    adsl %>%
      derive_disposition_dt(
        dataset_ds = ds,
        new_var = RFICDT,
        dtc = DSSTDTC,
        filter = DSCAT == "PROTOCOL MILESTONE" & DSDECOD == "INFORMED CONSENT OBTAINED",
        date_imputation = NULL),
    "deprecated",
    fixed = TRUE
  )
})


test_that("a warning is issued when using `derive_disposition_status()`", {
  dm <- tibble::tribble(
    ~STUDYID, ~USUBJID,
    "TEST01", "PAT01",
    "TEST01", "PAT02",
    "TEST01", "PAT03",
    "TEST01", "PAT04"
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
    "TEST01", "PAT02", "OTHER EVENT", "DEATH", "2022-04",
    "TEST01", "PAT03", "PROTOCOL MILESTONE", "INFORMED CONSENT OBTAINED", "2021-04-02",
    "TEST01", "PAT03", "PROTOCOL MILESTONE", "RANDOMIZATION", "2021-04-11",
    "TEST01", "PAT03", "DISPOSITION EVENT", "PROGRESSIVE DISEASE", "2021-05-01",
    "TEST01", "PAT03", "OTHER EVENT", "DEATH", "2022-04",
    "TEST01", "PAT04", "PROTOCOL MILESTONE", "INFORMED CONSENT OBTAINED", "2021-04-02",
    "TEST01", "PAT04", "PROTOCOL MILESTONE", "RANDOMIZATION", "2021-04-11")

  expect_warning(
    derive_disposition_status(
      dataset = dm,
      dataset_ds = ds,
      new_var = EOSSTT,
      status_var = DSDECOD,
      filter_ds = DSCAT == "DISPOSITION EVENT"
    ),
    "deprecated",
    fixed = TRUE
  )
})

test_that("a warning is issued when using `derive_extreme_flag()`", {
  input <- tibble::tribble(
    ~USUBJID, ~AVISITN, ~AVAL,
    1, 1, 12,
    1, 3, 9,
    2, 2, 42,
    3, 3, 14,
    3, 3, 10
  )

  expect_warning(
    derive_extreme_flag(
      input,
      by_vars = vars(USUBJID),
      order = vars(AVISITN, desc(AVAL)),
      new_var = firstfl,
      mode = "first"
    ),
    "deprecated",
    fixed = TRUE
  )
})

test_that("a warning is issued when using `derive_obs_number()`", {
  input <- tibble::tribble(
    ~USUBJID, ~AVISITN, ~AVAL,
    1, 1, 12,
    1, 3, 9,
    2, 2, 42,
    3, 3, 14,
    3, 3, 10
  )

  expect_warning(
    derive_obs_number(
      input,
      by_vars = vars(USUBJID),
      order = vars(AVISITN, AVAL),
    ),
    "deprecated",
    fixed = TRUE
  )
})

test_that("a warning is issued when using `derive_last_dose()`", {
  input_ae <- tibble::tribble(
    ~STUDYID, ~USUBJID, ~AESEQ, ~AESTDTC,
    "my_study", "subject1", 1, "2020-01-02",
    "my_study", "subject1", 2, "2020-08-31",
    "my_study", "subject1", 3, "2020-10-10",
    "my_study", "subject2", 1, "2019-05-15",
    "my_study", "subject2", 2, "2020-02-20",
    "my_study", "subject3", 1, "2020-03-02",
    "my_study", "subject4", 1, "2020-11-02"
  )

  input_ex <- tibble::tribble(
    ~STUDYID,   ~USUBJID,   ~EXSTDTC,     ~EXENDTC,    ~EXSEQ, ~EXDOSE, ~EXTRT,
    "my_study", "subject1", "2020-01-01", "2020-01-01", 1,     10,      "treatment",
    "my_study", "subject1", "2020-08-29", "2020-08-29", 2,     10,      "treatment",
    "my_study", "subject1", "2020-09-02", "2020-09-02", 3,     10,      "treatment",
    "my_study", "subject1", "2020-10-20", "2020-10-20", 4,     10,      "treatment",
    "my_study", "subject2", "2019-05-25", "2019-05-25", 1,      0,      "placebo",
    "my_study", "subject2", "2020-01-20", "2020-01-20", 2,      0,      "placebo",
    "my_study", "subject3", "2020-03-15", "2020-03-15", 1,     10,      "treatment"
  ) %>%
    mutate(EXSTDTC = as.Date(EXSTDTC), EXENDTC = as.Date(EXENDTC))

  expect_warning(
    derive_last_dose(
      input_ae,
      input_ex,
      filter_ex = (EXDOSE > 0) | (EXDOSE == 0 & EXTRT == "placebo"),
      by_vars = vars(STUDYID, USUBJID),
      dose_start = EXSTDTC,
      dose_end = EXENDTC,
      analysis_date = AESTDTC,
      dataset_seq_var = AESEQ,
      new_var = LDOSEDTM,
      output_datetime = TRUE,
      check_dates_only = FALSE,
      traceability_vars = NULL
    ),
    "deprecated",
    fixed = TRUE
  )
})

test_that("a warning is issued when using `derive_disposition_reason()`", {
  dm <- tibble::tribble(
    ~STUDYID, ~USUBJID,
    "TEST01", "PAT01",
    "TEST01", "PAT02",
    "TEST01", "PAT03",
    "TEST01", "PAT04"
  )

  ds <- tibble::tribble(
    ~STUDYID, ~USUBJID, ~DSCAT, ~DSDECOD, ~DSTERM, ~DSSTDTC,
    "TEST01", "PAT01", "PROTOCOL MILESTONE", "INFORMED CONSENT OBTAINED", "INFORMED CONSENT OBTAINED", "2021-04-01", # nolint
    "TEST01", "PAT01", "PROTOCOL MILESTONE", "RANDOMIZATION", "RANDOMIZED", "2021-04-11",
    "TEST01", "PAT01", "DISPOSITION EVENT", "ADVERSE EVENT", "ADVERSE EVENT", "2021-12-01",
    "TEST01", "PAT01", "OTHER EVENT", "DEATH", "DEATH", "2022-02-01",
    "TEST01", "PAT02", "PROTOCOL MILESTONE", "INFORMED CONSENT OBTAINED", "INFORMED CONSENT OBTAINED", "2021-04-02", # nolint
    "TEST01", "PAT02", "PROTOCOL MILESTONE", "RANDOMIZATION", "RANDOMIZED", "2021-04-11",
    "TEST01", "PAT02", "DISPOSITION EVENT", "COMPLETED", NA_character_, "2021-12-01",
    "TEST01", "PAT02", "OTHER EVENT", "DEATH", "DEATH", "2022-04",
    "TEST01", "PAT03", "PROTOCOL MILESTONE", "INFORMED CONSENT OBTAINED", "INFORMED CONSENT OBTAINED", "2021-04-02", # nolint
    "TEST01", "PAT03", "PROTOCOL MILESTONE", "RANDOMIZATION", "RANDOMIZED", "2021-04-11",
    "TEST01", "PAT03", "DISPOSITION EVENT", "PROGRESSIVE DISEASE", "DISEASE PROGRESSION", "2021-05-01", # nolint
    "TEST01", "PAT03", "OTHER EVENT", "DEATH", "DEATH", "2022-04",
    "TEST01", "PAT04", "PROTOCOL MILESTONE", "INFORMED CONSENT OBTAINED", "INFORMED CONSENT OBTAINED", "2021-04-02", # nolint
    "TEST01", "PAT04", "PROTOCOL MILESTONE", "RANDOMIZATION", "RANDOMIZED", "2021-04-11"
  )

  expect_warning(
    derive_disposition_reason(
      dataset = dm,
      dataset_ds = ds,
      new_var = DCSREAS,
      reason_var = DSDECOD,
      filter_ds = DSCAT == "DISPOSITION EVENT"),
    "deprecated",
    fixed = TRUE
  )
})

test_that("a warning is issued when using `derive_params_exposure()", {
  dataset <- tibble::tribble(
  ~USUBJID,      ~VISIT,     ~PARAMCD, ~AVAL, ~AVALC, ~EXSTDTC,     ~EXENDTC,
  "01-701-1015", "BASELINE", "DOSE",   80,    NA,     "2020-07-01", "2020-07-14",
  "01-701-1015", "WEEK 2",   "DOSE",   80,    NA,     "2020-07-15", "2020-09-23",
  "01-701-1015", "WEEK 12",  "DOSE",   65,    NA,     "2020-09-24", "2020-12-16",
  "01-701-1015", "WEEK 24",  "DOSE",   65,    NA,     "2020-12-17", "2021-06-02",
  "01-701-1015", "BASELINE", "ADJ",    NA,    NA,     "2020-07-01", "2020-07-14",
  "01-701-1015", "WEEK 2",   "ADJ",    NA,    "Y",    "2020-07-15", "2020-09-23",
  "01-701-1015", "WEEK 12",  "ADJ",    NA,    "Y",    "2020-09-24", "2020-12-16",
  "01-701-1015", "WEEK 24",  "ADJ",    NA,    NA,     "2020-12-17", "2021-06-02",
  "01-701-1281", "BASELINE", "DOSE",   80,    NA,     "2020-07-03", "2020-07-18",
  "01-701-1281", "WEEK 2",   "DOSE",   80,    NA,     "2020-07-19", "2020-10-01",
  "01-701-1281", "WEEK 12",  "DOSE",   82,    NA,     "2020-10-02", "2020-12-01",
  "01-701-1281", "BASELINE", "ADJ",    NA,    NA,     "2020-07-03", "2020-07-18",
  "01-701-1281", "WEEK 2",   "ADJ",    NA,    NA,     "2020-07-19", "2020-10-01",
  "01-701-1281", "WEEK 12",  "ADJ",    NA,    NA,     "2020-10-02", "2020-12-01"
  ) %>%
  mutate(
    ASTDTM = ymd_hms(paste(EXSTDTC, "T00:00:00")),
    ASTDT = date(ASTDTM),
    AENDTM = ymd_hms(paste(EXENDTC, "T00:00:00")),
    AENDT = date(AENDTM)
  )

  expect_warning(
    derive_params_exposure(
      dataset,
      by_vars = vars(USUBJID),
      input_code = "DOSE",
      analysis_var = AVAL,
      summary_fun = function(x) sum(x, na.rm = TRUE),
      filter = NULL,
      set_values_to = vars(PARAMCD = "TDOSE", PARCAT1 = "OVERALL")

    ),
    "deprecated",
    fixed = TRUE
  )
})

test_that("a warning is issued when using `derive_agegr_fda()`", {
  expect_warning(
    derive_agegr_fda(adsl, age_var = AGE, new_var = AGEGR2),
    "deprecated",
    fixed = TRUE
  )
})

test_that("a warning is issued when using `derive_agegr_ema()`", {
  expect_warning(
    derive_agegr_ema(adsl, age_var = AGE, new_var = AGEGR2),
    "deprecated",
    fixed = TRUE
  )
})

test_that("a warning is issued when using `derive_worst_flag()`", {
  input <- tibble::tribble(
    ~STUDYID, ~USUBJID, ~PARAMCD,  ~AVISIT,    ~ADT,                 ~AVAL,
    "TEST01", "PAT01",  "PARAM01", "BASELINE", as.Date("2021-04-27"), 15.0,
    "TEST01", "PAT01",  "PARAM01", "BASELINE", as.Date("2021-04-25"), 14.0,
    "TEST01", "PAT01",  "PARAM01", "BASELINE", as.Date("2021-04-23"), 15.0,
    "TEST01", "PAT01",  "PARAM01", "WEEK 1",   as.Date("2021-04-27"), 10.0,
    "TEST01", "PAT01",  "PARAM01", "WEEK 2",   as.Date("2021-04-30"), 12.0,

    "TEST01", "PAT02",  "PARAM01", "SCREENING",as.Date("2021-04-27"), 15.0,
    "TEST01", "PAT02",  "PARAM01", "BASELINE", as.Date("2021-04-25"), 14.0,
    "TEST01", "PAT02",  "PARAM01", "BASELINE", as.Date("2021-04-23"), 15.0,
    "TEST01", "PAT02",  "PARAM01", "WEEK 1",   as.Date("2021-04-27"), 10.0,
    "TEST01", "PAT02",  "PARAM01", "WEEK 2",   as.Date("2021-04-30"), 12.0,

    "TEST01", "PAT01",  "PARAM02", "SCREENING",as.Date("2021-04-27"), 15.0,
    "TEST01", "PAT01",  "PARAM02", "SCREENING",as.Date("2021-04-25"), 14.0,
    "TEST01", "PAT01",  "PARAM02", "SCREENING",as.Date("2021-04-23"), 15.0,
    "TEST01", "PAT01",  "PARAM02", "BASELINE", as.Date("2021-04-27"), 10.0,
    "TEST01", "PAT01",  "PARAM02", "WEEK 2",   as.Date("2021-04-30"), 12.0,

    "TEST01", "PAT02",  "PARAM02", "SCREENING",as.Date("2021-04-27"), 15.0,
    "TEST01", "PAT02",  "PARAM02", "BASELINE", as.Date("2021-04-25"), 14.0,
    "TEST01", "PAT02",  "PARAM02", "WEEK 1",   as.Date("2021-04-23"), 15.0,
    "TEST01", "PAT02",  "PARAM02", "WEEK 1",   as.Date("2021-04-27"), 10.0,
    "TEST01", "PAT02",  "PARAM02", "BASELINE", as.Date("2021-04-30"), 12.0,

    "TEST01", "PAT02",  "PARAM03", "SCREENING",as.Date("2021-04-27"), 15.0,
    "TEST01", "PAT02",  "PARAM03", "BASELINE", as.Date("2021-04-25"), 14.0,
    "TEST01", "PAT02",  "PARAM03", "WEEK 1",   as.Date("2021-04-23"), 15.0,
    "TEST01", "PAT02",  "PARAM03", "WEEK 1",   as.Date("2021-04-27"), 10.0,
    "TEST01", "PAT02",  "PARAM03", "BASELINE", as.Date("2021-04-30"), 12.0
  )

  expect_warning(
    derive_worst_flag(
      input,
      by_vars = vars(USUBJID, PARAMCD, AVISIT),
      order = vars(desc(ADT)),
      new_var = WORSTFL,
      param_var = PARAMCD,
      analysis_var = AVAL,
      worst_high = c("PARAM01", "PARAM03"),
      worst_low = "PARAM02"
    ),
    "deprecated",
    fixed = TRUE
  )
})
