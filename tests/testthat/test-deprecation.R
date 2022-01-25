library(admiral.test)

test_that("a warning is issued when specifying `derive_var_extreme_flag(flag_filter = )`", {
  data(advs)

  expect_warning(
    derive_var_extreme_flag(
      advs[1:100, ],
      by_vars = vars(USUBJID, PARAMCD),
      order = vars(ADT),
      new_var = ABLFL,
      mode = "last",
      flag_filter = AVISIT == "BASELINE"
    ),
    "deprecated",
    fixed = TRUE
  )
})

test_that("a warning is issued when specifying `dthcaus_source(date_var = )", {
  expect_warning(
    dthcaus_source(
      dataset_name = "ae",
      filter = AEOUT == "FATAL",
      date_var = AEDTHDTC,
      mode = "first",
      dthcaus = AEDECOD
    ),
    "deprecated",
    fixed = TRUE
  )
})

test_that("a warning is issued when specifying `dthcaus_source(traceabilty_vars = )", {
  expect_warning(
    dthcaus_source(
      dataset_name = "ae",
      filter = AEOUT == "FATAL",
      date = AEDTHDTC,
      mode = "first",
      dthcaus = AEDECOD,
      traceabilty_vars = vars(DTHDOM = "AE", DTHSEQ = AESEQ)
    ),
    "deprecated",
    fixed = TRUE
  )
})

test_that("a warning is issued when specifying `lstalvdt_source(date_var = )", {
  expect_warning(
    lstalvdt_source(dataset_name = "adsl", date_var = TRTEDT),
    "deprecated",
    fixed = TRUE
  )
})

test_that("a warning is issued when using `derive_suppqual_vars()", {
  data(ae)
  data(suppae)

  expect_warning(
    derive_suppqual_vars(ae[1:100, ], suppae[1:100, ]),
    "deprecated",
    fixed = TRUE
  )
})

test_that("a warning is issued when using `derive_query_vars()", {
  data(queries)
  adae <- tibble::tribble(
    ~USUBJID, ~ASTDTM, ~AETERM, ~AESEQ, ~AEDECOD, ~AELLT, ~AELLTCD,
    "01", "2020-06-02 23:59:59", "ALANINE AMINOTRANSFERASE ABNORMAL",
      3, "Alanine aminotransferase abnormal", NA_character_, NA_integer_,
    "02", "2020-06-05 23:59:59", "BASEDOW'S DISEASE",
      5, "Basedow's disease", NA_character_, 1L,
    "03", "2020-06-07 23:59:59", "SOME TERM",
      2, "Some query", "Some term", NA_integer_,
    "05", "2020-06-09 23:59:59", "ALVEOLAR PROTEINOSIS",
      7, "Alveolar proteinosis", NA_character_,  NA_integer_
  )

  expect_warning(
    derive_query_vars(adae, queries),
    "deprecated",
    fixed = TRUE
  )
})

test_that("a warning is issued when using `derive_duration()", {
  adsl <- tibble::tribble(
    ~BRTHDT, ~RANDDT,
    lubridate::ymd("1984-09-06"), lubridate::ymd("2020-02-24")
  )

  expect_warning(
    derive_duration(
      adsl,
      new_var = AAGE,
      new_var_unit = AAGEU,
      start_date = BRTHDT,
      end_date = RANDDT,
      out_unit = "years",
      add_one = FALSE,
      trunc_out = TRUE
    ),
    "deprecated",
    fixed = TRUE
  )
})

test_that("a warning is issued when using `derive_aage()", {
  adsl <- tibble::tribble(
    ~BRTHDT, ~RANDDT,
    lubridate::ymd("1984-09-06"), lubridate::ymd("2020-02-24")
  )

  expect_warning(
    derive_aage(adsl),
    "deprecated",
    fixed = TRUE
  )
})

test_that("a warning is issued when specifying `derive_var_ontrtfl(date = )", {
  data(advs)

  expect_warning(
    derive_var_ontrtfl(
      advs[1:100, ],
      date = ADT,
      ref_start_date = TRTSDT,
      ref_end_date = TRTEDT
    ),
    "deprecated",
    fixed = TRUE
  )
})

test_that("a warning is issued when specifying `derive_summary_records(filter_rows = )", {
  data(advs)

  expect_warning(
    derive_summary_records(
      advs[1:100, ],
      by_vars = vars(USUBJID, PARAM, AVISIT),
      analysis_var = AVAL,
      summary_fun = function(x) mean(x, na.rm = TRUE),
      filter_rows = dplyr::n() > 2,
      set_values_to = vars(DTYPE = "AVERAGE")
    ),
    "deprecated",
    fixed = TRUE
  )
})

test_that("an error is thrown when specifying `derive_summary_records(fns = )", {
  data(advs)

  expect_error(
    derive_summary_records(
      advs[1:100, ],
      by_vars = vars(USUBJID, PARAM, AVISIT),
      fns = AVAL ~ mean(., na.rm = TRUE),
      filter = dplyr::n() > 2,
      set_values_to = vars(DTYPE = "AVERAGE")
    ),
    "deprecated",
    fixed = TRUE
  )
})

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

test_that("a warning is issued when using `derive_var_lstalvdt()`", {
  adsl <- tibble::tribble(
    ~STUDYID, ~USUBJID, ~TRTEDTM, ~DTHDTC,
    "STUDY01",  "1", ymd_hms("2020-01-01T12:00:00"), NA_character_,
    "STUDY01",  "2", NA, "2020-06",
    "STUDY01",  "3", ymd_hms("2020-04-12T13:15:00"), NA_character_
  )

  adsl_trtdate <- lstalvdt_source(
    dataset_name = "adsl",
    date = TRTEDTM
  )

  adsl_dthdate <- lstalvdt_source(
    dataset_name = "adsl",
    date = DTHDTC,
    filter = nchar(DTHDTC) >= 10
  )

  expect_warning(
    derive_var_lstalvdt(
      adsl,
      source_datasets = list(adsl = adsl),
      ae_start,
      ae_end,
      adsl_trtdate,
      adsl_dthdate
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
