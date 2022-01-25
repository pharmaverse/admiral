adsl <- tibble::tribble(
  ~STUDYID, ~USUBJID, ~TRTEDTM, ~DTHDTC,
  "STUDY01",  "1", ymd_hms("2020-01-01T12:00:00"), NA_character_,
  "STUDY01",  "2", NA, "2020-06",
  "STUDY01",  "3", ymd_hms("2020-04-12T13:15:00"), NA_character_
)

ae <- tibble::tribble(
  ~STUDYID, ~USUBJID, ~AESTDTC, ~AEENDTC, ~AESEQ,
  "STUDY01",  "1", "2019-11", "2019-11-23", 1,
  "STUDY01",  "1", "2020-02", "2020-02", 2,
  "STUDY01",  "3", "2020-02-02", "2020-02-03", 1,
  "STUDY01",  "3", "2020-04-11", NA_character_, 2
)

test_that("LSTALVDT is derived", {
  ae_start <- date_source(
    dataset_name = "ae",
    date = AESTDTC,
    date_imputation = "first"
  )

  ae_end <- date_source(
    dataset_name = "ae",
    date = AEENDTC,
    date_imputation = "first"
  )

  adsl_trtdate <- date_source(
    dataset_name = "adsl",
    date = TRTEDTM
  )

  adsl_dthdate <- date_source(
    dataset_name = "adsl",
    date = DTHDTC,
    filter = nchar(DTHDTC) >= 10
  )

  expected_output <- adsl %>% mutate(LSTALVDT = c(ymd("2020-02-01"), NA, ymd("2020-04-12")))

  actual_output <- derive_vars_extreme_dt(
    adsl,
    new_var = LSTALVDT,
    source_datasets = list(ae = ae, adsl = adsl),
    ae_start, ae_end, adsl_trtdate, adsl_dthdate,
    mode = "last"
  )

  expect_dfs_equal(
    base = expected_output,
    compare = actual_output,
    keys = c("USUBJID")
  )
})

test_that("LSTALVDTM and traceability variables are derived", {
  ae_start <- date_source(
    dataset_name = "ae",
    date = AESTDTC,
    date_imputation = "first",
    time_imputation = "first",
    traceability_vars = vars(
      LALVDOM = "AE",
      LALVSEQ = AESEQ,
      LALVVAR = "AESTDTC"
    )
  )

  ae_end <- date_source(
    dataset_name = "ae",
    date = AEENDTC,
    date_imputation = "first",
    time_imputation = "first",
    traceability_vars = vars(
      LALVDOM = "AE",
      LALVSEQ = AESEQ,
      LALVVAR = "AEENDTC"
    )
  )

  adsl_trtdate <- date_source(
    dataset_name = "adsl",
    date = TRTEDTM,
    traceability_vars = vars(
      LALVDOM = "ADSL",
      LALVSEQ = NA_integer_,
      LALVVAR = "TRTEDTM"
    )
  )

  adsl_dthdate <- date_source(
    dataset_name = "adsl",
    date = DTHDTC,
    filter = nchar(DTHDTC) >= 10,
    traceability_vars = vars(
      LALVDOM = "ADSL",
      LALVSEQ = NA_integer_,
      LALVVAR = "DTHDTC"
    )
  )

  expected_output <- adsl %>%
    mutate(
      LSTALVDTM = c(ymd_hms("2020-02-01T00:00:00"), NA, ymd_hms("2020-04-12T00:00:00")),
      LALVDOM = c("AE", NA_character_, "ADSL"),
      LALVSEQ = c(2, NA_integer_, NA_integer_),
      LALVVAR = c("AEENDTC", NA_character_, "TRTEDTM")
    )

  actual_output <- derive_vars_extreme_dtm(
    adsl,
    new_var = LSTALVDTM,
    source_datasets = list(ae = ae, adsl = adsl),
    ae_start, ae_end, adsl_trtdate, adsl_dthdate,
    mode = "last"
  )

  expect_dfs_equal(
    base = expected_output,
    compare = actual_output,
    keys = c("USUBJID")
  )
})

test_that("LSTALVDT is derived for Date class as well", {
  adsl <- tibble::tribble(
    ~STUDYID,  ~USUBJID, ~TRTEDTM,
    "STUDY01", "1",      ymd_hms("2020-01-01T12:00:00"),
    "STUDY01", "2",      as.POSIXct(ymd("2020-02-03")),
    "STUDY01", "3",      ymd_hms("2020-04-12T13:15:00")
  ) %>%
    mutate(TRTEDTM = as.Date(TRTEDTM))

  adsl_trtdate <- date_source(
    dataset_name = "adsl",
    date = TRTEDTM
  )

  expected_output <- adsl %>%
    mutate(LSTALVDT = c(ymd("2020-01-01"), ymd("2020-02-03"), ymd("2020-04-12")))

  actual_output <- derive_vars_extreme_dt(
    adsl,
    new_var = LSTALVDT,
    source_datasets = list(adsl = adsl),
    adsl_trtdate,
    mode = "last"
  )

  expect_dfs_equal(
    base = expected_output,
    compare = actual_output,
    keys = c("USUBJID")
  )
})
