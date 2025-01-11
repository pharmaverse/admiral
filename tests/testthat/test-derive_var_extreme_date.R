adsl <- tibble::tribble(
  ~STUDYID,  ~USUBJID, ~TRTEDTM,                       ~DTHDTC,
  "STUDY01", "1",      ymd_hms("2020-01-01T12:00:00"), NA_character_,
  "STUDY01", "2",      NA,                             "2020-06",
  "STUDY01", "3",      ymd_hms("2020-04-12T13:15:00"), NA_character_
) %>%
  mutate(
    DTHDT = c(ymd(""), ymd("2020-06-01"), ymd(""))
  )

ae <- tibble::tribble(
  ~STUDYID,  ~USUBJID, ~AESTDTC,     ~AEENDTC,      ~AESEQ,
  "STUDY01", "1",      "2019-11-01", "2019-11-23",  1,
  "STUDY01", "1",      "2020-02-01", "2020-02-01",  2,
  "STUDY01", "3",      "2020-02-02", "2020-02-03",  1,
  "STUDY01", "3",      "2020-04-11", NA_character_, 2
) %>%
  mutate(
    AESTDT = ymd(AESTDTC),
    AEENDT = ymd(AEENDTC),
    AESTDTM = ymd_hms(paste(AESTDTC, "12:00:00")),
    AEENDTM = ymd_hms(if_else(is.na(AEENDTC), "", paste(AEENDTC, "12:00:00")))
  )

# derive_var_extreme_dt ----
## Test 1: deprecation message if function is called ----
test_that("derive_var_extreme_dt Test 1: deprecation message if function is called", {
  expect_snapshot({
    ae_start <- date_source(
      dataset_name = "ae",
      date = AESTDTM
    )

    ae_end <- date_source(
      dataset_name = "ae",
      date = AEENDTM
    )

    adsl_trtdate <- date_source(
      dataset_name = "adsl",
      date = TRTEDTM
    )

    adsl_dthdate <- date_source(
      dataset_name = "adsl",
      date = DTHDT,
      filter = nchar(DTHDTC) >= 10
    )


    derive_var_extreme_dt(
      adsl,
      new_var = LSTALVDT,
      source_datasets = list(ae = ae, adsl = adsl),
      ae_start, ae_end, adsl_trtdate, adsl_dthdate,
      mode = "last"
    )
  })
})

# derive_var_extreme_dt ----
## Test 2: LSTALVDT is derived ----
test_that("derive_var_extreme_dt Test 2: LSTALVDT is derived", {
  # Suppress lifecycle messages within the test environment
  withr::local_options(list(lifecycle_verbosity = "quiet"))

  ae_start <- date_source(
    dataset_name = "ae",
    date = AESTDTM
  )

  ae_end <- date_source(
    dataset_name = "ae",
    date = AEENDTM
  )

  adsl_trtdate <- date_source(
    dataset_name = "adsl",
    date = TRTEDTM
  )

  adsl_dthdate <- date_source(
    dataset_name = "adsl",
    date = DTHDT,
    filter = nchar(DTHDTC) >= 10
  )

  expected_output <- adsl %>% mutate(LSTALVDT = c(ymd("2020-02-01"), NA, ymd("2020-04-12")))

  actual_output <- derive_var_extreme_dt(
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

## Test 3: LSTALVDT is derived for Date class as well ----
test_that("derive_var_extreme_dt Test 3: LSTALVDT is derived for Date class as well", {
  # Suppress lifecycle messages within the test environment
  withr::local_options(list(lifecycle_verbosity = "quiet"))

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

  actual_output <- derive_var_extreme_dt(
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

## Test 4: `NA` dates are excluded ----
test_that("derive_var_extreme_dt Test 4: `NA` dates are excluded", {
  # Suppress lifecycle messages within the test environment
  withr::local_options(list(lifecycle_verbosity = "quiet"))

  ae_end <- date_source(
    dataset_name = "ae",
    date = AEENDTM
  )

  expected_output <- adsl %>% mutate(LSTALVDT = c(ymd("2020-02-01"), NA, ymd("2020-02-03")))

  actual_output <- derive_var_extreme_dt(
    adsl,
    new_var = LSTALVDT,
    source_datasets = list(ae = ae),
    ae_end,
    mode = "last"
  )

  expect_dfs_equal(
    base = expected_output,
    compare = actual_output,
    keys = c("USUBJID")
  )
})

# derive_var_extreme_dtm ----
## Test 5: Message sent to users ----
test_that("derive_var_extreme_dtm Test 5: Message sent to users", {
  expect_snapshot({
    ae_start <- date_source(
      dataset_name = "ae",
      date = convert_dtc_to_dtm(AESTDTC),
      set_values_to = exprs(
        LALVDOM = "AE",
        LALVSEQ = AESEQ,
        LALVVAR = "AESTDTC"
      )
    )

    ae_end <- date_source(
      dataset_name = "ae",
      date = AEENDTM,
      set_values_to = exprs(
        LALVDOM = "AE",
        LALVSEQ = AESEQ,
        LALVVAR = "AEENDTC"
      )
    )

    adsl_trtdate <- date_source(
      dataset_name = "adsl",
      date = TRTEDTM,
      set_values_to = exprs(
        LALVDOM = "ADSL",
        LALVSEQ = NA_integer_,
        LALVVAR = "TRTEDTM"
      )
    )

    adsl_dthdate <- date_source(
      dataset_name = "adsl",
      date = DTHDT,
      filter = nchar(DTHDTC) >= 10,
      set_values_to = exprs(
        LALVDOM = "ADSL",
        LALVSEQ = NA_integer_,
        LALVVAR = "DTHDTC"
      )
    )

    derive_var_extreme_dtm(
      adsl,
      new_var = LSTALVDTM,
      source_datasets = list(ae = ae, adsl = adsl),
      ae_start, ae_end, adsl_trtdate, adsl_dthdate,
      mode = "last"
    )
  })
})

# derive_var_extreme_dtm ----
## Test 6: `LSTALVDTM` and traceability variables are derived ----
test_that("derive_var_extreme_dtm Test 6: `LSTALVDTM` and traceability variables are derived", {
  # Suppress lifecycle messages within the test environment
  withr::local_options(list(lifecycle_verbosity = "quiet"))

  ae_start <- date_source(
    dataset_name = "ae",
    date = convert_dtc_to_dtm(AESTDTC),
    set_values_to = exprs(
      LALVDOM = "AE",
      LALVSEQ = AESEQ,
      LALVVAR = "AESTDTC"
    )
  )

  ae_end <- date_source(
    dataset_name = "ae",
    date = AEENDTM,
    set_values_to = exprs(
      LALVDOM = "AE",
      LALVSEQ = AESEQ,
      LALVVAR = "AEENDTC"
    )
  )

  adsl_trtdate <- date_source(
    dataset_name = "adsl",
    date = TRTEDTM,
    set_values_to = exprs(
      LALVDOM = "ADSL",
      LALVSEQ = NA_integer_,
      LALVVAR = "TRTEDTM"
    )
  )

  adsl_dthdate <- date_source(
    dataset_name = "adsl",
    date = DTHDT,
    filter = nchar(DTHDTC) >= 10,
    set_values_to = exprs(
      LALVDOM = "ADSL",
      LALVSEQ = NA_integer_,
      LALVVAR = "DTHDTC"
    )
  )

  expected_output <- adsl %>%
    mutate(
      LSTALVDTM = c(ymd_hms("2020-02-01T12:00:00"), NA, ymd_hms("2020-04-12T13:15:00")),
      LALVDOM = c("AE", NA_character_, "ADSL"),
      LALVSEQ = c(2, NA_integer_, NA_integer_),
      LALVVAR = c("AEENDTC", NA_character_, "TRTEDTM")
    )

  actual_output <- derive_var_extreme_dtm(
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

## Test 7: error is issued if `--DTC` variable is specified ----
test_that("derive_var_extreme_dtm Test 7: error is issued if `--DTC` variable is specified", {
  # Suppress lifecycle messages within the test environment
  withr::local_options(list(lifecycle_verbosity = "quiet"))

  ae_start <- date_source(
    dataset_name = "ae",
    date = AESTDTC,
    set_values_to = exprs(
      LALVDOM = "AE",
      LALVSEQ = AESEQ,
      LALVVAR = "AESTDTC"
    )
  )

  expect_error(
    derive_var_extreme_dtm(
      adsl,
      new_var = LSTALVDTM,
      source_datasets = list(ae = ae),
      ae_start,
      mode = "last"
    ),
    class = "assert_date_var"
  )
})

## Test 8: error if source dataset is not available ----
test_that("derive_var_extreme_dtm Test 8: error if source dataset is not available", {
  # Suppress lifecycle messages within the test environment
  withr::local_options(list(lifecycle_verbosity = "quiet"))

  expect_snapshot(
    {
      ae_start <- date_source(
        dataset_name = "ae",
        date = AESTDT,
        set_values_to = exprs(
          LALVDOM = "AE",
          LALVSEQ = AESEQ,
          LALVVAR = "AESTDTC"
        )
      )

      derive_var_extreme_dtm(
        adsl,
        new_var = LSTALVDTM,
        source_datasets = list(ea = ae),
        ae_start,
        mode = "last"
      )
    },
    error = TRUE
  )
})
