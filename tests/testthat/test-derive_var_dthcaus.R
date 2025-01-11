## Test 1: deprecation message if function is called ----
test_that("derive_var_dthcaus Test 1: deprecation message if function is called", {
  adsl <- tibble::tribble(
    ~STUDYID, ~USUBJID,
    "TEST01", "PAT01",
    "TEST01", "PAT02",
    "TEST01", "PAT03"
  )

  ae <- tibble::tribble(
    ~STUDYID, ~USUBJID, ~AESEQ, ~AEDECOD, ~AEOUT, ~AEDTHDTC,
    "TEST01", "PAT03", 12, "SUDDEN DEATH", "FATAL", "2021-04-04"
  ) %>%
    mutate(
      AEDTHDT = ymd(AEDTHDTC)
    )

  ds <- tibble::tribble(
    ~STUDYID, ~USUBJID, ~DSSEQ, ~DSDECOD, ~DSTERM, ~DSSTDTC,
    "TEST01", "PAT01", 1, "INFORMED CONSENT OBTAINED", "INFORMED CONSENT OBTAINED", "2021-04-01",
    "TEST01", "PAT01", 2, "RANDOMIZATION", "RANDOMIZATION", "2021-04-11",
    "TEST01", "PAT01", 3, "ADVERSE EVENT", "ADVERSE EVENT", "2021-12-01",
    "TEST01", "PAT01", 4, "DEATH", "DEATH DUE TO progression of disease", "2022-02-01",
    "TEST01", "PAT02", 1, "INFORMED CONSENT OBTAINED", "INFORMED CONSENT OBTAINED", "2021-04-02",
    "TEST01", "PAT02", 2, "RANDOMIZATION", "RANDOMIZATION", "2021-04-11",
    "TEST01", "PAT02", 3, "COMPLETED", "PROTOCOL COMPLETED", "2021-12-01",
    "TEST01", "PAT03", 1, "INFORMED CONSENT OBTAINED", "INFORMED CONSENT OBTAINED", "2021-04-03",
    "TEST01", "PAT03", 2, "RANDOMIZATION", "RANDOMIZATION", "2021-04-11",
    "TEST01", "PAT03", 3, "COMPLETED", "PROTOCOL COMPLETED", "2021-12-01"
  )

  expect_snapshot({
    src_ae <- dthcaus_source(
      dataset_name = "ae",
      filter = AEOUT == "FATAL",
      date = AEDTHDT,
      mode = "first",
      dthcaus = AEDECOD
    )

    src_ds <- dthcaus_source(
      dataset_name = "ds",
      filter = DSDECOD == "DEATH" & grepl("DEATH DUE TO", DSTERM),
      date = convert_dtc_to_dt(DSSTDTC),
      mode = "first",
      dthcaus = str_to_upper(DSTERM)
    )

    derive_var_dthcaus(
      adsl,
      source_datasets = list(ae = ae, ds = ds),
      src_ae, src_ds
    )
  })
})

## Test 2: error on invalid mode ----
test_that("derive_var_dthcaus Test 2: error on invalid mode", {
  # Suppress lifecycle messages within the test environment
  withr::local_options(list(lifecycle_verbosity = "quiet"))

  expect_error(dthcaus_source(
    dataset_name = "ae",
    filter = AEOUT == "FATAL",
    date = AEDTHDTC,
    mode = "blah",
    dthcaus = AEDECOD
  ))
})

# derive_var_dthcaus ----
## Test 3: DTHCAUS is added from AE and DS ----
test_that("derive_var_dthcaus Test 3: DTHCAUS is added from AE and DS", {
  # Suppress lifecycle messages within the test environment
  withr::local_options(list(lifecycle_verbosity = "quiet"))

  adsl <- tibble::tribble(
    ~STUDYID, ~USUBJID,
    "TEST01", "PAT01",
    "TEST01", "PAT02",
    "TEST01", "PAT03"
  )

  ae <- tibble::tribble(
    ~STUDYID, ~USUBJID, ~AESEQ, ~AEDECOD, ~AEOUT, ~AEDTHDTC,
    "TEST01", "PAT03", 12, "SUDDEN DEATH", "FATAL", "2021-04-04"
  ) %>%
    mutate(
      AEDTHDT = ymd(AEDTHDTC)
    )

  ds <- tibble::tribble(
    ~STUDYID, ~USUBJID, ~DSSEQ, ~DSDECOD, ~DSTERM, ~DSSTDTC,
    "TEST01", "PAT01", 1, "INFORMED CONSENT OBTAINED", "INFORMED CONSENT OBTAINED", "2021-04-01",
    "TEST01", "PAT01", 2, "RANDOMIZATION", "RANDOMIZATION", "2021-04-11",
    "TEST01", "PAT01", 3, "ADVERSE EVENT", "ADVERSE EVENT", "2021-12-01",
    "TEST01", "PAT01", 4, "DEATH", "DEATH DUE TO progression of disease", "2022-02-01",
    "TEST01", "PAT02", 1, "INFORMED CONSENT OBTAINED", "INFORMED CONSENT OBTAINED", "2021-04-02",
    "TEST01", "PAT02", 2, "RANDOMIZATION", "RANDOMIZATION", "2021-04-11",
    "TEST01", "PAT02", 3, "COMPLETED", "PROTOCOL COMPLETED", "2021-12-01",
    "TEST01", "PAT03", 1, "INFORMED CONSENT OBTAINED", "INFORMED CONSENT OBTAINED", "2021-04-03",
    "TEST01", "PAT03", 2, "RANDOMIZATION", "RANDOMIZATION", "2021-04-11",
    "TEST01", "PAT03", 3, "COMPLETED", "PROTOCOL COMPLETED", "2021-12-01"
  )

  src_ae <- dthcaus_source(
    dataset_name = "ae",
    filter = AEOUT == "FATAL",
    date = AEDTHDT,
    mode = "first",
    dthcaus = AEDECOD
  )

  src_ds <- dthcaus_source(
    dataset_name = "ds",
    filter = DSDECOD == "DEATH" & grepl("DEATH DUE TO", DSTERM),
    date = convert_dtc_to_dt(DSSTDTC),
    mode = "first",
    dthcaus = str_to_upper(DSTERM)
  )

  expected_output <- tibble::tribble(
    ~STUDYID, ~USUBJID, ~DTHCAUS,
    "TEST01", "PAT01", "DEATH DUE TO PROGRESSION OF DISEASE",
    "TEST01", "PAT02", NA,
    "TEST01", "PAT03", "SUDDEN DEATH"
  )

  actual_output <- derive_var_dthcaus(
    adsl,
    source_datasets = list(ae = ae, ds = ds),
    src_ae, src_ds
  )

  expect_dfs_equal(expected_output, actual_output, keys = "USUBJID")
})

## Test 4: `dthcaus` handles symbols and string literals correctly ----
test_that("derive_var_dthcaus Test 4: `dthcaus` handles symbols and string literals correctly", {
  # Suppress lifecycle messages within the test environment
  withr::local_options(list(lifecycle_verbosity = "quiet"))

  adsl <- tibble::tribble(
    ~STUDYID, ~USUBJID,
    "TEST01", "PAT01",
    "TEST01", "PAT02"
  )

  ae <- tibble::tribble(
    ~STUDYID, ~USUBJID, ~AESEQ, ~AEDECOD, ~AEOUT, ~AEDTHDTC,
    "TEST01", "PAT01", 12, "SUDDEN DEATH", "FATAL", "2021-04-04"
  ) %>%
    mutate(
      AEDTHDT = ymd(AEDTHDTC)
    )

  ds <- tibble::tribble(
    ~STUDYID, ~USUBJID, ~DSSEQ, ~DSDECOD, ~DSTERM, ~DSSTDTC,
    "TEST01", "PAT01", 1, "INFORMED CONSENT OBTAINED", "INFORMED CONSENT OBTAINED", "2021-04-02",
    "TEST01", "PAT01", 2, "RANDOMIZATION", "RANDOMIZATION", "2021-04-11",
    "TEST01", "PAT01", 3, "COMPLETED", "PROTOCOL COMPLETED", "2021-12-01",
    "TEST01", "PAT02", 1, "INFORMED CONSENT OBTAINED", "INFORMED CONSENT OBTAINED", "2021-04-01",
    "TEST01", "PAT02", 2, "RANDOMIZATION", "RANDOMIZATION", "2021-04-11",
    "TEST01", "PAT02", 3, "ADVERSE EVENT", "ADVERSE EVENT", "2021-12-01",
    "TEST01", "PAT02", 4, "DEATH", "DEATH DUE TO PROGRESSION OF DISEASE", "2022-02-01"
  ) %>%
    mutate(
      DSSTDT = ymd(DSSTDTC)
    )

  src_ae <- dthcaus_source(
    dataset_name = "ae",
    filter = AEOUT == "FATAL",
    date = AEDTHDT,
    mode = "first",
    dthcaus = "Adverse Event"
  )

  src_ds <- dthcaus_source(
    dataset_name = "ds",
    filter = DSDECOD == "DEATH" & grepl("DEATH DUE TO", DSTERM),
    date = DSSTDT,
    mode = "first",
    dthcaus = DSTERM
  )

  expected_output <- tibble::tribble(
    ~STUDYID, ~USUBJID, ~DTHCAUS,
    "TEST01", "PAT01", "Adverse Event",
    "TEST01", "PAT02", "DEATH DUE TO PROGRESSION OF DISEASE"
  )

  actual_output <- derive_var_dthcaus(
    adsl,
    source_datasets = list(ae = ae, ds = ds),
    src_ae, src_ds
  )

  expect_dfs_equal(expected_output, actual_output, keys = "USUBJID")
})

## Test 5: traceability variables are added from AE and DS ----
test_that("derive_var_dthcaus Test 5: traceability variables are added from AE and DS", {
  # Suppress lifecycle messages within the test environment
  withr::local_options(list(lifecycle_verbosity = "quiet"))

  adsl <- tibble::tribble(
    ~STUDYID, ~USUBJID,
    "TEST01", "PAT01",
    "TEST01", "PAT02",
    "TEST01", "PAT03"
  )

  ae <- tibble::tribble(
    ~STUDYID, ~USUBJID, ~AESEQ, ~AEDECOD, ~AEOUT, ~AEDTHDTC,
    "TEST01", "PAT03", 12, "SUDDEN DEATH", "FATAL", "2021-04-04"
  ) %>%
    mutate(
      AEDTHDT = ymd(AEDTHDTC)
    )

  ds <- tibble::tribble(
    ~STUDYID, ~USUBJID, ~DSSEQ, ~DSDECOD, ~DSTERM, ~DSSTDTC,
    "TEST01", "PAT01", 1, "INFORMED CONSENT OBTAINED", "INFORMED CONSENT OBTAINED", "2021-04-01",
    "TEST01", "PAT01", 2, "RANDOMIZATION", "RANDOMIZATION", "2021-04-11",
    "TEST01", "PAT01", 3, "ADVERSE EVENT", "ADVERSE EVENT", "2021-12-01",
    "TEST01", "PAT01", 4, "DEATH", "DEATH DUE TO PROGRESSION OF DISEASE", "2022-02-01",
    "TEST01", "PAT02", 1, "INFORMED CONSENT OBTAINED", "INFORMED CONSENT OBTAINED", "2021-04-02",
    "TEST01", "PAT02", 2, "RANDOMIZATION", "RANDOMIZATION", "2021-04-11",
    "TEST01", "PAT02", 3, "COMPLETED", "PROTOCOL COMPLETED", "2021-12-01",
    "TEST01", "PAT03", 1, "INFORMED CONSENT OBTAINED", "INFORMED CONSENT OBTAINED", "2021-04-03",
    "TEST01", "PAT03", 2, "RANDOMIZATION", "RANDOMIZATION", "2021-04-11",
    "TEST01", "PAT03", 3, "COMPLETED", "PROTOCOL COMPLETED", "2021-12-01"
  ) %>%
    mutate(
      DSSTDT = ymd(DSSTDTC)
    )

  src_ae <- dthcaus_source(
    dataset_name = "ae",
    filter = AEOUT == "FATAL",
    date = AEDTHDT,
    mode = "first",
    dthcaus = AEDECOD,
    set_values_to = exprs(DTHDOM = "AE", DTHSEQ = AESEQ)
  )

  src_ds <- dthcaus_source(
    dataset_name = "ds",
    filter = DSDECOD == "DEATH" & grepl("DEATH DUE TO", DSTERM),
    date = DSSTDT,
    mode = "first",
    dthcaus = DSTERM,
    set_values_to = exprs(DTHDOM = "DS", DTHSEQ = DSSEQ)
  )

  expected_output <- tibble::tribble(
    ~STUDYID, ~USUBJID, ~DTHCAUS, ~DTHDOM, ~DTHSEQ,
    "TEST01", "PAT01", "DEATH DUE TO PROGRESSION OF DISEASE", "DS", 4,
    "TEST01", "PAT02", NA, NA, NA,
    "TEST01", "PAT03", "SUDDEN DEATH", "AE", 12
  )

  actual_output <- derive_var_dthcaus(
    adsl,
    source_datasets = list(ae = ae, ds = ds),
    src_ae, src_ds
  )

  expect_dfs_equal(expected_output, actual_output, keys = "USUBJID")
})

## Test 6: DTHCAUS/traceabiity are added from 2 input datasets ----
test_that("derive_var_dthcaus Test 6: DTHCAUS/traceabiity are added from 2 input datasets", {
  # Suppress lifecycle messages within the test environment
  withr::local_options(list(lifecycle_verbosity = "quiet"))

  adsl <- tibble::tribble(
    ~STUDYID, ~USUBJID,
    "TEST01", "PAT01",
    "TEST01", "PAT02",
    "TEST01", "PAT03"
  )

  ae <- tibble::tribble(
    ~STUDYID, ~USUBJID, ~AESEQ, ~AEDECOD, ~AEOUT, ~AEDTHDTC,
    "TEST01", "PAT01", 14, "SUDDEN DEATH", "FATAL", "2021-04-04",
    "TEST01", "PAT03", 12, "SUDDEN DEATH", "FATAL", "2021-04-04"
  ) %>%
    mutate(
      AEDTHDT = ymd(AEDTHDTC)
    )

  ds <- tibble::tribble(
    ~STUDYID, ~USUBJID, ~DSSEQ, ~DSDECOD, ~DSTERM, ~DSSTDTC,
    "TEST01", "PAT01", 1, "INFORMED CONSENT OBTAINED", "INFORMED CONSENT OBTAINED", "2021-04-01",
    "TEST01", "PAT01", 2, "RANDOMIZATION", "RANDOMIZATION", "2021-04-11",
    "TEST01", "PAT01", 3, "ADVERSE EVENT", "ADVERSE EVENT", "2021-12-01",
    "TEST01", "PAT01", 4, "DEATH", "DEATH DUE TO PROGRESSION OF DISEASE", "2021-02-03",
    "TEST01", "PAT02", 1, "INFORMED CONSENT OBTAINED", "INFORMED CONSENT OBTAINED", "2021-04-02",
    "TEST01", "PAT02", 2, "RANDOMIZATION", "RANDOMIZATION", "2021-04-11",
    "TEST01", "PAT02", 3, "COMPLETED", "PROTOCOL COMPLETED", "2021-12-01",
    "TEST01", "PAT03", 1, "INFORMED CONSENT OBTAINED", "INFORMED CONSENT OBTAINED", "2021-04-03",
    "TEST01", "PAT03", 2, "RANDOMIZATION", "RANDOMIZATION", "2021-04-11",
    "TEST01", "PAT03", 3, "COMPLETED", "PROTOCOL COMPLETED", "2021-12-01"
  ) %>%
    mutate(
      DSSTDT = ymd(DSSTDTC)
    )

  src_ae <- dthcaus_source(
    dataset_name = "ae",
    filter = AEOUT == "FATAL",
    date = AEDTHDT,
    mode = "first",
    dthcaus = AEDECOD,
    set_values_to = exprs(DTHDOM = "AE", DTHSEQ = AESEQ)
  )

  src_ds <- dthcaus_source(
    dataset_name = "ds",
    filter = DSDECOD == "DEATH" & grepl("DEATH DUE TO", DSTERM),
    date = DSSTDT,
    mode = "first",
    dthcaus = DSTERM,
    set_values_to = exprs(DTHDOM = "DS", DTHSEQ = DSSEQ)
  )

  expected_output <- tibble::tribble(
    ~STUDYID, ~USUBJID, ~DTHCAUS, ~DTHDOM, ~DTHSEQ,
    "TEST01", "PAT01", "DEATH DUE TO PROGRESSION OF DISEASE", "DS", 4,
    "TEST01", "PAT02", NA, NA, NA,
    "TEST01", "PAT03", "SUDDEN DEATH", "AE", 12
  )

  actual_output <- derive_var_dthcaus(
    adsl,
    source_datasets = list(ae = ae, ds = ds),
    src_ae, src_ds
  )

  expect_dfs_equal(expected_output, actual_output, keys = "USUBJID")
})

## Test 7: DTHCAUS is added from AE and DS if filter is not specified ----
test_that("derive_var_dthcaus Test 7: DTHCAUS is added from AE and DS if filter is not specified", {
  # Suppress lifecycle messages within the test environment
  withr::local_options(list(lifecycle_verbosity = "quiet"))

  # test based on covr report - the case for unspecified filter has not been tested

  adsl <- tibble::tribble(
    ~STUDYID, ~USUBJID,
    "TEST01", "PAT01",
    "TEST01", "PAT02",
    "TEST01", "PAT03"
  )

  ae <- tibble::tribble(
    ~STUDYID, ~USUBJID, ~AESEQ, ~AEDECOD, ~AEOUT, ~AEDTHDTC,
    "TEST01", "PAT03", 12, "SUDDEN DEATH", "FATAL", "2021-04-04"
  ) %>%
    mutate(
      AEDTHDT = ymd(AEDTHDTC)
    )

  ds <- tibble::tribble(
    ~STUDYID, ~USUBJID, ~DSSEQ, ~DSDECOD, ~DSTERM, ~DSSTDTC,
    "TEST01", "PAT01", 1, "INFORMED CONSENT OBTAINED", "INFORMED CONSENT OBTAINED", "2021-04-01",
    "TEST01", "PAT01", 2, "RANDOMIZATION", "RANDOMIZATION", "2021-04-11",
    "TEST01", "PAT01", 3, "ADVERSE EVENT", "ADVERSE EVENT", "2021-12-01",
    "TEST01", "PAT01", 4, "DEATH", "DEATH DUE TO PROGRESSION OF DISEASE", "2022-02-01",
    "TEST01", "PAT02", 1, "INFORMED CONSENT OBTAINED", "INFORMED CONSENT OBTAINED", "2021-04-02",
    "TEST01", "PAT02", 2, "RANDOMIZATION", "RANDOMIZATION", "2021-04-11",
    "TEST01", "PAT02", 3, "COMPLETED", "PROTOCOL COMPLETED", "2021-12-01",
    "TEST01", "PAT03", 1, "INFORMED CONSENT OBTAINED", "INFORMED CONSENT OBTAINED", "2021-04-03",
    "TEST01", "PAT03", 2, "RANDOMIZATION", "RANDOMIZATION", "2021-04-11",
    "TEST01", "PAT03", 3, "COMPLETED", "PROTOCOL COMPLETED", "2021-12-01"
  ) %>%
    mutate(
      DSSTDT = ymd(DSSTDTC)
    )

  src_ae <- dthcaus_source(
    dataset_name = "ae",
    filter = AEOUT == "FATAL",
    date = AEDTHDT,
    mode = "first",
    dthcaus = AEDECOD
  )

  src_ds <- dthcaus_source(
    dataset_name = "ds",
    filter = NULL,
    date = DSSTDT,
    mode = "first",
    dthcaus = DSTERM
  )

  expected_output <- tibble::tribble(
    ~STUDYID, ~USUBJID, ~DTHCAUS,
    "TEST01", "PAT01", "INFORMED CONSENT OBTAINED",
    "TEST01", "PAT02", "INFORMED CONSENT OBTAINED",
    "TEST01", "PAT03", "INFORMED CONSENT OBTAINED"
  )

  actual_output <- derive_var_dthcaus(
    adsl,
    source_datasets = list(ae = ae, ds = ds),
    src_ae, src_ds
  )

  expect_dfs_equal(expected_output, actual_output, keys = "USUBJID")
})

## Test 8: error on a dthcaus_source object with invalid order ----
test_that("derive_var_dthcaus Test 8: error on a dthcaus_source object with invalid order", {
  # Suppress lifecycle messages within the test environment
  withr::local_options(list(lifecycle_verbosity = "quiet"))

  expect_error(dthcaus_source(
    dataset_name = "ae",
    filter = AEOUT == "FATAL",
    date = AEDTHDTC,
    order = c(AESEQ),
    mode = "first",
    dthcaus = AEDECOD
  ))
})

## Test 9: `dataset` is sorted using the `order` parameter ----
test_that("derive_var_dthcaus Test 9: `dataset` is sorted using the `order` parameter", {
  # Suppress lifecycle messages within the test environment
  withr::local_options(list(lifecycle_verbosity = "quiet"))

  adsl <- tibble::tribble(
    ~STUDYID, ~USUBJID,
    "TEST01", "PAT01",
    "TEST01", "PAT02"
  )

  ae <- tibble::tribble(
    ~STUDYID, ~USUBJID, ~AESEQ, ~AEDECOD, ~AEOUT, ~AEDTHDTC,
    "TEST01", "PAT01", 12, "SUDDEN DEATH", "FATAL", "2021-02-04"
  ) %>%
    mutate(
      AEDTHDT = ymd(AEDTHDTC)
    )

  ds <- tibble::tribble(
    ~STUDYID, ~USUBJID, ~DSSEQ, ~DSDECOD, ~DSTERM, ~DSSTDTC,
    "TEST01", "PAT01", 1, "INFORMED CONSENT OBTAINED", "INFORMED CONSENT OBTAINED", "2021-04-02",
    "TEST01", "PAT01", 2, "RANDOMIZATION", "RANDOMIZATION", "2021-04-11",
    "TEST01", "PAT01", 3, "COMPLETED", "PROTOCOL COMPLETED", "2021-12-01",
    "TEST01", "PAT02", 1, "INFORMED CONSENT OBTAINED", "INFORMED CONSENT OBTAINED", "2021-04-01",
    "TEST01", "PAT02", 2, "RANDOMIZATION", "RANDOMIZATION", "2021-04-11",
    "TEST01", "PAT02", 3, "DEATH", "DEATH DUE TO ADVERSE EVENT", "2022-02-02",
    "TEST01", "PAT02", 4, "DEATH", "DEATH DUE TO PROGRESSION OF DISEASE", "2022-02-02"
  ) %>%
    mutate(
      DSSTDT = ymd(DSSTDTC)
    )

  src_ae <- dthcaus_source(
    dataset_name = "ae",
    filter = AEOUT == "FATAL",
    date = AEDTHDT,
    mode = "first",
    dthcaus = "Adverse Event"
  )

  src_ds <- dthcaus_source(
    dataset_name = "ds",
    filter = DSDECOD == "DEATH" & grepl("DEATH DUE TO", DSTERM),
    date = DSSTDT,
    order = exprs(DSSEQ),
    mode = "last",
    dthcaus = DSTERM
  )

  expected_output <- tibble::tribble(
    ~STUDYID, ~USUBJID, ~DTHCAUS,
    "TEST01", "PAT01", "Adverse Event",
    "TEST01", "PAT02", "DEATH DUE TO PROGRESSION OF DISEASE"
  )

  actual_output <- derive_var_dthcaus(
    adsl,
    source_datasets = list(ae = ae, ds = ds),
    src_ae, src_ds
  )

  expect_dfs_equal(expected_output, actual_output, keys = "USUBJID")
})

## Test 10: multiple observations from different sources ----
test_that("derive_var_dthcaus Test 10: multiple observations from different sources", {
  # Suppress lifecycle messages within the test environment
  withr::local_options(list(lifecycle_verbosity = "quiet"))

  expected <- tibble::tribble(
    ~STUDYID, ~USUBJID, ~DTHCAUS,
    "TEST01", "PAT01",  "SUDDEN DEATH",
    "TEST01", "PAT02",  NA_character_,
    "TEST01", "PAT03",  "DEATH DUE TO progression of disease"
  )

  adsl <- select(expected, -DTHCAUS)

  ae <- tibble::tribble(
    ~STUDYID, ~USUBJID, ~AESEQ, ~AEDECOD, ~AEOUT, ~AEDTHDTC,
    "TEST01", "PAT01", 12, "SUDDEN DEATH", "FATAL", "2021-04-04"
  ) %>%
    mutate(
      AEDTHDT = ymd(AEDTHDTC)
    )

  ds <- tibble::tribble(
    ~STUDYID, ~USUBJID, ~DSSEQ, ~DSDECOD, ~DSTERM, ~DSSTDTC,
    "TEST01", "PAT01", 4, "DEATH", "DEATH DUE TO progression of disease", "2021-04-05",
    "TEST01", "PAT02", 1, "INFORMED CONSENT OBTAINED", "INFORMED CONSENT OBTAINED", "2021-04-02",
    "TEST01", "PAT02", 2, "RANDOMIZATION", "RANDOMIZATION", "2021-04-11",
    "TEST01", "PAT02", 3, "COMPLETED", "PROTOCOL COMPLETED", "2021-12-01",
    "TEST01", "PAT03", 1, "DEATH", "DEATH DUE TO progression of disease", "2021-04-07",
    "TEST01", "PAT03", 2, "RANDOMIZATION", "RANDOMIZATION", "2021-04-11",
    "TEST01", "PAT03", 3, "COMPLETED", "PROTOCOL COMPLETED", "2021-12-01"
  )

  # Derive `DTHCAUS` only - for on-study deaths only
  src_ae <- dthcaus_source(
    dataset_name = "ae",
    filter = AEOUT == "FATAL",
    date = convert_dtc_to_dt(AEDTHDTC),
    mode = "first",
    dthcaus = AEDECOD
  )

  src_ds <- dthcaus_source(
    dataset_name = "ds",
    filter = DSDECOD == "DEATH" & grepl("DEATH DUE TO", DSTERM),
    date = convert_dtc_to_dt(DSSTDTC),
    mode = "first",
    dthcaus = DSTERM
  )
  actual <- adsl %>%
    derive_var_dthcaus(src_ae, src_ds, source_datasets = list(ae = ae, ds = ds))

  expect_dfs_equal(
    base = expected,
    compare = actual,
    keys = c("USUBJID")
  )
})

## Test 11: multiple observations with same date ----
test_that("derive_var_dthcaus Test 11: multiple observations with same date", {
  # Suppress lifecycle messages within the test environment
  withr::local_options(list(lifecycle_verbosity = "quiet"))

  expected <- tibble::tribble(
    ~STUDYID, ~USUBJID, ~DTHCAUS,
    "TEST01", "PAT01",  "SUDDEN DEATH",
    "TEST01", "PAT02",  NA_character_,
    "TEST01", "PAT03",  "DEATH DUE TO progression of disease"
  )

  adsl <- select(expected, -DTHCAUS)

  ae <- tibble::tribble(
    ~STUDYID, ~USUBJID, ~AESEQ, ~AEDECOD, ~AEOUT, ~AEDTHDTC,
    "TEST01", "PAT01", 12, "SUDDEN DEATH", "FATAL", "2021-04-05"
  ) %>%
    mutate(
      AEDTHDT = ymd(AEDTHDTC)
    )

  ds <- tibble::tribble(
    ~STUDYID, ~USUBJID, ~DSSEQ, ~DSDECOD, ~DSTERM, ~DSSTDTC,
    "TEST01", "PAT01", 4, "DEATH", "DEATH DUE TO progression of disease", "2021-04-05",
    "TEST01", "PAT02", 1, "INFORMED CONSENT OBTAINED", "INFORMED CONSENT OBTAINED", "2021-04-02",
    "TEST01", "PAT02", 2, "RANDOMIZATION", "RANDOMIZATION", "2021-04-11",
    "TEST01", "PAT02", 3, "COMPLETED", "PROTOCOL COMPLETED", "2021-12-01",
    "TEST01", "PAT03", 1, "DEATH", "DEATH DUE TO progression of disease", "2021-04-07",
    "TEST01", "PAT03", 2, "RANDOMIZATION", "RANDOMIZATION", "2021-04-11",
    "TEST01", "PAT03", 3, "COMPLETED", "PROTOCOL COMPLETED", "2021-12-01"
  )

  # Derive `DTHCAUS` only - for on-study deaths only
  src_ae <- dthcaus_source(
    dataset_name = "ae",
    filter = AEOUT == "FATAL",
    date = convert_dtc_to_dt(AEDTHDTC),
    mode = "first",
    dthcaus = AEDECOD
  )

  src_ds <- dthcaus_source(
    dataset_name = "ds",
    filter = DSDECOD == "DEATH" & grepl("DEATH DUE TO", DSTERM),
    date = convert_dtc_to_dt(DSSTDTC),
    mode = "first",
    dthcaus = DSTERM
  )
  actual <- adsl %>%
    derive_var_dthcaus(src_ae, src_ds, source_datasets = list(ae = ae, ds = ds))

  expect_dfs_equal(
    base = expected,
    compare = actual,
    keys = c("USUBJID")
  )
})

## Test 12: error if source dataset is not available ----
test_that("derive_var_dthcaus Test 12: error if source dataset is not available", {
  # Suppress lifecycle messages within the test environment
  withr::local_options(list(lifecycle_verbosity = "quiet"))

  adsl <- tibble::tribble(
    ~STUDYID, ~USUBJID,
    "TEST01", "PAT01",
    "TEST01", "PAT02",
    "TEST01", "PAT03"
  )

  ae <- tibble::tribble(
    ~STUDYID, ~USUBJID, ~AESEQ, ~AEDECOD, ~AEOUT, ~AEDTHDTC,
    "TEST01", "PAT03", 12, "SUDDEN DEATH", "FATAL", "2021-04-04"
  ) %>%
    mutate(
      AEDTHDT = ymd(AEDTHDTC)
    )

  ds <- tibble::tribble(
    ~STUDYID, ~USUBJID, ~DSSEQ, ~DSDECOD, ~DSTERM, ~DSSTDTC,
    "TEST01", "PAT01", 1, "INFORMED CONSENT OBTAINED", "INFORMED CONSENT OBTAINED", "2021-04-01",
    "TEST01", "PAT01", 2, "RANDOMIZATION", "RANDOMIZATION", "2021-04-11",
    "TEST01", "PAT01", 3, "ADVERSE EVENT", "ADVERSE EVENT", "2021-12-01",
    "TEST01", "PAT01", 4, "DEATH", "DEATH DUE TO progression of disease", "2022-02-01",
    "TEST01", "PAT02", 1, "INFORMED CONSENT OBTAINED", "INFORMED CONSENT OBTAINED", "2021-04-02",
    "TEST01", "PAT02", 2, "RANDOMIZATION", "RANDOMIZATION", "2021-04-11",
    "TEST01", "PAT02", 3, "COMPLETED", "PROTOCOL COMPLETED", "2021-12-01",
    "TEST01", "PAT03", 1, "INFORMED CONSENT OBTAINED", "INFORMED CONSENT OBTAINED", "2021-04-03",
    "TEST01", "PAT03", 2, "RANDOMIZATION", "RANDOMIZATION", "2021-04-11",
    "TEST01", "PAT03", 3, "COMPLETED", "PROTOCOL COMPLETED", "2021-12-01"
  )

  src_ae <- dthcaus_source(
    dataset_name = "ae",
    filter = AEOUT == "FATAL",
    date = AEDTHDT,
    mode = "first",
    dthcaus = AEDECOD
  )

  src_ds <- dthcaus_source(
    dataset_name = "ds",
    filter = DSDECOD == "DEATH" & grepl("DEATH DUE TO", DSTERM),
    date = convert_dtc_to_dt(DSSTDTC),
    mode = "first",
    dthcaus = str_to_upper(DSTERM)
  )

  expect_snapshot(
    derive_var_dthcaus(
      adsl,
      source_datasets = list(ae = ae, dd = ds),
      src_ae, src_ds
    ),
    error = TRUE
  )
})
