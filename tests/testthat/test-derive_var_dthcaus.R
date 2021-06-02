context("test-derive_dthcaus")

test_that("error on a dthcaus_source object with invalid mode", {
  tmp <- dthcaus_source(
    dataset = ae,
    filter = AEOUT == "FATAL",
    date_var = AEDTHDTC,
    mode = "blah",
    dthcaus = AEDECOD
  )
  expect_error(validate_dthcaus_source(tmp))
})

test_that("DTHCAUS is added from AE and DS", {
  adsl <- tibble::tribble(
    ~STUDYID, ~USUBJID,
    "TEST01", "PAT01",
    "TEST01", "PAT02",
    "TEST01", "PAT03"
  )

  ae <- tibble::tribble(
    ~STUDYID, ~USUBJID, ~AESEQ, ~AEDECOD, ~AEOUT, ~AEDTHDTC,
    "TEST01", "PAT03", 12, "SUDDEN DEATH", "FATAL", "2021-04-04"
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
  )

  src_ae <- dthcaus_source(
    dataset = ae,
    filter = AEOUT == "FATAL",
    date_var = AEDTHDTC,
    mode = "first",
    dthcaus = AEDECOD
  )

  src_ds <- dthcaus_source(
    dataset = ds,
    filter = DSDECOD == "DEATH" & grepl("DEATH DUE TO", DSTERM),
    date_var = DSSTDTC,
    mode = "first",
    dthcaus = DSTERM
  )

  expected_output <- tibble::tribble(
    ~STUDYID, ~USUBJID, ~DTHCAUS,
    "TEST01", "PAT01", "DEATH DUE TO PROGRESSION OF DISEASE",
    "TEST01", "PAT02", NA,
    "TEST01", "PAT03", "SUDDEN DEATH"
  )

  actual_output <- derive_var_dthcaus(dataset = adsl, src_ae, src_ds)

  expect_dfs_equal(expected_output, actual_output, keys = "USUBJID")
})

test_that("DTHCAUS and traceabilty variables are added from AE and DS", {
  adsl <- tibble::tribble(
    ~STUDYID, ~USUBJID,
    "TEST01", "PAT01",
    "TEST01", "PAT02",
    "TEST01", "PAT03"
  )

  ae <- tibble::tribble(
    ~STUDYID, ~USUBJID, ~AESEQ, ~AEDECOD, ~AEOUT, ~AEDTHDTC,
    "TEST01", "PAT03", 12, "SUDDEN DEATH", "FATAL", "2021-04-04"
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
  )

  src_ae <- dthcaus_source(
    dataset = ae,
    filter = AEOUT == "FATAL",
    date_var = AEDTHDTC,
    mode = "first",
    dthcaus = AEDECOD,
    traceabilty_vars = vars(DTHDOM = "AE", DTHSEQ = AESEQ)
  )

  src_ds <- dthcaus_source(
    dataset = ds,
    filter = DSDECOD == "DEATH" & grepl("DEATH DUE TO", DSTERM),
    date = DSSTDTC,
    mode = "first",
    dthcaus = DSTERM,
    traceabilty_vars = vars(DTHDOM = "DS", DTHSEQ = DSSEQ)
  )

  expected_output <- tibble::tribble(
    ~STUDYID, ~USUBJID, ~DTHCAUS, ~DTHDOM, ~DTHSEQ,
    "TEST01", "PAT01", "DEATH DUE TO PROGRESSION OF DISEASE", "DS", 4,
    "TEST01", "PAT02", NA, NA, NA,
    "TEST01", "PAT03", "SUDDEN DEATH", "AE", 12
  )

  actual_output <- derive_var_dthcaus(dataset = adsl, src_ae, src_ds)

  expect_dfs_equal(expected_output, actual_output, keys = "USUBJID")
})

test_that("DTHCAUS/traceabiity are added from AE and DS, info available in 2 input datasets", {
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
  )

  src_ae <- dthcaus_source(
    dataset = ae,
    filter = AEOUT == "FATAL",
    date_var = AEDTHDTC,
    mode = "first",
    dthcaus = AEDECOD,
    traceabilty_vars = vars(DTHDOM = "AE", DTHSEQ = AESEQ)
  )

  src_ds <- dthcaus_source(
    dataset = ds,
    filter = DSDECOD == "DEATH" & grepl("DEATH DUE TO", DSTERM),
    date_var = DSSTDTC,
    mode = "first",
    dthcaus = DSTERM,
    traceabilty_vars = vars(DTHDOM = "DS", DTHSEQ = DSSEQ)
  )

  expected_output <- tibble::tribble(
    ~STUDYID, ~USUBJID, ~DTHCAUS, ~DTHDOM, ~DTHSEQ,
    "TEST01", "PAT01", "DEATH DUE TO PROGRESSION OF DISEASE", "DS", 4,
    "TEST01", "PAT02", NA, NA, NA,
    "TEST01", "PAT03", "SUDDEN DEATH", "AE", 12
  )

  actual_output <- derive_var_dthcaus(dataset = adsl, src_ae, src_ds)

  expect_dfs_equal(expected_output, actual_output, keys = "USUBJID")
})
test_that("DTHCAUS/traceabiity are added from AE and DS, info available in 2 input datasets, partial dates", { # nolint
  adsl <- tibble::tribble(
    ~STUDYID, ~USUBJID,
    "TEST01", "PAT01",
    "TEST01", "PAT02",
    "TEST01", "PAT03"
  )

  ae <- tibble::tribble(
    ~STUDYID, ~USUBJID, ~AESEQ, ~AEDECOD, ~AEOUT, ~AEDTHDTC,
    "TEST01", "PAT01", 14, "SUDDEN DEATH", "FATAL", "2021-05",
    "TEST01", "PAT03", 12, "SUDDEN DEATH", "FATAL", "2021-04-04"
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
  )

  src_ae <- dthcaus_source(
    dataset = ae,
    filter = AEOUT == "FATAL",
    date_var = AEDTHDTC,
    mode = "first",
    dthcaus = AEDECOD,
    traceabilty_vars = vars(DTHDOM = "AE", DTHSEQ = AESEQ)
  )

  src_ds <- dthcaus_source(
    dataset = ds,
    filter = DSDECOD == "DEATH" & grepl("DEATH DUE TO", DSTERM),
    date_var = DSSTDTC,
    mode = "first",
    dthcaus = DSTERM,
    traceabilty = vars(DTHDOM = "DS", DTHSEQ = DSSEQ)
  )

  expected_output <- tibble::tribble(
    ~STUDYID, ~USUBJID, ~DTHCAUS, ~DTHDOM, ~DTHSEQ,
    "TEST01", "PAT01", "DEATH DUE TO PROGRESSION OF DISEASE", "DS", 4,
    "TEST01", "PAT02", NA, NA, NA,
    "TEST01", "PAT03", "SUDDEN DEATH", "AE", 12
  )

  actual_output <- derive_var_dthcaus(dataset = adsl, src_ae, src_ds)

  expect_dfs_equal(expected_output, actual_output, keys = "USUBJID")
})
