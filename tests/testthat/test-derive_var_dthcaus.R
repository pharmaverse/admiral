library(tibble)
library(dplyr)
library(lubridate)

# dthcaus_source ----
## Test 1: error on invalid mode ----
test_that("dthcaus_source Test 1: error on invalid mode", {
  expect_error(dthcaus_source(
    dataset_name = "ae",
    filter = AEOUT == "FATAL",
    date = AEDTHDTC,
    mode = "blah",
    dthcaus = AEDECOD
  ))
})

# derive_var_dthcaus ----
## Test 2: DTHCAUS is added from AE and DS ----
test_that("derive_var_dthcaus Test 2: DTHCAUS is added from AE and DS", {
  adsl <- tribble(
    ~STUDYID, ~USUBJID,
    "TEST01", "PAT01",
    "TEST01", "PAT02",
    "TEST01", "PAT03"
  )

  ae <- tribble(
    ~STUDYID, ~USUBJID, ~AESEQ, ~AEDECOD, ~AEOUT, ~AEDTHDTC,
    "TEST01", "PAT03", 12, "SUDDEN DEATH", "FATAL", "2021-04-04"
  ) %>%
    mutate(
      AEDTHDT = ymd(AEDTHDTC)
    )

  ds <- tribble(
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
    filter = DSDECOD == "DEATH" & grepl("DEATH DUE TO", DSTERM),
    date = DSSTDT,
    mode = "first",
    dthcaus = DSTERM
  )

  expected_output <- tribble(
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

## Test 3: `dthcaus` handles symbols and string literals correctly ----
test_that("derive_var_dthcaus Test 3: `dthcaus` handles symbols and string literals correctly", {
  adsl <- tribble(
    ~STUDYID, ~USUBJID,
    "TEST01", "PAT01",
    "TEST01", "PAT02"
  )

  ae <- tribble(
    ~STUDYID, ~USUBJID, ~AESEQ, ~AEDECOD, ~AEOUT, ~AEDTHDTC,
    "TEST01", "PAT01", 12, "SUDDEN DEATH", "FATAL", "2021-04-04"
  ) %>%
    mutate(
      AEDTHDT = ymd(AEDTHDTC)
    )

  ds <- tribble(
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

  expected_output <- tribble(
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

## Test 4: DTHCAUS and traceability vars are added from AE and DS ----
test_that("derive_var_dthcaus Test 4: DTHCAUS and traceability vars are added from AE and DS", {
  adsl <- tribble(
    ~STUDYID, ~USUBJID,
    "TEST01", "PAT01",
    "TEST01", "PAT02",
    "TEST01", "PAT03"
  )

  ae <- tribble(
    ~STUDYID, ~USUBJID, ~AESEQ, ~AEDECOD, ~AEOUT, ~AEDTHDTC,
    "TEST01", "PAT03", 12, "SUDDEN DEATH", "FATAL", "2021-04-04"
  ) %>%
    mutate(
      AEDTHDT = ymd(AEDTHDTC)
    )

  ds <- tribble(
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
    traceability_vars = vars(DTHDOM = "AE", DTHSEQ = AESEQ)
  )

  src_ds <- dthcaus_source(
    dataset_name = "ds",
    filter = DSDECOD == "DEATH" & grepl("DEATH DUE TO", DSTERM),
    date = DSSTDT,
    mode = "first",
    dthcaus = DSTERM,
    traceability_vars = vars(DTHDOM = "DS", DTHSEQ = DSSEQ)
  )

  expected_output <- tribble(
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

## Test 5: DTHCAUS/traceabiity are added from 2 input datasets ----
test_that("derive_var_dthcaus Test 5: DTHCAUS/traceabiity are added from 2 input datasets", {
  adsl <- tribble(
    ~STUDYID, ~USUBJID,
    "TEST01", "PAT01",
    "TEST01", "PAT02",
    "TEST01", "PAT03"
  )

  ae <- tribble(
    ~STUDYID, ~USUBJID, ~AESEQ, ~AEDECOD, ~AEOUT, ~AEDTHDTC,
    "TEST01", "PAT01", 14, "SUDDEN DEATH", "FATAL", "2021-04-04",
    "TEST01", "PAT03", 12, "SUDDEN DEATH", "FATAL", "2021-04-04"
  ) %>%
    mutate(
      AEDTHDT = ymd(AEDTHDTC)
    )

  ds <- tribble(
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
    traceability_vars = vars(DTHDOM = "AE", DTHSEQ = AESEQ)
  )

  src_ds <- dthcaus_source(
    dataset_name = "ds",
    filter = DSDECOD == "DEATH" & grepl("DEATH DUE TO", DSTERM),
    date = DSSTDT,
    mode = "first",
    dthcaus = DSTERM,
    traceability_vars = vars(DTHDOM = "DS", DTHSEQ = DSSEQ)
  )

  expected_output <- tribble(
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

## Test 6: DTHCAUS is added from AE and DS if filter is not specified ----
test_that("derive_var_dthcaus Test 6: DTHCAUS is added from AE and DS if filter is not specified", {
  # test based on covr report - the case for unspecified filter has not been tested

  adsl <- tribble(
    ~STUDYID, ~USUBJID,
    "TEST01", "PAT01",
    "TEST01", "PAT02",
    "TEST01", "PAT03"
  )

  ae <- tribble(
    ~STUDYID, ~USUBJID, ~AESEQ, ~AEDECOD, ~AEOUT, ~AEDTHDTC,
    "TEST01", "PAT03", 12, "SUDDEN DEATH", "FATAL", "2021-04-04"
  ) %>%
    mutate(
      AEDTHDT = ymd(AEDTHDTC)
    )

  ds <- tribble(
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

  expected_output <- tribble(
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
