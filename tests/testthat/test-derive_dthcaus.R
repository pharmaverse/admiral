context("test-derive_dthcaus")

test_that("variable is added from the first observation in each group", {

  adsl <- tibble::tribble(
    ~STUDYID, ~USUBJID,
    "TEST01", "PAT01",
    "TEST01", "PAT02",
    "TEST01", "PAT03"
  )

  ae <- tibble::tribble(
    ~STUDYID, ~USUBJID, ~AEDECOD, ~AEOUT, ~AEDTHDTC,
    "TEST01", "PAT03", "SUDDEN DEATH" , "FATAL", "2021-04-04"
  )

  ds <- tibble::tribble(
    ~STUDYID, ~USUBJID, ~DSDECOD, ~DSTERM, ~DSSTDTC,
    "TEST01", "PAT01", "INFORMED CONSENT OBTAINED", "INFORMED CONSENT OBTAINED", "2021-04-01",
    "TEST01", "PAT01", "RANDOMIZATION", "RANDOMIZATION", "2021-04-11",
    "TEST01", "PAT01", "ADVERSE EVENT", "ADVERSE EVENT", "2021-12-01",
    "TEST01", "PAT01", "DEATH", "DEATH DUE TO PROGRESSION OF DISEASE", "2022-02-01",
    "TEST01", "PAT02", "INFORMED CONSENT OBTAINED", "INFORMED CONSENT OBTAINED", "2021-04-02",
    "TEST01", "PAT02", "RANDOMIZATION", "RANDOMIZATION", "2021-04-11",
    "TEST01", "PAT02", "COMPLETED", "PROTOCOL COMPLETED", "2021-12-01",
    "TEST01", "PAT03", "INFORMED CONSENT OBTAINED", "INFORMED CONSENT OBTAINED", "2021-04-03",
    "TEST01", "PAT03", "RANDOMIZATION", "RANDOMIZATION", "2021-04-11",
    "TEST01", "PAT03", "COMPLETED", "PROTOCOL COMPLETED", "2021-12-01"
  )

  actual_output <- derive_dthcaus(dataset = adsl, dataset_ae = ae, dataset_ds = ds)

})
