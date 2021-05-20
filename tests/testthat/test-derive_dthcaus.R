context("test-derive_dthcaus")

test_that("error on a dthcaus_source object with invalid mode", {
  tmp <- list(
    dataset = ae,
    filter = expr(AEOUT == "FATAL"),
    order = exprs(AEDTHDTC),
    mode = "blah",
    dthdom = "AE",
    dthcaus = expr(AEDECOD)
  )
  class(tmp) <- "dthcaus_source"
  expect_error(validate_dthcaus_source(tmp))
})

test_that("DTHCAUS and DTHDOM are added from AE and DS", {

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

  src_ae <- dthcaus_source(
    dataset = ae,
    filter = expr(AEOUT == "FATAL"),
    order = exprs(AEDTHDTC),
    mode = "first",
    dthdom = "AE",
    dthcaus = expr(AEDECOD)
  )

  src_ds <- dthcaus_source(
    dataset = ds,
    filter = expr(DSDECOD == "DEATH" & grepl("DEATH DUE TO", DSTERM)),
    order = exprs(DSSTDTC),
    mode = "first",
    dthdom = "DS",
    dthcaus = expr(DSTERM)
  )

  expected_output <- tibble::tribble(
    ~STUDYID, ~USUBJID, ~DTHDOM, ~DTHCAUS,
    "TEST01", "PAT01", "DS", "DEATH DUE TO PROGRESSION OF DISEASE",
    "TEST01", "PAT02", NA, NA,
    "TEST01", "PAT03", "AE", "SUDDEN DEATH"
  )

  actual_output <- derive_dthcaus(dataset = adsl, src_ae, src_ds)

  expect_dfs_equal(expected_output, actual_output, keys = "USUBJID")
})
