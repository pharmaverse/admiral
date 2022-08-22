test_that("`CHG` is calculated as `AVAL - BASE`", {
  library(tibble)
  input <- tribble(
    ~STUDYID, ~USUBJID, ~PARAMCD, ~AVAL, ~ABLFL, ~BASETYPE, ~BASE,
    "TEST01", "PAT01", "PARAM01", 10.12, "Y", "LAST", 10.12,
    "TEST01", "PAT01", "PARAM01", 9.7, "", "LAST", 10.12,
    "TEST01", "PAT01", "PARAM01", 15.01, "", "LAST", 10.12,
    "TEST01", "PAT01", "PARAM02", 8.35, "Y", "LAST", 8.35,
    "TEST01", "PAT01", "PARAM02", NA, "", "LAST", 8.35,
    "TEST01", "PAT01", "PARAM02", 8.35, "", "LAST", 8.35,
    "TEST01", "PAT02", "PARAM01", 29, "Y", "LAST", 29,
    "TEST01", "PAT02", "PARAM01", 19.7, "", "LAST", 29,
    "TEST01", "PAT02", "PARAM01", 18.01, "", "LAST", 29,
    "TEST01", "PAT02", "PARAM02", 8.9, "Y", "LAST", 8.9,
    "TEST01", "PAT02", "PARAM02", 9, "", "LAST", 8.9,
    "TEST01", "PAT02", "PARAM02", 5.35, "", "LAST", 8.9
  )
  expected_output <- tribble(
    ~STUDYID, ~USUBJID, ~PARAMCD, ~AVAL, ~ABLFL, ~BASETYPE, ~BASE, ~CHG,
    "TEST01", "PAT01", "PARAM01", 10.12, "Y", "LAST", 10.12, 0,
    "TEST01", "PAT01", "PARAM01", 9.7, "", "LAST", 10.12, -0.42,
    "TEST01", "PAT01", "PARAM01", 15.01, "", "LAST", 10.12, 4.89,
    "TEST01", "PAT01", "PARAM02", 8.35, "Y", "LAST", 8.35, 0,
    "TEST01", "PAT01", "PARAM02", NA, "", "LAST", 8.35, NA,
    "TEST01", "PAT01", "PARAM02", 8.35, "", "LAST", 8.35, 0,
    "TEST01", "PAT02", "PARAM01", 29, "Y", "LAST", 29, 0,
    "TEST01", "PAT02", "PARAM01", 19.7, "", "LAST", 29, -9.3,
    "TEST01", "PAT02", "PARAM01", 18.01, "", "LAST", 29, -10.99,
    "TEST01", "PAT02", "PARAM02", 8.9, "Y", "LAST", 8.9, 0,
    "TEST01", "PAT02", "PARAM02", 9, "", "LAST", 8.9, 0.1,
    "TEST01", "PAT02", "PARAM02", 5.35, "", "LAST", 8.9, -3.55
  )

  expect_equal(derive_var_chg(input)$CHG, expected_output$CHG)
})

test_that("`PCHG` is calculated as `(AVAL - BASE) / abs(BASE) * 100`", {
  library(tibble)
  input <- tribble(
    ~STUDYID, ~USUBJID, ~PARAMCD, ~AVAL, ~ABLFL, ~BASETYPE, ~BASE,
    "TEST01", "PAT01", "PARAM01", -10.12, "Y", "LAST", -10.12,
    "TEST01", "PAT01", "PARAM01", -9.7, "", "LAST", -10.12,
    "TEST01", "PAT01", "PARAM01", -15.01, "", "LAST", -10.12,
    "TEST01", "PAT01", "PARAM02", 8.35, "Y", "LAST", 8.35,
    "TEST01", "PAT01", "PARAM02", NA, "", "LAST", 8.35,
    "TEST01", "PAT01", "PARAM02", 8.35, "", "LAST", 8.35,
    "TEST01", "PAT02", "PARAM01", 29, "Y", "LAST", 29,
    "TEST01", "PAT02", "PARAM01", 19.7, "", "LAST", 29,
    "TEST01", "PAT02", "PARAM01", 18.01, "", "LAST", 29,
    "TEST01", "PAT02", "PARAM02", 8.9, "Y", "LAST", 8.9,
    "TEST01", "PAT02", "PARAM02", 9, "", "LAST", 8.9,
    "TEST01", "PAT02", "PARAM02", 5.35, "", "LAST", 8.9
  )

  expected_output <- tribble(
    ~STUDYID, ~USUBJID, ~PARAMCD, ~AVAL, ~ABLFL, ~BASETYPE, ~BASE, ~PCHG,
    "TEST01", "PAT01", "PARAM01", -10.12, "Y", "LAST", -10.12, 0,
    "TEST01", "PAT01", "PARAM01", -9.7, "", "LAST", -10.12, 4.150198,
    "TEST01", "PAT01", "PARAM01", -15.01, "", "LAST", -10.12, -48.32016,
    "TEST01", "PAT01", "PARAM02", 8.35, "Y", "LAST", 8.35, 0,
    "TEST01", "PAT01", "PARAM02", NA, "", "LAST", 8.35, NA,
    "TEST01", "PAT01", "PARAM02", 8.35, "", "LAST", 8.35, 0,
    "TEST01", "PAT02", "PARAM01", 29, "Y", "LAST", 29, 0,
    "TEST01", "PAT02", "PARAM01", 19.7, "", "LAST", 29, -32.06897,
    "TEST01", "PAT02", "PARAM01", 18.01, "", "LAST", 29, -37.89655,
    "TEST01", "PAT02", "PARAM02", 8.9, "Y", "LAST", 8.9, 0,
    "TEST01", "PAT02", "PARAM02", 9, "", "LAST", 8.9, 1.123596,
    "TEST01", "PAT02", "PARAM02", 5.35, "", "LAST", 8.9, -39.88764
  )

  expect_equal(derive_var_pchg(input)$PCHG, expected_output$PCHG, tolerance = 1e-5)
})

test_that("`PCHG` is set to `NA` if `BASE == 0`", {
  library(tibble)
  input <- tribble(
    ~STUDYID, ~USUBJID, ~PARAMCD, ~AVAL, ~ABLFL, ~BASETYPE, ~BASE, ~CHG,
    "TEST01", "PAT01", "PARAM01", 0, "Y", "LAST", 0, 0,
    "TEST01", "PAT01", "PARAM01", 1.7, "", "LAST", 0, 1.7,
    "TEST01", "PAT01", "PARAM01", 3.01, "", "LAST", 0, 3.01,
  )

  expected_output <- tribble(
    ~STUDYID, ~USUBJID, ~PARAMCD, ~AVAL, ~ABLFL, ~BASETYPE, ~BASE, ~CHG, ~PCHG,
    "TEST01", "PAT01", "PARAM01", 0, "Y", "LAST", 0, 0, NA_real_,
    "TEST01", "PAT01", "PARAM01", 1.7, "", "LAST", 0, 1.7, NA_real_,
    "TEST01", "PAT01", "PARAM01", 3.01, "", "LAST", 0, 3.01, NA_real_,
  )

  expect_equal(derive_var_pchg(input)$PCHG, expected_output$PCHG, tolerance = 1e-7)
})
