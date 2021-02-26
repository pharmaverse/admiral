context("test-derive_var_chg")

test_that("`CHG` is calculated as `AVAL - BASE`", {
  input <- tibble::tribble(
    ~STUDYID, ~USUBJID, ~PARAMCD,  ~AVAL, ~ABLFL, ~BASETYPE, ~BASE,
    "TEST01", "PAT01",  "PARAM01", 10.12, "Y",    "LAST",     10.12,
    "TEST01", "PAT01",  "PARAM01",  9.7,  "",     "LAST",     10.12,
    "TEST01", "PAT01",  "PARAM01", 15.01, "",     "LAST",     10.12,
    "TEST01", "PAT01",  "PARAM02",  8.35, "Y",    "LAST",      8.35,
    "TEST01", "PAT01",  "PARAM02", NA,    "",     "LAST",      8.35,
    "TEST01", "PAT01",  "PARAM02",  8.35, "",     "LAST",      8.35,

    "TEST01", "PAT02",  "PARAM01", 29,    "Y",    "LAST",     29,
    "TEST01", "PAT02",  "PARAM01", 19.7,  "",     "LAST",     29,
    "TEST01", "PAT02",  "PARAM01", 18.01, "",     "LAST",     29,
    "TEST01", "PAT02",  "PARAM02",  8.9,  "Y",    "LAST",      8.9,
    "TEST01", "PAT02",  "PARAM02",  9,    "",     "LAST",      8.9,
    "TEST01", "PAT02",  "PARAM02",  5.35, "",     "LAST",      8.9
  )
  expected_output <- tibble::tribble(
    ~STUDYID, ~USUBJID, ~PARAMCD,  ~AVAL, ~ABLFL, ~BASETYPE, ~BASE,  ~CHG,
    "TEST01", "PAT01",  "PARAM01", 10.12, "Y",    "LAST",     10.12,  0,
    "TEST01", "PAT01",  "PARAM01",  9.7,  "",     "LAST",     10.12, -0.42,
    "TEST01", "PAT01",  "PARAM01", 15.01, "",     "LAST",     10.12,  4.89,
    "TEST01", "PAT01",  "PARAM02",  8.35, "Y",    "LAST",      8.35,  0,
    "TEST01", "PAT01",  "PARAM02", NA,    "",     "LAST",      8.35, NA,
    "TEST01", "PAT01",  "PARAM02",  8.35, "",     "LAST",      8.35,  0,

    "TEST01", "PAT02",  "PARAM01", 29,    "Y",    "LAST",     29,     0,
    "TEST01", "PAT02",  "PARAM01", 19.7,  "",     "LAST",     29,    -9.3,
    "TEST01", "PAT02",  "PARAM01", 18.01, "",     "LAST",     29,   -10.99,
    "TEST01", "PAT02",  "PARAM02",  8.9,  "Y",    "LAST",      8.9,   0,
    "TEST01", "PAT02",  "PARAM02",  9,    "",     "LAST",      8.9,   0.1,
    "TEST01", "PAT02",  "PARAM02",  5.35, "",     "LAST",      8.9,  -3.55
  )

  expect_equal(derive_var_chg(input)$CHG, expected_output$CHG)
})

test_that("`PCHG` is calculated as `CHG / BASE`", {
  input <- tibble::tribble(
    ~STUDYID, ~USUBJID, ~PARAMCD,  ~AVAL, ~ABLFL, ~BASETYPE, ~BASE,  ~CHG,
    "TEST01", "PAT01",  "PARAM01", 10.12, "Y",    "LAST",     10.12,  0,
    "TEST01", "PAT01",  "PARAM01",  9.7,  "",     "LAST",     10.12, -0.42,
    "TEST01", "PAT01",  "PARAM01", 15.01, "",     "LAST",     10.12,  4.89,
    "TEST01", "PAT01",  "PARAM02",  8.35, "Y",    "LAST",      8.35,  0,
    "TEST01", "PAT01",  "PARAM02", NA,    "",     "LAST",      8.35, NA,
    "TEST01", "PAT01",  "PARAM02",  8.35, "",     "LAST",      8.35,  0,

    "TEST01", "PAT02",  "PARAM01", 29,    "Y",    "LAST",     29,     0,
    "TEST01", "PAT02",  "PARAM01", 19.7,  "",     "LAST",     29,    -9.3,
    "TEST01", "PAT02",  "PARAM01", 18.01, "",     "LAST",     29,   -10.99,
    "TEST01", "PAT02",  "PARAM02",  8.9,  "Y",    "LAST",      8.9,   0,
    "TEST01", "PAT02",  "PARAM02",  9,    "",     "LAST",      8.9,   0.1,
    "TEST01", "PAT02",  "PARAM02",  5.35, "",     "LAST",      8.9,  -3.55
  )

  expected_output <- tibble::tribble(
    ~STUDYID, ~USUBJID, ~PARAMCD,  ~AVAL, ~ABLFL, ~BASETYPE, ~BASE,  ~CHG,  ~PCHG,
    "TEST01", "PAT01",  "PARAM01", 10.12, "Y",    "LAST",     10.12,  0,     0,
    "TEST01", "PAT01",  "PARAM01",  9.7,  "",     "LAST",     10.12, -0.42, -0.04150198,
    "TEST01", "PAT01",  "PARAM01", 15.01, "",     "LAST",     10.12,  4.89,  0.4832016,
    "TEST01", "PAT01",  "PARAM02",  8.35, "Y",    "LAST",      8.35,  0,     0,
    "TEST01", "PAT01",  "PARAM02", NA,    "",     "LAST",      8.35, NA,    NA,
    "TEST01", "PAT01",  "PARAM02",  8.35, "",     "LAST",      8.35,  0,     0,

    "TEST01", "PAT02",  "PARAM01", 29,    "Y",    "LAST",     29,     0,     0,
    "TEST01", "PAT02",  "PARAM01", 19.7,  "",     "LAST",     29,    -9.3,  -0.3206897,
    "TEST01", "PAT02",  "PARAM01", 18.01, "",     "LAST",     29,   -10.99, -0.3789655,
    "TEST01", "PAT02",  "PARAM02",  8.9,  "Y",    "LAST",      8.9,   0,     0,
    "TEST01", "PAT02",  "PARAM02",  9,    "",     "LAST",      8.9,   0.1,   0.01123596,
    "TEST01", "PAT02",  "PARAM02",  5.35, "",     "LAST",      8.9,  -3.55, -0.3988764
  )

  expect_equal(derive_var_pchg(input)$PCHG, expected_output$PCHG, tolerance = 1e-7)
})
