context("test-derive_var_base")

test_that("`BASE` is set to `AVAL` where `ABLFL == 'Y'`", {
  input <- tibble::tribble(
    ~STUDYID, ~USUBJID, ~PARAMCD,  ~AVAL, ~ABLFL,
    "TEST01", "PAT01",  "PARAM01", 10.12, "Y",
    "TEST01", "PAT01",  "PARAM01",  9.7,  "",
    "TEST01", "PAT01",  "PARAM01", 15.01, "",
    "TEST01", "PAT01",  "PARAM02",  8.35, "Y",
    "TEST01", "PAT01",  "PARAM02", NA,    "",
    "TEST01", "PAT01",  "PARAM02",  8.35, "",

    "TEST01", "PAT02",  "PARAM01", 29,    "Y",
    "TEST01", "PAT02",  "PARAM01", 19.7,  "",
    "TEST01", "PAT02",  "PARAM01", 18.01, "",
    "TEST01", "PAT02",  "PARAM02",  8.9,  "Y",
    "TEST01", "PAT02",  "PARAM02",  9,    "",
    "TEST01", "PAT02",  "PARAM02",  5.35, "",
  )
  expected_output <- tibble::tribble(
    ~STUDYID, ~USUBJID, ~PARAMCD,  ~AVAL, ~ABLFL, ~BASE,
    "TEST01", "PAT01",  "PARAM01", 10.12, "Y",      10.12,
    "TEST01", "PAT01",  "PARAM01",  9.7,  "",       10.12,
    "TEST01", "PAT01",  "PARAM01", 15.01, "",       10.12,
    "TEST01", "PAT01",  "PARAM02",  8.35, "Y",       8.35,
    "TEST01", "PAT01",  "PARAM02", NA,    "",        8.35,
    "TEST01", "PAT01",  "PARAM02",  8.35, "",        8.35,

    "TEST01", "PAT02",  "PARAM01", 29,    "Y",      29,
    "TEST01", "PAT02",  "PARAM01", 19.7,  "",       29,
    "TEST01", "PAT02",  "PARAM01", 18.01, "",       29,
    "TEST01", "PAT02",  "PARAM02",  8.9,  "Y",       8.9,
    "TEST01", "PAT02",  "PARAM02",  9,    "",        8.9,
    "TEST01", "PAT02",  "PARAM02",  5.35, "",        8.9,
  )

  expect_identical(derive_var_base(input), expected_output)
})

test_that("`BASE` is set to `NA` if a baseline record is missing", {
  input <- tibble::tribble(
    ~STUDYID, ~USUBJID, ~PARAMCD,  ~AVAL, ~ABLFL,
    "TEST01", "PAT01",  "PARAM01", 10.12, "Y",
    "TEST01", "PAT01",  "PARAM01",  9.7,  "",
    "TEST01", "PAT01",  "PARAM01", 15.01, "",
    "TEST01", "PAT01",  "PARAM02",  4.9,  "",
    "TEST01", "PAT01",  "PARAM02",  7.1,  "",
    "TEST01", "PAT01",  "PARAM02",  8.35, "",
  )
  expected_output <- tibble::tribble(
    ~STUDYID, ~USUBJID, ~PARAMCD,  ~AVAL, ~ABLFL, ~BASE,
    "TEST01", "PAT01",  "PARAM01", 10.12, "Y",      10.12,
    "TEST01", "PAT01",  "PARAM01",  9.7,  "",       10.12,
    "TEST01", "PAT01",  "PARAM01", 15.01, "",       10.12,
    "TEST01", "PAT01",  "PARAM02",  4.9,  "",       NA,
    "TEST01", "PAT01",  "PARAM02",  7.1,  "",       NA,
    "TEST01", "PAT01",  "PARAM02",  8.35, "",       NA,
  )

  expect_identical(derive_var_base(input), expected_output)
})

test_that("`derive_var_base()` only adds the `BASE` variable", {
  input <- tibble::tribble(
    ~STUDYID, ~USUBJID, ~PARAMCD,  ~AVAL, ~ABLFL, ~ANL01FL,
    "TEST01", "PAT01",  "PARAM01", 10.12, "Y",      "Y",
    "TEST01", "PAT01",  "PARAM01",  9.7,  "",       "Y",
    "TEST01", "PAT01",  "PARAM01", 15.01, "",       "Y",
    "TEST01", "PAT01",  "PARAM02",  8.35, "Y",      "Y",
    "TEST01", "PAT01",  "PARAM02", NA,    "",       "Y",
    "TEST01", "PAT01",  "PARAM02",  8.35, "",       "Y",
  )
  expected_output <- tibble::tribble(
    ~STUDYID, ~USUBJID, ~PARAMCD,  ~AVAL, ~ABLFL, ~ANL01FL, ~BASE,
    "TEST01", "PAT01",  "PARAM01", 10.12, "Y",      "Y",      10.12,
    "TEST01", "PAT01",  "PARAM01",  9.7,  "",       "Y",      10.12,
    "TEST01", "PAT01",  "PARAM01", 15.01, "",       "Y",      10.12,
    "TEST01", "PAT01",  "PARAM02",  8.35, "Y",      "Y",       8.35,
    "TEST01", "PAT01",  "PARAM02", NA,    "",       "Y",       8.35,
    "TEST01", "PAT01",  "PARAM02",  8.35, "",       "Y",       8.35
  )

  expect_identical(derive_var_base(input), expected_output)
})

test_that("`derive_var_base()` fails when required variables are missing", {
  input <- tibble::tribble(
    ~STUDYID, ~USUBJID, ~PARAMCD,  ~AVAL, ~ABLFL,
    "TEST01", "PAT01",  "PARAM01", 10.12, "Y",
    "TEST01", "PAT01",  "PARAM01",  9.7,  "",
    "TEST01", "PAT01",  "PARAM01", 15.01, "",
    "TEST01", "PAT01",  "PARAM02",  8.35, "Y",
    "TEST01", "PAT01",  "PARAM02", NA,    "",
    "TEST01", "PAT01",  "PARAM02",  8.35, "",

    "TEST01", "PAT02",  "PARAM01", 29,    "Y",
    "TEST01", "PAT02",  "PARAM01", 19.7,  "",
    "TEST01", "PAT02",  "PARAM01", 18.01, "",
    "TEST01", "PAT02",  "PARAM02",  8.9,  "Y",
    "TEST01", "PAT02",  "PARAM02",  9,    "",
    "TEST01", "PAT02",  "PARAM02",  5.35, "",
  )

  expect_error(
    input %>% select(-STUDYID) %>% derive_var_base(),
    "Required variable `STUDYID` is missing."
  )
  expect_error(
    input %>% select(-USUBJID) %>% derive_var_base(),
    "Required variable `USUBJID` is missing."
  )
  expect_error(
    input %>% select(-PARAMCD) %>% derive_var_base(),
    "Required variable `PARAMCD` is missing."
  )
  expect_error(
    input %>% select(-AVAL) %>% derive_var_base(),
    "Required variable `AVAL` is missing."
  )
  expect_error(
    input %>% select(-ABLFL) %>% derive_var_base(),
    "Required variable `ABLFL` is missing."
  )
  expect_error(
    input %>% select(-c(AVAL, ABLFL)) %>% derive_var_base(),
    "Required variables `AVAL` and `ABLFL` are missing."
  )
})

test_that("`BASEC` is set to `AVALC` where `ABLFL == 'Y'`", {
  input <- tibble::tribble(
    ~STUDYID, ~USUBJID, ~PARAMCD,  ~AVALC,   ~ABLFL,
    "TEST01", "PAT01",  "PARAM01", "LOW",    "Y",
    "TEST01", "PAT01",  "PARAM01", "LOW",    "",
    "TEST01", "PAT01",  "PARAM01", "MEDIUM", "",
    "TEST01", "PAT01",  "PARAM02", "HIGH",   "Y",
    "TEST01", "PAT01",  "PARAM02", "HIGH",   "",
    "TEST01", "PAT01",  "PARAM02", "MEDIUM", "",

    "TEST01", "PAT02",  "PARAM01", "MEDIUM", "Y",
    "TEST01", "PAT02",  "PARAM01", "MEDIUM", "",
    "TEST01", "PAT02",  "PARAM01", "MEDIUM", "",
    "TEST01", "PAT02",  "PARAM02", "LOW",    "Y",
    "TEST01", "PAT02",  "PARAM02", "LOW",    "",
    "TEST01", "PAT02",  "PARAM02", "HIGH",   "",
  )
  expected_output <- tibble::tribble(
    ~STUDYID, ~USUBJID, ~PARAMCD,  ~AVALC,  ~ABLFL, ~BASEC,
    "TEST01", "PAT01",  "PARAM01", "LOW",    "Y",     "LOW",
    "TEST01", "PAT01",  "PARAM01", "LOW",    "",      "LOW",
    "TEST01", "PAT01",  "PARAM01", "MEDIUM", "",      "LOW",
    "TEST01", "PAT01",  "PARAM02", "HIGH",   "Y",     "HIGH",
    "TEST01", "PAT01",  "PARAM02", "HIGH",   "",      "HIGH",
    "TEST01", "PAT01",  "PARAM02", "MEDIUM", "",      "HIGH",

    "TEST01", "PAT02",  "PARAM01", "MEDIUM", "Y",     "MEDIUM",
    "TEST01", "PAT02",  "PARAM01", "MEDIUM", "",      "MEDIUM",
    "TEST01", "PAT02",  "PARAM01", "MEDIUM", "",      "MEDIUM",
    "TEST01", "PAT02",  "PARAM02", "LOW",    "Y",     "LOW",
    "TEST01", "PAT02",  "PARAM02", "LOW",    "",      "LOW",
    "TEST01", "PAT02",  "PARAM02", "HIGH",   "",      "LOW"
  )

  expect_identical(derive_var_basec(input), expected_output)
})

