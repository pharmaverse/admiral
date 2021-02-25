context("test-derive_var_base")

test_that("`BASE` is set to `AVAL` where `ABL01FL == 'Y'`", {
  input <- tibble::tribble(
    ~STUDYID, ~USUBJID, ~PARAMCD,  ~AVAL, ~ABL01FL,
    "TEST01", "PAT01",  "PARAM01", 10.12, "Y",
    "TEST01", "PAT01",  "PARAM01",  9.7,  "N",
    "TEST01", "PAT01",  "PARAM01", 15.01, "N",
    "TEST01", "PAT01",  "PARAM02",  8.35, "Y",
    "TEST01", "PAT01",  "PARAM02", NA,    "N",
    "TEST01", "PAT01",  "PARAM02",  8.35, "N",

    "TEST01", "PAT02",  "PARAM01", 29,    "Y",
    "TEST01", "PAT02",  "PARAM01", 19.7,  "N",
    "TEST01", "PAT02",  "PARAM01", 18.01, "N",
    "TEST01", "PAT02",  "PARAM02",  8.9,  "Y",
    "TEST01", "PAT02",  "PARAM02",  9,    "N",
    "TEST01", "PAT02",  "PARAM02",  5.35, "N",
  )
  expected_output <- tibble::tribble(
    ~STUDYID, ~USUBJID, ~PARAMCD,  ~AVAL, ~ABL01FL, ~BASE,
    "TEST01", "PAT01",  "PARAM01", 10.12, "Y",      10.12,
    "TEST01", "PAT01",  "PARAM01",  9.7,  "N",      10.12,
    "TEST01", "PAT01",  "PARAM01", 15.01, "N",      10.12,
    "TEST01", "PAT01",  "PARAM02",  8.35, "Y",       8.35,
    "TEST01", "PAT01",  "PARAM02", NA,    "N",       8.35,
    "TEST01", "PAT01",  "PARAM02",  8.35, "N",       8.35,

    "TEST01", "PAT02",  "PARAM01", 29,    "Y",      29,
    "TEST01", "PAT02",  "PARAM01", 19.7,  "N",      29,
    "TEST01", "PAT02",  "PARAM01", 18.01, "N",      29,
    "TEST01", "PAT02",  "PARAM02",  8.9,  "Y",       8.9,
    "TEST01", "PAT02",  "PARAM02",  9,    "N",       8.9,
    "TEST01", "PAT02",  "PARAM02",  5.35, "N",       8.9,
  )

  expect_identical(derive_var_base(input), expected_output)
})

test_that("`derive_var_base()` only adds the `BASE` variable", {
  input <- tibble::tribble(
    ~STUDYID, ~USUBJID, ~PARAMCD,  ~AVAL, ~ABL01FL, ~ANL01FL,
    "TEST01", "PAT01",  "PARAM01", 10.12, "Y",      "Y",
    "TEST01", "PAT01",  "PARAM01",  9.7,  "N",      "Y",
    "TEST01", "PAT01",  "PARAM01", 15.01, "N",      "Y",
    "TEST01", "PAT01",  "PARAM02",  8.35, "Y",      "Y",
    "TEST01", "PAT01",  "PARAM02", NA,    "N",      "Y",
    "TEST01", "PAT01",  "PARAM02",  8.35, "N",      "Y",
  )
  expected_output <- tibble::tribble(
    ~STUDYID, ~USUBJID, ~PARAMCD,  ~AVAL, ~ABL01FL, ~ANL01FL, ~BASE,
    "TEST01", "PAT01",  "PARAM01", 10.12, "Y",      "Y",      10.12,
    "TEST01", "PAT01",  "PARAM01",  9.7,  "N",      "Y",      10.12,
    "TEST01", "PAT01",  "PARAM01", 15.01, "N",      "Y",      10.12,
    "TEST01", "PAT01",  "PARAM02",  8.35, "Y",      "Y",       8.35,
    "TEST01", "PAT01",  "PARAM02", NA,    "N",      "Y",       8.35,
    "TEST01", "PAT01",  "PARAM02",  8.35, "N",      "Y",       8.35
  )

  expect_identical(derive_var_base(input), expected_output)
})

test_that("`derive_var_base()` fails when required variables are missing", {
  input <- tibble::tribble(
    ~STUDYID, ~USUBJID, ~PARAMCD,  ~AVAL, ~ABL01FL,
    "TEST01", "PAT01",  "PARAM01", 10.12, "Y",
    "TEST01", "PAT01",  "PARAM01",  9.7,  "N",
    "TEST01", "PAT01",  "PARAM01", 15.01, "N",
    "TEST01", "PAT01",  "PARAM02",  8.35, "Y",
    "TEST01", "PAT01",  "PARAM02", NA,    "N",
    "TEST01", "PAT01",  "PARAM02",  8.35, "N",

    "TEST01", "PAT02",  "PARAM01", 29,    "Y",
    "TEST01", "PAT02",  "PARAM01", 19.7,  "N",
    "TEST01", "PAT02",  "PARAM01", 18.01, "N",
    "TEST01", "PAT02",  "PARAM02",  8.9,  "Y",
    "TEST01", "PAT02",  "PARAM02",  9,    "N",
    "TEST01", "PAT02",  "PARAM02",  5.35, "N",
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
    input %>% select(-ABL01FL) %>% derive_var_base(),
    "Required variable `ABL01FL` is missing."
  )
  expect_error(
    input %>% select(-c(AVAL, ABL01FL)) %>% derive_var_base(),
    "Required variables `AVAL` and `ABL01FL` are missing."
  )
})
