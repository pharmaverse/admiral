## Test 1: `target` is set to `source` where `ABLFL == 'Y'` ----
test_that("derive_var_base Test 1: `target` is set to `source` where `ABLFL == 'Y'`", {
  input <- tibble::tribble(
    ~STUDYID, ~USUBJID, ~PARAMCD, ~ASEQ, ~AVAL, ~ABLFL, ~BASETYPE,
    "TEST01", "PAT01", "PARAM01", 1, 10.12, "Y", "LAST",
    "TEST01", "PAT01", "PARAM01", 2, 9.7, NA_character_, "LAST",
    "TEST01", "PAT01", "PARAM01", 3, 15.01, NA_character_, "LAST",
    "TEST01", "PAT01", "PARAM02", 1, 8.35, "Y", "LAST",
    "TEST01", "PAT01", "PARAM02", 2, NA, NA_character_, "LAST",
    "TEST01", "PAT01", "PARAM02", 3, 8.35, NA_character_, "LAST",
    "TEST01", "PAT02", "PARAM01", 1, 29, "Y", "LAST",
    "TEST01", "PAT02", "PARAM01", 2, 19.7, NA_character_, "LAST",
    "TEST01", "PAT02", "PARAM01", 3, 18.01, NA_character_, "LAST",
    "TEST01", "PAT02", "PARAM02", 1, 8.9, "Y", "LAST",
    "TEST01", "PAT02", "PARAM02", 2, 9, NA_character_, "LAST",
    "TEST01", "PAT02", "PARAM02", 3, 5.35, NA_character_, "LAST"
  )
  expected_output <- tibble::tribble(
    ~STUDYID, ~USUBJID, ~PARAMCD, ~ASEQ, ~AVAL, ~ABLFL, ~BASETYPE, ~BASE,
    "TEST01", "PAT01", "PARAM01", 1, 10.12, "Y", "LAST", 10.12,
    "TEST01", "PAT01", "PARAM01", 2, 9.7, NA_character_, "LAST", 10.12,
    "TEST01", "PAT01", "PARAM01", 3, 15.01, NA_character_, "LAST", 10.12,
    "TEST01", "PAT01", "PARAM02", 1, 8.35, "Y", "LAST", 8.35,
    "TEST01", "PAT01", "PARAM02", 2, NA, NA_character_, "LAST", 8.35,
    "TEST01", "PAT01", "PARAM02", 3, 8.35, NA_character_, "LAST", 8.35,
    "TEST01", "PAT02", "PARAM01", 1, 29, "Y", "LAST", 29,
    "TEST01", "PAT02", "PARAM01", 2, 19.7, NA_character_, "LAST", 29,
    "TEST01", "PAT02", "PARAM01", 3, 18.01, NA_character_, "LAST", 29,
    "TEST01", "PAT02", "PARAM02", 1, 8.9, "Y", "LAST", 8.9,
    "TEST01", "PAT02", "PARAM02", 2, 9, NA_character_, "LAST", 8.9,
    "TEST01", "PAT02", "PARAM02", 3, 5.35, NA_character_, "LAST", 8.9
  )
  actual_output <- derive_var_base(
    input,
    by_vars = exprs(USUBJID, PARAMCD, BASETYPE),
    source_var = AVAL,
    new_var = BASE
  )

  expect_dfs_equal(
    expected_output,
    actual_output,
    keys = c("STUDYID", "USUBJID", "PARAMCD", "ASEQ")
  )
})

## Test 2: `target` is set to `NA` if a baseline record is missing ----
test_that("derive_var_base Test 2: `target` is set to `NA` if a baseline record is missing", {
  input <- tibble::tribble(
    ~STUDYID, ~USUBJID, ~PARAMCD, ~ASEQ, ~AVAL, ~ABLFL,        ~BASETYPE,
    "TEST01", "PAT01", "PARAM01",     1, 10.12, "Y",           "LAST",
    "TEST01", "PAT01", "PARAM01",     2,   9.7, NA_character_, "LAST",
    "TEST01", "PAT01", "PARAM01",     3, 15.01, NA_character_, "LAST",
    "TEST01", "PAT01", "PARAM02",     1,   4.9, NA_character_, "LAST",
    "TEST01", "PAT01", "PARAM02",     2,   7.1, NA_character_, "LAST",
    "TEST01", "PAT01", "PARAM02",     3,  8.35, NA_character_, "LAST"
  )
  expected_output <- tibble::tribble(
    ~STUDYID, ~USUBJID, ~PARAMCD, ~ASEQ, ~AVAL, ~ABLFL,        ~BASETYPE, ~BASE,
    "TEST01", "PAT01", "PARAM01",     1, 10.12, "Y",           "LAST",    10.12,
    "TEST01", "PAT01", "PARAM01",     2,   9.7, NA_character_, "LAST",    10.12,
    "TEST01", "PAT01", "PARAM01",     3, 15.01, NA_character_, "LAST",    10.12,
    "TEST01", "PAT01", "PARAM02",     1,   4.9, NA_character_, "LAST",       NA,
    "TEST01", "PAT01", "PARAM02",     2,   7.1, NA_character_, "LAST",       NA,
    "TEST01", "PAT01", "PARAM02",     3,  8.35, NA_character_, "LAST",       NA
  )
  actual_output <- derive_var_base(
    input,
    by_vars = exprs(USUBJID, PARAMCD, BASETYPE),
    source_var = AVAL,
    new_var = BASE
  )

  expect_dfs_equal(
    actual_output,
    expected_output,
    keys = c("STUDYID", "USUBJID", "PARAMCD", "ASEQ")
  )
})

## Test 3: only the `target` variable is added to the input dataset ----
test_that("derive_var_base Test 3: only the `target` variable is added to the input dataset", {
  input <- tibble::tribble(
    ~STUDYID, ~USUBJID, ~PARAMCD, ~ASEQ, ~AVAL, ~ABLFL,     ~BASETYPE, ~ANL01FL,
    "TEST01", "PAT01", "PARAM01",     1, 10.12, "Y",           "LAST", "Y",
    "TEST01", "PAT01", "PARAM01",     2,   9.7, NA_character_, "LAST", "Y",
    "TEST01", "PAT01", "PARAM01",     3, 15.01, NA_character_, "LAST", "Y",
    "TEST01", "PAT01", "PARAM02",     1,  8.35, "Y",           "LAST", "Y",
    "TEST01", "PAT01", "PARAM02",     2,    NA, NA_character_, "LAST", "Y",
    "TEST01", "PAT01", "PARAM02",     3,  8.35, NA_character_, "LAST", "Y"
  )
  expected_output <- tibble::tribble(
    ~STUDYID, ~USUBJID, ~PARAMCD, ~ASEQ, ~AVAL, ~ABLFL,        ~BASETYPE, ~ANL01FL, ~BASE,
    "TEST01", "PAT01", "PARAM01",     1, 10.12, "Y",           "LAST",    "Y",      10.12,
    "TEST01", "PAT01", "PARAM01",     2,   9.7, NA_character_, "LAST",    "Y",      10.12,
    "TEST01", "PAT01", "PARAM01",     3, 15.01, NA_character_, "LAST",    "Y",      10.12,
    "TEST01", "PAT01", "PARAM02",     1,  8.35, "Y",           "LAST",    "Y",       8.35,
    "TEST01", "PAT01", "PARAM02",     2,    NA, NA_character_, "LAST",    "Y",       8.35,
    "TEST01", "PAT01", "PARAM02",     3,  8.35, NA_character_, "LAST",    "Y",       8.35
  )
  actual_output <- derive_var_base(
    input,
    by_vars = exprs(USUBJID, PARAMCD, BASETYPE),
    source_var = AVAL,
    new_var = BASE
  )

  expect_dfs_equal(
    actual_output,
    expected_output,
    keys = c("STUDYID", "USUBJID", "PARAMCD", "ASEQ")
  )
})

## Test 4: error if multiple baseline records ----
test_that("derive_var_base Test 4: error if multiple baseline records", {
  input <- tibble::tribble(
    ~STUDYID, ~USUBJID, ~PARAMCD, ~AVALC,   ~ABLFL,        ~BASETYPE,
    "TEST01", "PAT01", "PARAM01", "LOW",    "Y",           "LAST",
    "TEST01", "PAT01", "PARAM01", "MEDIUM", "Y",           "LAST",
    "TEST01", "PAT01", "PARAM01", "LOW",    NA_character_, "LAST",
    "TEST01", "PAT01", "PARAM01", "MEDIUM", NA_character_, "LAST",
    "TEST01", "PAT02", "PARAM02", "HIGH",   "Y",           "LAST",
    "TEST01", "PAT02", "PARAM02", "HIGH",   "Y",           "LAST",
    "TEST01", "PAT02", "PARAM02", "MEDIUM", NA_character_, "LAST",
  )

  expect_snapshot(
    derive_var_base(
      input,
      by_vars = exprs(USUBJID, PARAMCD, BASETYPE),
      source_var = AVALC,
      new_var = BASEC
    ),
    error = TRUE
  )
})
