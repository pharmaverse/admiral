context("test_derive_baseline")

test_that("`target` is set to `source` where `ABLFL == 'Y'`", {
  input <- tibble::tribble(
    ~STUDYID, ~USUBJID, ~PARAMCD,  ~ASEQ, ~AVAL, ~ABLFL, ~BASETYPE,
    "TEST01", "PAT01",  "PARAM01", 1,     10.12, "Y",    "LAST",
    "TEST01", "PAT01",  "PARAM01", 2,      9.7,  "",     "LAST",
    "TEST01", "PAT01",  "PARAM01", 3,     15.01, "",     "LAST",
    "TEST01", "PAT01",  "PARAM02", 1,      8.35, "Y",    "LAST",
    "TEST01", "PAT01",  "PARAM02", 2,     NA,    "",     "LAST",
    "TEST01", "PAT01",  "PARAM02", 3,      8.35, "",     "LAST",

    "TEST01", "PAT02",  "PARAM01", 1,     29,    "Y",    "LAST",
    "TEST01", "PAT02",  "PARAM01", 2,     19.7,  "",     "LAST",
    "TEST01", "PAT02",  "PARAM01", 3,     18.01, "",     "LAST",
    "TEST01", "PAT02",  "PARAM02", 1,      8.9,  "Y",    "LAST",
    "TEST01", "PAT02",  "PARAM02", 2,      9,    "",     "LAST",
    "TEST01", "PAT02",  "PARAM02", 3,      5.35, "",     "LAST"
  )
  expected_output <- tibble::tribble(
    ~STUDYID, ~USUBJID, ~PARAMCD,  ~ASEQ, ~AVAL, ~ABLFL, ~BASETYPE, ~BASE,
    "TEST01", "PAT01",  "PARAM01", 1,     10.12, "Y",    "LAST",     10.12,
    "TEST01", "PAT01",  "PARAM01", 2,      9.7,  "",     "LAST",     10.12,
    "TEST01", "PAT01",  "PARAM01", 3,     15.01, "",     "LAST",     10.12,
    "TEST01", "PAT01",  "PARAM02", 1,      8.35, "Y",    "LAST",      8.35,
    "TEST01", "PAT01",  "PARAM02", 2,     NA,    "",     "LAST",      8.35,
    "TEST01", "PAT01",  "PARAM02", 3,      8.35, "",     "LAST",      8.35,

    "TEST01", "PAT02",  "PARAM01", 1,     29,    "Y",    "LAST",     29,
    "TEST01", "PAT02",  "PARAM01", 2,     19.7,  "",     "LAST",     29,
    "TEST01", "PAT02",  "PARAM01", 3,     18.01, "",     "LAST",     29,
    "TEST01", "PAT02",  "PARAM02", 1,      8.9,  "Y",    "LAST",      8.9,
    "TEST01", "PAT02",  "PARAM02", 2,      9,    "",     "LAST",      8.9,
    "TEST01", "PAT02",  "PARAM02", 3,      5.35, "",     "LAST",      8.9
  )
  actual_output <- derive_baseline(
    input,
    by_vars = c("USUBJID", "PARAMCD", "BASETYPE"),
    source_var = AVAL,
    new_var = BASE
  )


  expect_dfs_equal(
    expected_output,
    actual_output,
    keys = c("STUDYID", "USUBJID", "PARAMCD", "ASEQ")
  )
})

test_that("`target` is set to `NA` if a baseline record is missing", {
  input <- tibble::tribble(
    ~STUDYID, ~USUBJID, ~PARAMCD,  ~ASEQ, ~AVAL, ~ABLFL, ~BASETYPE,
    "TEST01", "PAT01",  "PARAM01", 1,     10.12, "Y",    "LAST",
    "TEST01", "PAT01",  "PARAM01", 2,      9.7,  "",     "LAST",
    "TEST01", "PAT01",  "PARAM01", 3,     15.01, "",     "LAST",
    "TEST01", "PAT01",  "PARAM02", 1,      4.9,  "",     "LAST",
    "TEST01", "PAT01",  "PARAM02", 2,      7.1,  "",     "LAST",
    "TEST01", "PAT01",  "PARAM02", 3,      8.35, "",     "LAST"
  )
  expected_output <- tibble::tribble(
    ~STUDYID, ~USUBJID, ~PARAMCD,  ~ASEQ, ~AVAL, ~ABLFL, ~BASETYPE, ~BASE,
    "TEST01", "PAT01",  "PARAM01", 1,     10.12, "Y",    "LAST",    10.12,
    "TEST01", "PAT01",  "PARAM01", 2,      9.7,  "",     "LAST",    10.12,
    "TEST01", "PAT01",  "PARAM01", 3,     15.01, "",     "LAST",    10.12,
    "TEST01", "PAT01",  "PARAM02", 1,      4.9,  "",     "LAST",    NA,
    "TEST01", "PAT01",  "PARAM02", 2,      7.1,  "",     "LAST",    NA,
    "TEST01", "PAT01",  "PARAM02", 3,      8.35, "",     "LAST",    NA
  )
  actual_output <- derive_baseline(
    input,
    by_vars = c("USUBJID", "PARAMCD", "BASETYPE"),
    source_var = AVAL,
    new_var = BASE
  )

  expect_dfs_equal(
    actual_output,
    expected_output,
    keys = c("STUDYID", "USUBJID", "PARAMCD", "ASEQ")
  )
})

test_that("only the `target` variable is added to the input dataset", {
  input <- tibble::tribble(
    ~STUDYID, ~USUBJID, ~PARAMCD,  ~ASEQ, ~AVAL, ~ABLFL, ~BASETYPE, ~ANL01FL,
    "TEST01", "PAT01",  "PARAM01", 1,     10.12, "Y",    "LAST",    "Y",
    "TEST01", "PAT01",  "PARAM01", 2,      9.7,  "",     "LAST",    "Y",
    "TEST01", "PAT01",  "PARAM01", 3,     15.01, "",     "LAST",    "Y",
    "TEST01", "PAT01",  "PARAM02", 1,      8.35, "Y",    "LAST",    "Y",
    "TEST01", "PAT01",  "PARAM02", 2,     NA,    "",     "LAST",    "Y",
    "TEST01", "PAT01",  "PARAM02", 3,      8.35, "",     "LAST",    "Y"
  )
  expected_output <- tibble::tribble(
    ~STUDYID, ~USUBJID, ~PARAMCD,  ~ASEQ, ~AVAL, ~ABLFL, ~BASETYPE, ~ANL01FL, ~BASE,
    "TEST01", "PAT01",  "PARAM01", 1,     10.12, "Y",    "LAST",    "Y",      10.12,
    "TEST01", "PAT01",  "PARAM01", 2,      9.7,  "",     "LAST",    "Y",      10.12,
    "TEST01", "PAT01",  "PARAM01", 3,     15.01, "",     "LAST",    "Y",      10.12,
    "TEST01", "PAT01",  "PARAM02", 1,      8.35, "Y",    "LAST",    "Y",       8.35,
    "TEST01", "PAT01",  "PARAM02", 2,     NA,    "",     "LAST",    "Y",       8.35,
    "TEST01", "PAT01",  "PARAM02", 3,      8.35, "",     "LAST",    "Y",       8.35
  )
  actual_output <- derive_baseline(
    input,
    by_vars = c("USUBJID", "PARAMCD", "BASETYPE"),
    source_var = AVAL,
    new_var = BASE
  )

  expect_dfs_equal(
    actual_output,
    expected_output,
    keys = c("STUDYID", "USUBJID", "PARAMCD", "ASEQ")
  )
})

test_that("An error is thrown if a subject has multiple records per `PARAMCD` and `BASETYPE`", {
  input <- tibble::tribble(
    ~STUDYID, ~USUBJID, ~PARAMCD,  ~AVALC,   ~ABLFL, ~BASETYPE,
    "TEST01", "PAT01",  "PARAM01", "LOW",    "Y",    "LAST",
    "TEST01", "PAT01",  "PARAM01", "MEDIUM", "Y",    "LAST",
    "TEST01", "PAT01",  "PARAM01", "LOW",    "",     "LAST",
    "TEST01", "PAT01",  "PARAM01", "MEDIUM", "",     "LAST",
    "TEST01", "PAT02",  "PARAM02", "HIGH",   "Y",    "LAST",
    "TEST01", "PAT02",  "PARAM02", "HIGH",   "Y",    "LAST",
    "TEST01", "PAT02",  "PARAM02", "MEDIUM", "",     "LAST",
  )

  expect_error(
    derive_baseline(
      input,
      by_vars = c("USUBJID", "PARAMCD", "BASETYPE"),
      source_var = AVALC,
      new_var = BASEC
    ),
    "Dataset contains multiple baseline records."
  )
})

test_that("a `BASEC` column of type `character` is added to the input dataset", {
  input <- tibble::tribble(
    ~STUDYID, ~USUBJID, ~PARAMCD,  ~AVALC,   ~ABLFL, ~BASETYPE,
    "TEST01", "PAT01",  "PARAM01", "LOW",    "Y",    "LAST",
    "TEST01", "PAT01",  "PARAM01", "LOW",    "",     "LAST",
    "TEST01", "PAT01",  "PARAM01", "MEDIUM", "",     "LAST",
    "TEST01", "PAT01",  "PARAM02", "HIGH",   "Y",    "LAST",
    "TEST01", "PAT01",  "PARAM02", "HIGH",   "",     "LAST",
    "TEST01", "PAT01",  "PARAM02", "MEDIUM", "",     "LAST",
  )
  output <- derive_var_basec(input)

  expect_true("BASEC" %in% colnames(output))
  expect_true(is.character(output$BASEC))
})

test_that("a `BASE` column of type `numeric` is added to the input dataset", {
  input <- tibble::tribble(
    ~STUDYID, ~USUBJID, ~PARAMCD,  ~AVAL,   ~ABLFL, ~BASETYPE,
    "TEST01", "PAT01",  "PARAM01",  1,      "Y",    "LAST",
    "TEST01", "PAT01",  "PARAM01",  1,      "",     "LAST",
    "TEST01", "PAT01",  "PARAM01",  2.5,    "",     "LAST",
    "TEST01", "PAT01",  "PARAM02", 32.1,    "Y",    "LAST",
    "TEST01", "PAT01",  "PARAM02", 35.9,    "",     "LAST",
    "TEST01", "PAT01",  "PARAM02", 34.45,   "",     "LAST",

  )
  output <- derive_var_base(input)

  expect_true("BASE" %in% colnames(output))
  expect_true(is.numeric(output$BASE))
})
