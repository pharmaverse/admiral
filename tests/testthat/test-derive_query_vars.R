context("test-derive_query_vars")

test_that("Derive CQ and SMQ variables with two term levels", {

  queries <- tibble::tribble(
    ~VAR_PREFIX, ~QUERY_NAME, ~QUERY_ID, ~QUERY_SCOPE, ~QUERY_SCOPE_NUM, ~TERM_LEVEL, ~TERM_NAME,
    "CQ01", "Immune-Mediated Hepatitis (Diagnosis and Lab Abnormalities)", 20000008, "NARROW", 1, "AEDECOD", "ALANINE AMINOTRANSFERASE ABNORMAL",
    "CQ01", "Immune-Mediated Hepatitis (Diagnosis and Lab Abnormalities)", 20000008, "NARROW", 1, "AEDECOD", "AMMONIA ABNORMALL",
    "SMQ03", "Immune-Mediated Hypothyroidism", 20000161, "NARROW", 1, "AEDECOD", "BASEDOW'S DISEASE",
    "SMQ05", "Immune-Mediated Pneumonitis", NA_integer_, "NARROW", 1, "AEDECOD", "ALVEOLAR PROTEINOSIS",
    "CQ06", "Some query", 11111, NA_character_, NA_integer_, "AELLT", "SOME TERM"
  ) %>% dplyr::mutate(
    TERM_ID = as.integer(as.factor(.data$TERM_NAME))
  )

  adae <- tibble::tribble(
    ~USUBJID, ~ASTDTM, ~AETERM, ~AESEQ, ~AEDECOD, ~AELLT,
    "01", "2020-06-02 23:59:59", "ALANINE AMINOTRANSFERASE ABNORMAL", 3, "Alanine aminotransferase abnormal", NA_character_,
    "02", "2020-06-05 23:59:59", "BASEDOW'S DISEASE", 5, "Basedow's disease", NA_character_,
    "03", "2020-06-07 23:59:59", "SOME TERM", 2, "Some query", "Some term",
    "05", "2020-06-09 23:59:59", "ALVEOLAR PROTEINOSIS", 7, "Alveolar proteinosis", NA_character_
  )

  expected_output <- tibble::tribble(
    ~USUBJID, ~ASTDTM, ~AETERM, ~AESEQ, ~AEDECOD, ~AELLT, ~CQ01NAM, ~CQ01CD, ~CQ01SC, ~CQ01SCN, ~SMQ03NAM, ~SMQ03CD, ~SMQ03SC, ~SMQ03SCN, ~SMQ05NAM, ~SMQ05SC, ~SMQ05SCN, ~CQ06NAM, ~CQ06CD,
    "01", "2020-06-02 23:59:59", "ALANINE AMINOTRANSFERASE ABNORMAL", 3, "Alanine aminotransferase abnormal", NA_character_,
    "Immune-Mediated Hepatitis (Diagnosis and Lab Abnormalities)", 20000008, "NARROW", 1, # CQ01
    NA_character_, NA_integer_, NA_character_, NA_integer_, # SMQ03
    NA_character_, NA_character_, NA_integer_, # SMQ05
    NA_character_, NA_integer_, # CQ06
    "02", "2020-06-05 23:59:59", "BASEDOW'S DISEASE", 5, "Basedow's disease", NA_character_,
    NA_character_, NA_integer_, NA_character_, NA_integer_, # CQ01
    "Immune-Mediated Hypothyroidism", 20000161, "NARROW", 1, # SMQ03
    NA_character_, NA_character_, NA_integer_, # SMQ05
    NA_character_, NA_integer_, # CQ06
    "03", "2020-06-07 23:59:59", "SOME TERM", 2, "Some query", "Some term",
    NA_character_, NA_integer_, NA_character_, NA_integer_, # CQ01
    NA_character_, NA_integer_, NA_character_, NA_integer_, # SMQ03
    NA_character_, NA_character_, NA_integer_, # CQ06
    "Some query", 11111,
    "05", "2020-06-09 23:59:59", "ALVEOLAR PROTEINOSIS", 7, "Alveolar proteinosis", NA_character_,
    NA_character_, NA_integer_, NA_character_, NA_integer_, # CQ01
    NA_character_, NA_integer_, NA_character_, NA_integer_, # SMQ03
    "Immune-Mediated Pneumonitis", "NARROW", 1, # SMQ05
    NA_character_, NA_integer_ # CQ06
  )

  actual_output <- derive_query_vars(adae, queries)

  expect_dfs_equal(expected_output, actual_output, keys = "USUBJID")
})

test_that("Derive when dataset does not have a unique key when excluding `TERM_LEVEL` columns", {
  query <- tibble::tribble(
    ~VAR_PREFIX, ~QUERY_NAME, ~TERM_LEVEL, ~TERM_NAME, ~QUERY_ID, ~TERM_ID,
    "CQ42", "My Query", "AEDECOD", "PTSI", 1, NA_real_,
    "CQ42", "My Query", "AELLT", "LLTSI", 1, NA_real_
  )

  my_ae <- tibble::tribble(
    ~USUBJID, ~ASTDY, ~AEDECOD, ~AELLT,
    "1", 1, "PTSI", "other",
    "1", 2, "something", "LLTSI",
    "1", 2, "PTSI", "LLTSI",
    "1", 2, "something", "other"
  )

  expected_output <- tibble::tribble(
    ~USUBJID, ~ASTDY, ~AEDECOD, ~AELLT, ~CQ42NAM, ~CQ42CD,
    "1", 1, "PTSI", "other", "My Query", 1,
    "1", 2, "something", "LLTSI", "My Query", 1,
    "1", 2, "PTSI", "LLTSI", "My Query", 1,
    "1", 2, "something", "other", NA_character_, NA_integer_
  )

  actual_output <- derive_query_vars(my_ae, queries = query)

  expect_equal(expected_output, actual_output)
})

test_that("Derive when an adverse event is in multiple baskets", {
  query <- tibble::tribble(
    ~VAR_PREFIX, ~QUERY_NAME, ~TERM_LEVEL, ~TERM_NAME, ~QUERY_ID, ~TERM_ID,
    "CQ40", "My Query 1", "AEDECOD", "PTSI", 1, NA_real_,
    "CQ42", "My Query 2", "AELLT", "LLTSI", 2, NA_real_
  )

  my_ae <- tibble::tribble(
    ~USUBJID, ~ASTDY, ~AEDECOD, ~AELLT,
    "1", 1, "PTSI", "other",
    "1", 2, "something", "LLTSI",
    "1", 2, "PTSI", "LLTSI",
    "1", 2, "something", "other"
  )

  expected_output <- tibble::tribble(
    ~USUBJID, ~ASTDY, ~AEDECOD, ~AELLT, ~CQ40NAM, ~CQ42NAM, ~CQ40CD, ~CQ42CD,
    "1", 1, "PTSI", "other", "My Query 1", NA_character_, 1, NA_integer_,
    "1", 2, "something", "LLTSI", NA_character_, "My Query 2", NA_integer_, 2,
    "1", 2, "PTSI", "LLTSI", "My Query 1", "My Query 2", 1, 2,
    "1", 2, "something", "other", NA_character_, NA_character_, NA_integer_, NA_integer_
  )

  actual_output <- derive_query_vars(my_ae, queries = query)

  expect_equal(expected_output, actual_output)
})


test_that("Derive when query dataset does not have QUERY_ID or QUERY_SCOPE column", {

  query <- tibble::tribble(
    ~VAR_PREFIX, ~QUERY_NAME, ~TERM_LEVEL, ~TERM_NAME, ~TERM_ID,
    "CQ42", "My Query", "AEDECOD", "PTSI", NA_real_,
    "CQ42", "My Query", "AELLT", "LLTSI", NA_real_
  )

  my_ae <- tibble::tribble(
    ~USUBJID, ~ASTDY, ~AEDECOD, ~AELLT,
    "1", 1, "PTSI", "other",
    "1", 2, "something", "LLTSI",
    "1", 2, "PTSI", "LLTSI",
    "1", 2, "something", "other"
  )

  actual_output <- derive_query_vars(my_ae, queries = query)

  expected_output <- tibble::tribble(
    ~USUBJID, ~ASTDY, ~AEDECOD, ~AELLT, ~CQ42NAM,
    "1", 1, "PTSI", "other", "My Query",
    "1", 2, "something", "LLTSI", "My Query",
    "1", 2, "PTSI", "LLTSI", "My Query",
    "1", 2, "something", "other", NA_character_
  )

  expect_equal(expected_output, actual_output)
})

test_that("Derive decides between TERM_NAME and TERM_ID based on the type of the variable", {

  query <- tibble::tribble(
    ~VAR_PREFIX, ~QUERY_NAME, ~TERM_LEVEL, ~TERM_NAME, ~QUERY_ID, ~TERM_ID,
    "CQ40", "My Query 1", "AEDECOD", "PTSI", 1, NA,
    "CQ42", "My Query 2", "AELLTCD", NA_character_, 2, 1
  )

  my_ae <- tibble::tribble(
    ~USUBJID, ~ASTDY, ~AEDECOD, ~AELLT, ~AELLTCD,
    "1", 1, "PTSI", "other", NA,
    "1", 2, "PTSI", "LLTSI", NA,
    "1", 3, NA, NA, 1
  )

  actual_output <- derive_query_vars(my_ae, queries = query)

  expected_output <- tibble::tribble(
    ~USUBJID, ~ASTDY, ~AEDECOD, ~AELLT, ~AELLTCD, ~CQ40NAM, ~CQ42NAM, ~CQ40CD, ~CQ42CD,
    "1", 1, "PTSI", "other", NA, "My Query 1", NA, 1, NA,
    "1", 2, "PTSI", "LLTSI", NA, "My Query 1", NA, 1, NA,
    "1", 3, NA, NA, 1, NA, "My Query 2", NA, 2
  )

  expect_equal(expected_output, actual_output)

  expect_error(
    derive_query_vars(mutate(my_ae, AELLTCD = as.logical(AELLTCD)), query),
    regexp = ".* is of type logical, numeric or character is required"
  )

})


test_that("assert_valid_queries checks VAR_PREFIX values", {

  query <- tibble::tribble(
    ~VAR_PREFIX, ~QUERY_NAME, ~TERM_LEVEL, ~TERM_NAME, ~QUERY_ID, ~TERM_ID,
    "CQ40", "My Query 1", "AEDECOD", "PTSI", 1, NA,
    "CQ42", "My Query 2", "AELLTCD", NA_character_, 2, 1
  )

  expect_error(
    assert_valid_queries(
      mutate(query, VAR_PREFIX = c("30", "55")),
      "test"
    ),
    regexp = "`VAR_PREFIX` in `test` must start with 2-3 letters.. Problem with `30` and `55`."
  )

  expect_error(
    assert_valid_queries(
      mutate(query, VAR_PREFIX = c("AA", "BB")),
      "test"
    ),
    regexp = "`VAR_PREFIX` in `test` must end with 2-digit numbers. Issue with `AA` and `BB`."
  )

  expect_error(
    assert_valid_queries(
      mutate(query, QUERY_NAME = c("", "A")),
      "test"
    ),
    regexp = "`QUERY_NAME` in `test` cannot be empty string or NA."
  )

  expect_error(
    assert_valid_queries(
      mutate(query, QUERY_SCOPE = letters[1:2]),
      "test"
    ),
    regexp = "`QUERY_SCOPE` in `test` can only be 'BROAD', 'NARROW' or `NA`."
  )

  expect_error(
    assert_valid_queries(
      mutate(query, QUERY_SCOPE_NUM = 10:11),
      "test"
    ),
    regexp = "`QUERY_SCOPE_NUM` in `test` must be one of 1, 2, or NA. Issue with `10` and `11`."
  )

  expect_error(
    assert_valid_queries(
      mutate(query, TERM_NAME = c(NA, NA)),
      "test"
    ),
    regexp = paste(
      "Either `TERM_NAME` or `TERM_ID` need to be specified in `test`.",
      "They both cannot be NA or empty."
    )
  )

  expect_error(
    assert_valid_queries(
      mutate(query, VAR_PREFIX = c("CQ40", "CQ40")),
      "test"
    ),
    regexp = "In `test`, `QUERY_NAME` of 'CQ40' is not unique."
  )

  expect_error(
    assert_valid_queries(
      mutate(
        query,
        VAR_PREFIX = c("CQ40", "CQ40"),
        QUERY_NAME = c("My Query 1", "My Query 1")
      ),
      "test"
    ),
    regexp = "In `test`, `QUERY_ID` of 'CQ40' is not unique."
  )

})
