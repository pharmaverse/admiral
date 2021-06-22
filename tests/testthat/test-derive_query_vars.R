context("test-derive_query_vars")

test_that("Derive CQ and SMQ variables with two term levels", {
  queries <- tibble::tribble(
    ~VAR_PREFIX, ~QUERY_NAME, ~QUERY_ID, ~QUERY_SCOPE, ~QUERY_SCOPE_NUM, ~TERM_LEVEL, ~TERM_NAME,
    "CQ01", "Immune-Mediated Hepatitis (Diagnosis and Lab Abnormalities)", 20000008, "NARROW", 1, "AEDECOD", "ALANINE AMINOTRANSFERASE ABNORMAL",
    "CQ01", "Immune-Mediated Hepatitis (Diagnosis and Lab Abnormalities)", 20000008, "NARROW", 1, "AEDECOD", "AMMONIA ABNORMALL",
    "SMQ03", "Immune-Mediated Hypothyroidism", 20000161, "NARROW", 1, "AEDECOD", "BASEDOW'S DISEASE",
    "SMQ05", "Immune-Mediated Pneumonitis", NA_integer_, "NARROW", 1, "AEDECOD", "ALVEOLAR PROTEINOSIS",
    "CQ06", "Some query", 11111, NA_character_, NA_integer_, "AELLT", "SOME TERM"
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
  query <- tribble(
    ~VAR_PREFIX, ~QUERY_NAME, ~TERM_LEVEL, ~TERM_NAME, ~QUERY_ID,
    "CQ42", "My Query", "AEDECOD", "PTSI", 1,
    "CQ42", "My Query", "AELLT", "LLTSI", 1
  )

  my_ae <- tribble(
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
