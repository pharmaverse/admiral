context("test-derive_query_vars")

test_that("Derive CQ and SMQ variables with two term levels", {
  queries <- tibble::tribble(
    ~VAR_PREFIX, ~QUERY_NAME, ~QUERY_ID, ~QUERY_SCOPE, ~TERM_LEVEL, ~TERM_NAME,
    "CQ01", "Immune-Mediated Hepatitis (Diagnosis and Lab Abnormalities)", "20000008", "NARROW", "AEDECOD", "ALANINE AMINOTRANSFERASE ABNORMAL",
    "SMQ03", "Immune-Mediated Hypothyroidism", "20000161", "NARROW", "AEDECOD", "BASEDOW'S DISEASE",
    "CQ06", "Some query", "11111", NA_character_, "AELLT", "SOME TERM"
  )

  adae <- tibble::tribble(
    ~USUBJID, ~ASTDTM, ~AETERM, ~AESEQ, ~AEDECOD, ~AELLT,
    "01", "2020-06-02 23:59:59", "ALANINE AMINOTRANSFERASE ABNORMAL", 3, "Alanine aminotransferase abnormal", NA_character_,
    "02", "2020-06-05 23:59:59", "BASEDOW'S DISEASE", 5, "Basedow's disease", NA_character_,
    "03", "2020-06-07 23:59:59", "SOME TERM", 2, "Some query", "Some term"
  )

  expected_output <- tibble::tribble(
    ~USUBJID, ~ASTDTM, ~AETERM, ~AESEQ, ~AEDECOD, ~AELLT, ~CQ01NAM, ~SMQ03NAM, ~SMQ03CD, ~SMQ03SC, ~CQ06NAM,
    "01", "2020-06-02 23:59:59", "ALANINE AMINOTRANSFERASE ABNORMAL", 3, "Alanine aminotransferase abnormal", NA_character_, "Immune-Mediated Hepatitis (Diagnosis and Lab Abnormalities)", NA_character_, NA_character_, NA_character_, NA_character_,
    "02", "2020-06-05 23:59:59", "BASEDOW'S DISEASE", 5, "Basedow's disease", NA_character_, NA_character_, "Immune-Mediated Hypothyroidism", "20000161", "NARROW", NA_character_,
    "03", "2020-06-07 23:59:59", "SOME TERM", 2, NA_character_, "Some term",  NA_character_, NA_character_, NA_character_, NA_character_, "Some query"
  )

  actual_output <- derive_query_vars(adae, queries, c("USUBJID", "ASTDTM", "AETERM", "AESEQ"))

  expect_dfs_equal(expected_output, actual_output, keys = "USUBJID")
})
