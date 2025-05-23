# derive_vars_query ----
## Test 1: Derive CQ and SMQ variables with two term levels ----
test_that("derive_vars_query Test 1: Derive CQ and SMQ variables with two term levels", {
  # nolint start
  queries <- tibble::tribble(
    ~PREFIX, ~GRPNAME, ~GRPID, ~SCOPE, ~SCOPEN, ~SRCVAR, ~TERMCHAR,
    "CQ01", "Immune-Mediated Hepatitis (Diagnosis and Lab Abnormalities)", 20000008, "NARROW", 1, "AEDECOD", "Alanine AMINOTRANSFERASE ABNORMAL",
    "CQ01", "Immune-Mediated Hepatitis (Diagnosis and Lab Abnormalities)", 20000008, "NARROW", 1, "AEDECOD", "AMMONIA ABNORMALL",
    "SMQ03", "Immune-Mediated Hypothyroidism", 20000161, "NARROW", 1, "AEDECOD", "BASEDOW'S DISEASE",
    "SMQ05", "Immune-Mediated Pneumonitis", NA, "NARROW", 1, "AEDECOD", "ALVEOLAR PROTEINOSIS",
    "CQ06", "Some query", 11111, NA, NA, "AELLT", "SOME TERM"
  ) %>% dplyr::mutate(
    TERMNUM = as.integer(as.factor(.data$TERMCHAR))
  )

  adae <- tibble::tribble(
    ~USUBJID,           ~ASTDTM,                             ~AETERM, ~AESEQ,                            ~AEDECOD,      ~AELLT,
    "01", "2020-06-02 23:59:59", "ALANINE AMINOTRANSFERASE ABNORMAL",      3, "ALANINE aminotransferase abnormal",          NA,
    "02", "2020-06-05 23:59:59",                 "BASEDOW'S DISEASE",      5,                 "Basedow's disease",          NA,
    "03", "2020-06-07 23:59:59",                         "SOME TERM",      2,                        "Some query", "Some term",
    "05", "2020-06-09 23:59:59",              "ALVEOLAR PROTEINOSIS",      7,              "Alveolar proteinosis",          NA
  )

  expected_output <- tibble::tribble(
    ~USUBJID,           ~ASTDTM,                             ~AETERM, ~AESEQ,                            ~AEDECOD,      ~AELLT,                        ~SMQ03NAM, ~SMQ03CD, ~SMQ03SC, ~SMQ03SCN,                     ~SMQ05NAM, ~SMQ05SC, ~SMQ05SCN,                                                      ~CQ01NAM,  ~CQ01CD,  ~CQ01SC, ~CQ01SCN,     ~CQ06NAM, ~CQ06CD,
    "01", "2020-06-02 23:59:59", "ALANINE AMINOTRANSFERASE ABNORMAL",      3, "ALANINE aminotransferase abnormal",          NA,                               NA,       NA,       NA,        NA,                            NA,       NA,        NA, "Immune-Mediated Hepatitis (Diagnosis and Lab Abnormalities)", 20000008, "NARROW",        1,           NA,      NA,
    "02", "2020-06-05 23:59:59",                 "BASEDOW'S DISEASE",      5,                 "Basedow's disease",          NA, "Immune-Mediated Hypothyroidism", 20000161, "NARROW",         1,                            NA,       NA,        NA,                                                            NA,       NA,       NA,       NA,           NA,      NA,
    "03", "2020-06-07 23:59:59",                         "SOME TERM",      2,                        "Some query", "Some term",                               NA,       NA,       NA,        NA,                            NA,       NA,        NA,                                                            NA,       NA,       NA,       NA, "Some query",   11111,
    "05", "2020-06-09 23:59:59",              "ALVEOLAR PROTEINOSIS",      7,              "Alveolar proteinosis",          NA,                               NA,       NA,       NA,        NA, "Immune-Mediated Pneumonitis", "NARROW",         1,                                                            NA,       NA,       NA,       NA,           NA,      NA
  )
  # nolint end

  actual_output <- derive_vars_query(adae, queries)

  expect_dfs_equal(expected_output, actual_output, keys = "USUBJID")
})

## Test 2: Derive when no unique key excluding `SRCVAR` columns ----
test_that("derive_vars_query Test 2: Derive when no unique key excluding `SRCVAR` columns", {
  query <- tibble::tribble(
    ~PREFIX, ~GRPNAME, ~SRCVAR, ~TERMCHAR, ~GRPID, ~TERMNUM,
    "CQ42", "My Query", "AEDECOD", "PTSI", 1, NA_real_,
    "CQ42", "My Query", "AELLT", "LLTSI", 1, NA_real_
  )

  my_ae <- tibble::tribble(
    ~USUBJID, ~ASTDY,    ~AEDECOD,  ~AELLT,
    "1",           1,      "PTSI", "other",
    "1",           2, "something", "LLTSI",
    "1",           2,      "PTSI", "LLTSI",
    "1",           2, "something", "other"
  )

  expected_output <- tibble::tribble(
    ~USUBJID, ~ASTDY,    ~AEDECOD, ~AELLT,       ~CQ42NAM,     ~CQ42CD,
    "1",           1,      "PTSI", "other",    "My Query",           1,
    "1",           2, "something", "LLTSI",    "My Query",           1,
    "1",           2,      "PTSI", "LLTSI",    "My Query",           1,
    "1",           2, "something", "other", NA_character_, NA_integer_
  )

  actual_output <- derive_vars_query(my_ae, dataset_queries = query)

  expect_equal(expected_output, actual_output)
})

## Test 3: Derive when an adverse event is in multiple baskets ----
test_that("derive_vars_query Test 3: Derive when an adverse event is in multiple baskets", {
  query <- tibble::tribble(
    ~PREFIX, ~GRPNAME, ~SRCVAR, ~TERMCHAR, ~GRPID, ~TERMNUM,
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
    ~USUBJID, ~ASTDY, ~AEDECOD, ~AELLT, ~CQ40NAM, ~CQ40CD, ~CQ42NAM, ~CQ42CD,
    "1", 1, "PTSI", "other", "My Query 1", 1, NA_character_, NA_integer_,
    "1", 2, "something", "LLTSI", NA_character_, NA_integer_, "My Query 2", 2,
    "1", 2, "PTSI", "LLTSI", "My Query 1", 1, "My Query 2", 2,
    "1", 2, "something", "other", NA_character_, NA_integer_, NA_character_, NA_integer_
  )

  actual_output <- derive_vars_query(my_ae, dataset_queries = query)

  expect_equal(expected_output, actual_output)
})


## Test 4: Derive when no GRPID or SCOPE column ----
test_that("derive_vars_query Test 4: Derive when no GRPID or SCOPE column", {
  query <- tibble::tribble(
    ~PREFIX, ~GRPNAME, ~SRCVAR, ~TERMCHAR, ~TERMNUM,
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

  actual_output <- derive_vars_query(my_ae, dataset_queries = query)

  expected_output <- tibble::tribble(
    ~USUBJID, ~ASTDY, ~AEDECOD, ~AELLT, ~CQ42NAM,
    "1", 1, "PTSI", "other", "My Query",
    "1", 2, "something", "LLTSI", "My Query",
    "1", 2, "PTSI", "LLTSI", "My Query",
    "1", 2, "something", "other", NA_character_
  )

  expect_equal(expected_output, actual_output)
})

## Test 5: Derive decides between TERMCHAR and TERMNUM based on type ----
test_that("derive_vars_query Test 5: Derive decides between TERMCHAR and TERMNUM based on type", {
  query <- tibble::tribble(
    ~PREFIX, ~GRPNAME, ~SRCVAR, ~TERMCHAR, ~GRPID, ~TERMNUM,
    "CQ40", "My Query 1", "AEDECOD", "PTSI", 1, NA,
    "CQ42", "My Query 2", "AELLTCD", NA_character_, 2, 1
  )

  my_ae <- tibble::tribble(
    ~USUBJID, ~ASTDY, ~AEDECOD, ~AELLT, ~AELLTCD,
    "1", 1, "PTSI", "other", NA,
    "1", 2, "PTSI", "LLTSI", NA,
    "1", 3, NA, NA, 1
  )

  actual_output <- derive_vars_query(my_ae, dataset_queries = query)

  expected_output <- tibble::tribble(
    ~USUBJID, ~ASTDY, ~AEDECOD, ~AELLT, ~AELLTCD, ~CQ40NAM, ~CQ40CD, ~CQ42NAM, ~CQ42CD,
    "1", 1, "PTSI", "other", NA, "My Query 1", 1, NA, NA,
    "1", 2, "PTSI", "LLTSI", NA, "My Query 1", 1, NA, NA,
    "1", 3, NA, NA, 1, NA, NA, "My Query 2", 2
  )

  expect_equal(expected_output, actual_output)

  expect_snapshot(
    derive_vars_query(mutate(my_ae, AELLTCD = as.logical(AELLTCD)), query),
    error = TRUE
  )
})

## Test 6: Error is given if both TERMCHAR/TERMNUM are NA/empty ----
test_that("derive_vars_query Test 6: Error is given if both TERMCHAR/TERMNUM are NA/empty", {
  query <- tibble::tribble(
    ~PREFIX, ~GRPNAME, ~SRCVAR, ~TERMCHAR, ~GRPID, ~TERMNUM,
    "CQ40", "My Query 1", "AEDECOD", NA_character_, 1, NA,
    "CQ42", "My Query 2", "AELLTCD", NA_character_, 2, 1
  )

  my_ae <- tibble::tribble(
    ~USUBJID, ~ASTDY, ~AEDECOD, ~AELLT, ~AELLTCD,
    "1", 1, "PTSI", "other", NA,
    "1", 2, "PTSI", "LLTSI", NA,
    "1", 3, NA, NA, 1
  )

  snap <- expect_snapshot(
    derive_vars_query(my_ae, query),
    error = TRUE
  )
})

## Test 7: character SRCVAR and just TERMCHAR is provided ----
test_that("derive_vars_query Test 7: character SRCVAR and just TERMCHAR is provided", {
  query <- tibble::tribble(
    ~PREFIX, ~GRPNAME, ~SRCVAR, ~TERMCHAR, ~GRPID,
    "CQ40", "My Query 1", "AEDECOD", "PTSI", 1
  )

  my_ae <- tibble::tribble(
    ~USUBJID, ~ASTDY, ~AEDECOD, ~AELLT, ~AELLTCD,
    "1", 1, "PTSI", "other", NA,
    "1", 2, "PTSI", "LLTSI", NA,
    "1", 3, NA, NA, 1
  )
  actual_output <- derive_vars_query(my_ae, dataset_queries = query)

  expected_output <- tibble::tribble(
    ~USUBJID, ~ASTDY, ~AEDECOD, ~AELLT, ~AELLTCD, ~CQ40NAM, ~CQ40CD,
    "1", 1, "PTSI", "other", NA, "My Query 1", 1,
    "1", 2, "PTSI", "LLTSI", NA, "My Query 1", 1,
    "1", 3, NA, NA, 1, NA, NA
  )

  expect_equal(expected_output, actual_output)
})

## Test 8: numeric SRCVAR and just TERMNUM is provided ----
test_that("derive_vars_query Test 8: numeric SRCVAR and just TERMNUM is provided", {
  query <- tibble::tribble(
    ~PREFIX, ~GRPNAME, ~SRCVAR, ~GRPID, ~TERMNUM,
    "CQ42", "My Query 2", "AELLTCD", 2, 1
  )

  my_ae <- tibble::tribble(
    ~USUBJID, ~ASTDY, ~AEDECOD, ~AELLT, ~AELLTCD,
    "1", 1, "PTSI", "other", NA,
    "1", 2, "PTSI", "LLTSI", NA,
    "1", 3, NA, NA, 1
  )

  actual_output <- derive_vars_query(my_ae, dataset_queries = query)

  expected_output <- tibble::tribble(
    ~USUBJID, ~ASTDY, ~AEDECOD, ~AELLT, ~AELLTCD, ~CQ42NAM, ~CQ42CD,
    "1", 1, "PTSI", "other", NA, NA, NA,
    "1", 2, "PTSI", "LLTSI", NA, NA, NA,
    "1", 3, NA, NA, 1, "My Query 2", 2
  )

  expect_equal(expected_output, actual_output)
})

## Test 9: Error if requested variables already exist ----
test_that("derive_vars_query Test 9: Error if requested variables already exist", {
  # nolint start
  queries <- tibble::tribble(
    ~PREFIX, ~GRPNAME,                         ~GRPID,   ~SCOPE,   ~SCOPEN, ~SRCVAR,   ~TERMCHAR,
    "SMQ03", "Immune-Mediated Hypothyroidism", 20000161, "NARROW", 1,       "AEDECOD", "BASEDOW'S DISEASE",
  )

  adae <- tibble::tribble(
    ~USUBJID, ~AEDECOD,                            ~SMQ03NAM,                        ~SMQ03CD, ~SMQ03SC, ~SMQ03SCN,
    "01",     "ALANINE aminotransferase abnormal", NA_character_,                    NA,       NA,       NA,
    "02",     "Basedow's disease",                 "Immune-Mediated Hypothyroidism", 20000161, "NARROW", 1,
    "03",     "Some term",                         NA_character_,                    NA,       NA,       NA,
    "05",     "Alveolar proteinosis",              NA_character_,                    NA,       NA,       NA
  )
  # nolint end

  expect_snapshot(
    derive_vars_query(adae, queries),
    error = TRUE
  )
})


# assert_valid_queries ----
## Test 10: assert_valid_queries checks ----
test_that("assert_valid_queries Test 10: assert_valid_queries checks", {
  query <- tibble::tribble(
    ~PREFIX, ~GRPNAME,     ~SRCVAR,   ~TERMCHAR,     ~GRPID, ~TERMNUM,
    "CQ40",  "My Query 1", "AEDECOD", "PTSI",             1,       NA,
    "CQ42",  "My Query 2", "AELLTCD", NA_character_,      2,        1
  )

  expect_snapshot(
    assert_valid_queries(
      mutate(query, PREFIX = c("30", "55")),
      "test"
    ),
    error = TRUE
  )

  expect_snapshot(
    assert_valid_queries(
      mutate(query, PREFIX = c("AA", "BB")),
      "test"
    ),
    error = TRUE
  )

  expect_snapshot(
    assert_valid_queries(
      mutate(query, GRPNAME = c("", "A")),
      "test"
    ),
    error = TRUE
  )

  expect_snapshot(
    assert_valid_queries(
      mutate(query, GRPID = as.character(GRPID)),
      "test"
    ),
    error = TRUE
  )

  expect_snapshot(
    assert_valid_queries(
      mutate(query, SCOPE = letters[1:2]),
      "test"
    ),
    error = TRUE
  )

  expect_snapshot(
    assert_valid_queries(
      mutate(query, SCOPEN = 10:11),
      "test"
    ),
    error = TRUE
  )

  expect_snapshot(
    assert_valid_queries(
      mutate(query, PREFIX = c("CQ40", "CQ40")),
      "test"
    ),
    error = TRUE
  )

  expect_snapshot(
    assert_valid_queries(
      mutate(
        query,
        PREFIX = c("CQ40", "CQ40"),
        GRPNAME = c("My Query 1", "My Query 1")
      ),
      "test"
    ),
    error = TRUE
  )

  expect_error(
    assert_valid_queries(
      mutate(
        query,
        SCOPE = c("BROAD", "NARROW"),
        SCOPEN = c(1, 1)
      ),
      "test"
    ),
    class = "assert_one_to_one"
  )

  expect_snapshot(
    assert_valid_queries(
      queries = bind_rows(query, query),
      queries_name = "test"
    ),
    error = TRUE
  )
})
