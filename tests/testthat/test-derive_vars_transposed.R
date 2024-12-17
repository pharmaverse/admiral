dataset <- tibble::tribble(
  ~STUDYID,  ~USUBJID, ~VAR1,
  "STUDY01", "P01",        3,
  "STUDY01", "P02",       31,
  "STUDY01", "P03",       42
)
dataset_merge <- tibble::tribble(
  ~STUDYID,  ~USUBJID, ~TESTCD, ~VALUE,
  "STUDY01", "P01",    "T01",       31,
  "STUDY01", "P01",    "T02",        5,
  "STUDY01", "P02",    "T01",        3,
  "STUDY01", "P03",    "T02",        9
)

## Test 1: the merge dataset is transposed and merged correctly ----
test_that("derive_vars_transposed Test 1: the merge dataset is transposed and merged correctly", {
  expected_output <- tibble::tribble(
    ~STUDYID,  ~USUBJID, ~VAR1, ~T01, ~T02,
    "STUDY01", "P01",        3,   31,    5,
    "STUDY01", "P02",       31,    3,   NA,
    "STUDY01", "P03",       42,   NA,    9
  )
  actual_output <- derive_vars_transposed(
    dataset,
    dataset_merge,
    by_vars = get_admiral_option("subject_keys"),
    key_var = TESTCD,
    value_var = VALUE
  )

  expect_dfs_equal(expected_output, actual_output, keys = "USUBJID")
})


## Test 2: filtering the merge dataset works ----
test_that("derive_vars_transposed Test 2: filtering the merge dataset works", {
  expected_output <- tibble::tribble(
    ~STUDYID,  ~USUBJID, ~VAR1, ~T01,
    "STUDY01", "P01",        3,   31,
    "STUDY01", "P02",       31,    3,
    "STUDY01", "P03",       42,   NA
  )

  dataset_merge <- tibble::tribble(
    ~STUDYID,  ~USUBJID, ~TESTCD, ~VALUE,
    "STUDY01", "P01",    "T01",       31,
    "STUDY01", "P01",    "T02",        5,
    "STUDY01", "P02",    "T01",        3,
    "STUDY01", "P01",    "T02",        9
  )

  actual_output <- derive_vars_transposed(
    dataset,
    dataset_merge,
    by_vars = get_admiral_option("subject_keys"),
    key_var = TESTCD,
    value_var = VALUE,
    filter = TESTCD == "T01"
  )

  expect_dfs_equal(expected_output, actual_output, keys = "USUBJID")
})

## Test 3: filter merge dataset 'many-to-one' ----
test_that("derive_vars_transposed Test 3: filter merge dataset 'many-to-one'", {
  expect_snapshot(
    derive_vars_transposed(
      dataset,
      dataset_merge,
      by_vars = get_admiral_option("subject_keys"),
      key_var = TESTCD,
      value_var = VALUE,
      filter = TESTCD == "T01",
      relationship = "many-to-one"
    )
  )
})

## Test 4: error if `relationship` is unexpected ----
test_that("derive_vars_transposed Test 4: error if `relationship` is unexpected", {
  cm <- tibble::tribble(
    ~STUDYID,  ~USUBJID,       ~CMGRPID, ~CMREFID,  ~CMDECOD,
    "STUDY01", "BP40257-1001", "14",     "1192056", "PARACETAMOL",
    "STUDY01", "BP40257-1001", "18",     "2007001", "SOLUMEDROL",
    "STUDY01", "BP40257-1002", "19",     "2791596", "SPIRONOLACTONE"
  )
  facm <- tibble::tribble(
    ~STUDYID,  ~USUBJID,       ~FAGRPID, ~FAREFID,  ~FATESTCD,  ~FASTRESC,
    "STUDY01", "BP40257-1001", "1",      "1192056", "CMATC1CD", "N",
    "STUDY01", "BP40257-1001", "1",      "1192056", "CMATC2CD", "N02",
    "STUDY01", "BP40257-1001", "1",      "1192056", "CMATC3CD", "N02B",
    "STUDY01", "BP40257-1001", "1",      "1192056", "CMATC4CD", "N02BE",
    "STUDY01", "BP40257-1001", "1",      "2007001", "CMATC1CD", "D",
    "STUDY01", "BP40257-1001", "1",      "2007001", "CMATC2CD", "D10",
    "STUDY01", "BP40257-1001", "1",      "2007001", "CMATC3CD", "D10A",
    "STUDY01", "BP40257-1001", "1",      "2007001", "CMATC4CD", "D10AA",
    "STUDY01", "BP40257-1001", "2",      "2007001", "CMATC1CD", "D",
    "STUDY01", "BP40257-1001", "2",      "2007001", "CMATC2CD", "D07",
    "STUDY01", "BP40257-1001", "2",      "2007001", "CMATC3CD", "D07A",
    "STUDY01", "BP40257-1001", "2",      "2007001", "CMATC4CD", "D07AA",
    "STUDY01", "BP40257-1001", "3",      "2007001", "CMATC1CD", "H",
    "STUDY01", "BP40257-1001", "3",      "2007001", "CMATC2CD", "H02",
    "STUDY01", "BP40257-1001", "3",      "2007001", "CMATC3CD", "H02A",
    "STUDY01", "BP40257-1001", "3",      "2007001", "CMATC4CD", "H02AB",
    "STUDY01", "BP40257-1002", "1",      "2791596", "CMATC1CD", "C",
    "STUDY01", "BP40257-1002", "1",      "2791596", "CMATC2CD", "C03",
    "STUDY01", "BP40257-1002", "1",      "2791596", "CMATC3CD", "C03D",
    "STUDY01", "BP40257-1002", "1",      "2791596", "CMATC4CD", "C03DA"
  )

  expect_snapshot(
    cm %>%
      derive_vars_transposed(
        facm,
        by_vars = exprs(USUBJID, CMREFID = FAREFID),
        id_vars = exprs(FAGRPID),
        key_var = FATESTCD,
        value_var = FASTRESC,
        relationship = "one-to-one"
      ),
    error = TRUE
  )

  expect_snapshot(
    cm %>%
      derive_vars_transposed(
        facm,
        by_vars = exprs(USUBJID, CMREFID = FAREFID),
        id_vars = exprs(FAGRPID),
        key_var = FATESTCD,
        value_var = FASTRESC,
        relationship = "many-to-one"
      ),
    error = TRUE
  )
})
