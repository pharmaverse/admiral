dataset <- tibble::tribble(
  ~USUBJID, ~VAR1,
  "P01",    3,
  "P02",    31,
  "P03",    42
)
dataset_merge <- tibble::tribble(
  ~USUBJID, ~TESTCD, ~VALUE,
  "P01",    "T01",   31,
  "P01",    "T02",   5,
  "P02",    "T01",   3,
  "P03",    "T02",   9
)

test_that("the merge dataset is transposed and merged correctly", {
  expected_output <- tibble::tribble(
    ~USUBJID, ~VAR1, ~T01, ~T02,
    "P01",     3,    31,   5,
    "P02",    31,     3,   NA,
    "P03",    42,    NA,   9
  )
  actual_output <- derive_vars_transposed(
    dataset,
    dataset_merge,
    by_vars = vars(USUBJID),
    key_var = TESTCD,
    value_var = VALUE
  )

  expect_dfs_equal(expected_output, actual_output, keys = "USUBJID")
})


test_that("filtering the merge dataset works", {
  expected_output <- tibble::tribble(
    ~USUBJID, ~VAR1, ~T01,
    "P01",     3,    31,
    "P02",    31,     3,
    "P03",    42,    NA
  )
  actual_output <- derive_vars_transposed(
    dataset,
    dataset_merge,
    by_vars = vars(USUBJID),
    key_var = TESTCD,
    value_var = VALUE,
    filter = TESTCD == "T01"
  )

  expect_dfs_equal(expected_output, actual_output, keys = "USUBJID")
})

test_that("ATC variables are merged properly", {
  cm <- tibble::tribble(
          ~USUBJID, ~CMGRPID,  ~CMREFID,            ~CMDECOD,
    "BP40257-1001",     "14", "1192056",       "PARACETAMOL",
    "BP40257-1001",     "18", "2007001",       "SOLUMEDROL",
    "BP40257-1002",     "19", "2791596",    "SPIRONOLACTONE"
  )
  facm <- tibble::tribble(
          ~USUBJID, ~FAGRPID,  ~FAREFID,   ~FATESTCD, ~FASTRESC,
    "BP40257-1001",      "1", "1192056",  "CMATC1CD",       "N",
    "BP40257-1001",      "1", "1192056",  "CMATC2CD",     "N02",
    "BP40257-1001",      "1", "1192056",  "CMATC3CD",    "N02B",
    "BP40257-1001",      "1", "1192056",  "CMATC4CD",   "N02BE",

    "BP40257-1001",      "1", "2007001",  "CMATC1CD",       "D",
    "BP40257-1001",      "1", "2007001",  "CMATC2CD",     "D10",
    "BP40257-1001",      "1", "2007001",  "CMATC3CD",    "D10A",
    "BP40257-1001",      "1", "2007001",  "CMATC4CD",   "D10AA",
    "BP40257-1001",      "2", "2007001",  "CMATC1CD",       "D",
    "BP40257-1001",      "2", "2007001",  "CMATC2CD",     "D07",
    "BP40257-1001",      "2", "2007001",  "CMATC3CD",    "D07A",
    "BP40257-1001",      "2", "2007001",  "CMATC4CD",   "D07AA",
    "BP40257-1001",      "3", "2007001",  "CMATC1CD",       "H",
    "BP40257-1001",      "3", "2007001",  "CMATC2CD",     "H02",
    "BP40257-1001",      "3", "2007001",  "CMATC3CD",    "H02A",
    "BP40257-1001",      "3", "2007001",  "CMATC4CD",   "H02AB",

    "BP40257-1002",      "1", "2791596",  "CMATC1CD",       "C",
    "BP40257-1002",      "1", "2791596",  "CMATC2CD",     "C03",
    "BP40257-1002",      "1", "2791596",  "CMATC3CD",    "C03D",
    "BP40257-1002",      "1", "2791596",  "CMATC4CD",   "C03DA"
  )
  # nolint start
  expected_output <- tibble::tribble(
          ~USUBJID, ~CMGRPID,  ~CMREFID,            ~CMDECOD, ~ATC1CD, ~ATC2CD, ~ATC3CD, ~ATC4CD,
    "BP40257-1001",     "14", "1192056",       "PARACETAMOL",     "N",   "N02",  "N02B", "N02BE",
    "BP40257-1001",     "18", "2007001",        "SOLUMEDROL",     "D",   "D07",  "D07A", "D07AA",
    "BP40257-1001",     "18", "2007001",        "SOLUMEDROL",     "D",   "D10",  "D10A", "D10AA",
    "BP40257-1001",     "18", "2007001",        "SOLUMEDROL",     "H",   "H02",  "H02A", "H02AB",
    "BP40257-1002",     "19", "2791596",    "SPIRONOLACTONE",     "C",   "C03",  "C03D", "C03DA"
  )
  # nolint end
  actual_output <- derive_vars_atc(cm, facm)

  expect_dfs_equal(expected_output, actual_output, keys = c("USUBJID", "CMDECOD", "ATC4CD"))
})
