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
          ~USUBJID, ~CMSEQ, ~CMGRPID,  ~CMREFID,            ~CMDECOD,
    "BP40257-1001",     1L,     "14", "1192056",       "PARACETAMOL",
    "BP40257-1001",     2L,      "9", "1192057",       "PARACETAMOL",
    "BP40257-1002",     1L,     "19", "2791596",    "SPIRONOLACTONE",
    "BP40257-1002",     2L,     "12", "1265064", "ENOXAPARIN SODIUM"
  )
  facm <- tibble::tribble(
          ~USUBJID, ~FAGRPID,  ~FAREFID, ~FATESTCD,                            ~FAORRES,
    "BP40257-1001",      "1", "1192056",  "CMATC1",                    "NERVOUS SYSTEM",
    "BP40257-1001",      "1", "1192056",  "CMATC2",                        "ANALGESICS",
    "BP40257-1001",      "1", "1192056",  "CMATC3", "OTHER ANALGESICS AND ANTIPYRETICS",
    "BP40257-1001",      "1", "1192056",  "CMATC4",                          "ANILIDES",
    "BP40257-1001",      "1", "1192057",  "CMATC1",                    "NERVOUS SYSTEM",
    "BP40257-1001",      "1", "1192057",  "CMATC2",                        "ANALGESICS",
    "BP40257-1001",      "1", "1192057",  "CMATC3", "OTHER ANALGESICS AND ANTIPYRETICS",
    "BP40257-1001",      "1", "1192057",  "CMATC4",                          "ANILIDES",
    "BP40257-1002",      "1", "1265064",  "CMATC1",    "BLOOD AND BLOOD FORMING ORGANS",
    "BP40257-1002",      "1", "1265064",  "CMATC2",             "ANTITHROMBOTIC AGENTS",
    "BP40257-1002",      "1", "1265064",  "CMATC3",             "ANTITHROMBOTIC AGENTS",
    "BP40257-1002",      "1", "1265064",  "CMATC4",                     "HEPARIN GROUP",
    "BP40257-1002",      "1", "2791596",  "CMATC1",             "CARDIOVASCULAR SYSTEM",
    "BP40257-1002",      "1", "2791596",  "CMATC2",                         "DIURETICS",
    "BP40257-1002",      "1", "2791596",  "CMATC3",          "POTASSIUM-SPARING AGENTS",
    "BP40257-1002",      "1", "2791596",  "CMATC4",           "ALDOSTERONE ANTAGONISTS"
  )
  # nolint start
  expected_output <- tibble::tribble(
          ~USUBJID, ~CMSEQ, ~CMGRPID,  ~CMREFID,            ~CMDECOD,                                ~ATC1,                      ~ATC2,                               ~ATC3,                     ~ATC4,
    "BP40257-1001",     1L,     "14", "1192056",       "PARACETAMOL",                     "NERVOUS SYSTEM",               "ANALGESICS", "OTHER ANALGESICS AND ANTIPYRETICS", "ANILIDES",
    "BP40257-1001",     2L,      "9", "1192057",       "PARACETAMOL",                     "NERVOUS SYSTEM",               "ANALGESICS", "OTHER ANALGESICS AND ANTIPYRETICS", "ANILIDES",
    "BP40257-1002",     1L,     "19", "2791596",    "SPIRONOLACTONE",              "CARDIOVASCULAR SYSTEM",                "DIURETICS",          "POTASSIUM-SPARING AGENTS", "ALDOSTERONE ANTAGONISTS",
    "BP40257-1002",     2L,     "12", "1265064", "ENOXAPARIN SODIUM",     "BLOOD AND BLOOD FORMING ORGANS",    "ANTITHROMBOTIC AGENTS",             "ANTITHROMBOTIC AGENTS", "HEPARIN GROUP"
  )
  # nolint end
  actual_output <- derive_vars_atc(cm, facm)

  expect_dfs_equal(expected_output, actual_output, keys = c("USUBJID", "CMSEQ"))
})

