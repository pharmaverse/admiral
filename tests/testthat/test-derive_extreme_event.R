# derive_extreme_records ----
## Test 1: `mode` = first ----
test_that("derive_extreme_records Test 1: `mode` = first", {
  input <- tibble::tribble(
    ~USUBJID, ~PARAMCD,       ~AVALC,        ~ADY,
    "1",      "NO SLEEP",     "N",              1,
    "1",      "WAKE UP",      "N",              2,
    "1",      "FALL ASLEEP",  "N",              3,
    "2",      "NO SLEEP",     "N",              1,
    "2",      "WAKE UP",      "Y",              2,
    "2",      "WAKE UP",      "Y",              3,
    "2",      "FALL ASLEEP",  "N",              4,
    "3",      "NO SLEEP",     NA_character_,    1
  )

  expected_output <- bind_rows(
    input,
    tibble::tribble(
      ~USUBJID, ~PARAMCD, ~AVALC,                            ~AVAL, ~ADY,
      "1",      "WSP",    "No sleeping problems",                4,    1,
      "2",      "WSP",    "Waking up more than three times",     2,    2,
      "3",      "WSP",    "Missing",                            99,    1
    )
  )

  actual_output <- derive_extreme_event(
    input,
    by_vars = exprs(USUBJID),
    events = list(
      event(
        condition = PARAMCD == "NO SLEEP" & AVALC == "Y",
        set_values_to = exprs(AVALC = "No sleep", AVAL = 1)
      ),
      event(
        condition = PARAMCD == "WAKE UP" & AVALC == "Y",
        set_values_to = exprs(AVALC = "Waking up more than three times", AVAL = 2)
      ),
      event(
        condition = PARAMCD == "FALL ASLEEP" & AVALC == "Y",
        set_values_to = exprs(AVALC = "More than 30 mins to fall asleep", AVAL = 3)
      ),
      event(
        condition = all(AVALC == "N"),
        set_values_to = exprs(
          AVALC = "No sleeping problems", AVAL = 4
        )
      ),
      event(
        condition = TRUE,
        set_values_to = exprs(AVALC = "Missing", AVAL = 99)
      )
    ),
    order = exprs(ADY),
    mode = "first",
    set_values_to = exprs(
      PARAMCD = "WSP"
    ),
    check_type = "none"
  )

  expect_dfs_equal(
    base = expected_output,
    compare = actual_output,
    keys = c("USUBJID", "PARAMCD", "ADY")
  )
})

## Test 2: `mode` = last ----
test_that("derive_extreme_records Test 2: `mode` = last", {
  input <- tibble::tribble(
    ~USUBJID, ~PARAMCD,       ~AVALC,        ~ADY,
    "1",      "NO SLEEP",     "N",              1,
    "1",      "WAKE UP",      "N",              2,
    "1",      "FALL ASLEEP",  "N",              3,
    "2",      "NO SLEEP",     "N",              1,
    "2",      "WAKE UP",      "Y",              2,
    "2",      "WAKE UP",      "Y",              3,
    "2",      "FALL ASLEEP",  "N",              4,
    "3",      "NO SLEEP",     NA_character_,    1
  )

  expected_output <- bind_rows(
    input,
    tibble::tribble(
      ~USUBJID, ~PARAMCD, ~AVALC,                            ~AVAL, ~ADY,
      "1",      "WSP",    "No sleeping problems",                4,    3,
      "2",      "WSP",    "Waking up more than three times",     2,    3,
      "3",      "WSP",    "Missing",                            99,    1
    )
  )

  actual_output <- derive_extreme_event(
    input,
    by_vars = exprs(USUBJID),
    events = list(
      event(
        condition = PARAMCD == "NO SLEEP" & AVALC == "Y",
        set_values_to = exprs(AVALC = "No sleep", AVAL = 1)
      ),
      event(
        condition = PARAMCD == "WAKE UP" & AVALC == "Y",
        set_values_to = exprs(AVALC = "Waking up more than three times", AVAL = 2)
      ),
      event(
        condition = PARAMCD == "FALL ASLEEP" & AVALC == "Y",
        set_values_to = exprs(AVALC = "More than 30 mins to fall asleep", AVAL = 3)
      ),
      event(
        condition = all(AVALC == "N"),
        set_values_to = exprs(
          AVALC = "No sleeping problems", AVAL = 4
        )
      ),
      event(
        condition = TRUE,
        set_values_to = exprs(AVALC = "Missing", AVAL = 99)
      )
    ),
    order = exprs(ADY),
    mode = "last",
    set_values_to = exprs(
      PARAMCD = "WSP"
    ),
    check_type = "none"
  )

  expect_dfs_equal(
    base = expected_output,
    compare = actual_output,
    keys = c("USUBJID", "PARAMCD", "ADY")
  )
})

## Test 3: `source_datasets` works ----
test_that("derive_extreme_records Test 3: `source_datasets` works", {
  adsl <- tibble::tribble(
    ~USUBJID, ~TRTSDTC,      ~CHECKKEPTCOL, ~CHECKNOTKEPTCOL,
    "1",      "2020-01-01",  "001",         "991",
    "2",      "2019-12-12",  "002",         "992",
    "3",      "2019-11-11",  "003",         "993",
    "4",      "2019-12-30",  "004",         "994",
    "5",      "2020-01-01",  "005",         "995",
    "6",      "2020-02-02",  "006",         "996",
    "7",      "2020-02-02",  "007",         "997",
    "8",      "2020-04-01",  "008",         "999"
  ) %>%
    mutate(
      TRTSDT = lubridate::ymd(TRTSDTC),
      STUDYID = "XX1234"
    )

  adrs <- tibble::tribble(
    ~USUBJID, ~ADTC, ~AVALC, ~CHECKKEPTCOL,
    "1", "2020-01-01", "PR", "001",
    "1", "2020-02-01", "CR", "001",
    "1", "2020-02-16", "NE", "001",
    "1", "2020-03-01", "CR", "001",
    "1", "2020-04-01", "SD", "001",
    "2", "2020-01-01", "SD", "002",
    "2", "2020-02-01", "PR", "002",
    "2", "2020-03-01", "SD", "002",
    "2", "2020-03-13", "CR", "002",
    "3", "2019-11-12", "CR", "003",
    "3", "2019-12-02", "CR", "003",
    "3", "2020-01-01", "SD", "003",
    "4", "2020-01-01", "PR", "004",
    "4", "2020-03-01", "SD", "004",
    "4", "2020-04-01", "SD", "004",
    "4", "2020-05-01", "PR", "004",
    "4", "2020-05-15", "NON-CR/NON-PD", "004",
    "5", "2020-01-01", "PR", "005",
    "5", "2020-01-10", "SD", "005",
    "5", "2020-01-20", "PR", "005",
    "5", "2020-05-15", "NON-CR/NON-PD", "005",
    "6", "2020-02-06", "PR", "006",
    "6", "2020-02-16", "CR", "006",
    "6", "2020-03-30", "PR", "006",
    "7", "2020-02-06", "PR", "007",
    "7", "2020-02-16", "CR", "007",
    "7", "2020-04-01", "NE", "007"
  ) %>%
    mutate(
      PARAMCD = "OVR",
      ADT = lubridate::ymd(ADTC),
      STUDYID = "XX1234"
    ) %>%
    select(-ADTC) %>%
    derive_vars_merged(
      dataset_add = adsl,
      by_vars = exprs(STUDYID, USUBJID),
      new_vars = exprs(TRTSDT)
    )
  expected <- bind_rows(
    adrs,
    tibble::tribble(
      ~USUBJID, ~ADTC, ~AVALC, ~AVAL, ~TRTSDTC, ~CHECKKEPTCOL,
      "1", "2020-02-01", "CR", 11, "2020-01-01", "001",
      "2", "2020-03-13", "CR", 11, "2019-12-12", "002",
      "3", "2019-11-12", "CR", 11, "2019-11-11", "003",
      "4", "2020-01-01", "PR", 22, "2019-12-30", "004",
      "5", "2020-01-01", "PR", 22, "2020-01-01", "005",
      "6", "2020-02-16", "CR", 11, "2020-02-02", "006",
      "7", "2020-02-16", "CR", 11, "2020-02-02", "007",
      "8", "", "MISSING", 77, "2020-04-01", "008"
    ) %>%
      mutate(
        ADT = lubridate::ymd(ADTC),
        TRTSDT = lubridate::ymd(TRTSDTC),
        STUDYID = "XX1234",
        PARAMCD = "BOR",
        PARAM = "Best Overall Response"
      ) %>%
      select(-ADTC, -TRTSDTC)
  )

  actual <- derive_extreme_event(
    dataset = adrs,
    by_vars = exprs(STUDYID, USUBJID),
    order = exprs(ADT),
    mode = "first",
    filter_source = PARAMCD == "OVR",
    source_datasets = list(adsl = adsl),
    events = list(
      event(
        condition = AVALC == "CR",
        set_values_to = exprs(
          AVALC = "CR"
        )    )
    ),
    reference_date = TRTSDT,
    ref_start_window = 28,
    set_values_to = exprs(
      AVAL = {{ aval_fun_pass }}(AVALC),
      PARAMCD = "BOR",
      PARAM = "Best Overall Response"
    )
  )

  expect_dfs_equal(
    base    = expected,
    compare = actual,
    keys    = c("USUBJID", "PARAMCD", "ADT")
  )
})
