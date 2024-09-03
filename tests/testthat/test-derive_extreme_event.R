# derive_extreme_event ----
## Test 1: `mode` = first ----
test_that("derive_extreme_event Test 1: `mode` = first", {
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
    tmp_event_nr_var = event_nr,
    order = exprs(event_nr, ADY),
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
test_that("derive_extreme_event Test 2: `mode` = last", {
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
    tmp_event_nr_var = event_nr,
    order = exprs(desc(event_nr), ADY),
    mode = "last",
    set_values_to = exprs(
      PARAMCD = "WSP"
    ),
    check_type = "warning"
  )

  expect_dfs_equal(
    base = expected_output,
    compare = actual_output,
    keys = c("USUBJID", "PARAMCD", "ADY")
  )
})

## Test 3: `source_datasets` works ----
test_that("derive_extreme_event Test 3: `source_datasets` works", {
  adsl <- tibble::tribble(
    ~USUBJID, ~TRTSDTC,
    "1",      "2020-01-01",
    "2",      "2019-12-12",
    "3",      "2019-11-11",
    "4",      "2019-12-30",
    "5",      "2020-01-01",
    "6",      "2020-02-02",
    "7",      "2020-02-02",
    "8",      "2020-04-01"
  ) %>%
    mutate(
      TRTSDT = lubridate::ymd(TRTSDTC),
      STUDYID = "XX1234"
    )

  adrs <- tibble::tribble(
    ~USUBJID, ~ADTC, ~AVALC,
    "1", "2020-01-01", "PR",
    "1", "2020-02-01", "CR",
    "1", "2020-02-16", "NE",
    "1", "2020-03-01", "CR",
    "1", "2020-04-01", "SD",
    "2", "2020-01-01", "SD",
    "2", "2020-02-01", "PR",
    "2", "2020-03-01", "SD",
    "2", "2020-03-13", "CR",
    "3", "2019-11-12", "CR",
    "3", "2019-12-02", "CR",
    "3", "2020-01-01", "SD",
    "4", "2020-01-01", "PR",
    "4", "2020-03-01", "SD",
    "4", "2020-04-01", "SD",
    "4", "2020-05-01", "PR",
    "4", "2020-05-15", "NON-CR/NON-PD",
    "5", "2020-01-01", "PR",
    "5", "2020-01-10", "SD",
    "5", "2020-01-20", "PR",
    "5", "2020-05-15", "NON-CR/NON-PD",
    "6", "2020-02-06", "PR",
    "6", "2020-02-16", "CR",
    "6", "2020-03-30", "PR",
    "7", "2020-02-06", "PR",
    "7", "2020-02-16", "CR",
    "7", "2020-04-01", "NE"
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
      ~USUBJID, ~ADTC,        ~AVALC,    ~TRTSDTC,
      "1",      "2020-02-01", "CR",      "2020-01-01",
      "2",      "2020-03-13", "CR",      "2019-12-12",
      "3",      "2019-11-12", "CR",      "2019-11-11",
      "4",      "2020-01-01", "PR",      "2019-12-30",
      "5",      "2020-01-01", "PR",      "2020-01-01",
      "6",      "2020-02-16", "CR",      "2020-02-02",
      "7",      "2020-02-16", "CR",      "2020-02-02",
      "8",      "",           "MISSING", "2020-04-01"
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
    tmp_event_nr_var = event_nr,
    order = exprs(event_nr, ADT),
    mode = "first",
    source_datasets = list(adsl = adsl),
    events = list(
      event(
        condition = AVALC == "CR",
        set_values_to = exprs(
          AVALC = "CR"
        )
      ),
      event(
        condition = AVALC == "PR",
        set_values_to = exprs(
          AVALC = "PR"
        )
      ),
      event(
        condition = AVALC == "SD" & ADT >= TRTSDT + 28,
        set_values_to = exprs(
          AVALC = "SD"
        )
      ),
      event(
        condition = AVALC == "NON-CR/NON-PD" & ADT >= TRTSDT + 28,
        set_values_to = exprs(
          AVALC = "NON-CR/NON-PD"
        )
      ),
      event(
        condition = AVALC == "PD",
        set_values_to = exprs(
          AVALC = "PD"
        )
      ),
      event(
        condition = AVALC %in% c("SD", "NON-CR/NON-PD"),
        set_values_to = exprs(
          AVALC = "NE"
        )
      ),
      event(
        dataset_name = "adsl",
        condition = TRUE,
        set_values_to = exprs(
          AVALC = "MISSING"
        ),
        keep_source_vars = exprs(TRTSDT)
      )
    ),
    set_values_to = exprs(
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

## Test 4: event-specific mode ----
test_that("derive_extreme_event Test 4: event-specific mode", {
  adhy <- tibble::tribble(
    ~USUBJID, ~AVISITN, ~CRIT1FL,
    "1",             1, "Y",
    "1",             2, "Y",
    "2",             1, "Y",
    "2",             2, NA_character_,
    "2",             3, "Y",
    "2",             4, NA_character_
  ) %>%
    mutate(
      PARAMCD = "ALKPH",
      PARAM = "Alkaline Phosphatase (U/L)"
    )

  actual <- derive_extreme_event(
    adhy,
    by_vars = exprs(USUBJID),
    events = list(
      event(
        condition = is.na(CRIT1FL),
        set_values_to = exprs(AVALC = "N")
      ),
      event(
        condition = CRIT1FL == "Y",
        mode = "last",
        set_values_to = exprs(AVALC = "Y")
      )
    ),
    tmp_event_nr_var = event_nr,
    order = exprs(event_nr, AVISITN),
    mode = "first",
    keep_source_vars = exprs(AVISITN),
    set_values_to = exprs(
      PARAMCD = "ALK2",
      PARAM = "ALKPH <= 2 times ULN"
    )
  )

  expected <- bind_rows(
    adhy,
    tribble(
      ~USUBJID, ~AVISITN, ~AVALC,
      "1",             2, "Y",
      "2",             2, "N"
    ) %>%
      mutate(
        PARAMCD = "ALK2",
        PARAM = "ALKPH <= 2 times ULN"
      )
  )

  expect_dfs_equal(
    base = expected,
    compare = actual,
    keys = c("USUBJID", "AVISITN", "PARAMCD")
  )
})


## Test 5: event_joined() is handled correctly ----
test_that("derive_extreme_event Test 5: event_joined() is handled correctly", {
  adsl <- tibble::tribble(
    ~USUBJID, ~TRTSDTC,
    "1",      "2020-01-01",
    "2",      "2019-12-12",
    "3",      "2019-11-11",
    "4",      "2019-12-30",
    "5",      "2020-01-01",
    "6",      "2020-02-02",
    "7",      "2020-02-02",
    "8",      "2020-04-01",
    "9",      "2020-02-01"
  ) %>%
    mutate(
      TRTSDT = lubridate::ymd(TRTSDTC),
      STUDYID = "XX1234"
    )

  adrs <- tibble::tribble(
    ~USUBJID, ~ADTC,        ~AVALC,
    "1",      "2020-01-01", "PR",
    "1",      "2020-02-01", "CR",
    "1",      "2020-02-16", "NE",
    "1",      "2020-03-01", "CR",
    "1",      "2020-04-01", "SD",
    "2",      "2020-01-01", "SD",
    "2",      "2020-02-01", "PR",
    "2",      "2020-03-01", "SD",
    "2",      "2020-03-13", "CR",
    "3",      "2019-11-12", "CR",
    "3",      "2019-12-02", "CR",
    "3",      "2020-01-01", "SD",
    "4",      "2020-01-01", "PR",
    "4",      "2020-03-01", "SD",
    "4",      "2020-04-01", "SD",
    "4",      "2020-05-01", "PR",
    "4",      "2020-05-15", "NON-CR/NON-PD",
    "5",      "2020-01-01", "PR",
    "5",      "2020-01-10", "SD",
    "5",      "2020-01-20", "PR",
    "5",      "2020-05-15", "NON-CR/NON-PD",
    "6",      "2020-02-06", "PR",
    "6",      "2020-02-16", "CR",
    "6",      "2020-03-30", "PR",
    "7",      "2020-02-06", "PR",
    "7",      "2020-02-16", "CR",
    "7",      "2020-04-01", "NE",
    "9",      "2020-02-16", "PD"
  ) %>%
    mutate(
      PARAMCD = "OVR",
      ADT = lubridate::ymd(ADTC),
      STUDYID = "XX1234"
    ) %>%
    derive_vars_merged(
      dataset_add = adsl,
      by_vars = exprs(STUDYID, USUBJID),
      new_vars = exprs(TRTSDT)
    )

  actual <-
    derive_extreme_event(
      adrs,
      by_vars = exprs(STUDYID, USUBJID),
      tmp_event_nr_var = event_nr,
      order = exprs(event_nr, ADT),
      mode = "first",
      source_datasets = list(adsl = adsl),
      events = list(
        event_joined(
          join_vars = exprs(AVALC, ADT),
          join_type = "after",
          first_cond_upper = AVALC.join == "CR" &
            ADT.join >= ADT + 28,
          condition = AVALC == "CR" &
            all(AVALC.join %in% c("CR", "NE")) &
            count_vals(var = AVALC.join, val = "NE") <= 1,
          set_values_to = exprs(
            AVALC = "CR"
          )
        ),
        event_joined(
          join_vars = exprs(AVALC, ADT),
          join_type = "after",
          first_cond_upper = AVALC.join %in% c("CR", "PR") &
            ADT.join >= ADT + 28,
          condition = AVALC == "PR" &
            all(AVALC.join %in% c("CR", "PR", "NE")) &
            count_vals(var = AVALC.join, val = "NE") <= 1 &
            (
              min_cond(
                var = ADT.join,
                cond = AVALC.join == "CR"
              ) > max_cond(var = ADT.join, cond = AVALC.join == "PR") |
                count_vals(var = AVALC.join, val = "CR") == 0 |
                count_vals(var = AVALC.join, val = "PR") == 0
            ),
          set_values_to = exprs(
            AVALC = "PR"
          )
        ),
        event(
          condition = AVALC %in% c("CR", "PR", "SD") & ADT >= TRTSDT + 28,
          set_values_to = exprs(
            AVALC = "SD"
          )
        ),
        event(
          condition = AVALC == "NON-CR/NON-PD" & ADT >= TRTSDT + 28,
          set_values_to = exprs(
            AVALC = "NON-CR/NON-PD"
          )
        ),
        event(
          condition = AVALC == "PD",
          set_values_to = exprs(
            AVALC = "PD"
          )
        ),
        event(
          condition = AVALC %in% c("CR", "PR", "SD", "NON-CR/NON-PD", "NE"),
          set_values_to = exprs(
            AVALC = "NE"
          )
        ),
        event(
          dataset_name = "adsl",
          condition = TRUE,
          set_values_to = exprs(
            AVALC = "MISSING"
          ),
          keep_source_vars = exprs(TRTSDT)
        )
      ),
      set_values_to = exprs(
        PARAMCD = "CBOR",
        PARAM = "Best Confirmed Overall Response by Investigator"
      )
    )

  expected <- bind_rows(
    adrs,
    tibble::tribble(
      ~USUBJID, ~ADTC,         ~AVALC,
      "1",      "2020-02-01",  "CR",
      "2",      "2020-02-01",  "SD",
      "3",      "2020-01-01",  "SD",
      "4",      "2020-03-01",  "SD",
      "5",      "2020-05-15",  "NON-CR/NON-PD",
      "6",      "2020-03-30",  "SD",
      "7",      "2020-02-06",  "NE",
      "8",      NA_character_, "MISSING",
      "9",      "2020-02-16",  "PD"
    ) %>%
      mutate(
        ADT = lubridate::ymd(ADTC),
        STUDYID = "XX1234",
        PARAMCD = "CBOR",
        PARAM = "Best Confirmed Overall Response by Investigator"
      ) %>%
      derive_vars_merged(
        dataset_add = adsl,
        by_vars = exprs(STUDYID, USUBJID),
        new_vars = exprs(TRTSDT)
      )
  )

  expect_dfs_equal(
    base = expected,
    compare = actual,
    keys = c("USUBJID", "PARAMCD", "ADT")
  )
})

## Test 6: no tmp_event_nr_var ----
test_that("derive_extreme_event Test 6: no tmp_event_nr_var", {
  adrs <- tibble::tribble(
    ~USUBJID, ~AVISITN, ~AVALC,
    "1",             1, "PR",
    "1",             2, "CR",
    "1",             3, "CR"
  ) %>%
    mutate(PARAMCD = "OVR")

  actual <- derive_extreme_event(
    adrs,
    by_vars = exprs(USUBJID),
    order = exprs(AVISITN),
    mode = "first",
    events = list(
      event_joined(
        join_vars = exprs(AVALC),
        join_type = "after",
        first_cond_upper = AVALC.join == "CR",
        condition = AVALC == "CR",
        set_values_to = exprs(AVALC = "Y")
      ),
      event_joined(
        join_vars = exprs(AVALC),
        join_type = "after",
        first_cond_upper = AVALC.join %in% c("CR", "PR"),
        condition = AVALC == "PR",
        set_values_to = exprs(AVALC = "Y")
      )
    ),
    set_values_to = exprs(
      PARAMCD = "CRSP"
    )
  )

  expected <- bind_rows(
    adrs,
    tibble::tribble(
      ~USUBJID, ~AVISITN, ~AVALC, ~PARAMCD,
      "1",             1, "Y",    "CRSP"
    )
  )

  expect_dfs_equal(
    base = expected,
    compare = actual,
    keys = c("USUBJID", "PARAMCD", "AVISITN")
  )
})

## Test 7: mode and condition used in event() ----
test_that("derive_extreme_event Test 7: mode and condition used in event()", {
  mydata <- tibble::tribble(
    ~USUBJID, ~CRIT1FL, ~ADY,
    "1",      "Y",         1,
    "1",      "Y",         2,
    "2",      "N",         1
  )

  expected <- tibble::tribble(
    ~USUBJID, ~CRIT1FL, ~ADY,
    "1",      "Y",         1
  )

  expect_dfs_equal(
    base = expected,
    compare = derive_extreme_event(
      source_datasets = list(mydata = mydata),
      by_vars = exprs(USUBJID),
      order = exprs(ADY),
      mode = "first",
      events = list(
        event(
          dataset_name = "mydata",
          condition = CRIT1FL == "Y",
          mode = "first"
        )
      )
    ),
    keys = "USUBJID"
  )
})

## Test 8: error if source dataset not available ----
test_that("derive_extreme_event Test 8: error if source dataset not available", {
  adhy <- tibble::tribble(
    ~USUBJID, ~AVISITN, ~CRIT1FL,
    "1",             1, "Y",
    "1",             2, "Y",
    "2",             1, "Y",
    "2",             2, NA_character_,
    "2",             3, "Y",
    "2",             4, NA_character_
  ) %>%
    mutate(
      PARAMCD = "ALKPH",
      PARAM = "Alkaline Phosphatase (U/L)"
    )

  expect_snapshot(
    derive_extreme_event(
      adhy,
      by_vars = exprs(USUBJID),
      events = list(
        event(
          dataset_name = "adyh",
          condition = is.na(CRIT1FL),
          set_values_to = exprs(AVALC = "N")
        ),
        event(
          condition = CRIT1FL == "Y",
          mode = "last",
          set_values_to = exprs(AVALC = "Y")
        )
      ),
      source_datasets = list(adhy = adhy),
      tmp_event_nr_var = event_nr,
      order = exprs(event_nr, AVISITN),
      mode = "first",
      keep_source_vars = exprs(AVISITN),
      set_values_to = exprs(
        PARAMCD = "ALK2",
        PARAM = "ALKPH <= 2 times ULN"
      )
    ),
    error = TRUE
  )
})

## Test 9: test for duplicates: one warning ----
test_that("derive_extreme_event Test 11: test for duplicates: one warning", {
  ad1 <- tribble(
    ~USUBJID, ~AVALC, ~ADY, ~ASEQ,
    "1",      "Y",       3,     1,
    "2",      "Y",       5,     1,
    "3",      "N",       2,     1,
    "4",      "N",       4,     1
  )

  ad2 <- tribble(
    ~USUBJID, ~AVALC, ~ADY, ~ASEQ,
    "1",      "Y",       3,     1,
    "2",      "Y",       3,     1,
    "3",      "Y",       2,     1
  )


  expect_warning(
    derive_extreme_event(
      by_vars = exprs(USUBJID),
      source_datasets = list(ad1 = ad1, ad2 = ad2),
      order = exprs(event_nr, ADY),
      mode = "first",
      check_type = "warning",
      tmp_event_nr_var = event_nr,
      events = list(
        event(
          dataset_name = "ad1",
          condition = AVALC == "Y",
          mode = "first",
          set_values_to = exprs(
            event_nr = 1,
            AVALC = "Y"
          )
        ),
        event(
          dataset_name = "ad2",
          condition = AVALC == "Y",
          mode = "first",
          set_values_to = exprs(
            event_nr = 1,
            AVALC = "Y"
          )
        ),
        event(
          dataset_name = "ad1",
          mode = "last",
          set_values_to = exprs(
            event_nr = 2,
            AVALC = "N"
          )
        )
      )
    ),
    "Check duplicates*"
  )
})


## Test 10: test for duplicates: with error ----
test_that("derive_extreme_event Test 12: test for duplicates: with error", {
  ad1 <- tribble(
    ~USUBJID, ~AVALC, ~ADY, ~ASEQ,
    "1",      "Y",       3,     1,
    "1",      "Y",       3,     2,
    "2",      "Y",       5,     1,
    "3",      "N",       2,     1,
    "4",      "N",       4,     1
  )

  ad2 <- tribble(
    ~USUBJID, ~AVALC, ~ADY, ~ASEQ,
    "1",      "Y",       3,     1,
    "2",      "Y",       3,     1,
    "3",      "Y",       2,     1
  )


  expect_error(
    derive_extreme_event(
      by_vars = exprs(USUBJID),
      source_datasets = list(ad1 = ad1, ad2 = ad2),
      order = exprs(event_nr, ADY),
      mode = "first",
      check_type = "error",
      tmp_event_nr_var = event_nr,
      events = list(
        event(
          dataset_name = "ad1",
          condition = AVALC == "Y",
          mode = "first",
          set_values_to = exprs(
            event_nr = 1,
            AVALC = "Y"
          )
        ),
        event(
          dataset_name = "ad2",
          condition = AVALC == "Y",
          mode = "first",
          set_values_to = exprs(
            event_nr = 1,
            AVALC = "Y"
          )
        ),
        event(
          dataset_name = "ad1",
          mode = "last",
          set_values_to = exprs(
            event_nr = 2,
            AVALC = "N"
          )
        )
      )
    ),
    NULL
  )
})
