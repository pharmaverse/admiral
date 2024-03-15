# derive_param_extreme_record ----
## Test 1: Analysis date are derived correctly ----
test_that("derive_param_extreme_record Test 1: Analysis date are derived correctly", {
  aevent <- tibble::tribble(
    ~STUDYID, ~USUBJID,     ~LBSTDTC, ~PARAMCD, ~PARAM,
    "1001",        "1", "2023-01-01",    "TST", "TEST",
    "1001",        "2", "2023-01-01",    "TST", "TEST",
    "1001",        "3", "2023-01-01",    "TST", "TEST"
  )

  cm <- tibble::tribble(
    ~STUDYID, ~USUBJID, ~CMDECOD,     ~CMSTDTC,
    "1001",        "1",    "ACT", "2020-12-25"
  )

  pr <- tibble::tribble(
    ~STUDYID, ~USUBJID, ~PRDECOD,     ~PRSTDTC,
    "1001",        "1",    "ACS", "2020-12-27",
    "1001",        "2",    "ACS", "2021-12-25",
    "1001",        "3",    "ACS", "2022-12-25",
  )
  expected_output <- tibble::tribble(
    ~STUDYID, ~USUBJID,     ~LBSTDTC,   ~PARAMCD,                      ~PARAM,                         ~ADT, ~AVALC, # nolint
    "1001",        "1", "2023-01-01",      "TST",                      "TEST",                           NA,     NA, # nolint
    "1001",        "2", "2023-01-01",      "TST",                      "TEST",                           NA,     NA, # nolint
    "1001",        "3", "2023-01-01",      "TST",                      "TEST",                           NA,     NA, # nolint
    "1001",        "1",           NA, "FIRSTACT", "First Anti-Cancer Therapy", lubridate::ymd("2020-12-25"),  "ACT", # nolint
    "1001",        "2",           NA, "FIRSTACT", "First Anti-Cancer Therapy", lubridate::ymd("2021-12-25"),  "ACS", # nolint
    "1001",        "3",           NA, "FIRSTACT", "First Anti-Cancer Therapy", lubridate::ymd("2022-12-25"),  "ACS" # nolint
  )
  actual_output <- derive_param_extreme_record(
    dataset = aevent,
    sources = list(
      records_source(
        dataset_name = "cm",
        filter = CMDECOD == "ACT",
        new_vars = exprs(
          ADT = convert_dtc_to_dt(CMSTDTC),
          AVALC = CMDECOD
        )
      ),
      records_source(
        dataset_name = "pr",
        filter = PRDECOD == "ACS",
        new_vars = exprs(
          ADT = convert_dtc_to_dt(PRSTDTC),
          AVALC = PRDECOD
        )
      )
    ),
    source_datasets = list(cm = cm, pr = pr),
    by_vars = exprs(STUDYID, USUBJID),
    order = exprs(ADT),
    mode = "first",
    set_values_to = exprs(
      PARAMCD = "FIRSTACT",
      PARAM = "First Anti-Cancer Therapy"
    )
  )

  expect_dfs_equal(expected_output, actual_output, keys = c("USUBJID", "PARAMCD", "PARAM", "ADT", "AVALC")) # nolint
})

## Test 2: Error given when order variable is not inside source datasets ----
test_that("derive_param_extreme_record Test 2: Error given when order variable is not inside source datasets", { # nolint
  aevent <- tibble::tribble(
    ~STUDYID, ~USUBJID,     ~LBSTDTC, ~PARAMCD, ~PARAM,
    "1001",        "1", "2023-01-01",    "TST", "TEST",
    "1001",        "2", "2023-01-01",    "TST", "TEST",
    "1001",        "3", "2023-01-01",    "TST", "TEST"
  )

  cm <- tibble::tribble(
    ~STUDYID, ~USUBJID, ~CMDECOD,     ~CMSTDTC,
    "1001",        "1",    "ACT", "2020-12-25"
  )

  pr <- tibble::tribble(
    ~STUDYID, ~USUBJID, ~PRDECOD,     ~PRSTDTC,
    "1001",        "1",    "ACS", "2020-12-27",
    "1001",        "2",    "ACS", "2021-12-25",
    "1001",        "3",    "ACS", "2022-12-25",
  )
  expect_error(
    derive_param_extreme_record(
      dataset = aevent,
      sources = list(
        records_source(
          dataset_name = "cm",
          filter = CMDECOD == "ACT",
          new_vars = exprs(
            ADT = convert_dtc_to_dt(CMSTDTC),
            AVALC = CMDECOD
          )
        ),
        records_source(
          dataset_name = "pr",
          filter = PRDECOD == "ACS",
          new_vars = exprs(
            ADT = convert_dtc_to_dt(PRSTDTC),
            AVALC = PRDECOD
          )
        )
      ),
      source_datasets = list(cm = cm, pr = pr),
      by_vars = exprs(STUDYID, USUBJID),
      order = exprs(ADT2),
      mode = "first",
      set_values_to = exprs(
        PARAMCD = "FIRSTACT",
        PARAM = "First Anti-Cancer Therapy"
      )
    ),
    regexp = "Required variable `ADT2` is missing"
  )
})

## Test 3: Error given when sources is not in proper list format ----
test_that("derive_param_extreme_record Test 3: Error given when sources is not in proper list format", { # nolint
  aevent <- tibble::tribble(
    ~STUDYID, ~USUBJID,     ~LBSTDTC, ~PARAMCD, ~PARAM,
    "1001",        "1", "2023-01-01",    "TST", "TEST",
    "1001",        "2", "2023-01-01",    "TST", "TEST",
    "1001",        "3", "2023-01-01",    "TST", "TEST"
  )

  pr <- tibble::tribble(
    ~STUDYID, ~USUBJID, ~PRDECOD,     ~PRSTDTC,
    "1001",        "1",    "ACS", "2020-12-27",
    "1001",        "2",    "ACS", "2021-12-25",
    "1001",        "3",    "ACS", "2022-12-25",
  )
  expect_error(
    derive_param_extreme_record(
      dataset = aevent,
      sources = records_source(
        dataset_name = "pr",
        filter = PRDECOD == "ACS",
        new_vars = exprs(
          ADT = convert_dtc_to_dt(PRSTDTC),
          AVALC = PRDECOD
        )
      ),
      source_datasets = list(pr = pr),
      by_vars = exprs(STUDYID, USUBJID),
      order = exprs(ADT),
      mode = "first",
      set_values_to = exprs(
        PARAMCD = "FIRSTACT",
        PARAM = "First Anti-Cancer Therapy"
      )
    ),
    class = "assert_list_of"
  )
})

## Test 4: Non-existent/missing values are accounted for ----
test_that("derive_param_extreme_record Test 4: Non-existent/missing values are accounted for", {
  aevent <- tibble::tribble(
    ~STUDYID, ~USUBJID,     ~LBSTDTC, ~PARAMCD, ~PARAM,
    "1001",        "1", "2023-01-01",    "TST", "TEST",
    "1001",        "2", "2023-01-01",    "TST", "TEST",
    "1001",        "3", "2023-01-01",    "TST", "TEST"
  )

  pr <- tibble::tribble(
    ~STUDYID, ~USUBJID, ~PRDECOD,     ~PRSTDTC,
    "1001",        "2",    "ACS", "2021-12-25",
    "1001",        "3",    "ACS", "2022-12-25",
  )
  expected_output <- tibble::tribble(
    ~STUDYID, ~USUBJID,     ~LBSTDTC,   ~PARAMCD,                      ~PARAM,                         ~ADT, ~AVALC, # nolint
    "1001",        "1", "2023-01-01",      "TST",                      "TEST",                           NA,     NA, # nolint
    "1001",        "2", "2023-01-01",      "TST",                      "TEST",                           NA,     NA, # nolint
    "1001",        "3", "2023-01-01",      "TST",                      "TEST",                           NA,     NA, # nolint
    "1001",        "2",           NA, "FIRSTACT", "First Anti-Cancer Therapy", lubridate::ymd("2021-12-25"),  "ACS", # nolint
    "1001",        "3",           NA, "FIRSTACT", "First Anti-Cancer Therapy", lubridate::ymd("2022-12-25"),  "ACS" # nolint
  )
  actual_output <- derive_param_extreme_record(
    dataset = aevent,
    sources = list(
      records_source(
        dataset_name = "pr",
        filter = PRDECOD == "ACS",
        new_vars = exprs(
          ADT = convert_dtc_to_dt(PRSTDTC),
          AVALC = PRDECOD
        )
      )
    ),
    source_datasets = list(pr = pr),
    by_vars = exprs(STUDYID, USUBJID),
    order = exprs(ADT),
    mode = "first",
    set_values_to = exprs(
      PARAMCD = "FIRSTACT",
      PARAM = "First Anti-Cancer Therapy"
    )
  )

  expect_dfs_equal(expected_output, actual_output, keys = c("USUBJID", "PARAMCD", "PARAM", "ADT", "AVALC")) # nolint
})
