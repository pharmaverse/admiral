# derive_extreme_records ----
## Test 1: add last observation for each group ----
test_that("derive_extreme_records Test 1: add last observation for each group", {
  input <- tibble::tribble(
    ~USUBJID, ~AVISITN, ~AVAL, ~LBSEQ,
    "1",             1,    12,      1,
    "1",             3,     9,      2,
    "2",             2,    42,      1,
    "3",             3,    14,      1,
    "3",             3,    10,      2
  )

  expected_output <- bind_rows(
    input,
    tibble::tribble(
      ~USUBJID, ~AVISITN, ~AVAL, ~LBSEQ,
      "1",             3,     9,      2,
      "2",             2,    42,      1,
      "3",             3,    10,      2
    ) %>%
      mutate(DTYPE = "LOV")
  )

  actual_output <- derive_extreme_records(
    input,
    order = exprs(AVISITN, LBSEQ),
    by_vars = exprs(USUBJID),
    mode = "last",
    set_values_to = exprs(DTYPE = "LOV")
  )

  expect_dfs_equal(
    base = expected_output,
    compare = actual_output,
    keys = c("USUBJID", "AVISITN", "LBSEQ", "DTYPE")
  )
})

## Test 2: derive first PD date ----
test_that("derive_extreme_records Test 2: derive first PD date", {
  adsl <- tibble::tribble(
    ~USUBJID,
    "1",
    "2",
    "3"
  ) %>%
    mutate(STUDYID = "XX1234")

  adrs <- tibble::tribble(
    ~USUBJID, ~ADTC,        ~AVALC, ~PARAMCD,
    "1",      "2020-01-02", "PR",   "OVR",
    "1",      "2020-02-01", "CR",   "OVR",
    "1",      "2020-03-01", "CR",   "OVR",
    "1",      "2020-04-01", "SD",   "OVR",
    "2",      "2021-06-15", "SD",   "OVR",
    "2",      "2021-07-16", "PD",   "OVR",
    "2",      "2021-09-14", "PD",   "OVR",
    "1",      "2020-01-02", "PR",   "OVRF",
    "1",      "2020-02-01", "CR",   "OVRF",
    "1",      "2020-03-01", "CR",   "OVRF",
    "1",      "2020-04-01", "PD",   "OVRF",
    "2",      "2021-06-15", "SD",   "OVRF",
    "2",      "2021-07-16", "PD",   "OVRF",
    "2",      "2021-09-14", "PD",   "OVRF"
  ) %>%
    mutate(
      STUDYID = "XX1234",
      ADT = ymd(ADTC)
    ) %>%
    select(-ADTC)

  actual <- derive_extreme_records(
    adrs,
    dataset_ref = adsl,
    dataset_add = adrs,
    by_vars = exprs(USUBJID),
    filter_add = PARAMCD == "OVR" & AVALC == "PD",
    exist_flag = AVALC,
    order = exprs(ADT),
    mode = "first",
    set_values_to = exprs(
      PARAMCD = "PD",
      ANL01FL = "Y",
      ADT = ADT
    )
  )

  expected <- bind_rows(
    adrs,
    tibble::tribble(
      ~USUBJID, ~ADT,              ~AVALC,
      "1",      ymd(""),           "N",
      "2",      ymd("2021-07-16"), "Y",
      "3",      ymd(""),           "N"
    ) %>%
      mutate(
        STUDYID = "XX1234",
        PARAMCD = "PD",
        ANL01FL = "Y"
      )
  )

  expect_dfs_equal(
    base = expected,
    comp = actual,
    keys = c("USUBJID", "PARAMCD", "ADT")
  )
})

## Test 3: derive death date parameter ----
test_that("derive_extreme_records Test 3: derive death date parameter", {
  adsl <- tibble::tribble(
    ~USUBJID, ~DTHDT,
    "1",      ymd("2022-05-13"),
    "2",      ymd(""),
    "3",      ymd("")
  ) %>%
    mutate(STUDYID = "XX1234")

  actual <- derive_extreme_records(
    dataset_ref = adsl,
    dataset_add = adsl,
    by_vars = exprs(STUDYID, USUBJID),
    filter_add = !is.na(DTHDT),
    exist_flag = AVAL,
    true_value = 1,
    false_value = 0,
    mode = "first",
    set_values_to = exprs(
      PARAMCD = "DEATH",
      ANL01FL = "Y",
      ADT = DTHDT
    )
  )

  expected <- tibble::tribble(
    ~USUBJID, ~ADT,              ~AVAL, ~DTHDT,
    "1",      ymd("2022-05-13"), 1,     ymd("2022-05-13"),
    "2",      ymd(""),           0,     ymd(""),
    "3",      ymd(""),           0,     ymd("")
  ) %>%
    mutate(
      STUDYID = "XX1234",
      PARAMCD = "DEATH",
      ANL01FL = "Y"
    )

  expect_dfs_equal(
    base = expected,
    comp = actual,
    keys = c("USUBJID", "PARAMCD", "ADT")
  )
})

## Test 4: latest evaluable tumor assessment date parameter ----
test_that("derive_extreme_records Test 4: latest evaluable tumor assessment date parameter", {
  adsl <- tibble::tribble(
    ~USUBJID,
    "1",
    "2",
    "3"
  ) %>%
    mutate(STUDYID = "XX1234")

  adrs <- tibble::tribble(
    ~USUBJID, ~ADTC,        ~AVALC, ~PARAMCD,
    "1",      "2020-01-02", "PR",   "OVR",
    "1",      "2020-02-01", "CR",   "OVR",
    "1",      "2020-03-01", "NE",   "OVR",
    "1",      "2020-04-01", "SD",   "OVR",
    "2",      "2021-06-15", "SD",   "OVR",
    "2",      "2021-07-16", "SD",   "OVR",
    "2",      "2021-09-14", "NE",   "OVR",
    "3",      "2021-08-03", "NE",   "OVR",
    "1",      "2020-01-02", "PR",   "OVRF",
    "1",      "2020-02-01", "CR",   "OVRF",
    "1",      "2020-03-01", "NE",   "OVRF",
    "1",      "2020-04-01", "SD",   "OVRF",
    "2",      "2021-06-15", "SD",   "OVRF",
    "2",      "2021-07-16", "SD",   "OVRF",
    "2",      "2021-09-14", "NE",   "OVRF",
    "3",      "2021-08-03", "NE",   "OVRF"
  ) %>%
    mutate(
      STUDYID = "XX1234",
      ADT = ymd(ADTC)
    ) %>%
    select(-ADTC)

  actual <- derive_extreme_records(
    dataset = adrs,
    dataset_ref = adsl,
    dataset_add = adrs,
    by_vars = exprs(STUDYID, USUBJID),
    filter_add = PARAMCD == "OVR" & AVALC != "NE",
    order = exprs(ADT),
    exist_flag = AVALC,
    true_value = "Y",
    false_value = "N",
    mode = "last",
    set_values_to = exprs(
      PARAMCD = "LSTEVLDT",
      ANL01FL = "Y",
      ADT = ADT
    )
  )

  expected <- bind_rows(
    adrs,
    tibble::tribble(
      ~USUBJID, ~ADT,              ~AVALC,
      "1",      ymd("2020-04-01"), "Y",
      "2",      ymd("2021-07-16"), "Y",
      "3",      ymd(""),           "N"
    ) %>%
      mutate(
        STUDYID = "XX1234",
        PARAMCD = "LSTEVLDT",
        ANL01FL = "Y"
      )
  )

  expect_dfs_equal(
    base = expected,
    comp = actual,
    keys = c("USUBJID", "PARAMCD", "ADT")
  )
})

## Test 5: latest evaluable tumor assessment date parameter without overwriting existing result ----
test_that("derive_extreme_records Test 5: latest evaluable tumor assessment date parameter without overwriting existing result", { # nolint
  adsl <- tibble::tribble(
    ~USUBJID,
    "1",
    "2",
    "3"
  ) %>%
    mutate(STUDYID = "XX1234")

  adrs <- tibble::tribble(
    ~USUBJID, ~ADTC,        ~AVALC, ~PARAMCD,
    "1",      "2020-01-02", "PR",   "OVR",
    "1",      "2020-02-01", "CR",   "OVR",
    "1",      "2020-03-01", "NE",   "OVR",
    "1",      "2020-04-01", "SD",   "OVR",
    "2",      "2021-06-15", "SD",   "OVR",
    "2",      "2021-07-16", "SD",   "OVR",
    "2",      "2021-09-14", "NE",   "OVR",
    "3",      "2021-08-03", "NE",   "OVR",
    "1",      "2020-01-02", "PR",   "OVRF",
    "1",      "2020-02-01", "CR",   "OVRF",
    "1",      "2020-03-01", "NE",   "OVRF",
    "1",      "2020-04-01", "SD",   "OVRF",
    "2",      "2021-06-15", "SD",   "OVRF",
    "2",      "2021-07-16", "SD",   "OVRF",
    "2",      "2021-09-14", "NE",   "OVRF",
    "3",      "2021-08-03", "NE",   "OVRF"
  ) %>%
    mutate(
      STUDYID = "XX1234",
      ADT = ymd(ADTC)
    ) %>%
    select(-ADTC)

  actual <- derive_extreme_records(
    dataset = adrs,
    dataset_ref = adsl,
    dataset_add = adrs,
    by_vars = exprs(STUDYID, USUBJID),
    filter_add = PARAMCD == "OVR" & AVALC != "NE",
    order = exprs(ADT),
    mode = "last",
    set_values_to = exprs(
      PARAMCD = "LSTEVLDT",
      ANL01FL = "Y",
      ADT = ADT
    )
  )

  expected <- bind_rows(
    adrs,
    tibble::tribble(
      ~USUBJID, ~ADT,              ~AVALC,
      "1",      ymd("2020-04-01"), "SD",
      "2",      ymd("2021-07-16"), "SD",
      "3",      ymd(""),           NA
    ) %>%
      mutate(
        STUDYID = "XX1234",
        PARAMCD = "LSTEVLDT",
        ANL01FL = "Y"
      )
  )

  expect_dfs_equal(
    base = expected,
    comp = actual,
    keys = c("USUBJID", "PARAMCD", "ADT")
  )
})

## Test 6: warning if filter argument is used ----
test_that("derive_extreme_records Test 6: warning if filter argument is used", {
  adsl <- tibble::tribble(
    ~USUBJID,
    "1",
    "2",
    "3"
  ) %>%
    mutate(STUDYID = "XX1234")

  adrs <- tibble::tribble(
    ~USUBJID, ~ADTC,        ~AVALC, ~PARAMCD,
    "1",      "2020-01-02", "PR",   "OVR",
    "1",      "2020-02-01", "CR",   "OVR",
    "1",      "2020-03-01", "NE",   "OVR",
    "1",      "2020-04-01", "SD",   "OVR",
    "2",      "2021-06-15", "SD",   "OVR",
    "2",      "2021-07-16", "SD",   "OVR",
    "2",      "2021-09-14", "NE",   "OVR",
    "3",      "2021-08-03", "NE",   "OVR",
  ) %>%
    mutate(
      STUDYID = "XX1234",
      ADT = ymd(ADTC)
    ) %>%
    select(-ADTC)

  actual <- derive_extreme_records(
    adrs,
    dataset_ref = adsl,
    dataset_add = adrs,
    by_vars = exprs(USUBJID),
    filter_add = PARAMCD == "OVR" & AVALC == "PD",
    exist_flag = AVALC,
    order = exprs(ADT),
    mode = "first",
    set_values_to = exprs(
      PARAMCD = "PD",
      ANL01FL = "Y",
      ADT = ADT
    )
  )

  expect_error(
    derive_extreme_records(
      adrs,
      dataset_ref = adsl,
      dataset_add = adrs,
      by_vars = exprs(USUBJID),
      filter = PARAMCD == "OVR" & AVALC == "PD",
      exist_flag = AVALC,
      order = exprs(ADT),
      mode = "first",
      set_values_to = exprs(
        PARAMCD = "PD",
        ANL01FL = "Y",
        ADT = ADT
      )
    ),
    class = "lifecycle_error_deprecated"
  )
})

## Test 7: error if no input data ----
test_that("derive_extreme_records Test 7: error if no input data", {
  expect_error(
    derive_extreme_records(
      set_values_to = exprs(PARAMCD = "HELLO")
    ),
    regexp = paste(
      "Neither `dataset` nor `dataset_add` is specified.",
      "At least one of them must be specified.",
      sep = "\n"
    ),
    fixed = TRUE
  )
})

## Test 8: keep vars in `keep_source_vars` in the new records ----
test_that("derive_extreme_records Test 8: keep vars in `keep_source_vars` in the new records", {
  input <- tibble::tribble(
    ~USUBJID, ~AVISITN, ~AVAL, ~LBSEQ,
    1, 1, 12, 1,
    1, 3, 9, 2,
    2, 2, 42, 1,
    3, 3, 14, 1,
    3, 3, 10, 2
  )

  expected_output <- bind_rows(
    input,
    tibble::tribble(
      ~USUBJID, ~AVISITN, ~AVAL, ~LBSEQ,
      1, 3, 9, 2,
      2, 2, 42, 1,
      3, 3, 10, 2
    ) %>%
      select(USUBJID, AVISITN, AVAL) %>%
      mutate(DTYPE = "LOV")
  )

  actual_output <- derive_extreme_records(
    input,
    order = exprs(AVISITN, LBSEQ),
    by_vars = exprs(USUBJID),
    mode = "last",
    keep_source_vars = exprs(AVISITN, AVAL),
    set_values_to = exprs(DTYPE = "LOV")
  )

  expect_dfs_equal(
    base = expected_output,
    compare = actual_output,
    keys = c("USUBJID", "AVISITN", "LBSEQ", "DTYPE")
  )
})

## Test 9: keep all vars in the new records when `keep_source_vars` is 'exprs(everything())' ----
test_that("derive_extreme_records Test 9: keep all vars in the new records when `keep_source_vars` is 'exprs(everything())'", { # nolint
  input <- tibble::tribble(
    ~USUBJID, ~AVISITN, ~AVAL, ~LBSEQ,
    1, 1, 12, 1,
    1, 3, 9, 2,
    2, 2, 42, 1,
    3, 3, 14, 1,
    3, 3, 10, 2
  )

  expected_output <- bind_rows(
    input,
    tibble::tribble(
      ~USUBJID, ~AVISITN, ~AVAL, ~LBSEQ,
      1, 3, 9, 2,
      2, 2, 42, 1,
      3, 3, 10, 2
    ) %>%
      mutate(DTYPE = "LOV")
  )

  actual_output <- derive_extreme_records(
    input,
    order = exprs(AVISITN, LBSEQ),
    by_vars = exprs(USUBJID),
    mode = "last",
    keep_source_vars = exprs(everything()),
    set_values_to = exprs(DTYPE = "LOV")
  )

  expect_dfs_equal(
    base = expected_output,
    compare = actual_output,
    keys = c("USUBJID", "AVISITN", "LBSEQ", "DTYPE")
  )
})
