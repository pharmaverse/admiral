adsl <- tibble::tribble(
  ~USUBJID, ~DTHDT,
  "1",      ymd("2022-05-13"),
  "2",      ymd(""),
  "3",      ymd("")
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

# derive_param_extreme_event ----
## Test 1: derive first PD date ----
test_that("derive_param_extreme_event Test 1: derive first PD date", {
  actual <- derive_param_extreme_event(
    adrs,
    dataset_adsl = adsl,
    dataset_source = adrs,
    filter_source = PARAMCD == "OVR" & AVALC == "PD",
    order = vars(ADT),
    set_values_to = vars(
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

## Test 2: derive death date parameter ----
test_that("derive_param_extreme_event Test 2: derive death date parameter", {
  actual <- derive_param_extreme_event(
    dataset_adsl = adsl,
    dataset_source = adsl,
    filter_source = !is.na(DTHDT),
    new_var = AVAL,
    true_value = 1,
    false_value = 0,
    mode = "first",
    set_values_to = vars(
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

## Test 3: derive latest evaluable tumor assessment date parameter ----
test_that("derive_param_extreme_event Test 3: latest evaluable tumor assessment date parameter", {
  actual <- derive_param_extreme_event(
    dataset = adrs,
    dataset_adsl = adsl,
    dataset_source = adrs,
    filter_source = PARAMCD == "OVR" & AVALC != "NE",
    order = vars(ADT),
    new_var = AVALC,
    true_value = "Y",
    false_value = "N",
    mode = "last",
    set_values_to = vars(
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
