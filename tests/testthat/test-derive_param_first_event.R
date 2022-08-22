library(tibble)
adsl <- tribble(
  ~USUBJID, ~DTHDT,
  "1",      ymd("2022-05-13"),
  "2",      ymd(""),
  "3",      ymd("")
) %>%
  mutate(STUDYID = "XX1234")

adrs <- tribble(
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

# derive_param_first_event ----
## derive_param_first_event Test 1: derive first PD date ----
test_that("derive_param_first_event Test 1: derive first PD date", {
  library(tibble)
  actual <- derive_param_first_event(
    adrs,
    dataset_adsl = adsl,
    dataset_source = adrs,
    filter_source = PARAMCD == "OVR" & AVALC == "PD",
    date_var = ADT,
    set_values_to = vars(
      PARAMCD = "PD",
      ANL01FL = "Y"
    )
  )

  expected <- bind_rows(
    adrs,
    tribble(
      ~USUBJID, ~ADT,              ~AVALC, ~AVAL,
      "1",      ymd(""),           "N",    0,
      "2",      ymd("2021-07-16"), "Y",    1,
      "3",      ymd(""),           "N",    0
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

## derive_param_first_event Test 2: derive death date parameter ----
test_that("derive_param_first_event Test 2: derive death date parameter", {
  library(tibble)
  actual <- derive_param_first_event(
    dataset = adrs,
    dataset_adsl = adsl,
    dataset_source = adsl,
    filter_source = !is.na(DTHDT),
    date_var = DTHDT,
    set_values_to = vars(
      PARAMCD = "DEATH",
      ANL01FL = "Y"
    )
  )

  expected <- bind_rows(
    adrs,
    tribble(
      ~USUBJID, ~ADT,              ~AVALC, ~AVAL, ~DTHDT,
      "1",      ymd("2022-05-13"), "Y",    1,     ymd("2022-05-13"),
      "2",      ymd(""),           "N",    0,     ymd(""),
      "3",      ymd(""),           "N",    0,     ymd("")
    ) %>%
      mutate(
        STUDYID = "XX1234",
        PARAMCD = "DEATH",
        ANL01FL = "Y"
      )
  )

  expect_dfs_equal(
    base = expected,
    comp = actual,
    keys = c("USUBJID", "PARAMCD", "ADT")
  )
})
