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
## Test 1: deprecation error if function is called ----
test_that("derive_param_extreme_event Test 1: deprecation error if function is called", {
  expect_error(
    derive_param_extreme_event(
      adrs,
      dataset_adsl = adsl,
      dataset_source = adrs,
      filter_source = PARAMCD == "OVR" & AVALC == "PD",
      new_var = AVALC,
      order = exprs(ADT),
      set_values_to = exprs(
        PARAMCD = "PD",
        ANL01FL = "Y",
        ADT = ADT
      )
    ),
    class = "lifecycle_error_deprecated"
  )
})
