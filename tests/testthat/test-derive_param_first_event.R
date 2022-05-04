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

# derive_param_first_event ----
## derive_param_first_event Test 1: derive first PD date ----
test_that("derive_param_first_event Test 1: derive first PD date", {
  actual <- derive_param_first_event(
    adrs,
    dataset_adsl = adsl,
    source_param = "OVR",
    condition = AVALC == "PD",
    set_values_to = vars(
      PARAMCD = "PD",
      ANL01FL = "Y"
    )
  )

  expected <- bind_rows(
    adrs,
    tibble::tribble(
      ~USUBJID, ~ADT,              ~AVALC, ~AVALN,
      "1",      ymd(""),           "N",    0,
      "2",      ymd("2021-07-16"), "Y",    1,
      "3",      ymd(""),           "N",    0) %>%
      mutate(
        STUDYID = "XX1234",
        PARAMCD = "PD",
        ANL01FL = "Y")
    )

  expect_dfs_equal(base = expected,
                   comp = actual,
                   keys = c("USUBJID", "PARAMCD", "ADT"))
})
