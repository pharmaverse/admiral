adsl <- tibble::tribble(
  ~USUBJID,
  "1",
  "2",
  "3"
) %>%
  mutate(STUDYID = "XX1234")

adrs <- tibble::tribble(
  ~USUBJID, ~AVALC, ~PARAMCD,
  "1",      "PR",   "OVR",
  "1",      "CR",   "OVR",
  "1",      "CR",   "OVR",
  "1",      "SD",   "OVR",
  "2",      "SD",   "OVR",
  "2",      "PD",   "OVR",
  "2",      "PD",   "OVR",
  "1",      "PR",   "OVRF",
  "1",      "CR",   "OVRF",
  "1",      "CR",   "OVRF",
  "1",      "PD",   "OVRF",
  "2",      "SD",   "OVRF",
  "2",      "PD",   "OVRF",
  "2",      "PD",   "OVRF"
) %>%
  mutate(
    STUDYID = "XX1234",
    ANL01FL = "Y"
  )

# derive_param_merged_exist_flag ----
## Test 1: derive parameter indicating PD ----
test_that("derive_param_merged_exist_flag Test 1: derive parameter indicating PD", {
  actual <- derive_param_exist_flag(
    dataset_ref = adsl,
    dataset_add = adrs,
    filter_add = PARAMCD == "OVR",
    condition = AVALC == "PD",
    false_value = "N",
    set_values_to = exprs(
      AVAL = yn_to_numeric(AVALC),
      PARAMCD = "PD",
      ANL01FL = "Y"
    )
  )

  expected <- tibble::tribble(
    ~USUBJID, ~AVALC,        ~AVAL,
    "1",      "N",           0,
    "2",      "Y",           1,
    "3",      NA_character_, NA_real_
  ) %>%
    mutate(
      STUDYID = "XX1234",
      PARAMCD = "PD",
      ANL01FL = "Y"
    )

  expect_dfs_equal(
    base = expected,
    comp = actual,
    keys = c("USUBJID", "PARAMCD")
  )
})



## Test 2: error is issued if parameter already exists in dataset ----
test_that("derive_param_merged_exist_flag Test 2: error is issued if parameter already exists in dataset", { # nolint
  expect_error(
    derive_param_exist_flag(
      dataset = adrs,
      dataset_ref = adsl,
      dataset_add = adrs,
      filter_add = PARAMCD == "OVR",
      condition = AVALC == "PD",
      false_value = "N",
      set_values_to = exprs(
        PARAMCD = "OVR",
        ANL01FL = "Y"
      )
    ),
    class = "assert_param_does_not_exist"
  )
})
