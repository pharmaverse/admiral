context("test-derive_param_map")


test_that("new observations for MAP based on DIABP and SYSBP are derived correctly", {
  input <- tibble::tribble(
    ~USUBJID,      ~PARAMCD, ~PARAM,                            ~AVAL, ~AVALU, ~VISIT,
    "01-701-1015", "DIABP",  "Diastolic Blood Pressure (mmHg)",  51,   "mmHg", "BASELINE",
    "01-701-1015", "DIABP",  "Diastolic Blood Pressure (mmHg)",  50,   "mmHg", "WEEK 2",
    "01-701-1015", "SYSBP",  "Systolic Blood Pressure (mmHg)",  121,   "mmHg", "BASELINE",
    "01-701-1015", "SYSBP",  "Systolic Blood Pressure (mmHg)",  121,   "mmHg", "WEEK 2",
    "01-701-1028", "DIABP",  "Diastolic Blood Pressure (mmHg)",  79,   "mmHg", "BASELINE",
    "01-701-1028", "DIABP",  "Diastolic Blood Pressure (mmHg)",  80,   "mmHg", "WEEK 2",
    "01-701-1028", "SYSBP",  "Systolic Blood Pressure (mmHg)",  130,   "mmHg", "BASELINE",
    "01-701-1028", "SYSBP",  "Systolic Blood Pressure (mmHg)",  132,   "mmHg", "WEEK 2"
  )
  new_obs <-
    inner_join(input %>% filter(PARAMCD == "DIABP") %>% select(USUBJID, VISIT, AVAL, AVALU),
               input %>% filter(PARAMCD == "SYSBP") %>% select(USUBJID, VISIT, AVAL),
               by = c("USUBJID", "VISIT"),
               suffix = c(".DIABP", ".SYSBP")) %>%
    mutate(AVAL = (2 * AVAL.DIABP + AVAL.SYSBP) / 3,
           PARAMCD = "MAP") %>%
    select(-AVAL.DIABP, -AVAL.SYSBP)
  expected_output <- bind_rows(input, new_obs)

  expect_dfs_equal(
    derive_param_map(
      input,
      unit_var = AVALU,
      by_vars = vars(USUBJID, VISIT),
      set_values_to = vars(PARAMCD = "MAP")
    ),
    expected_output,
    keys = c("USUBJID", "PARAMCD", "VISIT")
  )
})

test_that("new observations for MAP based on DIABP and SYSBP are derived correctly without unit", {
  input <- tibble::tribble(
    ~USUBJID,      ~PARAMCD, ~PARAM,                            ~AVAL, ~VISIT,
    "01-701-1015", "DIABP",  "Diastolic Blood Pressure (mmHg)",  51,   "BASELINE",
    "01-701-1015", "DIABP",  "Diastolic Blood Pressure (mmHg)",  50,   "WEEK 2",
    "01-701-1015", "SYSBP",  "Systolic Blood Pressure (mmHg)",  121,   "BASELINE",
    "01-701-1015", "SYSBP",  "Systolic Blood Pressure (mmHg)",  121,   "WEEK 2",
    "01-701-1028", "DIABP",  "Diastolic Blood Pressure (mmHg)",  79,   "BASELINE",
    "01-701-1028", "DIABP",  "Diastolic Blood Pressure (mmHg)",  80,   "WEEK 2",
    "01-701-1028", "SYSBP",  "Systolic Blood Pressure (mmHg)",  130,   "BASELINE",
    "01-701-1028", "SYSBP",  "Systolic Blood Pressure (mmHg)",  132,   "WEEK 2"
  )
  new_obs <-
    inner_join(input %>% filter(PARAMCD == "DIABP") %>% select(USUBJID, VISIT, AVAL),
               input %>% filter(PARAMCD == "SYSBP") %>% select(USUBJID, VISIT, AVAL),
               by = c("USUBJID", "VISIT"),
               suffix = c(".DIABP", ".SYSBP")) %>%
    mutate(AVAL = (2 * AVAL.DIABP + AVAL.SYSBP) / 3,
           PARAMCD = "MAP") %>%
    select(-AVAL.DIABP, -AVAL.SYSBP)
  expected_output <- bind_rows(input, new_obs)

  expect_dfs_equal(
    derive_param_map(
      input,
      by_vars = vars(USUBJID, VISIT),
      set_values_to = vars(PARAMCD = "MAP")
    ),
    expected_output,
    keys = c("USUBJID", "PARAMCD", "VISIT")
  )
})

test_that("new observations for MAP based on DIABP, SYSBP, and HR are derived correctly", {
  input <- tibble::tribble(
    ~USUBJID,      ~PARAMCD, ~PARAM,                            ~AVAL, ~AVALU, ~VISIT,
    "01-701-1015", "PULSE",  "Pulse (beats/min)",                59,   "beats/min", "BASELINE",
    "01-701-1015", "PULSE",  "Pulse (beats/min)",                61,   "beats/min", "WEEK 2",
    "01-701-1015", "DIABP",  "Diastolic Blood Pressure (mmHg)",  51,   "mmHg",      "BASELINE",
    "01-701-1015", "DIABP",  "Diastolic Blood Pressure (mmHg)",  50,   "mmHg",      "WEEK 2",
    "01-701-1015", "SYSBP",  "Systolic Blood Pressure (mmHg)",  121,   "mmHg",      "BASELINE",
    "01-701-1015", "SYSBP",  "Systolic Blood Pressure (mmHg)",  121,   "mmHg",      "WEEK 2",
    "01-701-1028", "PULSE",  "Pulse (beats/min)",                62,   "beats/min", "BASELINE",
    "01-701-1028", "PULSE",  "Pulse (beats/min)",                77,   "beats/min", "WEEK 2",
    "01-701-1028", "DIABP",  "Diastolic Blood Pressure (mmHg)",  79,   "mmHg",      "BASELINE",
    "01-701-1028", "DIABP",  "Diastolic Blood Pressure (mmHg)",  80,   "mmHg",      "WEEK 2",
    "01-701-1028", "SYSBP",  "Systolic Blood Pressure (mmHg)",  130,   "mmHg",      "BASELINE",
    "01-701-1028", "SYSBP",  "Systolic Blood Pressure (mmHg)",  132,   "mmHg",      "WEEK 2"
  )
  new_obs <-
    inner_join(input %>% filter(PARAMCD == "DIABP") %>% select(USUBJID, VISIT, AVAL, AVALU),
               input %>% filter(PARAMCD == "SYSBP") %>% select(USUBJID, VISIT, AVAL),
               by = c("USUBJID", "VISIT"),
               suffix = c(".DIABP", ".SYSBP")) %>%
    inner_join(input %>% filter(PARAMCD == "PULSE") %>% select(USUBJID, VISIT, AVAL.PULSE = AVAL),
               by = c("USUBJID", "VISIT")) %>%
    mutate(AVAL = AVAL.DIABP + 0.01 * exp(4.14 - 40.74 / AVAL.PULSE) * (AVAL.SYSBP - AVAL.DIABP),
           PARAMCD = "MAP") %>%
    select(-AVAL.DIABP, -AVAL.SYSBP, -AVAL.PULSE)
  expected_output <- bind_rows(input, new_obs)

  expect_dfs_equal(
    derive_param_map(
      input,
      hr_code = "PULSE",
      unit_var = AVALU,
      by_vars = vars(USUBJID, VISIT),
      set_values_to = (vars(PARAMCD = "MAP"))
    ),
    expected_output,
    keys = c("USUBJID", "PARAMCD", "VISIT")
  )
})

test_that("an error is issued if PARAMCD is not set", {
  input <- tibble::tribble(
    ~USUBJID,      ~PARAMCD, ~PARAM,                            ~AVAL, ~AVALU, ~VISIT,
    "01-701-1015", "DIABP",  "Diastolic Blood Pressure (mmHg)",  51,   "mmHg", "BASELINE",
    "01-701-1015", "DIABP",  "Diastolic Blood Pressure (mmHg)",  50,   "mmHg", "WEEK 2",
    "01-701-1015", "SYSBP",  "Systolic Blood Pressure (mmHg)",  121,   "mmHg", "BASELINE",
    "01-701-1015", "SYSBP",  "Systolic Blood Pressure (mmHg)",  121,   "mmHg", "WEEK 2",
    "01-701-1028", "DIABP",  "Diastolic Blood Pressure (mmHg)",  79,   "mmHg", "BASELINE",
    "01-701-1028", "DIABP",  "Diastolic Blood Pressure (mmHg)",  80,   "mmHg", "WEEK 2",
    "01-701-1028", "SYSBP",  "Systolic Blood Pressure (mmHg)",  130,   "mmHg", "BASELINE",
    "01-701-1028", "SYSBP",  "Systolic Blood Pressure (mmHg)",  132,   "mmHg", "WEEK 2"
  )

  expect_error(
    derive_param_map(
      input,
      unit_var = AVALU,
      by_vars = vars(USUBJID, VISIT),
      set_values_to = vars(PARAM = "Mean Arterial Pressure")
    ),
    "The following required elements are missing in `set_values_to`: 'PARAMCD'"
  )
})
