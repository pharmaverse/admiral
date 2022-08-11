test_that("new observations are derived correctly", {
  input <- tibble::tribble(
    ~USUBJID, ~PARAMCD, ~PARAM, ~AVAL, ~AVALU, ~VISIT,
    "01-701-1015", "DIABP", "Diastolic Blood Pressure (mmHg)", 51, "mmHg", "BASELINE",
    "01-701-1015", "DIABP", "Diastolic Blood Pressure (mmHg)", 50, "mmHg", "WEEK 2",
    "01-701-1015", "SYSBP", "Systolic Blood Pressure (mmHg)", 121, "mmHg", "BASELINE",
    "01-701-1015", "SYSBP", "Systolic Blood Pressure (mmHg)", 121, "mmHg", "WEEK 2",
    "01-701-1028", "DIABP", "Diastolic Blood Pressure (mmHg)", 79, "mmHg", "BASELINE",
    "01-701-1028", "DIABP", "Diastolic Blood Pressure (mmHg)", 80, "mmHg", "WEEK 2",
    "01-701-1028", "SYSBP", "Systolic Blood Pressure (mmHg)", 130, "mmHg", "BASELINE",
    "01-701-1028", "SYSBP", "Systolic Blood Pressure (mmHg)", 132, "mmHg", "WEEK 2"
  )

  new_obs <-
    inner_join(input %>% filter(PARAMCD == "DIABP") %>% select(USUBJID, VISIT, AVAL),
      input %>% filter(PARAMCD == "SYSBP") %>% select(USUBJID, VISIT, AVAL),
      by = c("USUBJID", "VISIT"),
      suffix = c(".DIABP", ".SYSBP")
    ) %>%
    mutate(
      AVAL = (2 * AVAL.DIABP + AVAL.SYSBP) / 3,
      PARAMCD = "MAP",
      PARAM = "Mean arterial pressure (mmHg)",
      AVALU = "mmHg"
    ) %>%
    select(-AVAL.DIABP, -AVAL.SYSBP)
  expected_output <- bind_rows(input, new_obs)

  expect_dfs_equal(
    derive_param_computed(
      input,
      parameters = c("SYSBP", "DIABP"),
      by_vars = vars(USUBJID, VISIT),
      analysis_value = (AVAL.SYSBP + 2 * AVAL.DIABP) / 3,
      set_values_to = vars(
        PARAMCD = "MAP",
        PARAM = "Mean arterial pressure (mmHg)",
        AVALU = "mmHg"
      )
    ),
    expected_output,
    keys = c("USUBJID", "PARAMCD", "VISIT")
  )
})

test_that("new observations are derived correctly with constant parameters", {
  input <- tibble::tribble(
    ~USUBJID, ~PARAMCD, ~PARAM, ~AVAL, ~AVALU, ~VISIT,
    "01-701-1015", "HEIGHT", "Height (cm)", 147, "cm", "SCREENING",
    "01-701-1015", "WEIGHT", "Weight (kg)", 54.0, "kg", "SCREENING",
    "01-701-1015", "WEIGHT", "Weight (kg)", 54.4, "kg", "BASELINE",
    "01-701-1015", "WEIGHT", "Weight (kg)", 53.1, "kg", "WEEK 2",
    "01-701-1028", "HEIGHT", "Height (cm)", 163, "cm", "SCREENING",
    "01-701-1028", "WEIGHT", "Weight (kg)", 78.5, "kg", "SCREENING",
    "01-701-1028", "WEIGHT", "Weight (kg)", 80.3, "kg", "BASELINE",
    "01-701-1028", "WEIGHT", "Weight (kg)", 80.7, "kg", "WEEK 2"
  )

  new_obs <-
    inner_join(input %>% filter(PARAMCD == "HEIGHT") %>% select(USUBJID, AVAL),
      input %>% filter(PARAMCD == "WEIGHT") %>% select(USUBJID, VISIT, AVAL),
      by = c("USUBJID"),
      suffix = c(".HEIGHT", ".WEIGHT")
    ) %>%
    mutate(
      AVAL = AVAL.WEIGHT / (AVAL.HEIGHT / 100)^2,
      PARAMCD = "BMI",
      PARAM = "Body Mass Index (kg/m2)",
      AVALU = "kg/m2"
    ) %>%
    select(-AVAL.HEIGHT, -AVAL.WEIGHT)
  expected_output <- bind_rows(input, new_obs)

  expect_dfs_equal(
    derive_param_computed(
      input,
      parameters = c("WEIGHT"),
      by_vars = vars(USUBJID, VISIT),
      constant_parameters = c("HEIGHT"),
      constant_by_vars = vars(USUBJID),
      analysis_value = AVAL.WEIGHT / (AVAL.HEIGHT / 100)^2,
      set_values_to = vars(
        PARAMCD = "BMI",
        PARAM = "Body Mass Index (kg/m2)",
        AVALU = "kg/m2"
      )
    ),
    expected_output,
    keys = c("USUBJID", "PARAMCD", "VISIT")
  )
})

test_that("no new observations are added if filtered dataset is empty", {
  input <- tibble::tribble(
    ~USUBJID, ~PARAMCD, ~PARAM, ~AVAL, ~AVALU, ~VISIT,
    "01-701-1015", "DIABP", "Diastolic Blood Pressure (mmHg)", 51, "mmHg", "BASELINE",
    "01-701-1015", "DIABP", "Diastolic Blood Pressure (mmHg)", 50, "mmHg", "WEEK 2",
    "01-701-1015", "SYSBP", "Systolic Blood Pressure (mmHg)", 121, "mmHg", "BASELINE",
    "01-701-1015", "SYSBP", "Systolic Blood Pressure (mmHg)", 121, "mmHg", "WEEK 2",
    "01-701-1028", "DIABP", "Diastolic Blood Pressure (mmHg)", 79, "mmHg", "BASELINE",
    "01-701-1028", "DIABP", "Diastolic Blood Pressure (mmHg)", 80, "mmHg", "WEEK 2",
    "01-701-1028", "SYSBP", "Systolic Blood Pressure (mmHg)", 130, "mmHg", "BASELINE",
    "01-701-1028", "SYSBP", "Systolic Blood Pressure (mmHg)", 132, "mmHg", "WEEK 2"
  )

  expect_warning(
    derive_param_computed(
      input,
      filter = VISIT == "WEEK 24",
      parameters = c("SYSBP", "DIABP"),
      by_vars = vars(USUBJID, VISIT),
      analysis_value = (AVAL.SYSBP + 2 * AVAL.DIABP) / 3,
      set_values_to = vars(
        PARAMCD = "MAP",
        PARAM = "Mean arterial pressure (mmHg)",
        AVALU = "mmHg"
      )
    ),
    "The input dataset does not contain any observations fullfiling the filter condition .*"
  ) %>%
    expect_dfs_equal(input,
      keys = c("USUBJID", "PARAMCD", "VISIT")
    )
})

test_that("no new observations are added if a parameter is missing", {
  input <- tibble::tribble(
    ~USUBJID, ~PARAMCD, ~PARAM, ~AVAL, ~AVALU, ~VISIT,
    "01-701-1015", "DIABP", "Diastolic Blood Pressure (mmHg)", 51, "mmHg", "BASELINE",
    "01-701-1015", "DIABP", "Diastolic Blood Pressure (mmHg)", 50, "mmHg", "WEEK 2",
    "01-701-1015", "SYSBP", "Systolic Blood Pressure (mmHg)", 121, "mmHg", "BASELINE",
    "01-701-1015", "SYSBP", "Systolic Blood Pressure (mmHg)", 121, "mmHg", "WEEK 2",
    "01-701-1028", "DIABP", "Diastolic Blood Pressure (mmHg)", 79, "mmHg", "BASELINE",
    "01-701-1028", "DIABP", "Diastolic Blood Pressure (mmHg)", 80, "mmHg", "WEEK 2",
    "01-701-1028", "SYSBP", "Systolic Blood Pressure (mmHg)", 130, "mmHg", "BASELINE",
    "01-701-1028", "SYSBP", "Systolic Blood Pressure (mmHg)", 132, "mmHg", "WEEK 2"
  )

  expect_warning(
    derive_param_computed(
      input,
      filter = PARAMCD == "DIABP",
      parameters = c("SYSBP", "DIABP"),
      by_vars = vars(USUBJID, VISIT),
      analysis_value = (AVAL.SYSBP + 2 * AVAL.DIABP) / 3,
      set_values_to = vars(
        PARAMCD = "MAP",
        PARAM = "Mean arterial pressure (mmHg)",
        AVALU = "mmHg"
      )
    ),
    "The input dataset does not contain any observations fullfiling the filter condition .*"
  ) %>%
    expect_dfs_equal(input,
      keys = c("USUBJID", "PARAMCD", "VISIT")
    )
})
