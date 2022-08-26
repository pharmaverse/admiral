input <- tibble::tribble(
  # nolint start
  ~USUBJID, ~PARAMCD, ~PARAM, ~AVAL, ~AVALU, ~VISIT, ~AGE, ~SEX, ~SMOKEFL, ~DIABETFL, ~TRTHYPFL,
  "01-701-1015", "SYSBP", "Systolic Blood Pressure (mmHg)", 121, "mmHg", "BASELINE", 44, "F", "N", "N", "N",
  "01-701-1015", "SYSBP", "Systolic Blood Pressure (mmHg)", 115, "mmHg", "WEEK 2", 44, "F", "N", "N", "Y",
  "01-701-1015", "CHOL", "Total Cholesterol (mg/dL)", 216.16, "mg/dL", "BASELINE", 44, "F", "N", "N", "N",
  "01-701-1015", "CHOL", "Total Cholesterol (mg/dL)", 210.78, "mg/dL", "WEEK 2", 44, "F", "N", "N", "Y",
  "01-701-1015", "CHOLHDL", "Cholesteral/HDL-Cholesterol (mg/dL)", 54.91, "mg/dL", "BASELINE", 44, "F", "N", "N", "N",
  "01-701-1015", "CHOLHDL", "Cholesteral/HDL-Cholesterol (mg/dL)", 26.72, "mg/dL", "WEEK 2", 44, "F", "N", "N", "Y",
  "01-701-1028", "SYSBP", "Systolic Blood Pressure (mmHg)", 119, "mmHg", "BASELINE", 55, "M", "Y", "Y", "Y",
  "01-701-1028", "SYSBP", "Systolic Blood Pressure (mmHg)", 101, "mmHg", "WEEK 2", 55, "M", "Y", "Y", "Y",
  "01-701-1028", "CHOL", "Total Cholesterol (mg/dL)", 292.01, "mg/dL", "BASELINE", 55, "M", "Y", "Y", "Y",
  "01-701-1028", "CHOL", "Total Cholesterol (mg/dL)", 246.73, "mg/dL", "WEEK 2", 55, "M", "Y", "Y", "Y",
  "01-701-1028", "CHOLHDL", "Cholesteral/HDL-Cholesterol (mg/dL)", 65.55, "mg/dL", "BASELINE", 55, "M", "Y", "Y", "Y",
  "01-701-1028", "CHOLHDL", "Cholesteral/HDL-Cholesterol (mg/dL)", 44.62, "mg/dL", "WEEK 2", 55, "M", "Y", "Y", "Y"
  # nolint end
)

test_that("derive_param_framingham Test 1: New observations are derived correctly", {
  new_obs <- select(
    filter(input, PARAMCD == "SYSBP"), USUBJID, VISIT, AVAL, AVALU, AGE, SEX, SMOKEFL, DIABETFL, TRTHYPFL # nolint
  ) %>%
    rename(AVAL_SYSBP = AVAL) %>%
    left_join(select(
      filter(input, PARAMCD == "CHOL"), USUBJID, VISIT, AVAL
    ), by = c("USUBJID", "VISIT")) %>%
    rename(AVAL_CHOL = AVAL) %>%
    left_join(select(
      filter(input, PARAMCD == "CHOLHDL"), USUBJID, VISIT, AVAL
    ), by = c("USUBJID", "VISIT")) %>%
    rename(AVAL_CHOLHDL = AVAL) %>%
    mutate(
      AVAL = compute_framingham(
        sysbp = AVAL_SYSBP,
        chol = AVAL_CHOL,
        cholhdl = AVAL_CHOLHDL,
        age = AGE,
        sex = SEX,
        smokefl = SMOKEFL,
        diabetfl = DIABETFL,
        trthypfl = TRTHYPFL
      ),
      PARAMCD = "FCVD101",
      PARAM = "FCVD1-Framingham CVD 10-Year Risk Score (%)",
      AVALU = NA_character_
    ) %>%
    select(-AVAL_CHOLHDL, -AVAL_CHOL, -AVAL_SYSBP)
  expected_output <- bind_rows(input, new_obs)

  expect_dfs_equal(
    derive_param_framingham(
      input,
      by_vars = vars(USUBJID, VISIT),
      set_values_to = vars(
        PARAMCD = "FCVD101",
        PARAM = "FCVD1-Framingham CVD 10-Year Risk Score (%)"
      ),
      get_unit_expr = AVALU
    ),
    expected_output,
    keys = c("USUBJID", "PARAMCD", "VISIT")
  )
})
