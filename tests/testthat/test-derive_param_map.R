# compute_map ----

## compute_map: DBP & SBP ----
# ((2 x DBP) + SBP) / 3

## Test 1: MAP based on diastolic & systolic BP - single values ----
test_that("compute_map Test 1: MAP based on diastolic & systolic BP - single values", {
  expect_equal(round(compute_map(diabp = 51, sysbp = 121), 3L), 74.333)
})

## Test 2: MAP based on diastolic & systolic BP - vectors ----
test_that("compute_map Test 2: MAP based on diastolic & systolic BP - vectors", {
  expect_equal(
    round(compute_map(diabp = c(51, 61), sysbp = c(121, 141)), 3L), c(74.333, 87.667)
  )
})

## Test 3: MAP based on diastolic & systolic BP with missing values ----
test_that("compute_map Test 3: MAP based on diastolic & systolic BP with missing values", {
  expect_equal(
    compute_map(diabp = c(NA, 61), sysbp = c(121, NA)), c(NA_real_, NA_real_)
  )
})

## compute_map: DBP, SBP & HR ----
# DBP + 0.01 x exp(4.14 - 40.74 / PULSE) x (SBP - DBP)

## Test 4: MAP based on DBP & SBP & heart rate - single values ----
test_that("compute_map Test 4: MAP based on DBP & SBP & heart rate - single values", {
  expect_equal(
    round(compute_map(diabp = 51, sysbp = 121, hr = 59), 3L), 73.039
  )
})

## Test 5: MAP based on diastolic, systolic BP & heart rate - vectors ----
test_that("compute_map Test 5: MAP based on diastolic, systolic BP & heart rate - vectors", {
  expect_equal(
    round(compute_map(diabp = c(51, 91), sysbp = c(121, 101), hr = c(59, 62)), 3L),
    c(73.039, 94.255)
  )
})

## Test 6: MAP based on DBP, SBP & heart rate - with missing values ----
test_that("compute_map Test 6: MAP based on DBP, SBP & heart rate - with missing values", {
  expect_equal(
    compute_map(diabp = c(NA, 61, 51), sysbp = c(121, NA, 121), hr = c(59, 62, NA)),
    c(NA_real_, NA_real_, NA_real_)
  )
})

# derive_param_map ----

## derive_param_map: Error checks ----

## Test 7: MAP parameter NOT added - wrong DIABP unit ----
test_that("derive_param_map Test 7: MAP parameter NOT added - wrong DIABP unit", {
  input <- tibble::tribble(
    ~USUBJID, ~PARAMCD, ~PARAM, ~AVAL, ~VISIT,
    "01-701-1015", "DIABP", "Diastolic Blood Pressure (mHg)", 51, "BASELINE",
    "01-701-1015", "SYSBP", "Systolic Blood Pressure (mmHg)", 121, "BASELINE",
  )

  expect_error(
    derive_param_map(
      input,
      by_vars = exprs(USUBJID, VISIT),
      get_unit_expr = extract_unit(PARAM)
    ),
    class = "assert_unit"
  )
})

## Test 8: MAP parameter NOT added - wrong SYSBP unit ----
test_that("derive_param_map Test 8: MAP parameter NOT added - wrong SYSBP unit", {
  input <- tibble::tribble(
    ~USUBJID,      ~PARAMCD, ~PARAM,                            ~AVAL, ~VISIT,
    "01-701-1015", "DIABP",  "Diastolic Blood Pressure (mmHg)",    51, "BASELINE",
    "01-701-1015", "SYSBP",  "Systolic Blood Pressure (mHg)",     121, "BASELINE",
  )

  expect_error(
    derive_param_map(
      input,
      by_vars = exprs(USUBJID, VISIT),
      get_unit_expr = extract_unit(PARAM)
    ),
    class = "assert_unit"
  )
})

## Test 9: MAP parameter NOT added - wrong PULSE unit ----
test_that("derive_param_map Test 9: MAP parameter NOT added - wrong PULSE unit", {
  input <- tibble::tribble(
    ~USUBJID, ~PARAMCD, ~PARAM, ~AVAL, ~VISIT,
    "01-701-1015", "DIABP", "Diastolic Blood Pressure (mmHg)", 51, "BASELINE",
    "01-701-1015", "SYSBP", "Systolic Blood Pressure (mmHg)", 121, "BASELINE",
    "01-701-1015", "PULSE", "Pulse (beats/m)", 65, "BASELINE",
  )

  expect_error(
    derive_param_map(
      input,
      by_vars = exprs(USUBJID, VISIT),
      hr_code = "PULSE",
      get_unit_expr = extract_unit(PARAM)
    ),
    class = "assert_unit"
  )
})

## Test 10: MAP parameter NOT added - PARAMCD not set ----
test_that("derive_param_map Test 10: MAP parameter NOT added - PARAMCD not set", {
  input <- tibble::tribble(
    ~USUBJID, ~PARAMCD, ~PARAM, ~AVAL, ~VISIT,
    "01-701-1015", "DIABP", "Diastolic Blood Pressure (mmHg)", 51, "BASELINE",
    "01-701-1015", "SYSBP", "Systolic Blood Pressure (mmHg)", 121, "BASELINE",
  )

  expect_error(
    derive_param_map(
      input,
      by_vars = exprs(USUBJID, VISIT),
      set_values_to = exprs(PARAM = "Mean Arterial Pressure"),
      get_unit_expr = extract_unit(PARAM)
    ),
    class = "assert_varval_list"
  )
})

## derive_param_map: No obs added ----

## Test 11: MAP parameter NOT added ----
test_that("derive_param_map Test 11: MAP parameter NOT added", {
  expected_output <- tibble::tribble(
    ~USUBJID,      ~PARAMCD, ~PARAM,                            ~AVAL, ~VISIT,
    "01-701-1015", "DIABP",  "Diastolic Blood Pressure (mmHg)",    NA, "BASELINE",
    "01-701-1015", "DIABP",  "Diastolic Blood Pressure (mmHg)",    50, "WEEK 2",
    "01-701-1015", "DIABP",  "Diastolic Blood Pressure (mmHg)",    55, "WEEK 4",
    "01-701-1015", "SYSBP",  "Systolic Blood Pressure (mmHg)",    121, "BASELINE",
    "01-701-1015", "SYSBP",  "Systolic Blood Pressure (mmHg)",     NA, "WEEK 2",
    "01-701-1015", "SYSBP",  "Systolic Blood Pressure (mmHg)",    127, "WEEK 4",
    "01-701-1015", "PULSE",  "Pulse (beats/min)",                  65, "BASELINE",
    "01-701-1015", "PULSE",  "Pulse (beats/min)",                  68, "WEEK 2",
    "01-701-1015", "PULSE",  "Pulse (beats/min)",                  NA, "WEEK 4",
  )

  input <- expected_output

  expect_snapshot(
    result <- derive_param_map(
      input,
      by_vars = exprs(USUBJID, VISIT),
      hr_code = "PULSE",
      get_unit_expr = extract_unit(PARAM)
    )
  )

  expect_dfs_equal(
    result,
    expected_output,
    keys = c("USUBJID", "PARAMCD", "VISIT")
  )
})

## derive_param_map: Obs created ----

maphr <- function(sbp, dbp, hr) {
  dbp + 0.01 * exp(4.14 - 40.74 / hr) * (sbp - dbp)
}

## Test 12: MAP parameter (DBP/SBP/PULSE) is correctly added ----
test_that("derive_param_map Test 12: MAP parameter (DBP/SBP/PULSE) is correctly added", {
  expected_output <- tibble::tribble(
    ~USUBJID, ~PARAMCD, ~PARAM, ~VISIT, ~AVAL,
    "01-701-1015", "PULSE", "Pulse (beats/min)", "BASELINE", 59,
    "01-701-1015", "DIABP", "Diastolic Blood Pressure (mmHg)", "BASELINE", 51,
    "01-701-1015", "SYSBP", "Systolic Blood Pressure (mmHg)", "BASELINE", 121,
    # New row added for MAP for SUBJID="01-701-1015" and VISIT="BASELINE"
    # PULSE = 59 DIABP = 51 and SYSBP = 121
    "01-701-1015", "MAP", NA, "BASELINE", maphr(121, 51, 59),
    "01-701-1028", "PULSE", "Pulse (beats/min)", "WEEK 2", 61,
    "01-701-1028", "DIABP", "Diastolic Blood Pressure (mmHg)", "WEEK 2", 50,
    "01-701-1028", "SYSBP", "Systolic Blood Pressure (mmHg)", "WEEK 2", 125,
    # New row added for MAP for SUBJID="01-701-1028" and VISIT="WEEK 2"
    # PULSE = 61 DIABP = 50 and SYSBP = 125
    "01-701-1028", "MAP", NA, "WEEK 2", maphr(125, 50, 61),
  )

  input <- expected_output %>% filter(PARAMCD != "MAP")

  expect_dfs_equal(
    derive_param_map(
      input,
      by_vars = exprs(USUBJID, VISIT),
      hr_code = "PULSE",
      get_unit_expr = extract_unit(PARAM)
    ),
    expected_output,
    keys = c("USUBJID", "PARAMCD", "VISIT")
  )
})

map <- function(sbp, dbp) {
  (2 * dbp + sbp) / 3
}

## Test 13: MAP parameter (DBP/SBP) is correctly added ----
test_that("derive_param_map Test 13: MAP parameter (DBP/SBP) is correctly added", {
  expected_output <- tibble::tribble(
    ~USUBJID, ~PARAMCD, ~PARAM, ~VISIT, ~AVAL,
    "01-701-1015", "DIABP", "Diastolic Blood Pressure (mmHg)", "BASELINE", 51,
    "01-701-1015", "SYSBP", "Systolic Blood Pressure (mmHg)", "BASELINE", 121,
    # New row added for MAP for SUBJID="01-701-1015" and VISIT="BASELINE"
    # DIABP = 51 and SYSBP = 121
    "01-701-1015", "MAP", NA, "BASELINE", map(121, 51),
    "01-701-1028", "DIABP", "Diastolic Blood Pressure (mmHg)", "WEEK 2", 50,
    "01-701-1028", "SYSBP", "Systolic Blood Pressure (mmHg)", "WEEK 2", 125,
    # New row added for MAP for SUBJID="01-701-1028" and VISIT="WEEK 2"
    # DIABP = 50 and SYSBP = 125
    "01-701-1028", "MAP", NA, "WEEK 2", map(125, 50),
  )

  input <- expected_output %>% filter(PARAMCD != "MAP")

  expect_dfs_equal(
    derive_param_map(input,
      by_vars = exprs(USUBJID, VISIT),
      get_unit_expr = extract_unit(PARAM)
    ),
    expected_output,
    keys = c("USUBJID", "PARAMCD", "VISIT")
  )
})
