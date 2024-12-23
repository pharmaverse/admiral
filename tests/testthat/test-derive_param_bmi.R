# compute_bmi ----

## Test 1: BMI calculation - single height & weight values ----
test_that("compute_bmi Test 1: BMI calculation - single height & weight values", {
  # Expected values are taken from the Center of Disease Control and Prevention's
  # (CDC) 'Adult BMI Calculator' at
  # https://cdc.gov/healthyweight/assessing/bmi/adult_bmi/metric_bmi_calculator/bmi_calculator.html
  expect_equal(round(compute_bmi(height = 180, weight = 75), 3L), 23.148)
})

## Test 2: compute_bmi BMI calculation - height & weight vectors ----
test_that("compute_bmi Test 2: compute_bmi BMI calculation - height & weight vectors", {
  expect_equal(
    round(compute_bmi(height = c(180, 200), weight = c(75, 100)), 3L),
    c(23.148, 25)
  )
})

## Test 3: BMI height & weight vectors - missing values ----
test_that("compute_bmi Test 3: BMI height & weight vectors - missing values", {
  expect_equal(
    compute_bmi(height = c(NA, 200, 0), weight = c(75, NA, 75)),
    c(NA_real_, NA_real_, NA_real_)
  )
})

# derive_param_bmi ----

## derive_param_bmi: Error checks ----

## Test 4: BMI parameter NOT added - wrong hgt unit ----
test_that("derive_param_bmi Test 4: BMI parameter NOT added - wrong hgt unit", {
  input <- tibble::tribble(
    ~USUBJID, ~PARAMCD, ~PARAM, ~VISIT, ~VSSTRESU, ~AVAL,
    # Wrong unit for HEIGHT should be cm
    "01-701-1015", "HEIGHT", "Height (cm)", "BASELINE", "m", 170,
    "01-701-1015", "WEIGHT", "Weight (kg)", "BASELINE", "kg", 85,
  )

  expect_error(
    derive_param_bmi(input, by_vars = exprs(USUBJID, VISIT), get_unit_expr = VSSTRESU),
    class = "assert_unit"
  )
})

## Test 5: BMI parameter NOT added - wrong wgt unit ----
test_that("derive_param_bmi Test 5: BMI parameter NOT added - wrong wgt unit", {
  input <- tibble::tribble(
    ~USUBJID,      ~PARAMCD, ~PARAM,        ~VISIT,     ~VSSTRESU, ~AVAL,
    "01-701-1015", "HEIGHT", "Height (cm)", "BASELINE",      "cm",   170,
    # Wrong unit for WEIGHT should be kg
    "01-701-1015", "WEIGHT", "Weight (kg)", "BASELINE",       "g",    85,
  )

  expect_error(
    derive_param_bmi(input, by_vars = exprs(USUBJID, VISIT), get_unit_expr = VSSTRESU),
    class = "assert_unit"
  )
})

## Test 6: BMI parameter NOT added - multiple unit for wgt ----
test_that("derive_param_bmi Test 6: BMI parameter NOT added - multiple unit for wgt", {
  input <- tibble::tribble(
    ~USUBJID,      ~PARAMCD, ~PARAM,        ~VISIT,     ~VSSTRESU, ~AVAL,
    "01-701-1015", "HEIGHT", "Height (cm)", "BASELINE",      "cm",   170,
    # Multiple units for WEIGHT
    "01-701-1015", "WEIGHT", "Weight (kg)", "BASELINE",      "kg",    85,
    "01-701-1016", "WEIGHT", "Weight (kg)", "BASELINE",       "g",  8500,
  )

  expect_error(
    derive_param_bmi(input, by_vars = exprs(USUBJID, VISIT), get_unit_expr = VSSTRESU),
    class = "assert_unit"
  )
})

## Test 7: BMI parameter NOT added - PARAMCD not set ----
test_that("derive_param_bmi Test 7: BMI parameter NOT added - PARAMCD not set", {
  input <- tibble::tribble(
    ~USUBJID,      ~PARAMCD, ~PARAM,        ~VISIT,     ~VSSTRESU, ~AVAL,
    "01-701-1015", "HEIGHT", "Height (cm)", "BASELINE",      "cm",   170,
    "01-701-1015", "WEIGHT", "Weight (kg)", "BASELINE",      "kg",    85,
  )

  expect_error(
    derive_param_bmi(
      input,
      by_vars = exprs(USUBJID, VISIT),
      set_values_to = exprs(PARAM = "Body Mass Index"),
      get_unit_expr = VSSTRESU
    ),
    class = "assert_varval_list"
  )
})

## derive_param_bmi: No obs added  ----

## Test 8: BMI parameter NOT added ----
test_that("derive_param_bmi Test 8: BMI parameter NOT added", {
  expected_output <- tibble::tribble(
    ~USUBJID, ~PARAMCD, ~PARAM, ~VISIT, ~VSSTRESU, ~AVAL,
    "01-701-1015", "HEIGHT", "Height (cm)", "BASELINE", "cm", 170,
    # WEIGHT set to NA - so BMI not calculated
    "01-701-1015", "WEIGHT", "Weight (kg)", "BASELINE", "kg", NA,
    "01-701-1015", "WEIGHT", "Weight (kg)", "MONTH 1", "kg", 78,
    # HEIGHT set to NA - so BMI not calculated
    "01-701-1028", "HEIGHT", "Height (cm)", "BASELINE", "cm", NA,
    "01-701-1028", "WEIGHT", "Weight (kg)", "BASELINE", "kg", 90,
    "01-701-1028", "HEIGHT", "Height (cm)", "MONTH 1", "cm", 88,
  )

  input <- expected_output

  expect_snapshot(
    result <- derive_param_bmi(
      input,
      by_vars = exprs(USUBJID, VISIT),
      get_unit_expr = VSSTRESU
    )
  )

  expect_dfs_equal(
    result,
    expected_output,
    keys = c("USUBJID", "PARAMCD", "VISIT")
  )
})

## derive_param_bmi: Obs created ----

bmi <- function(hgt, wgt) {
  wgt / (hgt / 100)^2
}

## Test 9: BMI parameter is correctly added ----
test_that("derive_param_bmi Test 9: BMI parameter is correctly added", {
  expected_output <- tibble::tribble(
    ~USUBJID, ~PARAMCD, ~PARAM, ~VISIT, ~VSSTRESU, ~AVAL,
    "01-701-1015", "HEIGHT", "Height (cm)", "BASELINE", "cm", 170,
    "01-701-1015", "WEIGHT", "Weight (kg)", "BASELINE", "kg", 75,
    # New row added for BMI for SUBJID="01-701-1015" and VISIT="BASELINE"
    # WEIGHT = 75 and HEIGHT = 170
    "01-701-1015", "BMI", NA, "BASELINE", NA, bmi(170, 75),
    "01-701-1015", "WEIGHT", "Weight (kg)", "MONTH 1", "kg", 78,
    "01-701-1028", "HEIGHT", "Height (cm)", "BASELINE", "cm", 185,
    "01-701-1028", "WEIGHT", "Weight (kg)", "BASELINE", "kg", 90,
    # New row added for BMI for SUBJID="01-701-1028" and VISIT='BASELINE'
    # WEIGHT = 90 and HEIGHT = 185
    "01-701-1028", "BMI", NA, "BASELINE", NA, bmi(185, 90),
    "01-701-1028", "WEIGHT", "Weight (kg)", "MONTH 1", "kg", 88,
  )

  input <- expected_output %>% filter(PARAMCD != "BMI")

  expect_dfs_equal(
    derive_param_bmi(input, by_vars = exprs(USUBJID, VISIT), get_unit_expr = VSSTRESU),
    expected_output,
    keys = c("USUBJID", "PARAMCD", "VISIT")
  )
})


# Derive BMI where height is measured only once
## Test 10: Derive BMI where height is measured only once ----
test_that("derive_param_bmi Test 10: Derive BMI where height is measured only once", {
  input <- tibble::tribble(
    ~USUBJID,      ~PARAMCD, ~PARAM,        ~AVAL, ~AVALU, ~VISIT,
    "01-701-1015", "HEIGHT", "Height (cm)", 147.0, "cm",   "SCREENING",
    "01-701-1015", "WEIGHT", "Weight (kg)",  54.0, "kg",   "SCREENING",
    "01-701-1015", "WEIGHT", "Weight (kg)",  54.4, "kg",   "BASELINE",
    "01-701-1015", "WEIGHT", "Weight (kg)",  53.1, "kg",   "WEEK 2",
    "01-701-1028", "HEIGHT", "Height (cm)", 163.0, "cm",   "SCREENING",
    "01-701-1028", "WEIGHT", "Weight (kg)",  78.5, "kg",   "SCREENING",
    "01-701-1028", "WEIGHT", "Weight (kg)",  80.3, "kg",   "BASELINE",
    "01-701-1028", "WEIGHT", "Weight (kg)",  80.7, "kg",   "WEEK 2"
  )

  expected_output <- derive_param_computed(
    input,
    by_vars = exprs(USUBJID, VISIT),
    parameters = "WEIGHT",
    set_values_to = exprs(
      AVAL = AVAL.WEIGHT / (AVAL.HEIGHT / 100)^2,
      PARAMCD = "BMI",
      PARAM = "Body Mass Index (kg/m^2)",
      AVALU = "kg/m^2"
    ),
    constant_parameters = c("HEIGHT"),
    constant_by_vars = exprs(USUBJID)
  )

  expect_dfs_equal(
    expected_output,
    derive_param_bmi(
      input,
      by_vars = exprs(USUBJID, VISIT),
      weight_code = "WEIGHT",
      height_code = "HEIGHT",
      set_values_to = exprs(
        PARAMCD = "BMI",
        PARAM = "Body Mass Index (kg/m^2)",
        AVALU = "kg/m^2"
      ),
      get_unit_expr = extract_unit(PARAM),
      constant_by_vars = exprs(USUBJID)
    ),
    keys = c("USUBJID", "PARAMCD", "VISIT")
  )
})
