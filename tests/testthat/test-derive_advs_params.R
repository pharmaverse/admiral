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

# compute_bsa ----

## compute_bsa: Mosteller method ----
# sqrt (Height x Weight / 3600)

## Test 4: Mosteller method - single height & weight values ----
test_that("compute_bsa Test 4: Mosteller method - single height & weight values", {
  expect_equal(
    round(compute_bsa(height = 170, weight = 75, method = "Mosteller"), 3L),
    1.882
  )
})

## Test 5: Mosteller method - height & weight vectors ----
test_that("compute_bsa Test 5: Mosteller method - height & weight vectors", {
  expect_equal(
    round(compute_bsa(height = c(170, 185), weight = c(75, 90), method = "Mosteller"), 3L),
    c(1.882, 2.151)
  )
})

## Test 6: Mosteller method - height & weight vectors - missing values ----
test_that("compute_bsa Test 6: Mosteller method - height & weight vectors - missing values", {
  expect_equal(
    compute_bsa(height = c(NA, 185), weight = c(75, NA), method = "Mosteller"),
    c(NA_real_, NA_real_)
  )
})

## compute_bsa: DuBois-DuBois method ----
#  FORMULA : 0.20247 x (HGT/100)^0.725 x WGT^0.425

## Test 7: DuBois-DuBois method - single height & weight values ----
test_that("compute_bsa Test 7: DuBois-DuBois method - single height & weight values", {
  expect_equal(
    round(compute_bsa(height = 170, weight = 75, method = "DuBois-DuBois"), 3L),
    1.864
  )
})

## Test 8: DuBois-DuBois method - height & weight vectors ----
test_that("compute_bsa Test 8: DuBois-DuBois method - height & weight vectors", {
  expect_equal(
    round(compute_bsa(height = c(170, 185), weight = c(75, 90), method = "DuBois-DuBois"), 3L),
    c(1.864, 2.141)
  )
})

## Test 9: DuBois-DuBois method - hgt and wgt vectors - missing values ----
test_that("compute_bsa Test 9: DuBois-DuBois method - hgt and wgt vectors - missing values", {
  expect_equal(
    compute_bsa(height = c(NA, 185), weight = c(75, NA), method = "DuBois-DuBois"),
    c(NA_real_, NA_real_)
  )
})

## compute_bsa: Haycock method (Test 03.xx) ----
# 0.024265 x HGT^0.3964 x WGT^0.5378

## Test 10: Haycock method - single height & weight values ----
test_that("compute_bsa Test 10: Haycock method - single height & weight values", {
  expect_equal(
    round(compute_bsa(height = 170, weight = 75, method = "Haycock"), 3L),
    1.895
  )
})

## Test 11: Haycock method - height & weight vectors ----
test_that("compute_bsa Test 11: Haycock method - height & weight vectors", {
  expect_equal(
    round(compute_bsa(height = c(170, 185), weight = c(75, 90), method = "Haycock"), 3L),
    c(1.895, 2.161)
  )
})

## Test 12: Haycock method - height & weight vectors - missing values ----
test_that("compute_bsa Test 12: Haycock method - height & weight vectors - missing values", {
  expect_equal(
    compute_bsa(height = c(NA, 185), weight = c(75, NA), method = "Haycock"),
    c(NA_real_, NA_real_)
  )
})

## compute_bsa: Gehan-George method ----
# 0.0235 x HGT^0.42246 x WGT^0.51456

## Test 13: Gehan-George method - single height & weight values ----
test_that("compute_bsa Test 13: Gehan-George method - single height & weight values", {
  expect_equal(
    round(compute_bsa(height = 170, weight = 75, method = "Gehan-George"), 3L),
    1.897
  )
})

## Test 14: Gehan-George method - height & weight vectors ----
test_that("compute_bsa Test 14: Gehan-George method - height & weight vectors", {
  expect_equal(
    round(compute_bsa(height = c(170, 185), weight = c(75, 90), method = "Gehan-George"), 3L),
    c(1.897, 2.16)
  )
})

## Test 15: Gehan-George method - height & weight vectors - missing values ----
test_that("compute_bsa Test 15: Gehan-George method - height & weight vectors - missing values", {
  expect_equal(
    compute_bsa(height = c(NA, 185), weight = c(75, NA), method = "Gehan-George"),
    c(NA_real_, NA_real_)
  )
})

## compute_bsa: Boyd method ----
# 0.0003207 x (HGT^0.3) x (1000 x WGT)^(0.7285 - (0.0188 x log10(1000 x WGT)))

## Test 16: Boyd method - single height & weight values ----
test_that("compute_bsa Test 16: Boyd method - single height & weight values", {
  expect_equal(
    round(compute_bsa(height = 170, weight = 75, method = "Boyd"), 3L),
    1.905
  )
})

## Test 17: Boyd method - height & weight vectors ----
test_that("compute_bsa Test 17: Boyd method - height & weight vectors", {
  expect_equal(
    round(compute_bsa(height = c(170, 185), weight = c(75, 90), method = "Boyd"), 3L),
    c(1.905, 2.158)
  )
})

## Test 18: Boyd method - height & weight vectors - missing values ----
test_that("compute_bsa Test 18: Boyd method - height & weight vectors - missing values", {
  expect_equal(
    compute_bsa(height = c(NA, 185), weight = c(75, NA), method = "Boyd"),
    c(NA_real_, NA_real_)
  )
})

## compute_bsa: Fujimoto method ----
# 0.008883 x HGT^0.663 x WGT^0.444

## Test 19: Fujimoto method - single height & weight values ----
test_that("compute_bsa Test 19: Fujimoto method - single height & weight values", {
  expect_equal(
    round(compute_bsa(height = 170, weight = 75, method = "Fujimoto"), 3L),
    1.819
  )
})

## Test 20: Fujimoto method - height & weight vectors ----
test_that("compute_bsa Test 20: Fujimoto method - height & weight vectors", {
  expect_equal(
    round(compute_bsa(height = c(170, 185), weight = c(75, 90), method = "Fujimoto"), 3L),
    c(1.819, 2.086)
  )
})

## Test 21: Fujimoto method - height & weight vectors - missing values ----
test_that("compute_bsa Test 21: Fujimoto method - height & weight vectors - missing values", {
  expect_equal(
    compute_bsa(height = c(NA, 185), weight = c(75, NA), method = "Fujimoto"),
    c(NA_real_, NA_real_)
  )
})

## compute_bsa: Takahira method ----
# 0.007241 x HGT^0.725 x WGT^0.425

## Test 22: Takahira method - single height & weight values ----
test_that("compute_bsa Test 22: Takahira method - single height & weight values", {
  expect_equal(
    round(compute_bsa(height = 170, weight = 75, method = "Takahira"), 3L),
    1.878
  )
})

## Test 23: Takahira method - height & weight vectors ----
test_that("compute_bsa Test 23: Takahira method - height & weight vectors", {
  expect_equal(
    round(compute_bsa(height = c(170, 185), weight = c(75, 90), method = "Takahira"), 3L),
    c(1.878, 2.158)
  )
})

## Test 24: Takahira method - height & weight vectors - missing values ----
test_that("compute_bsa Test 24: Takahira method - height & weight vectors - missing values", {
  expect_equal(
    compute_bsa(height = c(NA, 185), weight = c(75, NA), method = "Takahira"),
    c(NA_real_, NA_real_)
  )
})

## compute_bsa: Check error messages ----

## Test 25: an error is issued if an invalid method is specified ----
test_that("compute_bsa Test 25: an error is issued if an invalid method is specified", {
  expect_error(
    compute_bsa(height = c(170, 185), weight = c(75, 90), method = "unknown-method"),
    class = "assert_character_scalar"
  )
})

# compute_map ----

## compute_map: DBP & SBP ----
# ((2 x DBP) + SBP) / 3

## Test 26: MAP based on diastolic & systolic BP - single values ----
test_that("compute_map Test 26: MAP based on diastolic & systolic BP - single values", {
  expect_equal(round(compute_map(diabp = 51, sysbp = 121), 3L), 74.333)
})

## Test 27: MAP based on diastolic & systolic BP - vectors ----
test_that("compute_map Test 27: MAP based on diastolic & systolic BP - vectors", {
  expect_equal(
    round(compute_map(diabp = c(51, 61), sysbp = c(121, 141)), 3L), c(74.333, 87.667)
  )
})

## Test 28: MAP based on diastolic & systolic BP with missing values ----
test_that("compute_map Test 28: MAP based on diastolic & systolic BP with missing values", {
  expect_equal(
    compute_map(diabp = c(NA, 61), sysbp = c(121, NA)), c(NA_real_, NA_real_)
  )
})

## compute_map: DBP, SBP & HR ----
# DBP + 0.01 x exp(4.14 - 40.74 / PULSE) x (SBP - DBP)

## Test 29: MAP based on DBP & SBP & heart rate - single values ----
test_that("compute_map Test 29: MAP based on DBP & SBP & heart rate - single values", {
  expect_equal(
    round(compute_map(diabp = 51, sysbp = 121, hr = 59), 3L), 73.039
  )
})

## Test 30: MAP based on diastolic, systolic BP & heart rate - vectors ----
test_that("compute_map Test 30: MAP based on diastolic, systolic BP & heart rate - vectors", {
  expect_equal(
    round(compute_map(diabp = c(51, 91), sysbp = c(121, 101), hr = c(59, 62)), 3L),
    c(73.039, 94.255)
  )
})

## Test 31: MAP based on DBP, SBP & heart rate - with missing values ----
test_that("compute_map Test 31: MAP based on DBP, SBP & heart rate - with missing values", {
  expect_equal(
    compute_map(diabp = c(NA, 61, 51), sysbp = c(121, NA, 121), hr = c(59, 62, NA)),
    c(NA_real_, NA_real_, NA_real_)
  )
})

# derive_param_bmi ----

## derive_param_bmi: Error checks ----

## Test 32: BMI parameter NOT added - wrong hgt unit ----
test_that("derive_param_bmi Test 32: BMI parameter NOT added - wrong hgt unit", {
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

## Test 33: BMI parameter NOT added - wrong wgt unit ----
test_that("derive_param_bmi Test 33: BMI parameter NOT added - wrong wgt unit", {
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

## Test 34: BMI parameter NOT added - multiple unit for wgt ----
test_that("derive_param_bmi Test 34: BMI parameter NOT added - multiple unit for wgt", {
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

## Test 35: BMI parameter NOT added - PARAMCD not set ----
test_that("derive_param_bmi Test 35: BMI parameter NOT added - PARAMCD not set", {
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

## Test 36: BMI parameter NOT added ----
test_that("derive_param_bmi Test 36: BMI parameter NOT added", {
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

  expect_dfs_equal(
    derive_param_bmi(input, by_vars = exprs(USUBJID, VISIT), get_unit_expr = VSSTRESU),
    expected_output,
    keys = c("USUBJID", "PARAMCD", "VISIT")
  )
})

## derive_param_bmi: Obs created ----

bmi <- function(hgt, wgt) {
  wgt / (hgt / 100)^2
}

## Test 37: BMI parameter is correctly added ----
test_that("derive_param_bmi Test 37: BMI parameter is correctly added", {
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
## Test 38: Derive BMI where height is measured only once ----
test_that("derive_param_bmi Test 38: Derive BMI where height is measured only once", {
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

# derive_param_bsa ----

## derive_param_bsa: Error checks ----

## Test 39: BSA parameter NOT added - wrong unit for height ----
test_that("derive_param_bsa Test 39: BSA parameter NOT added - wrong unit for height", {
  input <- tibble::tribble(
    ~USUBJID, ~PARAMCD, ~PARAM, ~VISIT, ~VSSTRESU, ~AVAL,
    # Wrong unit for HEIGHT should be cm
    "01-701-1015", "HEIGHT", "Height (cm)", "BASELINE", "m", 170,
    "01-701-1015", "WEIGHT", "Weight (kg)", "BASELINE", "kg", 85,
  )

  expect_error(
    derive_param_bsa(
      input,
      by_vars = exprs(USUBJID, VISIT),
      method = "Mosteller",
      get_unit_expr = VSSTRESU
    ),
    class = "assert_unit"
  )
})

## Test 40: BSA parameter NOT added - wrong unit for weight ----
test_that("derive_param_bsa Test 40: BSA parameter NOT added - wrong unit for weight", {
  input <- tibble::tribble(
    ~USUBJID,      ~PARAMCD, ~PARAM,        ~VISIT,     ~VSSTRESU, ~AVAL,
    "01-701-1015", "HEIGHT", "Height (cm)", "BASELINE",      "cm",   170,
    # Wrong unit for WEIGHT should be kg
    "01-701-1015", "WEIGHT", "Weight (kg)", "BASELINE",       "g",    85,
  )

  expect_error(
    derive_param_bsa(
      input,
      by_vars = exprs(USUBJID, VISIT),
      method = "Mosteller",
      get_unit_expr = VSSTRESU
    ),
    class = "assert_unit"
  )
})

## Test 41: BSA parameter NOT added - multiple unit for weight ----
test_that("derive_param_bsa Test 41: BSA parameter NOT added - multiple unit for weight", {
  input <- tibble::tribble(
    ~USUBJID,      ~PARAMCD, ~PARAM,        ~VISIT,     ~VSSTRESU, ~AVAL,
    "01-701-1015", "HEIGHT", "Height (cm)", "BASELINE",      "cm",   170,
    # Multiple units for WEIGHT
    "01-701-1015", "WEIGHT", "Weight (kg)", "BASELINE",      "kg",    85,
    "01-701-1016", "WEIGHT", "Weight (kg)", "BASELINE",       "g",  8500,
  )

  expect_error(
    derive_param_bsa(
      input,
      by_vars = exprs(USUBJID, VISIT),
      method = "Mosteller",
      get_unit_expr = VSSTRESU
    ),
    class = "assert_unit"
  )
})

## Test 42: BSA parameter NOT added - PARAMCD not set ----
test_that("derive_param_bsa Test 42: BSA parameter NOT added - PARAMCD not set", {
  input <- tibble::tribble(
    ~USUBJID,      ~PARAMCD, ~PARAM,        ~VISIT,     ~VSSTRESU, ~AVAL,
    "01-701-1015", "HEIGHT", "Height (cm)", "BASELINE",      "cm",   170,
    "01-701-1015", "WEIGHT", "Weight (kg)", "BASELINE",      "kg",    85,
  )

  expect_error(
    derive_param_bsa(
      input,
      by_vars = exprs(USUBJID, VISIT),
      method = "Mosteller",
      set_values_to = exprs(PARAM = "Body Surface Area"),
      get_unit_expr = VSSTRESU
    ),
    class = "assert_varval_list"
  )
})

## derive_param_bsa: No obs added ----

## Test 43: BSA parameter NOT added ----
test_that("derive_param_bsa Test 43: BSA parameter NOT added", {
  expected_output <- tibble::tribble(
    ~USUBJID, ~PARAMCD, ~PARAM, ~VISIT, ~VSSTRESU, ~AVAL,
    "01-701-1015", "HEIGHT", "Height (cm)", "BASELINE", "cm", 170,
    # WEIGHT set to NA - so BSA not calculated
    "01-701-1015", "WEIGHT", "Weight (kg)", "BASELINE", "kg", NA,
    "01-701-1015", "WEIGHT", "Weight (kg)", "MONTH 1", "kg", 78,
    # HEIGHT set to NA - so BSA not calculated
    "01-701-1028", "HEIGHT", "Height (cm)", "BASELINE", "cm", NA,
    "01-701-1028", "WEIGHT", "Weight (kg)", "BASELINE", "kg", 90,
    "01-701-1028", "HEIGHT", "Height (cm)", "MONTH 1", "cm", 88,
  )

  input <- expected_output

  expect_dfs_equal(
    derive_param_bsa(
      input,
      by_vars = exprs(USUBJID, VISIT),
      method = "Mosteller",
      get_unit_expr = VSSTRESU
    ),
    expected_output,
    keys = c("USUBJID", "PARAMCD", "VISIT")
  )
})

## derive_param_bsa: Obs created ----

mosteller <- function(hgt, wgt) {
  sqrt(hgt * wgt / 3600)
}

## Test 44: BSA parameter (Mosteller Method) is correctly added ----
test_that("derive_param_bsa Test 44: BSA parameter (Mosteller Method) is correctly added", {
  expected_output <- tibble::tribble(
    ~USUBJID, ~PARAMCD, ~PARAM, ~VISIT, ~VSSTRESU, ~AVAL,
    "01-701-1015", "HEIGHT", "Height (cm)", "BASELINE", "cm", 170,
    "01-701-1015", "WEIGHT", "Weight (kg)", "BASELINE", "kg", 75,
    # New row added for BMI for SUBJID="01-701-1015" and VISIT="BASELINE"
    # WEIGHT = 75 and HEIGHT = 170
    "01-701-1015", "BSA", NA, "BASELINE", NA, mosteller(170, 75),
    "01-701-1015", "WEIGHT", "Weight (kg)", "MONTH 1", "kg", 78,
    "01-701-1028", "HEIGHT", "Height (cm)", "BASELINE", "cm", 185,
    "01-701-1028", "WEIGHT", "Weight (kg)", "BASELINE", "kg", 90,
    # New row added for BMI for SUBJID="01-701-1028" and VISIT='BASELINE'
    # WEIGHT = 90 and HEIGHT = 185
    "01-701-1028", "BSA", NA, "BASELINE", NA, mosteller(185, 90),
    "01-701-1028", "WEIGHT", "Weight (kg)", "MONTH 1", "kg", 88,
  )

  input <- expected_output %>% filter(PARAMCD != "BSA")

  expect_dfs_equal(
    derive_param_bsa(input,
      by_vars = exprs(USUBJID, VISIT),
      method = "Mosteller",
      get_unit_expr = VSSTRESU
    ),
    expected_output,
    keys = c("USUBJID", "PARAMCD", "VISIT")
  )
})

dubois <- function(hgt, wgt) {
  0.20247 * (hgt / 100)^0.725 * wgt^0.425
}

## Test 45: BSA parameter (DuBois-DuBois method) is correctly added ----
test_that("derive_param_bsa Test 45: BSA parameter (DuBois-DuBois method) is correctly added", {
  expected_output <- tibble::tribble(
    ~USUBJID, ~PARAMCD, ~PARAM, ~VISIT, ~VSSTRESU, ~AVAL,
    "01-701-1015", "HEIGHT", "Height (cm)", "BASELINE", "cm", 170,
    "01-701-1015", "WEIGHT", "Weight (kg)", "BASELINE", "kg", 75,
    # New row added for BMI for SUBJID="01-701-1015" and VISIT="BASELINE"
    # WEIGHT = 75 and HEIGHT = 170
    "01-701-1015", "BSA", NA, "BASELINE", NA, dubois(170, 75),
    "01-701-1015", "WEIGHT", "Weight (kg)", "MONTH 1", "kg", 78,
    "01-701-1028", "HEIGHT", "Height (cm)", "BASELINE", "cm", 185,
    "01-701-1028", "WEIGHT", "Weight (kg)", "BASELINE", "kg", 90,
    # New row added for BMI for SUBJID="01-701-1028" and VISIT='BASELINE'
    # WEIGHT = 90 and HEIGHT = 185
    "01-701-1028", "BSA", NA, "BASELINE", NA, dubois(185, 90),
    "01-701-1028", "WEIGHT", "Weight (kg)", "MONTH 1", "kg", 88,
  )

  input <- expected_output %>% filter(PARAMCD != "BSA")

  expect_dfs_equal(
    derive_param_bsa(
      input,
      by_vars = exprs(USUBJID, VISIT),
      method = "DuBois-DuBois",
      get_unit_expr = VSSTRESU
    ),
    expected_output,
    keys = c("USUBJID", "PARAMCD", "VISIT")
  )
})

haycock <- function(hgt, wgt) {
  0.024265 * hgt^0.3964 * wgt^0.5378
}

## Test 46: BSA parameter (Haycock method) is correctly added ----
test_that("derive_param_bsa Test 46: BSA parameter (Haycock method) is correctly added", {
  expected_output <- tibble::tribble(
    ~USUBJID, ~PARAMCD, ~PARAM, ~VISIT, ~VSSTRESU, ~AVAL,
    "01-701-1015", "HEIGHT", "Height (cm)", "BASELINE", "cm", 170,
    "01-701-1015", "WEIGHT", "Weight (kg)", "BASELINE", "kg", 75,
    # New row added for BMI for SUBJID="01-701-1015" and VISIT="BASELINE"
    # WEIGHT = 75 and HEIGHT = 170
    "01-701-1015", "BSA", NA, "BASELINE", NA, haycock(170, 75),
    "01-701-1015", "WEIGHT", "Weight (kg)", "MONTH 1", "kg", 78,
    "01-701-1028", "HEIGHT", "Height (cm)", "BASELINE", "cm", 185,
    "01-701-1028", "WEIGHT", "Weight (kg)", "BASELINE", "kg", 90,
    # New row added for BMI for SUBJID="01-701-1028" and VISIT='BASELINE'
    # WEIGHT = 90 and HEIGHT = 185
    "01-701-1028", "BSA", NA, "BASELINE", NA, haycock(185, 90),
    "01-701-1028", "WEIGHT", "Weight (kg)", "MONTH 1", "kg", 88,
  )

  input <- expected_output %>% filter(PARAMCD != "BSA")

  expect_dfs_equal(
    derive_param_bsa(input,
      by_vars = exprs(USUBJID, VISIT),
      method = "Haycock",
      get_unit_expr = VSSTRESU
    ),
    expected_output,
    keys = c("USUBJID", "PARAMCD", "VISIT")
  )
})

gehan <- function(hgt, wgt) {
  0.0235 * hgt^0.42246 * wgt^0.51456
}

## Test 47: BSA parameter (Gehan-George method) is correctly added ----
test_that("derive_param_bsa Test 47: BSA parameter (Gehan-George method) is correctly added", {
  expected_output <- tibble::tribble(
    ~USUBJID, ~PARAMCD, ~PARAM, ~VISIT, ~VSSTRESU, ~AVAL,
    "01-701-1015", "HEIGHT", "Height (cm)", "BASELINE", "cm", 170,
    "01-701-1015", "WEIGHT", "Weight (kg)", "BASELINE", "kg", 75,
    # New row added for BMI for SUBJID="01-701-1015" and VISIT="BASELINE"
    # WEIGHT = 75 and HEIGHT = 170
    "01-701-1015", "BSA", NA, "BASELINE", NA, gehan(170, 75),
    "01-701-1015", "WEIGHT", "Weight (kg)", "MONTH 1", "kg", 78,
    "01-701-1028", "HEIGHT", "Height (cm)", "BASELINE", "cm", 185,
    "01-701-1028", "WEIGHT", "Weight (kg)", "BASELINE", "kg", 90,
    # New row added for BMI for SUBJID="01-701-1028" and VISIT='BASELINE'
    # WEIGHT = 90 and HEIGHT = 185
    "01-701-1028", "BSA", NA, "BASELINE", NA, gehan(185, 90),
    "01-701-1028", "WEIGHT", "Weight (kg)", "MONTH 1", "kg", 88,
  )

  input <- expected_output %>% filter(PARAMCD != "BSA")

  expect_dfs_equal(
    derive_param_bsa(
      input,
      by_vars = exprs(USUBJID, VISIT),
      method = "Gehan-George",
      get_unit_expr = VSSTRESU
    ),
    expected_output,
    keys = c("USUBJID", "PARAMCD", "VISIT")
  )
})

boyd <- function(hgt, wgt) {
  0.0003207 * (hgt^0.3) * (1000 * wgt)^(0.7285 - (0.0188 * log10(1000 * wgt))) # nolint
}

## Test 48: BSA parameter (Boyd method) is correctly added ----
test_that("derive_param_bsa Test 48: BSA parameter (Boyd method) is correctly added", {
  expected_output <- tibble::tribble(
    ~USUBJID, ~PARAMCD, ~PARAM, ~VISIT, ~VSSTRESU, ~AVAL,
    "01-701-1015", "HEIGHT", "Height (cm)", "BASELINE", "cm", 170,
    "01-701-1015", "WEIGHT", "Weight (kg)", "BASELINE", "kg", 75,
    # New row added for BMI for SUBJID="01-701-1015" and VISIT="BASELINE"
    # WEIGHT = 75 and HEIGHT = 170
    "01-701-1015", "BSA", NA, "BASELINE", NA, boyd(170, 75),
    "01-701-1015", "WEIGHT", "Weight (kg)", "MONTH 1", "kg", 78,
    "01-701-1028", "HEIGHT", "Height (cm)", "BASELINE", "cm", 185,
    "01-701-1028", "WEIGHT", "Weight (kg)", "BASELINE", "kg", 90,
    # New row added for BMI for SUBJID="01-701-1028" and VISIT='BASELINE'
    # WEIGHT = 90 and HEIGHT = 185
    "01-701-1028", "BSA", NA, "BASELINE", NA, boyd(185, 90),
    "01-701-1028", "WEIGHT", "Weight (kg)", "MONTH 1", "kg", 88,
  )

  input <- expected_output %>% filter(PARAMCD != "BSA")

  expect_dfs_equal(
    derive_param_bsa(input,
      by_vars = exprs(USUBJID, VISIT),
      method = "Boyd",
      get_unit_expr = VSSTRESU
    ),
    expected_output,
    keys = c("USUBJID", "PARAMCD", "VISIT")
  )
})

fujimoto <- function(hgt, wgt) {
  0.008883 * hgt^0.663 * wgt^0.444
}

## Test 49: BSA parameter (Fujimoto method) is correctly added ----
test_that("derive_param_bsa Test 49: BSA parameter (Fujimoto method) is correctly added", {
  expected_output <- tibble::tribble(
    ~USUBJID, ~PARAMCD, ~PARAM, ~VISIT, ~VSSTRESU, ~AVAL,
    "01-701-1015", "HEIGHT", "Height (cm)", "BASELINE", "cm", 170,
    "01-701-1015", "WEIGHT", "Weight (kg)", "BASELINE", "kg", 75,
    # New row added for BMI for SUBJID="01-701-1015" and VISIT="BASELINE"
    # WEIGHT = 75 and HEIGHT = 170
    "01-701-1015", "BSA", NA, "BASELINE", NA, fujimoto(170, 75),
    "01-701-1015", "WEIGHT", "Weight (kg)", "MONTH 1", "kg", 78,
    "01-701-1028", "HEIGHT", "Height (cm)", "BASELINE", "cm", 185,
    "01-701-1028", "WEIGHT", "Weight (kg)", "BASELINE", "kg", 90,
    # New row added for BMI for SUBJID="01-701-1028" and VISIT='BASELINE'
    # WEIGHT = 90 and HEIGHT = 185
    "01-701-1028", "BSA", NA, "BASELINE", NA, fujimoto(185, 90),
    "01-701-1028", "WEIGHT", "Weight (kg)", "MONTH 1", "kg", 88,
  )

  input <- expected_output %>% filter(PARAMCD != "BSA")

  expect_dfs_equal(
    derive_param_bsa(
      input,
      by_vars = exprs(USUBJID, VISIT),
      method = "Fujimoto",
      get_unit_expr = VSSTRESU
    ),
    expected_output,
    keys = c("USUBJID", "PARAMCD", "VISIT")
  )
})

takahira <- function(hgt, wgt) {
  0.007241 * hgt^0.725 * wgt^0.425
}

## Test 50: BSA parameter (Takahira method) is correctly added ----
test_that("derive_param_bsa Test 50: BSA parameter (Takahira method) is correctly added", {
  expected_output <- tibble::tribble(
    ~USUBJID, ~PARAMCD, ~PARAM, ~VISIT, ~VSSTRESU, ~AVAL,
    "01-701-1015", "HEIGHT", "Height (cm)", "BASELINE", "cm", 170,
    "01-701-1015", "WEIGHT", "Weight (kg)", "BASELINE", "kg", 75,
    # New row added for BMI for SUBJID="01-701-1015" and VISIT="BASELINE"
    # WEIGHT = 75 and HEIGHT = 170
    "01-701-1015", "BSA", NA, "BASELINE", NA, takahira(170, 75),
    "01-701-1015", "WEIGHT", "Weight (kg)", "MONTH 1", "kg", 78,
    "01-701-1028", "HEIGHT", "Height (cm)", "BASELINE", "cm", 185,
    "01-701-1028", "WEIGHT", "Weight (kg)", "BASELINE", "kg", 90,
    # New row added for BMI for SUBJID="01-701-1028" and VISIT='BASELINE'
    # WEIGHT = 90 and HEIGHT = 185
    "01-701-1028", "BSA", NA, "BASELINE", NA, takahira(185, 90),
    "01-701-1028", "WEIGHT", "Weight (kg)", "MONTH 1", "kg", 88,
  )

  input <- expected_output %>% filter(PARAMCD != "BSA")

  expect_dfs_equal(
    derive_param_bsa(
      input,
      by_vars = exprs(USUBJID, VISIT),
      method = "Takahira",
      get_unit_expr = VSSTRESU
    ),
    expected_output,
    keys = c("USUBJID", "PARAMCD", "VISIT")
  )
})

## Test 51: Derive BSA where height is measured only once ----
test_that("derive_param_bsa Test 51: Derive BSA where height is measured only once", {
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
      AVAL = compute_bsa(
        height = AVAL.HEIGHT, weight = AVAL.WEIGHT,
        method = "Mosteller"
      ),
      PARAMCD = "BSA",
      PARAM = "Body Surface Area (m^2)",
      AVALU = "m^2"
    ),
    constant_parameters = c("HEIGHT"),
    constant_by_vars = exprs(USUBJID)
  )

  expect_dfs_equal(
    expected_output,
    derive_param_bsa(
      input,
      by_vars = exprs(USUBJID, VISIT),
      method = "Mosteller",
      set_values_to = exprs(
        PARAMCD = "BSA",
        PARAM = "Body Surface Area (m^2)",
        AVALU = "m^2"
      ),
      get_unit_expr = extract_unit(PARAM),
      constant_by_vars = exprs(USUBJID)
    ),
    keys = c("USUBJID", "PARAMCD", "VISIT")
  )
})


# derive_param_map ----

## derive_param_map: Error checks ----

## Test 52: MAP parameter NOT added - wrong DIABP unit ----
test_that("derive_param_map Test 52: MAP parameter NOT added - wrong DIABP unit", {
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

## Test 53: MAP parameter NOT added - wrong SYSBP unit ----
test_that("derive_param_map Test 53: MAP parameter NOT added - wrong SYSBP unit", {
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

## Test 54: MAP parameter NOT added - wrong PULSE unit ----
test_that("derive_param_map Test 54: MAP parameter NOT added - wrong PULSE unit", {
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

## Test 55: MAP parameter NOT added - PARAMCD not set ----
test_that("derive_param_map Test 55: MAP parameter NOT added - PARAMCD not set", {
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

## Test 56: MAP parameter NOT added ----
test_that("derive_param_map Test 56: MAP parameter NOT added", {
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

## derive_param_map: Obs created ----

maphr <- function(sbp, dbp, hr) {
  dbp + 0.01 * exp(4.14 - 40.74 / hr) * (sbp - dbp)
}

## Test 57: MAP parameter (DBP/SBP/PULSE) is correctly added ----
test_that("derive_param_map Test 57: MAP parameter (DBP/SBP/PULSE) is correctly added", {
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

## Test 58: MAP parameter (DBP/SBP) is correctly added ----
test_that("derive_param_map Test 58: MAP parameter (DBP/SBP) is correctly added", {
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
