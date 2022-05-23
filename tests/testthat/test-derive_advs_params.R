# compute_bmi: (Test 01.xx) ----

test_that("compute_bmi Test 01.01: BMI calculation - single height and weight values", {
  # Expected values are taken from the Center of Disease Control and Prevention's
  # (CDC) 'Adult BMI Calculator' at
  # https://cdc.gov/healthyweight/assessing/bmi/adult_bmi/metric_bmi_calculator/bmi_calculator.html
  expect_equal(round(compute_bmi(height = 180, weight = 75), 3L), 23.148)
})

test_that("compute_bmi Test 01.02: BMI calculation - height and weight vectors", {
  expect_equal(
    round(compute_bmi(height = c(180, 200), weight = c(75, 100)), 3L),
    c(23.148, 25)
  )
})

test_that("compute_bmi Test 01.03: BMI calculation - height and weight vectors - missing values", {
  expect_equal(
    compute_bmi(height = c(NA, 200, 0), weight = c(75, NA, 75)),
    c(NA_real_, NA_real_, NA_real_)
  )
})

# compute_bsa ----

## compute_bsa: Mosteller method (Test 01.xx) ----
# sqrt (Height x Weight / 3600)

test_that("compute_bsa Test 01.01: Mosteller method - single height and weight values", {
  expect_equal(
    round(compute_bsa(height = 170, weight = 75, method = "Mosteller"), 3L),
    1.882
  )
})

test_that("compute_bsa Test 01.02: Mosteller method - height and weight vectors", {
  expect_equal(
    round(compute_bsa(height = c(170, 185), weight = c(75, 90), method = "Mosteller"), 3L),
    c(1.882, 2.151)
  )
})

test_that("compute_bsa Test 01.03: Mosteller method - height and weight vectors - missing values", {
  expect_equal(
    compute_bsa(height = c(NA, 185), weight = c(75, NA), method = "Mosteller"),
    c(NA_real_, NA_real_)
  )
})

## compute_bsa: DuBois-DuBois method (Test 02.xx) ----
#  FORMULA : 0.20247 x (HGT/100)^0.725 x WGT^0.425

test_that("compute_bsa Test 02.01: DuBois-DuBois method - single height and weight values", {
  expect_equal(
    round(compute_bsa(height = 170, weight = 75, method = "DuBois-DuBois"), 3L),
    1.864
  )
})

test_that("compute_bsa Test 02.02: DuBois-DuBois method - height and weight vectors", {
  expect_equal(
    round(compute_bsa(height = c(170, 185), weight = c(75, 90), method = "DuBois-DuBois"), 3L),
    c(1.864, 2.141)
  )
})

test_that("compute_bsa Test 02.03: DuBois-DuBois method - hgt and wgt vectors - missing values", {
  expect_equal(
    compute_bsa(height = c(NA, 185), weight = c(75, NA), method = "DuBois-DuBois"),
    c(NA_real_, NA_real_)
  )
})

## compute_bsa: Haycock method (Test 03.xx) ----
# 0.024265 x HGT^0.3964 x WGT^0.5378

test_that("compute_bsa Test 03.01: Haycock method - single height and weight values", {
  expect_equal(
    round(compute_bsa(height = 170, weight = 75, method = "Haycock"), 3L),
    1.895
  )
})

test_that("compute_bsa Test 03.02: Haycock method - height and weight vectors", {
  expect_equal(
    round(compute_bsa(height = c(170, 185), weight = c(75, 90), method = "Haycock"), 3L),
    c(1.895, 2.161)
  )
})

test_that("compute_bsa Test 03.03: Haycock method - height and weight vectors - missing values", {
  expect_equal(
    compute_bsa(height = c(NA, 185), weight = c(75, NA), method = "Haycock"),
    c(NA_real_, NA_real_)
  )
})

## compute_bsa: Gehan-George method (Test 04.xx) ----
# 0.0235 x HGT^0.42246 x WGT^0.51456

test_that("compute_bsa Test 04.01: Gehan-George method - single height and weight values", {
  expect_equal(
    round(compute_bsa(height = 170, weight = 75, method = "Gehan-George"), 3L),
    1.897
  )
})

test_that("compute_bsa Test 04.02: Gehan-George method - height and weight vectors", {
  expect_equal(
    round(compute_bsa(height = c(170, 185), weight = c(75, 90), method = "Gehan-George"), 3L),
    c(1.897, 2.16)
  )
})

test_that(paste(
  "compute_bsa Test 04.03: Gehan-George method - height and",
  "weight vectors - missing values"
), {
  expect_equal(
    compute_bsa(height = c(NA, 185), weight = c(75, NA), method = "Gehan-George"),
    c(NA_real_, NA_real_)
  )
})

## compute_bsa: Boyd method  (Test 05.xx) ----
# 0.0003207 x (HGT^0.3) x (1000 x WGT)^(0.7285 - (0.0188 x log10(1000 x WGT)))

test_that("compute_bsa Test 05.01: Boyd method - single height and weight values", {
  expect_equal(
    round(compute_bsa(height = 170, weight = 75, method = "Boyd"), 3L),
    1.905
  )
})

test_that("compute_bsa Test 05.02: Boyd method - height and weight vectors", {
  expect_equal(
    round(compute_bsa(height = c(170, 185), weight = c(75, 90), method = "Boyd"), 3L),
    c(1.905, 2.158)
  )
})

test_that("compute_bsa Test 05.03: Boyd method - height and weight vectors - missing values", {
  expect_equal(
    compute_bsa(height = c(NA, 185), weight = c(75, NA), method = "Boyd"),
    c(NA_real_, NA_real_)
  )
})

## compute_bsa: Fujimoto method (Test 06.xx) ----
# 0.008883 x HGT^0.663 x WGT^0.444

test_that("compute_bsa Test 06.01: Fujimoto method - single height and weight values", {
  expect_equal(
    round(compute_bsa(height = 170, weight = 75, method = "Fujimoto"), 3L),
    1.819
  )
})

test_that("compute_bsa Test 06.02: Fujimoto method - height and weight vectors", {
  expect_equal(
    round(compute_bsa(height = c(170, 185), weight = c(75, 90), method = "Fujimoto"), 3L),
    c(1.819, 2.086)
  )
})

test_that("compute_bsa Test 06.03: Fujimoto method - height and weight vectors - missing values", {
  expect_equal(
    compute_bsa(height = c(NA, 185), weight = c(75, NA), method = "Fujimoto"),
    c(NA_real_, NA_real_)
  )
})

## compute_bsa: Takahira method (Test 07.xx) ----
# 0.007241 x HGT^0.725 x WGT^0.425

test_that("compute_bsa Test 07.01: Takahira method - single height and weight values", {
  expect_equal(
    round(compute_bsa(height = 170, weight = 75, method = "Takahira"), 3L),
    1.878
  )
})

test_that("compute_bsa Test 07.02: Takahira method - height and weight vectors", {
  expect_equal(
    round(compute_bsa(height = c(170, 185), weight = c(75, 90), method = "Takahira"), 3L),
    c(1.878, 2.158)
  )
})

test_that("compute_bsa Test 07.03: Takahira method - height and weight vectors - missing values", {
  expect_equal(
    compute_bsa(height = c(NA, 185), weight = c(75, NA), method = "Takahira"),
    c(NA_real_, NA_real_)
  )
})

## compute_bsa: Check error messages (Test 08.xx) ----

test_that("compute_bsa Test 08.01: an error is issued if an invalid method is specified", {
  expect_error(
    compute_bsa(height = c(170, 185), weight = c(75, 90), method = "unknown-method"),
    paste(
      "`method` must be one of 'Mosteller', 'DuBois-DuBois', 'Haycock', 'Gehan-George',",
      "'Boyd', 'Fujimoto' or 'Takahira' but is 'unknown-method'"
    )
  )
})

# compute_map ----

## compute_map: DBP & SBP (Test 01.xx) ----
# ((2 x DBP) + SBP) / 3

test_that(paste(
  "compute_map Test 01.01: Mean Arterial Pressure based on diastolic",
  "& systolic BP - single values"
), {
  expect_equal(round(compute_map(diabp = 51, sysbp = 121), 3L), 74.333)
})

test_that(paste(
  "compute_map Test 01.02: Mean Arterial Pressure based on diastolic",
  "& systolic BP - vectors"
), {
  expect_equal(
    round(compute_map(diabp = c(51, 61), sysbp = c(121, 141)), 3L), c(74.333, 87.667)
  )
})

test_that(paste(
  "compute_map Test 01.03: Mean Arterial Pressure based on diastolic",
  "& systolic BP - vectors with missing values"
), {
  expect_equal(
    compute_map(diabp = c(NA, 61), sysbp = c(121, NA)), c(NA_real_, NA_real_)
  )
})

## compute_map: DBP, SBP & HR (Test 02.xx) ----
# DBP + 0.01 x exp(4.14 - 40.74 / PULSE) x (SBP - DBP)

test_that(paste(
  "compute_map Test 02.01: Mean Arterial Pressure based on diastolic,",
  "systolic BP & heart rate - single values"
), {
  expect_equal(
    round(compute_map(diabp = 51, sysbp = 121, hr = 59), 3L), 73.039
  )
})

test_that(paste(
  "compute_map Test 02.02: Mean Arterial Pressure based on diastolic,",
  "systolic BP & heart rate - vectors"
), {
  expect_equal(
    round(compute_map(diabp = c(51, 91), sysbp = c(121, 101), hr = c(59, 62)), 3L),
    c(73.039, 94.255)
  )
})

test_that(paste(
  "compute_map Test 02.03: Mean Arterial Pressure based on diastolic,",
  "systolic blood BP & heart rate - vectors with missing values"
), {
  expect_equal(
    compute_map(diabp = c(NA, 61, 51), sysbp = c(121, NA, 121), hr = c(59, 62, NA)),
    c(NA_real_, NA_real_, NA_real_)
  )
})

# derive_param_bmi  ----

## derive_param_bmi: Error checks (Test 01.xx) ----

test_that(paste(
  "derive_param_bmi Test 01.01: BMI parameter NOT added to input dataset",
  "- wrong unit for hgt"
), {
  input <- tibble::tribble(
    ~USUBJID, ~PARAMCD, ~PARAM, ~VISIT, ~VSSTRESU, ~AVAL,
    # Wrong unit for HEIGHT should be cm
    "01-701-1015", "HEIGHT", "Height (cm)", "BASELINE", "m", 170,
    "01-701-1015", "WEIGHT", "Weight (kg)", "BASELINE", "kg", 85,
  )

  expect_error(
    derive_param_bmi(input, by_vars = vars(USUBJID, VISIT), get_unit_expr = VSSTRESU),
    paste(
      "It is expected that 'HEIGHT' is measured in 'cm'.\nIn the",
      "input dataset it is measured in 'm'."
    )
  )
})

test_that(paste(
  "derive_param_bmi Test 01.02: BMI parameter NOT added to input dataset",
  "- wrong unit for wgt"
), {
  input <- tibble::tribble(
    ~USUBJID,      ~PARAMCD, ~PARAM,        ~VISIT,     ~VSSTRESU, ~AVAL,
    "01-701-1015", "HEIGHT", "Height (cm)", "BASELINE",      "cm",   170,
    # Wrong unit for WEIGHT should be kg
    "01-701-1015", "WEIGHT", "Weight (kg)", "BASELINE",       "g",    85,
  )

  expect_error(
    derive_param_bmi(input, by_vars = vars(USUBJID, VISIT), get_unit_expr = VSSTRESU),
    paste(
      "It is expected that 'WEIGHT' is measured in 'kg'.\nIn the",
      "input dataset it is measured in 'g'."
    )
  )
})

test_that(paste(
  "derive_param_bmi Test 01.03: BMI parameter NOT added to input dataset",
  "- multiple unit for wgt"
), {
  input <- tibble::tribble(
    ~USUBJID,      ~PARAMCD, ~PARAM,        ~VISIT,     ~VSSTRESU, ~AVAL,
    "01-701-1015", "HEIGHT", "Height (cm)", "BASELINE",      "cm",   170,
    # Multiple units for WEIGHT
    "01-701-1015", "WEIGHT", "Weight (kg)", "BASELINE",      "kg",    85,
    "01-701-1016", "WEIGHT", "Weight (kg)", "BASELINE",       "g",  8500,
  )

  expect_error(
    derive_param_bmi(input, by_vars = vars(USUBJID, VISIT), get_unit_expr = VSSTRESU),
    paste0(
      "Multiple units 'kg' and 'g' found for 'WEIGHT'.",
      "\nPlease review and update the units."
    )
  )
})

test_that(paste(
  "derive_param_bmi Test 01.04: BMI parameter NOT added to input dataset",
  "- PARAMCD not set"
), {
  input <- tibble::tribble(
    ~USUBJID,      ~PARAMCD, ~PARAM,        ~VISIT,     ~VSSTRESU, ~AVAL,
    "01-701-1015", "HEIGHT", "Height (cm)", "BASELINE",      "cm",   170,
    "01-701-1015", "WEIGHT", "Weight (kg)", "BASELINE",      "kg",    85,
  )

  expect_error(
    derive_param_bmi(
      input,
      by_vars = vars(USUBJID, VISIT),
      set_values_to = vars(PARAM = "Body Mass Index"),
      get_unit_expr = VSSTRESU
    ),
    "The following required elements are missing in `set_values_to`: 'PARAMCD'"
  )
})

## derive_param_bmi: No obs added (Test 02.xx)  ----

test_that("derive_param_bmi Test 02.01: BMI parameter NOT added to input dataset", {
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
    derive_param_bmi(input, by_vars = vars(USUBJID, VISIT), get_unit_expr = VSSTRESU),
    expected_output,
    keys = c("USUBJID", "PARAMCD", "VISIT")
  )
})

## derive_param_bmi: Obs created (Test 03.xx)  ----

bmi <- function(hgt, wgt) {
  wgt / (hgt / 100)^2
}

test_that("derive_param_bmi Test 03.01: BMI parameter is correctly added to input dataset", {
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
    derive_param_bmi(input, by_vars = vars(USUBJID, VISIT), get_unit_expr = VSSTRESU),
    expected_output,
    keys = c("USUBJID", "PARAMCD", "VISIT")
  )
})

# derive_param_bsa ----

## derive_param_bsa: Error checks (Test 01.xx) ----

test_that(paste(
  "derive_param_bsa Test 01.01: BSA parameter NOT added to input dataset",
  "- wrong unit for height"
), {
  input <- tibble::tribble(
    ~USUBJID, ~PARAMCD, ~PARAM, ~VISIT, ~VSSTRESU, ~AVAL,
    # Wrong unit for HEIGHT should be cm
    "01-701-1015", "HEIGHT", "Height (cm)", "BASELINE", "m", 170,
    "01-701-1015", "WEIGHT", "Weight (kg)", "BASELINE", "kg", 85,
  )

  expect_error(
    derive_param_bsa(
      input,
      by_vars = vars(USUBJID, VISIT),
      method = "Mosteller",
      get_unit_expr = VSSTRESU
    ),
    paste(
      "It is expected that 'HEIGHT' is measured in 'cm'.\nIn the",
      "input dataset it is measured in 'm'."
    )
  )
})

test_that(paste(
  "derive_param_bsa Test 01.02: BSA parameter NOT added to input dataset",
  "- wrong unit for weight"
), {
  input <- tibble::tribble(
    ~USUBJID,      ~PARAMCD, ~PARAM,        ~VISIT,     ~VSSTRESU, ~AVAL,
    "01-701-1015", "HEIGHT", "Height (cm)", "BASELINE",      "cm",   170,
    # Wrong unit for WEIGHT should be kg
    "01-701-1015", "WEIGHT", "Weight (kg)", "BASELINE",       "g",    85,
  )

  expect_error(
    derive_param_bsa(
      input,
      by_vars = vars(USUBJID, VISIT),
      method = "Mosteller",
      get_unit_expr = VSSTRESU
    ),
    paste(
      "It is expected that 'WEIGHT' is measured in 'kg'.\nIn the",
      "input dataset it is measured in 'g'."
    )
  )
})

test_that(paste(
  "derive_param_bsa Test 01.03: BSA parameter NOT added to input dataset",
  "- multiple unit for weight"
), {
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
      by_vars = vars(USUBJID, VISIT),
      method = "Mosteller",
      get_unit_expr = VSSTRESU
    ),
    paste0(
      "Multiple units 'kg' and 'g' found for 'WEIGHT'.",
      "\nPlease review and update the units."
    )
  )
})

test_that(paste(
  "derive_param_bsa Test 01.04: BSA parameter NOT added to input dataset",
  "- PARAMCD not set"
), {
  input <- tibble::tribble(
    ~USUBJID,      ~PARAMCD, ~PARAM,        ~VISIT,     ~VSSTRESU, ~AVAL,
    "01-701-1015", "HEIGHT", "Height (cm)", "BASELINE",      "cm",   170,
    "01-701-1015", "WEIGHT", "Weight (kg)", "BASELINE",      "kg",    85,
  )

  expect_error(
    derive_param_bsa(
      input,
      by_vars = vars(USUBJID, VISIT),
      method = "Mosteller",
      set_values_to = vars(PARAM = "Body Surface Area"),
      get_unit_expr = VSSTRESU
    ),
    paste(
      "The following required elements are missing in",
      "`set_values_to`: 'PARAMCD'"
    )
  )
})

## derive_param_bsa: No obs added (Test 02.xx) ----

test_that("derive_param_bsa Test 02.01: BSA parameter NOT added to input dataset", {
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
      by_vars = vars(USUBJID, VISIT),
      method = "Mosteller",
      get_unit_expr = VSSTRESU
    ),
    expected_output,
    keys = c("USUBJID", "PARAMCD", "VISIT")
  )
})

## derive_param_bsa: Obs created (Test 03.xx) ----

mosteller <- function(hgt, wgt) {
  sqrt(hgt * wgt / 3600)
}

test_that(paste(
  "derive_param_bsa Test 03.01: BSA parameter (Mosteller method) is",
  "correctly added to input dataset"
), {
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
      by_vars = vars(USUBJID, VISIT),
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

test_that(paste(
  "derive_param_bsa Test 03.02: BSA parameter (DuBois-DuBois method)",
  "is correctly added to input dataset"
), {
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
      by_vars = vars(USUBJID, VISIT),
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

test_that(paste(
  "derive_param_bsa Test 03.03: BSA parameter (Haycock method) is",
  "correctly added to input dataset"
), {
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
      by_vars = vars(USUBJID, VISIT),
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

test_that(paste(
  "derive_param_bsa Test 03.04: BSA parameter (Gehan-George method)",
  "is correctly added to input dataset"
), {
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
      by_vars = vars(USUBJID, VISIT),
      method = "Gehan-George",
      get_unit_expr = VSSTRESU
    ),
    expected_output,
    keys = c("USUBJID", "PARAMCD", "VISIT")
  )
})

boyd <- function(hgt, wgt) {
  0.0003207 * (hgt^0.3) * (1000 * wgt)^(0.7285 - (0.0188 * log10(1000 * wgt)))
}

test_that(paste(
  "derive_param_bsa Test 03.05: BSA parameter (Boyd method) is ",
  "correctly added to input dataset"
), {
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
      by_vars = vars(USUBJID, VISIT),
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

test_that(paste(
  "derive_param_bsa Test 03.06: BSA parameter (Fujimoto method) is",
  "correctly added to input dataset"
), {
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
      by_vars = vars(USUBJID, VISIT),
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
test_that(paste(
  "derive_param_bsa Test 03.07: BSA parameter (Takahira method) is",
  "correctly added to input dataset"
), {
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
      by_vars = vars(USUBJID, VISIT),
      method = "Takahira",
      get_unit_expr = VSSTRESU
    ),
    expected_output,
    keys = c("USUBJID", "PARAMCD", "VISIT")
  )
})

# derive_param_map ----

## derive_param_map: Error checks (Test 01.xx) ----

test_that(paste(
  "derive_param_map Test 01.01: MAP parameter NOT added to input dataset",
  "- wrong unit for DIABP"
), {
  input <- tibble::tribble(
    ~USUBJID, ~PARAMCD, ~PARAM, ~AVAL, ~VISIT,
    "01-701-1015", "DIABP", "Diastolic Blood Pressure (mHg)", 51, "BASELINE",
    "01-701-1015", "SYSBP", "Systolic Blood Pressure (mmHg)", 121, "BASELINE",
  )

  expect_error(
    derive_param_map(
      input,
      by_vars = vars(USUBJID, VISIT),
      get_unit_expr = extract_unit(PARAM)
    ),
    paste(
      "It is expected that 'DIABP' is measured in 'mmHg'.\nIn the",
      "input dataset it is measured in 'mHg'."
    )
  )
})

test_that(paste(
  "derive_param_map Test 01.02: MAP parameter NOT added to input dataset",
  "- wrong unit for SYSBP"
), {
  input <- tibble::tribble(
    ~USUBJID,      ~PARAMCD, ~PARAM,                            ~AVAL, ~VISIT,
    "01-701-1015", "DIABP",  "Diastolic Blood Pressure (mmHg)",    51, "BASELINE",
    "01-701-1015", "SYSBP",  "Systolic Blood Pressure (mHg)",     121, "BASELINE",
  )

  expect_error(
    derive_param_map(
      input,
      by_vars = vars(USUBJID, VISIT),
      get_unit_expr = extract_unit(PARAM)
    ),
    paste(
      "It is expected that 'SYSBP' is measured in 'mmHg'.\nIn the",
      "input dataset it is measured in 'mHg'."
    )
  )
})

test_that(paste(
  "derive_param_map Test 01.03: MAP parameter NOT added to input dataset",
  "- wrong unit for PULSE"
), {
  input <- tibble::tribble(
    ~USUBJID, ~PARAMCD, ~PARAM, ~AVAL, ~VISIT,
    "01-701-1015", "DIABP", "Diastolic Blood Pressure (mmHg)", 51, "BASELINE",
    "01-701-1015", "SYSBP", "Systolic Blood Pressure (mmHg)", 121, "BASELINE",
    "01-701-1015", "PULSE", "Pulse (beats/m)", 65, "BASELINE",
  )

  expect_error(
    derive_param_map(
      input,
      by_vars = vars(USUBJID, VISIT),
      hr_code = "PULSE",
      get_unit_expr = extract_unit(PARAM)
    ),
    paste(
      "It is expected that 'PULSE' is measured in 'beats/min'.\nIn the",
      "input dataset it is measured in 'beats/m'."
    )
  )
})

test_that(paste(
  "derive_param_map Test 01.04: MAP parameter NOT added to input dataset",
  "- PARAMCD not set"
), {
  input <- tibble::tribble(
    ~USUBJID, ~PARAMCD, ~PARAM, ~AVAL, ~VISIT,
    "01-701-1015", "DIABP", "Diastolic Blood Pressure (mmHg)", 51, "BASELINE",
    "01-701-1015", "SYSBP", "Systolic Blood Pressure (mmHg)", 121, "BASELINE",
  )

  expect_error(
    derive_param_map(
      input,
      by_vars = vars(USUBJID, VISIT),
      set_values_to = vars(PARAM = "Mean Arterial Pressure"),
      get_unit_expr = extract_unit(PARAM)
    ),
    paste(
      "The following required elements are missing in",
      "`set_values_to`: 'PARAMCD'"
    )
  )
})

## derive_param_map: No obs added (Test 02.xx) ----

test_that("derive_param_map Test 02.01: MAP parameter NOT added to input dataset", {
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
      by_vars = vars(USUBJID, VISIT),
      hr_code = "PULSE",
      get_unit_expr = extract_unit(PARAM)
    ),
    expected_output,
    keys = c("USUBJID", "PARAMCD", "VISIT")
  )
})

## derive_param_map: Obs created (Test 03.xx) ----

maphr <- function(sbp, dbp, hr) {
  dbp + 0.01 * exp(4.14 - 40.74 / hr) * (sbp - dbp)
}

test_that(paste(
  "derive_param_map Test 03.01: MAP parameter (DBP/SBP/PULSE) is correctly",
  "added to input dataset"
), {
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
      by_vars = vars(USUBJID, VISIT),
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

test_that(paste(
  "derive_param_map Test 03.02: MAP parameter (DBP/SBP) is correctly",
  "added to input dataset"
), {
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
      by_vars = vars(USUBJID, VISIT),
      get_unit_expr = extract_unit(PARAM)
    ),
    expected_output,
    keys = c("USUBJID", "PARAMCD", "VISIT")
  )
})
