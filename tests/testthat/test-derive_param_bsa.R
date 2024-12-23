# compute_bsa ----

## compute_bsa: Mosteller method ----
# sqrt (Height x Weight / 3600)

## Test 1: Mosteller method - single height & weight values ----
test_that("compute_bsa Test 1: Mosteller method - single height & weight values", {
  expect_equal(
    round(compute_bsa(height = 170, weight = 75, method = "Mosteller"), 3L),
    1.882
  )
})

## Test 2: Mosteller method - height & weight vectors ----
test_that("compute_bsa Test 2: Mosteller method - height & weight vectors", {
  expect_equal(
    round(compute_bsa(height = c(170, 185), weight = c(75, 90), method = "Mosteller"), 3L),
    c(1.882, 2.151)
  )
})

## Test 3: Mosteller method - height & weight vectors - missing values ----
test_that("compute_bsa Test 3: Mosteller method - height & weight vectors - missing values", {
  expect_equal(
    compute_bsa(height = c(NA, 185), weight = c(75, NA), method = "Mosteller"),
    c(NA_real_, NA_real_)
  )
})

## compute_bsa: DuBois-DuBois method ----
#  FORMULA : 0.007184 x (HGT)^0.725 x WGT^0.425

## Test 4: DuBois-DuBois method - single height & weight values ----
test_that("compute_bsa Test 4: DuBois-DuBois method - single height & weight values", {
  expect_equal(
    round(compute_bsa(height = 170, weight = 75, method = "DuBois-DuBois"), 3L),
    1.864
  )
})

## Test 5: DuBois-DuBois method - height & weight vectors ----
test_that("compute_bsa Test 5: DuBois-DuBois method - height & weight vectors", {
  expect_equal(
    round(compute_bsa(height = c(170, 185), weight = c(75, 90), method = "DuBois-DuBois"), 3L),
    c(1.864, 2.141)
  )
})

## Test 6: DuBois-DuBois method - hgt and wgt vectors - missing values ----
test_that("compute_bsa Test 6: DuBois-DuBois method - hgt and wgt vectors - missing values", {
  expect_equal(
    compute_bsa(height = c(NA, 185), weight = c(75, NA), method = "DuBois-DuBois"),
    c(NA_real_, NA_real_)
  )
})

## compute_bsa: Haycock method (Test 03.xx) ----
# 0.024265 x HGT^0.3964 x WGT^0.5378

## Test 7: Haycock method - single height & weight values ----
test_that("compute_bsa Test 7: Haycock method - single height & weight values", {
  expect_equal(
    round(compute_bsa(height = 170, weight = 75, method = "Haycock"), 3L),
    1.895
  )
})

## Test 8: Haycock method - height & weight vectors ----
test_that("compute_bsa Test 8: Haycock method - height & weight vectors", {
  expect_equal(
    round(compute_bsa(height = c(170, 185), weight = c(75, 90), method = "Haycock"), 3L),
    c(1.895, 2.161)
  )
})

## Test 9: Haycock method - height & weight vectors - missing values ----
test_that("compute_bsa Test 9: Haycock method - height & weight vectors - missing values", {
  expect_equal(
    compute_bsa(height = c(NA, 185), weight = c(75, NA), method = "Haycock"),
    c(NA_real_, NA_real_)
  )
})

## compute_bsa: Gehan-George method ----
# 0.0235 x HGT^0.42246 x WGT^0.51456

## Test 10: Gehan-George method - single height & weight values ----
test_that("compute_bsa Test 10: Gehan-George method - single height & weight values", {
  expect_equal(
    round(compute_bsa(height = 170, weight = 75, method = "Gehan-George"), 3L),
    1.897
  )
})

## Test 11: Gehan-George method - height & weight vectors ----
test_that("compute_bsa Test 11: Gehan-George method - height & weight vectors", {
  expect_equal(
    round(compute_bsa(height = c(170, 185), weight = c(75, 90), method = "Gehan-George"), 3L),
    c(1.897, 2.16)
  )
})

## Test 12: Gehan-George method - height & weight vectors - missing values ----
test_that("compute_bsa Test 12: Gehan-George method - height & weight vectors - missing values", {
  expect_equal(
    compute_bsa(height = c(NA, 185), weight = c(75, NA), method = "Gehan-George"),
    c(NA_real_, NA_real_)
  )
})

## compute_bsa: Boyd method ----
# 0.0003207 x (HGT^0.3) x (1000 x WGT)^(0.7285 - (0.0188 x log10(1000 x WGT)))

## Test 13: Boyd method - single height & weight values ----
test_that("compute_bsa Test 13: Boyd method - single height & weight values", {
  expect_equal(
    round(compute_bsa(height = 170, weight = 75, method = "Boyd"), 3L),
    1.905
  )
})

## Test 14: Boyd method - height & weight vectors ----
test_that("compute_bsa Test 14: Boyd method - height & weight vectors", {
  expect_equal(
    round(compute_bsa(height = c(170, 185), weight = c(75, 90), method = "Boyd"), 3L),
    c(1.905, 2.158)
  )
})

## Test 15: Boyd method - height & weight vectors - missing values ----
test_that("compute_bsa Test 15: Boyd method - height & weight vectors - missing values", {
  expect_equal(
    compute_bsa(height = c(NA, 185), weight = c(75, NA), method = "Boyd"),
    c(NA_real_, NA_real_)
  )
})

## compute_bsa: Fujimoto method ----
# 0.008883 x HGT^0.663 x WGT^0.444

## Test 16: Fujimoto method - single height & weight values ----
test_that("compute_bsa Test 16: Fujimoto method - single height & weight values", {
  expect_equal(
    round(compute_bsa(height = 170, weight = 75, method = "Fujimoto"), 3L),
    1.819
  )
})

## Test 17: Fujimoto method - height & weight vectors ----
test_that("compute_bsa Test 17: Fujimoto method - height & weight vectors", {
  expect_equal(
    round(compute_bsa(height = c(170, 185), weight = c(75, 90), method = "Fujimoto"), 3L),
    c(1.819, 2.086)
  )
})

## Test 18: Fujimoto method - height & weight vectors - missing values ----
test_that("compute_bsa Test 18: Fujimoto method - height & weight vectors - missing values", {
  expect_equal(
    compute_bsa(height = c(NA, 185), weight = c(75, NA), method = "Fujimoto"),
    c(NA_real_, NA_real_)
  )
})

## compute_bsa: Takahira method ----
# 0.007241 x HGT^0.725 x WGT^0.425

## Test 19: Takahira method - single height & weight values ----
test_that("compute_bsa Test 19: Takahira method - single height & weight values", {
  expect_equal(
    round(compute_bsa(height = 170, weight = 75, method = "Takahira"), 3L),
    1.878
  )
})

## Test 20: Takahira method - height & weight vectors ----
test_that("compute_bsa Test 20: Takahira method - height & weight vectors", {
  expect_equal(
    round(compute_bsa(height = c(170, 185), weight = c(75, 90), method = "Takahira"), 3L),
    c(1.878, 2.158)
  )
})

## Test 21: Takahira method - height & weight vectors - missing values ----
test_that("compute_bsa Test 21: Takahira method - height & weight vectors - missing values", {
  expect_equal(
    compute_bsa(height = c(NA, 185), weight = c(75, NA), method = "Takahira"),
    c(NA_real_, NA_real_)
  )
})

## compute_bsa: Check error messages ----

## Test 22: an error is issued if an invalid method is specified ----
test_that("compute_bsa Test 22: an error is issued if an invalid method is specified", {
  expect_error(
    compute_bsa(height = c(170, 185), weight = c(75, 90), method = "unknown-method"),
    class = "assert_character_scalar"
  )
})

# derive_param_bsa ----

## derive_param_bsa: Error checks ----

## Test 23: BSA parameter NOT added - wrong unit for height ----
test_that("derive_param_bsa Test 23: BSA parameter NOT added - wrong unit for height", {
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

## Test 24: BSA parameter NOT added - wrong unit for weight ----
test_that("derive_param_bsa Test 24: BSA parameter NOT added - wrong unit for weight", {
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

## Test 25: BSA parameter NOT added - multiple unit for weight ----
test_that("derive_param_bsa Test 25: BSA parameter NOT added - multiple unit for weight", {
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

## Test 26: BSA parameter NOT added - PARAMCD not set ----
test_that("derive_param_bsa Test 26: BSA parameter NOT added - PARAMCD not set", {
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

## Test 27: BSA parameter NOT added ----
test_that("derive_param_bsa Test 27: BSA parameter NOT added", {
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

  expect_snapshot(
    result <- derive_param_bsa(
      input,
      by_vars = exprs(USUBJID, VISIT),
      method = "Mosteller",
      get_unit_expr = VSSTRESU
    )
  )

  expect_dfs_equal(
    result,
    expected_output,
    keys = c("USUBJID", "PARAMCD", "VISIT")
  )
})

## derive_param_bsa: Obs created ----

mosteller <- function(hgt, wgt) {
  sqrt(hgt * wgt / 3600)
}

## Test 28: BSA parameter (Mosteller Method) is correctly added ----
test_that("derive_param_bsa Test 28: BSA parameter (Mosteller Method) is correctly added", {
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
  0.007184 * hgt^0.725 * wgt^0.425
}

## Test 29: BSA parameter (DuBois-DuBois method) is correctly added ----
test_that("derive_param_bsa Test 29: BSA parameter (DuBois-DuBois method) is correctly added", {
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

## Test 30: BSA parameter (Haycock method) is correctly added ----
test_that("derive_param_bsa Test 30: BSA parameter (Haycock method) is correctly added", {
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

## Test 31: BSA parameter (Gehan-George method) is correctly added ----
test_that("derive_param_bsa Test 31: BSA parameter (Gehan-George method) is correctly added", {
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

## Test 32: BSA parameter (Boyd method) is correctly added ----
test_that("derive_param_bsa Test 32: BSA parameter (Boyd method) is correctly added", {
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

## Test 33: BSA parameter (Fujimoto method) is correctly added ----
test_that("derive_param_bsa Test 33: BSA parameter (Fujimoto method) is correctly added", {
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

## Test 34: BSA parameter (Takahira method) is correctly added ----
test_that("derive_param_bsa Test 34: BSA parameter (Takahira method) is correctly added", {
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

## Test 35: Derive BSA where height is measured only once ----
test_that("derive_param_bsa Test 35: Derive BSA where height is measured only once", {
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
