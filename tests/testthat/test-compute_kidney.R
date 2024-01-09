# compute_egfr ----

## Test 1: EGFR CKD-EPI calculation ----
test_that("compute_egfr Test 1: EGFR CKD-EPI calculation", {
  # Expected values are taken from the National Kidney Foundation's
  # CKD-EPI Creatinine Equation (2021) calculator at
  # https://www.kidney.org/professionals/kdoqi/gfr_calculator
  expect_equal(round(compute_egfr(
    creat = 1.09, creatu = "mg/dL", age = 55, sex = "M", method = "CKD-EPI"
  ), 0L), 80)
})

## Test 2: EGFR CKD-EPI calculation ----
test_that("compute_egfr Test 2: EGFR CKD-EPI calculation", {
  # Expected values are taken from the National Kidney Foundation's
  # CKD-EPI Creatinine Equation (2021) calculator at
  # https://www.kidney.org/professionals/kdoqi/gfr_calculator
  expect_equal(round(compute_egfr(
    creat = 85, creatu = "umol/L", age = 65, sex = "F", method = "CKD-EPI"
  ), 0L), 66)
})


## Test 3: CRCL calculation ----
test_that("compute_egfr Test 3: CRCL calculation", {
  # Expected values are taken from the National Kidney Foundation's
  # CRCL Cockcroft and Gault (1973) calculator at
  # https://www.kidney.org/professionals/kdoqi/gfr_calculatorCoc
  expect_equal(round(compute_egfr(
    creat = 1.09, creatu = "mg/dL", age = 55, sex = "M", weight = 90, method = "CRCL"
  ), 0L), 97)
})


## Test 4: CRCL calculation ----
test_that("compute_egfr Test 4: CRCL calculation", {
  # Expected values are taken from the National Kidney Foundation's
  # CRCL Cockcroft and Gault (1973) calculator at
  # https://www.kidney.org/professionals/kdoqi/gfr_calculatorCoc
  expect_equal(round(compute_egfr(
    creat = 85, creatu = "umol/L", age = 65, sex = "F", weight = 60, method = "CRCL"
  ), 0L), 55)
})


## Test 5: EGFR MDRD calculation ----
test_that("compute_egfr Test 5: EGFR MDRD calculation", {
  # Expected values are taken from the MD Calc
  # MDRD GFR calculator at
  # https://www.mdcalc.com/calc/76/mdrd-gfr-equation
  expect_equal(round(compute_egfr(
    creat = 1.09, creatu = "mg/dL", age = 55, sex = "M", weight = 90, race = "WHITE",
    method = "MDRD"
  ), 1L), 70.2)
})

## Test 6: EGFR MDRD calculation ----
test_that("compute_egfr Test 6: EGFR MDRD calculation", {
  # Expected values are taken from the MD Calc
  # MDRD GFR calculator at
  # https://www.mdcalc.com/calc/76/mdrd-gfr-equation
  expect_equal(round(compute_egfr(
    creat = 85, creatu = "umol/L", age = 65, sex = "F",
    race = "BLACK OR AFRICAN AMERICAN", method = "MDRD"
  ), 1L), 70.6)
})


## Test 7: CKD-EPI calculated on input data ----
test_that("compute_egfr Test 7: CKD-EPI calculated on input data", {
  input <- tibble::tribble(
    ~STUDYID, ~USUBJID, ~AGE, ~SEX, ~RACE, ~WTBL, ~CREATBL, ~CREATBLU,
    "P01", "P01-1001", 55, "M", "WHITE", 90.7, 96.3, "umol/L",
    "P01", "P01-1002", 52, "F", "BLACK OR AFRICAN AMERICAN", 68, 70, "umol/L",
    "P01", "P01-1003", 67, "M", "BLACK OR AFRICAN AMERICAN", 85, 77, "umol/L",
    "P01", "P01-1004", 76, "F", "ASIAN", 60, 65, "umol/L"
  )

  expected_output <- tibble::tribble(
    ~STUDYID, ~USUBJID, ~AGE, ~SEX, ~RACE, ~WTBL, ~CREATBL, ~CREATBLU, ~EGFR,
    "P01", "P01-1001", 55, "M", "WHITE", 90.7, 96.3, "umol/L", 80.2293,
    "P01", "P01-1002", 52, "F", "BLACK OR AFRICAN AMERICAN", 68, 70, "umol/L", 89.7175,
    "P01", "P01-1003", 67, "M", "BLACK OR AFRICAN AMERICAN", 85, 77, "umol/L", 94.5453,
    "P01", "P01-1004", 76, "F", "ASIAN", 60, 65, "umol/L", 84.4646,
  )

  egfr <- input %>%
    dplyr::mutate(
      EGFR = compute_egfr(
        creat = CREATBL, creatu = CREATBLU, age = AGE, weight = WTBL, sex = SEX,
        method = "CKD-EPI"
      ),
      EGFR = round(EGFR, 4L)
    )

  expect_dfs_equal(
    egfr,
    expected_output,
    keys = c("USUBJID")
  )
})
