# derive_var_nfrlt ----

## Test 1: basic single dose study ----
test_that("derive_var_nfrlt Test 1: basic single dose study", {
  input <- tibble::tribble(
    ~USUBJID, ~VISITDY, ~PCTPT,
    "001",    1,        "Pre-dose",
    "001",    1,        "1H Post-dose",
    "001",    1,        "2H Post-dose",
    "001",    1,        "4H Post-dose",
    "001",    1,        "24H Post-dose"
  )

  result <- derive_var_nfrlt(
    input,
    new_var = NFRLT,
    tpt_var = PCTPT,
    visit_day = VISITDY
  )

  expect_equal(
    result$NFRLT,
    c(0, 1, 2, 4, 24)
  )
})

## Test 2: multiple dose study ----
test_that("derive_var_nfrlt Test 2: multiple dose study", {
  input <- tibble::tribble(
    ~USUBJID, ~VISITDY, ~PCTPT,
    "001",    1,        "Pre-dose",
    "001",    1,        "2H Post-dose",
    "001",    8,        "Pre-dose",
    "001",    8,        "2H Post-dose",
    "001",    15,       "Pre-dose",
    "001",    15,       "2H Post-dose"
  )

  result <- derive_var_nfrlt(
    input,
    new_var = NFRLT,
    tpt_var = PCTPT,
    visit_day = VISITDY
  )

  expected_nfrlt <- c(
    0, # Day 1, Pre-dose
    2, # Day 1, 2H Post-dose
    168, # Day 8, Pre-dose (7 * 24 = 168)
    170, # Day 8, 2H Post-dose
    336, # Day 15, Pre-dose (14 * 24 = 336)
    338 # Day 15, 2H Post-dose
  )

  expect_equal(result$NFRLT, expected_nfrlt)
})

## Test 3: custom treatment duration (IV infusion) ----
test_that("derive_var_nfrlt Test 3: custom treatment duration", {
  input <- tibble::tribble(
    ~USUBJID, ~VISITDY, ~PCTPT,
    "001",    1,        "Pre-dose",
    "001",    1,        "EOI",
    "001",    1,        "1H Post EOI",
    "001",    1,        "10MIN PRE EOI"
  )

  result <- derive_var_nfrlt(
    input,
    new_var = NFRLT,
    tpt_var = PCTPT,
    visit_day = VISITDY,
    treatment_duration = 2
  )

  expect_equal(
    result$NFRLT,
    c(0, 2, 3, 2 - 10 / 60)
  )
})

## Test 4: custom first dose day ----
test_that("derive_var_nfrlt Test 4: custom first dose day", {
  input <- tibble::tribble(
    ~USUBJID, ~VISITDY, ~PCTPT,
    "001",    3,        "Pre-dose",
    "001",    3,        "2H Post-dose",
    "001",    3,        "8H Post-dose",
    "001",    3,        "24H Post-dose"
  )

  result <- derive_var_nfrlt(
    input,
    new_var = NFRLT,
    tpt_var = PCTPT,
    visit_day = VISITDY,
    first_dose_day = 3
  )

  expect_equal(
    result$NFRLT,
    c(0, 2, 8, 24) # Day 4 is 24 hours after Day 3
  )
})

## Test 5: handles NA in timepoint ----
test_that("derive_var_nfrlt Test 5: handles NA in timepoint", {
  input <- tibble::tribble(
    ~USUBJID, ~VISITDY, ~PCTPT,
    "001",    1,        "Pre-dose",
    "001",    1,        NA_character_,
    "001",    1,        "2H Post-dose"
  )

  result <- derive_var_nfrlt(
    input,
    new_var = NFRLT,
    tpt_var = PCTPT,
    visit_day = VISITDY
  )

  expect_equal(
    result$NFRLT,
    c(0, NA_real_, 2)
  )
})

## Test 6: handles NA in visit day ----
test_that("derive_var_nfrlt Test 6: handles NA in visit day", {
  input <- tibble::tribble(
    ~USUBJID, ~VISITDY,  ~PCTPT,
    "001",    1,         "Pre-dose",
    "001",    NA_real_,  "2H Post-dose",
    "001",    1,         "4H Post-dose"
  )

  result <- derive_var_nfrlt(
    input,
    new_var = NFRLT,
    tpt_var = PCTPT,
    visit_day = VISITDY
  )

  expect_equal(
    result$NFRLT,
    c(0, NA_real_, 4)
  )
})

## Test 7: custom new variable name ----
test_that("derive_var_nfrlt Test 7: custom new variable name", {
  input <- tibble::tribble(
    ~USUBJID, ~VISITDY, ~PCTPT,
    "001",    1,        "Pre-dose",
    "001",    1,        "1H Post-dose"
  )

  result <- derive_var_nfrlt(
    input,
    new_var = NOMTIME,
    tpt_var = PCTPT,
    visit_day = VISITDY
  )

  expect_true("NOMTIME" %in% names(result))
  expect_equal(result$NOMTIME, c(0, 1))
})

## Test 8: preserves existing columns ----
test_that("derive_var_nfrlt Test 8: preserves existing columns", {
  input <- tibble::tribble(
    ~USUBJID, ~VISITDY, ~PCTPT, ~AVAL,
    "001", 1, "Pre-dose", 10.5,
    "001", 1, "1H Post-dose", 25.3
  )

  result <- derive_var_nfrlt(
    input,
    new_var = NFRLT,
    tpt_var = PCTPT,
    visit_day = VISITDY
  )

  expect_true(all(c("USUBJID", "VISITDY", "PCTPT", "AVAL") %in% names(result)))
  expect_equal(result$AVAL, c(10.5, 25.3))
})

## Test 9: does not error if tpt_var missing from dataset ----
test_that("derive_var_nfrlt Test 9: handles missing tpt_var from dataset gracefully", {
  input <- tibble::tribble(
    ~USUBJID, ~VISITDY,
    "001",    1,
    "001",    8
  )

  # Should not error - treats as if tpt_var was not provided
  expect_no_error(
    result <- derive_var_nfrlt(
      input,
      new_var = NFRLT,
      tpt_var = PCTPT, # PCTPT doesn't exist in dataset
      visit_day = VISITDY
    )
  )

  # NFRLT should be calculated with timepoint contribution = 0
  result <- derive_var_nfrlt(
    input,
    new_var = NFRLT,
    tpt_var = PCTPT,
    visit_day = VISITDY
  )

  expect_equal(
    result$NFRLT,
    c(0, 168) # Just day offsets, no timepoint contribution
  )
})

## Test 10: error if visit_day missing ----
test_that("derive_var_nfrlt Test 10: error if visit_day missing", {
  input <- tibble::tribble(
    ~USUBJID, ~PCTPT,
    "001",    "Pre-dose"
  )

  expect_error(
    derive_var_nfrlt(
      input,
      new_var = NFRLT,
      tpt_var = PCTPT,
      visit_day = VISITDY
    ),
    class = "assert_data_frame"
  )
})

## Test 11: error if first_dose_day not positive ----
test_that("derive_var_nfrlt Test 11: error if first_dose_day not positive", {
  input <- tibble::tribble(
    ~USUBJID, ~VISITDY, ~PCTPT,
    "001",    1,        "Pre-dose"
  )

  expect_error(
    derive_var_nfrlt(
      input,
      new_var = NFRLT,
      tpt_var = PCTPT,
      visit_day = VISITDY,
      first_dose_day = 0
    ),
    regexp = "must be positive"
  )

  expect_error(
    derive_var_nfrlt(
      input,
      new_var = NFRLT,
      tpt_var = PCTPT,
      visit_day = VISITDY,
      first_dose_day = -1
    ),
    regexp = "must be positive"
  )
})

## Test 12: error if treatment_duration is negative ----
test_that("derive_var_nfrlt Test 12: error if treatment_duration is negative", {
  input <- tibble::tribble(
    ~USUBJID, ~VISITDY, ~PCTPT,
    "001",    1,        "EOI"
  )

  expect_error(
    derive_var_nfrlt(
      input,
      new_var = NFRLT,
      tpt_var = PCTPT,
      visit_day = VISITDY,
      treatment_duration = -2
    ),
    regexp = "must be non-negative"
  )

  # Zero is allowed
  expect_no_error(
    derive_var_nfrlt(
      input,
      new_var = NFRLT,
      tpt_var = PCTPT,
      visit_day = VISITDY,
      treatment_duration = 0
    )
  )
})

## Test 13: works with different timepoint patterns ----
test_that("derive_var_nfrlt Test 13: works with different timepoint patterns", {
  input <- tibble::tribble(
    ~USUBJID, ~VISITDY, ~PCTPT,
    "001",    1,        "5 MIN PREDOSE",
    "001",    1,        "5 MIN BEFORE",
    "001",    1,        "30MIN AFTER END OF INFUSION",
    "001",    1,        "24 HR POST INF",
    "001",    8,        "1H30M"
  )

  result <- derive_var_nfrlt(
    input,
    new_var = NFRLT,
    tpt_var = PCTPT,
    visit_day = VISITDY,
    treatment_duration = 1
  )

  expect_equal(
    result$NFRLT,
    c(
      -5 / 60,
      -5 / 60,
      1.5,
      25,
      168 + 1.5
    )
  )
})

## Test 14: handles decimal visit days ----
test_that("derive_var_nfrlt Test 14: handles decimal visit days", {
  input <- tibble::tribble(
    ~USUBJID, ~VISITDY, ~PCTPT,
    "001",    1.5,      "Pre-dose",
    "001",    1.5,      "2H Post-dose"
  )

  result <- derive_var_nfrlt(
    input,
    new_var = NFRLT,
    tpt_var = PCTPT,
    visit_day = VISITDY
  )

  expect_equal(
    result$NFRLT,
    c(
      0.5 * 24,
      0.5 * 24 + 2
    )
  )
})

## Test 15: multiple subjects ----
test_that("derive_var_nfrlt Test 15: multiple subjects", {
  input <- tibble::tribble(
    ~USUBJID, ~VISITDY, ~PCTPT,
    "001",    1,        "Pre-dose",
    "001",    1,        "2H Post-dose",
    "002",    1,        "Pre-dose",
    "002",    1,        "4H Post-dose"
  )

  result <- derive_var_nfrlt(
    input,
    new_var = NFRLT,
    tpt_var = PCTPT,
    visit_day = VISITDY
  )

  expect_equal(
    result$NFRLT,
    c(0, 2, 0, 4)
  )
})

## Test 16: range_method parameter works ----
test_that("derive_var_nfrlt Test 16: range_method parameter works", {
  input <- tibble::tribble(
    ~USUBJID, ~VISITDY, ~PCTPT,
    "001",    1,        "0-6h Post-dose"
  )

  # Default midpoint
  result_midpoint <- derive_var_nfrlt(
    input,
    new_var = NFRLT,
    tpt_var = PCTPT,
    visit_day = VISITDY
  )
  expect_equal(result_midpoint$NFRLT, 3)

  # Start
  result_start <- derive_var_nfrlt(
    input,
    new_var = NFRLT,
    tpt_var = PCTPT,
    visit_day = VISITDY,
    range_method = "start"
  )
  expect_equal(result_start$NFRLT, 0)

  # End
  result_end <- derive_var_nfrlt(
    input,
    new_var = NFRLT,
    tpt_var = PCTPT,
    visit_day = VISITDY,
    range_method = "end"
  )
  expect_equal(result_end$NFRLT, 6)
})

## Test 17: using "Before" and "After" terminology ----
test_that("derive_var_nfrlt Test 17: works with Before and After terminology", {
  input <- tibble::tribble(
    ~USUBJID, ~VISITDY, ~PCTPT,
    "001",    1,        "Before",
    "001",    1,        "1H After",
    "001",    1,        "2H After"
  )

  result <- derive_var_nfrlt(
    input,
    new_var = NFRLT,
    tpt_var = PCTPT,
    visit_day = VISITDY
  )

  expect_equal(
    result$NFRLT,
    c(0, 1, 2)
  )
})

## Test 18: EOI and EOT treated identically ----
test_that("derive_var_nfrlt Test 18: EOI and EOT produce same results", {
  input_eoi <- tibble::tribble(
    ~USUBJID, ~VISITDY, ~PCTPT,
    "001",    1,        "EOI",
    "001",    1,        "1H POST EOI",
    "001",    1,        "10MIN PRE EOI"
  )

  input_eot <- tibble::tribble(
    ~USUBJID, ~VISITDY, ~PCTPT,
    "001",    1,        "EOT",
    "001",    1,        "1H POST EOT",
    "001",    1,        "10MIN PRE EOT"
  )

  result_eoi <- derive_var_nfrlt(
    input_eoi,
    new_var = NFRLT,
    tpt_var = PCTPT,
    visit_day = VISITDY,
    treatment_duration = 2
  )

  result_eot <- derive_var_nfrlt(
    input_eot,
    new_var = NFRLT,
    tpt_var = PCTPT,
    visit_day = VISITDY,
    treatment_duration = 2
  )

  expect_equal(result_eoi$NFRLT, result_eot$NFRLT)
})

## Test 19: default treatment_duration is 0 (oral) ----
test_that("derive_var_nfrlt Test 19: default treatment_duration is 0", {
  input <- tibble::tribble(
    ~USUBJID, ~VISITDY, ~PCTPT,
    "001",    1,        "Pre-dose",
    "001",    1,        "EOT",
    "001",    1,        "1H POST EOT"
  )

  result <- derive_var_nfrlt(
    input,
    new_var = NFRLT,
    tpt_var = PCTPT,
    visit_day = VISITDY
  )

  expect_equal(
    result$NFRLT,
    c(0, 0, 1)
  )
})

## Test 20: treatment_duration with PRE EOI/EOT patterns ----
test_that("derive_var_nfrlt Test 20: PRE EOI/EOT correctly uses treatment_duration", {
  input <- tibble::tribble(
    ~USUBJID, ~VISITDY, ~PCTPT,
    "001",    1,        "10MIN PRE EOI",
    "001",    1,        "10MIN BEFORE EOT"
  )

  # With 0 duration
  result_0 <- derive_var_nfrlt(
    input,
    new_var = NFRLT,
    tpt_var = PCTPT,
    visit_day = VISITDY,
    treatment_duration = 0
  )
  expect_equal(result_0$NFRLT, c(-10 / 60, -10 / 60))

  # With 2 hour duration
  result_2 <- derive_var_nfrlt(
    input,
    new_var = NFRLT,
    tpt_var = PCTPT,
    visit_day = VISITDY,
    treatment_duration = 2
  )
  expect_equal(result_2$NFRLT, c(2 - 10 / 60, 2 - 10 / 60))
})

## Test 21: works without tpt_var (for EX records) ----
test_that("derive_var_nfrlt Test 21: works without tpt_var", {
  input <- tibble::tribble(
    ~USUBJID, ~VISITDY,
    "001",    1,
    "001",    8,
    "001",    15
  )

  result <- derive_var_nfrlt(
    input,
    new_var = NFRLT,
    visit_day = VISITDY
  )

  # Without tpt_var, NFRLT should just be day offset * 24
  expect_equal(
    result$NFRLT,
    c(0, 168, 336)
  )
})

## Test 22: handles unscheduled visits with VISIT variable ----
test_that("derive_var_nfrlt Test 22: handles unscheduled visits", {
  input <- tibble::tribble(
    ~USUBJID, ~VISITDY, ~VISIT,        ~PCTPT,
    "001",    1,        "VISIT 1",     "Pre-dose",
    "001",    1,        "VISIT 1",     "2H Post-dose",
    "001",    NA_real_, "UNSCHEDULED", "Pre-dose",
    "001",    NA_real_, "UNSCHEDULED", "2H Post-dose"
  )

  result <- derive_var_nfrlt(
    input,
    new_var = NFRLT,
    tpt_var = PCTPT,
    visit_day = VISITDY,
    set_values_to_na = VISIT == "UNSCHEDULED"
  )

  expect_equal(
    result$NFRLT,
    c(0, 2, NA_real_, NA_real_)
  )
})

## Test 23: handles early discontinuation visits ----
test_that("derive_var_nfrlt Test 23: handles early discontinuation", {
  input <- tibble::tribble(
    ~USUBJID, ~VISITDY, ~VISIT,                              ~PCTPT,
    "001",    1,        "VISIT 1",                           "Pre-dose",
    "001",    1,        "VISIT 1",                           "2H Post-dose",
    "001",    NA_real_, "STUDY DRUG EARLY DISCONTINUATION",  "Pre-dose"
  )

  result <- derive_var_nfrlt(
    input,
    new_var = NFRLT,
    tpt_var = PCTPT,
    visit_day = VISITDY,
    set_values_to_na = VISIT == "STUDY DRUG EARLY DISCONTINUATION"
  )

  expect_equal(
    result$NFRLT,
    c(0, 2, NA_real_)
  )
})

## Test 24: handles multiple exclusion criteria ----
test_that("derive_var_nfrlt Test 24: handles multiple exclusion criteria", {
  input <- tibble::tribble(
    ~USUBJID, ~VISITDY, ~VISIT,                              ~PCTPT,
    "001",    1,        "VISIT 1",                           "Pre-dose",
    "001",    8,        "VISIT 2",                           "Pre-dose",
    "001",    NA_real_, "UNSCHEDULED",                       "Pre-dose",
    "001",    NA_real_, "STUDY DRUG EARLY DISCONTINUATION",  "2H Post-dose"
  )

  result <- derive_var_nfrlt(
    input,
    new_var = NFRLT,
    tpt_var = PCTPT,
    visit_day = VISITDY,
    set_values_to_na = VISIT %in% c("UNSCHEDULED", "STUDY DRUG EARLY DISCONTINUATION")
  )

  expect_equal(
    result$NFRLT,
    c(0, 168, NA_real_, NA_real_)
  )
})

## Test 25: works without tpt_var and with set_values_to_na ----
test_that("derive_var_nfrlt Test 25: works without tpt_var and with set_values_to_na", {
  input <- tibble::tribble(
    ~USUBJID, ~VISITDY, ~VISIT,
    "001",    1,        "VISIT 1",
    "001",    8,        "VISIT 2",
    "001",    NA_real_, "UNSCHEDULED"
  )

  result <- derive_var_nfrlt(
    input,
    new_var = NFRLT,
    visit_day = VISITDY,
    set_values_to_na = VISIT == "UNSCHEDULED"
  )

  expect_equal(
    result$NFRLT,
    c(0, 168, NA_real_)
  )
})

## Test 26: tpt_var provided but doesn't exist in dataset ----
test_that("derive_var_nfrlt Test 26: handles missing tpt_var gracefully", {
  input <- tibble::tribble(
    ~USUBJID, ~VISITDY,
    "001",    1,
    "001",    8
  )

  # If tpt_var doesn't exist in dataset, should behave as if not provided
  result <- derive_var_nfrlt(
    input,
    new_var = NFRLT,
    tpt_var = PCTPT, # PCTPT doesn't exist
    visit_day = VISITDY
  )

  expect_equal(
    result$NFRLT,
    c(0, 168)
  )
})

## Test 27: NA in VISITDY results in NA for NFRLT ----
test_that("derive_var_nfrlt Test 27: NA in VISITDY results in NA", {
  input <- tibble::tribble(
    ~USUBJID, ~VISITDY, ~PCTPT,
    "001",    1,        "Pre-dose",
    "001",    NA_real_, "Pre-dose"
  )

  result <- derive_var_nfrlt(
    input,
    new_var = NFRLT,
    tpt_var = PCTPT,
    visit_day = VISITDY
  )

  expect_equal(
    result$NFRLT,
    c(0, NA_real_)
  )
})

## Test 28: set_values_to_na with complex condition ----
test_that("derive_var_nfrlt Test 28: set_values_to_na with complex condition", {
  input <- tibble::tribble(
    ~USUBJID, ~VISITDY, ~VISIT,        ~PCTPT,
    "001",    1,        "VISIT 1",     "Pre-dose",
    "001",    NA_real_, "UNSCHEDULED", "Pre-dose",
    "002",    1,        "VISIT 1",     "Pre-dose",
    "002",    NA_real_, "VISIT 2",     "Pre-dose"
  )

  # Set to NA if VISITDY is NA OR VISIT is UNSCHEDULED
  result <- derive_var_nfrlt(
    input,
    new_var = NFRLT,
    tpt_var = PCTPT,
    visit_day = VISITDY,
    set_values_to_na = is.na(VISITDY) | VISIT == "UNSCHEDULED"
  )

  expect_equal(
    result$NFRLT,
    c(0, NA_real_, 0, NA_real_)
  )
})
## Test 29: treatment_duration as variable ----
test_that("derive_var_nfrlt Test 29: treatment_duration as variable", {
  input <- tibble::tribble(
    ~USUBJID, ~VISITDY, ~PCTPT,           ~EXDUR,
    "001",    1,        "Pre-dose",       1,
    "001",    1,        "EOI",            1,
    "001",    1,        "1H POST EOI",    1,
    "002",    1,        "Pre-dose",       2,
    "002",    1,        "EOI",            2,
    "002",    1,        "1H POST EOI",    2
  )

  result <- derive_var_nfrlt(
    input,
    new_var = NFRLT,
    tpt_var = PCTPT,
    visit_day = VISITDY,
    treatment_duration = EXDUR # Variable!
  )

  expect_equal(
    result$NFRLT,
    c(0, 1, 2, 0, 2, 3)
  )
})

## Test 30: treatment_duration as variable with NA ----
test_that("derive_var_nfrlt Test 30: handles NA in treatment_duration variable", {
  input <- tibble::tribble(
    ~USUBJID, ~VISITDY, ~PCTPT,     ~EXDUR,
    "001",    1,        "Pre-dose", 1,
    "001",    1,        "EOI",      NA_real_,
    "001",    1,        "Pre-dose", 2
  )

  result <- derive_var_nfrlt(
    input,
    new_var = NFRLT,
    tpt_var = PCTPT,
    visit_day = VISITDY,
    treatment_duration = EXDUR
  )

  expect_equal(
    result$NFRLT,
    c(0, NA_real_, 0)
  )
})

## Test 31: handles negative visit days (no Day 0) ----
test_that("derive_var_nfrlt Test 31: handles negative visit days correctly", {
  input <- tibble::tribble(
    ~USUBJID, ~VISITDY, ~PCTPT,
    "001",    -14,      "Pre-dose",
    "001",    -7,       "Pre-dose",
    "001",    -2,       "Pre-dose",
    "001",    -1,       "Pre-dose",
    "001",    1,        "Pre-dose",
    "001",    2,        "Pre-dose"
  )

  result <- derive_var_nfrlt(
    input,
    new_var = NFRLT,
    tpt_var = PCTPT,
    visit_day = VISITDY,
    first_dose_day = 1
  )

  expect_equal(
    result$NFRLT,
    c(-336, -168, -48, -24, 0, 24)
  )
})

## Test 32: negative days with timepoints ----
test_that("derive_var_nfrlt Test 32: negative days with timepoints", {
  input <- tibble::tribble(
    ~USUBJID, ~VISITDY, ~PCTPT,
    "001",    -1,       "Pre-dose",
    "001",    -1,       "2H Post-dose",
    "001",    1,        "Pre-dose",
    "001",    1,        "2H Post-dose"
  )

  result <- derive_var_nfrlt(
    input,
    new_var = NFRLT,
    tpt_var = PCTPT,
    visit_day = VISITDY,
    first_dose_day = 1
  )

  expect_equal(
    result$NFRLT,
    c(-24, -22, 0, 2)
  )
})

## Test 33: first_dose_day = 7 with negative days ----
test_that("derive_var_nfrlt Test 33: first_dose_day = 7 with negative days", {
  input <- tibble::tribble(
    ~USUBJID, ~VISITDY, ~PCTPT,
    "001",    -1,       "Pre-dose",
    "001",    1,        "Pre-dose",
    "001",    6,        "Pre-dose",
    "001",    7,        "Pre-dose",
    "001",    8,        "Pre-dose"
  )

  result <- derive_var_nfrlt(
    input,
    new_var = NFRLT,
    tpt_var = PCTPT,
    visit_day = VISITDY,
    first_dose_day = 7
  )

  expect_equal(
    result$NFRLT,
    c(-168, -144, -24, 0, 24)
  )
})

## Test 34: all positive days (no adjustment needed) ----
test_that("derive_var_nfrlt Test 34: all positive days require no adjustment", {
  input <- tibble::tribble(
    ~USUBJID, ~VISITDY, ~PCTPT,
    "001",    1,        "Pre-dose",
    "001",    2,        "Pre-dose",
    "001",    8,        "Pre-dose",
    "001",    15,       "Pre-dose"
  )

  result <- derive_var_nfrlt(
    input,
    new_var = NFRLT,
    tpt_var = PCTPT,
    visit_day = VISITDY,
    first_dose_day = 1
  )

  # All positive days, normal calculation
  expect_equal(
    result$NFRLT,
    c(0, 24, 168, 336)
  )
})

## Test 35: edge case - Day -1 to Day 1 boundary ----
test_that("derive_var_nfrlt Test 35: edge case Day -1 to Day 1", {
  input <- tibble::tribble(
    ~USUBJID, ~VISITDY, ~PCTPT,
    "001",    -1,       "12H Post-dose",
    "001",    1,        "12H Post-dose"
  )

  result <- derive_var_nfrlt(
    input,
    new_var = NFRLT,
    tpt_var = PCTPT,
    visit_day = VISITDY,
    first_dose_day = 1
  )

  expect_equal(
    result$NFRLT,
    c(-12, 12)
  )
})

## Test 36: negative days without tpt_var ----
test_that("derive_var_nfrlt Test 36: negative days without tpt_var", {
  input <- tibble::tribble(
    ~USUBJID, ~VISITDY,
    "001",    -7,
    "001",    -1,
    "001",    1,
    "001",    8
  )

  result <- derive_var_nfrlt(
    input,
    new_var = NFRLT,
    visit_day = VISITDY,
    first_dose_day = 1
  )

  expect_equal(
    result$NFRLT,
    c(-168, -24, 0, 168)
  )
})

## Test 37: first_dose_day = 8, mix of negative and positive ----
test_that("derive_var_nfrlt Test 37: first_dose_day = 8 with mixed days", {
  input <- tibble::tribble(
    ~USUBJID, ~VISITDY, ~PCTPT,
    "001",    -14,      "Pre-dose",
    "001",    -7,       "Pre-dose",
    "001",    -1,       "Pre-dose",
    "001",    1,        "Pre-dose",
    "001",    7,        "Pre-dose",
    "001",    8,        "Pre-dose",
    "001",    15,       "Pre-dose"
  )

  result <- derive_var_nfrlt(
    input,
    new_var = NFRLT,
    tpt_var = PCTPT,
    visit_day = VISITDY,
    first_dose_day = 8
  )

  expect_equal(
    result$NFRLT,
    c(-504, -336, -192, -168, -24, 0, 168)
  )
})

## Test 38: different output units ----
test_that("derive_var_nfrlt: different output units work correctly", {
  adpc <- tribble(
    ~USUBJID, ~VISITDY, ~PCTPT,
    "001",    1,        "Pre-dose",
    "001",    1,        "12H Post-dose",
    "001",    8,        "Pre-dose"
  )

  # Hours
  result_hours <- derive_var_nfrlt(
    adpc,
    new_var = NFRLT,
    out_unit = "hours",
    tpt_var = PCTPT,
    visit_day = VISITDY
  )
  expect_equal(result_hours$NFRLT, c(0, 12, 168))

  # Days
  result_days <- derive_var_nfrlt(
    adpc,
    new_var = NFRLTDY,
    out_unit = "days",
    tpt_var = PCTPT,
    visit_day = VISITDY
  )
  expect_equal(result_days$NFRLTDY, c(0, 0.5, 7))

  # Weeks
  result_weeks <- derive_var_nfrlt(
    adpc,
    new_var = NFRLTWK,
    out_unit = "weeks",
    tpt_var = PCTPT,
    visit_day = VISITDY
  )
  expect_equal(result_weeks$NFRLTWK, c(0, 12 / 168, 1))

  # Minutes
  result_min <- derive_var_nfrlt(
    adpc,
    new_var = NFRLTMIN,
    out_unit = "minutes",
    tpt_var = PCTPT,
    visit_day = VISITDY
  )
  expect_equal(result_min$NFRLTMIN, c(0, 720, 10080))
})

## Test 39: invalid out_unit ----
test_that("derive_var_nfrlt: invalid out_unit throws error", {
  adpc <- tribble(
    ~USUBJID, ~VISITDY, ~PCTPT,
    "001",    1,        "Pre-dose"
  )

  expect_error(
    derive_var_nfrlt(
      adpc,
      new_var = NFRLT,
      out_unit = "seconds",
      tpt_var = PCTPT,
      visit_day = VISITDY
    ),
    class = "assert_character_scalar"
  )
})

## Test 40: weeks output with weekly dosing ----
test_that("derive_var_nfrlt: weeks output for weekly dosing study", {
  adpc_weekly <- tribble(
    ~USUBJID, ~VISITDY, ~PCTPT,
    "001",    1,        "Pre-dose",
    "001",    8,        "Pre-dose",
    "001",    15,       "Pre-dose",
    "001",    22,       "Pre-dose"
  )

  result <- derive_var_nfrlt(
    adpc_weekly,
    new_var = NFRLTWK,
    out_unit = "weeks",
    tpt_var = PCTPT,
    visit_day = VISITDY
  )

  expect_equal(result$NFRLTWK, c(0, 1, 2, 3))
})

## Test 41: minutes output for short-term PK ----
test_that("derive_var_nfrlt: minutes output for short-term PK study", {
  adpc_short <- tribble(
    ~USUBJID, ~VISITDY, ~PCTPT,
    "001",    1,        "Pre-dose",
    "001",    1,        "5 MIN POST",
    "001",    1,        "15 MIN POST",
    "001",    1,        "30 MIN POST"
  )

  result <- derive_var_nfrlt(
    adpc_short,
    new_var = NFRLTMIN,
    out_unit = "minutes",
    tpt_var = PCTPT,
    visit_day = VISITDY
  )

  expect_equal(result$NFRLTMIN, c(0, 5, 15, 30))
})

## Test 42: unit variable creation ----
test_that("derive_var_nfrlt: unit variable is created correctly", {
  adpc <- tribble(
    ~USUBJID, ~VISITDY, ~PCTPT,
    "001",    1,        "Pre-dose",
    "001",    1,        "12H Post-dose",
    "001",    8,        "Pre-dose"
  )

  # Hours with unit variable
  result_hours <- derive_var_nfrlt(
    adpc,
    new_var = NFRLT,
    new_var_unit = FRLTU,
    out_unit = "hours",
    tpt_var = PCTPT,
    visit_day = VISITDY
  )
  expect_equal(result_hours$NFRLT, c(0, 12, 168))
  expect_equal(result_hours$FRLTU, c("hours", "hours", "hours"))

  # Days with unit variable
  result_days <- derive_var_nfrlt(
    adpc,
    new_var = NFRLTDY,
    new_var_unit = FRLTDYU,
    out_unit = "days",
    tpt_var = PCTPT,
    visit_day = VISITDY
  )
  expect_equal(result_days$NFRLTDY, c(0, 0.5, 7))
  expect_equal(result_days$FRLTDYU, c("days", "days", "days"))

  # Weeks with unit variable
  result_weeks <- derive_var_nfrlt(
    adpc,
    new_var = NFRLTWK,
    new_var_unit = FRLTU,
    out_unit = "weeks",
    tpt_var = PCTPT,
    visit_day = VISITDY
  )
  expect_equal(result_weeks$FRLTU, c("weeks", "weeks", "weeks"))

  # Minutes with unit variable
  result_min <- derive_var_nfrlt(
    adpc,
    new_var = NFRLTMIN,
    new_var_unit = FRLTU,
    out_unit = "minutes",
    tpt_var = PCTPT,
    visit_day = VISITDY
  )
  expect_equal(result_min$FRLTU, c("minutes", "minutes", "minutes"))
})

## Test 43: unit variable with NA values ----
test_that("derive_var_nfrlt Test 43: unit variable handles NA correctly", {
  adpc <- tribble(
    ~USUBJID, ~VISITDY, ~PCTPT,
    "001",    1,        "Pre-dose",
    "001",    NA_real_, "Pre-dose"
  )

  result <- derive_var_nfrlt(
    adpc,
    new_var = NFRLT,
    new_var_unit = FRLTU,
    tpt_var = PCTPT,
    visit_day = VISITDY
  )

  expect_equal(result$NFRLT, c(0, NA_real_))
  expect_equal(result$FRLTU, c("HOURS", NA_character_))
})

## Test 44: unit variable with set_values_to_na ----
test_that("derive_var_nfrlt Test 44: unit variable respects set_values_to_na", {
  adpc <- tribble(
    ~USUBJID, ~VISITDY, ~VISIT,        ~PCTPT,
    "001",    1,        "VISIT 1",     "Pre-dose",
    "001",    NA_real_, "UNSCHEDULED", "Pre-dose"
  )

  result <- derive_var_nfrlt(
    adpc,
    new_var = NFRLT,
    new_var_unit = FRLTU,
    tpt_var = PCTPT,
    visit_day = VISITDY,
    set_values_to_na = VISIT == "UNSCHEDULED"
  )

  expect_equal(result$NFRLT, c(0, NA_real_))
  expect_equal(result$FRLTU, c("HOURS", NA_character_))
})

## Test 45: no unit variable when not requested ----
test_that("derive_var_nfrlt Test 45: no unit variable created when not requested", {
  adpc <- tribble(
    ~USUBJID, ~VISITDY, ~PCTPT,
    "001",    1,        "Pre-dose"
  )

  result <- derive_var_nfrlt(
    adpc,
    new_var = NFRLT,
    tpt_var = PCTPT,
    visit_day = VISITDY
  )

  expect_false("FRLTU" %in% names(result))
})

## Test 46: unit variable preserves case ----
test_that("derive_var_nfrlt Test 46: unit variable preserves user's case", {
  adpc <- tribble(
    ~USUBJID, ~VISITDY, ~PCTPT,
    "001",    1,        "Pre-dose"
  )

  # Lowercase
  result_lower <- derive_var_nfrlt(
    adpc,
    new_var = NFRLT,
    new_var_unit = FRLTU,
    out_unit = "hours",
    tpt_var = PCTPT,
    visit_day = VISITDY
  )
  expect_equal(result_lower$FRLTU, "hours")

  # Uppercase
  result_upper <- derive_var_nfrlt(
    adpc,
    new_var = NFRLT,
    new_var_unit = FRLTU,
    out_unit = "HOURS",
    tpt_var = PCTPT,
    visit_day = VISITDY
  )
  expect_equal(result_upper$FRLTU, "HOURS")

  # Mixed case
  result_mixed <- derive_var_nfrlt(
    adpc,
    new_var = NFRLT,
    new_var_unit = FRLTU,
    out_unit = "Hours",
    tpt_var = PCTPT,
    visit_day = VISITDY
  )
  expect_equal(result_mixed$FRLTU, "Hours")
})

## Test 47: deriving multiple time variables ----
test_that("derive_var_nfrlt Test 47: deriving multiple time variables", {
  adpc <- tribble(
    ~USUBJID, ~VISITDY, ~PCTPT,
    "001",    1,        "Pre-dose",
    "001",    1,        "12H Post-dose",
    "001",    8,        "Pre-dose"
  )

  result <- adpc %>%
    derive_var_nfrlt(
      new_var = NFRLT,
      new_var_unit = FRLTU,
      tpt_var = PCTPT,
      visit_day = VISITDY
    ) %>%
    derive_var_nfrlt(
      new_var = NFRLTDY,
      new_var_unit = FRLTDYU,
      out_unit = "DAYS",
      tpt_var = PCTPT,
      visit_day = VISITDY
    )

  expect_equal(result$NFRLT, c(0, 12, 168))
  expect_equal(result$FRLTU, c("HOURS", "HOURS", "HOURS"))
  expect_equal(result$NFRLTDY, c(0, 0.5, 7))
  expect_equal(result$FRLTDYU, c("DAYS", "DAYS", "DAYS"))
})

## Test 48: case-insensitive out_unit validation ----
test_that("derive_var_nfrlt Test 48: out_unit is case-insensitive", {
  adpc <- tribble(
    ~USUBJID, ~VISITDY, ~PCTPT,
    "001",    1,        "Pre-dose"
  )

  # Should all work (case-insensitive validation)
  expect_no_error(
    derive_var_nfrlt(adpc,
      new_var = NFRLT, out_unit = "HOURS",
      tpt_var = PCTPT, visit_day = VISITDY
    )
  )
  expect_no_error(
    derive_var_nfrlt(adpc,
      new_var = NFRLT, out_unit = "Hours",
      tpt_var = PCTPT, visit_day = VISITDY
    )
  )
  expect_no_error(
    derive_var_nfrlt(adpc,
      new_var = NFRLT, out_unit = "DAYS",
      tpt_var = PCTPT, visit_day = VISITDY
    )
  )
  expect_no_error(
    derive_var_nfrlt(adpc,
      new_var = NFRLT, out_unit = "Weeks",
      tpt_var = PCTPT, visit_day = VISITDY
    )
  )
  expect_no_error(
    derive_var_nfrlt(adpc,
      new_var = NFRLT, out_unit = "MINUTES",
      tpt_var = PCTPT, visit_day = VISITDY
    )
  )
})

## Test 49: accepts all unit variations ----
test_that("derive_var_nfrlt Test 49: accepts all unit variations", {
  adpc <- tribble(
    ~USUBJID, ~VISITDY, ~PCTPT,
    "001",    1,        "Pre-dose",
    "001",    1,        "12H Post-dose"
  )

  # Days variations
  expect_no_error(derive_var_nfrlt(adpc,
    new_var = NFRLT, out_unit = "day",
    tpt_var = PCTPT, visit_day = VISITDY
  ))
  expect_no_error(derive_var_nfrlt(adpc,
    new_var = NFRLT, out_unit = "days",
    tpt_var = PCTPT, visit_day = VISITDY
  ))
  expect_no_error(derive_var_nfrlt(adpc,
    new_var = NFRLT, out_unit = "d",
    tpt_var = PCTPT, visit_day = VISITDY
  ))

  # Hours variations
  expect_no_error(derive_var_nfrlt(adpc,
    new_var = NFRLT, out_unit = "hour",
    tpt_var = PCTPT, visit_day = VISITDY
  ))
  expect_no_error(derive_var_nfrlt(adpc,
    new_var = NFRLT, out_unit = "hours",
    tpt_var = PCTPT, visit_day = VISITDY
  ))
  expect_no_error(derive_var_nfrlt(adpc,
    new_var = NFRLT, out_unit = "hr",
    tpt_var = PCTPT, visit_day = VISITDY
  ))
  expect_no_error(derive_var_nfrlt(adpc,
    new_var = NFRLT, out_unit = "hrs",
    tpt_var = PCTPT, visit_day = VISITDY
  ))
  expect_no_error(derive_var_nfrlt(adpc,
    new_var = NFRLT, out_unit = "h",
    tpt_var = PCTPT, visit_day = VISITDY
  ))

  # Minutes variations
  expect_no_error(derive_var_nfrlt(adpc,
    new_var = NFRLT, out_unit = "minute",
    tpt_var = PCTPT, visit_day = VISITDY
  ))
  expect_no_error(derive_var_nfrlt(adpc,
    new_var = NFRLT, out_unit = "minutes",
    tpt_var = PCTPT, visit_day = VISITDY
  ))
  expect_no_error(derive_var_nfrlt(adpc,
    new_var = NFRLT, out_unit = "min",
    tpt_var = PCTPT, visit_day = VISITDY
  ))
  expect_no_error(derive_var_nfrlt(adpc,
    new_var = NFRLT, out_unit = "mins",
    tpt_var = PCTPT, visit_day = VISITDY
  ))

  # Weeks variations
  expect_no_error(derive_var_nfrlt(adpc,
    new_var = NFRLT, out_unit = "week",
    tpt_var = PCTPT, visit_day = VISITDY
  ))
  expect_no_error(derive_var_nfrlt(adpc,
    new_var = NFRLT, out_unit = "weeks",
    tpt_var = PCTPT, visit_day = VISITDY
  ))
  expect_no_error(derive_var_nfrlt(adpc,
    new_var = NFRLT, out_unit = "wk",
    tpt_var = PCTPT, visit_day = VISITDY
  ))
  expect_no_error(derive_var_nfrlt(adpc,
    new_var = NFRLT, out_unit = "wks",
    tpt_var = PCTPT, visit_day = VISITDY
  ))
  expect_no_error(derive_var_nfrlt(adpc,
    new_var = NFRLT, out_unit = "w",
    tpt_var = PCTPT, visit_day = VISITDY
  ))

  # Case insensitive
  expect_no_error(derive_var_nfrlt(adpc,
    new_var = NFRLT, out_unit = "HOURS",
    tpt_var = PCTPT, visit_day = VISITDY
  ))
  expect_no_error(derive_var_nfrlt(adpc,
    new_var = NFRLT, out_unit = "Days",
    tpt_var = PCTPT, visit_day = VISITDY
  ))
})

## Test 50: unit variations produce correct results ----
test_that("derive_var_nfrlt Test 50: unit variations produce correct results", {
  adpc <- tribble(
    ~USUBJID, ~VISITDY, ~PCTPT,
    "001",    8,        "Pre-dose"
  )

  # All day variations should give same result
  result_day <- derive_var_nfrlt(adpc,
    new_var = NFRLT, out_unit = "day",
    tpt_var = PCTPT, visit_day = VISITDY
  )
  result_days <- derive_var_nfrlt(adpc,
    new_var = NFRLT, out_unit = "days",
    tpt_var = PCTPT, visit_day = VISITDY
  )
  result_d <- derive_var_nfrlt(adpc,
    new_var = NFRLT, out_unit = "d",
    tpt_var = PCTPT, visit_day = VISITDY
  )

  expect_equal(result_day$NFRLT, 7)
  expect_equal(result_days$NFRLT, 7)
  expect_equal(result_d$NFRLT, 7)

  # All hour variations should give same result
  result_hour <- derive_var_nfrlt(adpc,
    new_var = NFRLT, out_unit = "hour",
    tpt_var = PCTPT, visit_day = VISITDY
  )
  result_hr <- derive_var_nfrlt(adpc,
    new_var = NFRLT, out_unit = "hr",
    tpt_var = PCTPT, visit_day = VISITDY
  )
  result_h <- derive_var_nfrlt(adpc,
    new_var = NFRLT, out_unit = "h",
    tpt_var = PCTPT, visit_day = VISITDY
  )

  expect_equal(result_hour$NFRLT, 168)
  expect_equal(result_hr$NFRLT, 168)
  expect_equal(result_h$NFRLT, 168)
})
