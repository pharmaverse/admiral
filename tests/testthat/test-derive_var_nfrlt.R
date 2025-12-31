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
