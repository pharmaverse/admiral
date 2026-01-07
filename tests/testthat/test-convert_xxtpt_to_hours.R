# convert_xxtpt_to_hours ----

## Test 1: basic patterns ----
test_that("convert_xxtpt_to_hours Test 1: basic patterns work correctly", {
  expect_equal(
    convert_xxtpt_to_hours(c(
      "Screening",
      "Pre-dose",
      "Predose",
      "Pre-treatment",
      "Before",
      "30M",
      "1H",
      "2H POSTDOSE"
    )),
    c(0, 0, 0, 0, 0, 0.5, 1, 2)
  )
})

## Test 1b: treatment_duration parameter works ----
test_that("convert_xxtpt_to_hours Test 1b: treatment_duration parameter works", {
  expect_equal(
    convert_xxtpt_to_hours(c("EOI", "EOT", "End of Infusion"), treatment_duration = 2),
    c(2, 2, 2)
  )

  expect_equal(
    convert_xxtpt_to_hours(c("EOI", "End of Treatment"), treatment_duration = 0.5),
    c(0.5, 0.5)
  )

  expect_equal(
    convert_xxtpt_to_hours(c("EOI", "End of Infusion"), treatment_duration = 4),
    c(4, 4)
  )

  # Default is 0
  expect_equal(
    convert_xxtpt_to_hours(c("EOI", "EOT")),
    c(0, 0)
  )
})

## Test 1c: treatment_duration must be non-negative ----
test_that("convert_xxtpt_to_hours Test 1c: treatment_duration must be non-negative", {
  # Negative value should error
  expect_error(
    convert_xxtpt_to_hours("EOI", treatment_duration = -1),
    regexp = "must be non-negative"
  )

  # Zero is allowed
  expect_no_error(
    convert_xxtpt_to_hours("EOI", treatment_duration = 0)
  )

  # Must be numeric
  expect_error(
    convert_xxtpt_to_hours("EOI", treatment_duration = "1"),
    class = "assert_numeric_vector"
  )

  # Negative in vector should also error
  expect_error(
    convert_xxtpt_to_hours(c("EOI", "Pre-dose"), treatment_duration = c(1, -1)),
    regexp = "must be non-negative"
  )
})

## Test 2: returns expected values for days ----
test_that("convert_xxtpt_to_hours Test 2: returns expected values for days", {
  expect_equal(
    convert_xxtpt_to_hours(c(
      "Day 1",
      "2D",
      "2 days",
      "1.5 days",
      "30 DAYS AFTER LAST"
    )),
    c(24, 48, 48, 36, 720)
  )
})

## Test 3: returns expected values for hours+minutes combinations ----
test_that("convert_xxtpt_to_hours Test 3: returns expected values for hours+minutes", {
  expect_equal(
    convert_xxtpt_to_hours(c(
      "1H30M",
      "1 hour 30 min",
      "2HR15MIN",
      "1H30M POST",
      "1H30M AFTER"
    )),
    c(1.5, 1.5, 2.25, 1.5, 1.5)
  )
})

## Test 4: returns expected values for time ranges with default midpoint ----
test_that("convert_xxtpt_to_hours Test 4: returns expected values for time ranges", {
  expect_equal(
    convert_xxtpt_to_hours(c(
      "0-6h",
      "6-12h Post-dose",
      "0.5 - 6.5h",
      "6-12h After-dose"
    )),
    c(3, 9, 3.5, 9)
  )
})

## Test 4b: range_method parameter works ----
test_that("convert_xxtpt_to_hours Test 4b: range_method parameter works", {
  # Start method
  expect_equal(
    convert_xxtpt_to_hours("0-6h", range_method = "start"),
    0
  )

  # End method
  expect_equal(
    convert_xxtpt_to_hours("0-6h", range_method = "end"),
    6
  )

  # Midpoint method (default)
  expect_equal(
    convert_xxtpt_to_hours("0-6h", range_method = "midpoint"),
    3
  )

  # Multiple ranges
  expect_equal(
    convert_xxtpt_to_hours(
      c("0-6h", "6-12h", "12-24h"),
      range_method = "end"
    ),
    c(6, 12, 24)
  )
})

## Test 4c: range_method validation ----
test_that("convert_xxtpt_to_hours Test 4c: range_method must be valid", {
  expect_error(
    convert_xxtpt_to_hours("0-6h", range_method = "average"),
    regexp = "Argument `range_method` must be <character> with"
  )

  expect_error(
    convert_xxtpt_to_hours("0-6h", range_method = c("start", "end")),
    regexp = "EXPR must be a length 1 vector"
  )

  expect_error(
    convert_xxtpt_to_hours("0-6h", range_method = 123),
    class = "assert_character_vector"
  )
})

## Test 5: returns expected values for hours only ----
test_that("convert_xxtpt_to_hours Test 5: returns expected values for hours only", {
  expect_equal(
    convert_xxtpt_to_hours(c(
      "1H",
      "2 hours",
      "0.5HR",
      "4hr Post-dose",
      "1 HOUR POST",
      "4hr After-dose",
      "1 HOUR AFTER"
    )),
    c(1, 2, 0.5, 4, 1, 4, 1)
  )
})

## Test 6: returns expected values for minutes only ----
test_that("convert_xxtpt_to_hours Test 6: returns expected values for minutes only", {
  expect_equal(
    convert_xxtpt_to_hours(c(
      "30M",
      "45 min",
      "2.5 Min",
      "30 MIN POST",
      "30 MIN AFTER"
    )),
    c(0.5, 0.75, 2.5 / 60, 0.5, 0.5)
  )
})

## Test 7: returns expected values for predose/before patterns (negative time) ----
test_that("convert_xxtpt_to_hours Test 7: returns expected values for predose/before patterns", {
  expect_equal(
    convert_xxtpt_to_hours(c(
      "5 MIN PREDOSE",
      "1 HOUR PREDOSE",
      "5 MIN BEFORE",
      "1 HOUR BEFORE"
    )),
    c(-5 / 60, -1, -5 / 60, -1)
  )
})

## Test 8: returns expected values for post/after EOI/EOT patterns ----
test_that("convert_xxtpt_to_hours Test 8: returns expected values for post/after EOI/EOT", {
  # With default treatment_duration = 0
  expect_equal(
    convert_xxtpt_to_hours(c(
      "1 HOUR POST EOI",
      "30 MIN POST EOI",
      "1 HOUR AFTER EOT",
      "30 MIN AFTER EOT"
    )),
    c(1, 0.5, 1, 0.5)
  )

  # With treatment_duration = 1
  expect_equal(
    convert_xxtpt_to_hours(
      c(
        "1 HOUR POST EOI",
        "30 MIN POST EOI",
        "24 HR POST EOI",
        "24 HR POST EOT",
        "1 HOUR POST EOT",
        "1 HOUR AFTER EOI"
      ),
      treatment_duration = 1
    ),
    c(2, 1.5, 25, 25, 2, 2)
  )

  # With treatment_duration = 2
  expect_equal(
    convert_xxtpt_to_hours(
      c(
        "1 HOUR POST EOI",
        "30 MIN AFTER EOT"
      ),
      treatment_duration = 2
    ),
    c(3, 2.5)
  )
})

## Test 9: returns expected values for after end patterns ----
test_that("convert_xxtpt_to_hours Test 9: returns expected values for after end patterns", {
  # With default treatment_duration = 0
  expect_equal(
    convert_xxtpt_to_hours(c(
      "30MIN AFTER END OF INFUSION",
      "30MIN AFTER END OF TREATMENT"
    )),
    c(0.5, 0.5)
  )

  # With treatment_duration = 1
  expect_equal(
    convert_xxtpt_to_hours(
      c(
        "30MIN AFTER END OF INFUSION",
        "30MIN AFTER END OF TREATMENT"
      ),
      treatment_duration = 1
    ),
    c(1.5, 1.5)
  )

  # With treatment_duration = 2
  expect_equal(
    convert_xxtpt_to_hours(
      c(
        "30MIN AFTER END OF INFUSION",
        "1 HOUR AFTER END OF TREATMENT"
      ),
      treatment_duration = 2
    ),
    c(2.5, 3)
  )
})

## Test 10: returns expected values for HR POST INF patterns ----
test_that("convert_xxtpt_to_hours Test 10: returns expected values for HR POST INF", {
  # With default treatment_duration = 0
  expect_equal(
    convert_xxtpt_to_hours(c(
      "24 HR POST INF",
      "24 HR POST EOI",
      "24 HR POST EOT"
    )),
    c(24, 24, 24)
  )

  # With treatment_duration = 1
  expect_equal(
    convert_xxtpt_to_hours(
      c(
        "24 HR POST INF",
        "24 HR POST EOI",
        "24 HR POST EOT"
      ),
      treatment_duration = 1
    ),
    c(25, 25, 25)
  )

  # With treatment_duration = 0.5
  expect_equal(
    convert_xxtpt_to_hours(
      "24 HR POST INF",
      treatment_duration = 0.5
    ),
    24.5
  )
})

## Test 11: returns expected values for start of infusion/treatment patterns ----
test_that("convert_xxtpt_to_hours Test 11: returns expected values for start patterns", {
  expect_equal(
    convert_xxtpt_to_hours(c(
      "8H PRIOR START OF INFUSION",
      "8H POST START OF INFUSION",
      "8H BEFORE START OF TREATMENT",
      "8H AFTER START OF TREATMENT"
    )),
    c(-8, 8, -8, 8)
  )
})

## Test 12: returns expected values for MIN AFTER START INF ----
test_that("convert_xxtpt_to_hours Test 12: returns expected values for MIN AFTER START INF", {
  expect_equal(
    convert_xxtpt_to_hours("60 MIN AFTER START INF"),
    1
  )
})

## Test 13: returns expected values for MIN PRE/BEFORE EOI/EOT ----
test_that("convert_xxtpt_to_hours Test 13: returns expected values for MIN PRE/BEFORE EOI/EOT", {
  # With default treatment_duration = 0
  expect_equal(
    convert_xxtpt_to_hours(c(
      "10MIN PRE EOI",
      "10MIN BEFORE EOT"
    )),
    c(-10 / 60, -10 / 60)
  )

  # With treatment_duration = 2 (2 hours)
  # 10 minutes before end of 2-hour treatment = 2 - 10/60 = 1.833...
  expect_equal(
    convert_xxtpt_to_hours(
      c(
        "10MIN PRE EOI",
        "10MIN BEFORE EOT"
      ),
      treatment_duration = 2
    ),
    c(2 - 10 / 60, 2 - 10 / 60)
  )

  # With treatment_duration = 1
  expect_equal(
    convert_xxtpt_to_hours("10MIN PRE EOI", treatment_duration = 1),
    1 - 10 / 60
  )
})

## Test 14: handles NA values correctly ----
test_that("convert_xxtpt_to_hours Test 14: handles NA values correctly", {
  expect_equal(
    convert_xxtpt_to_hours(c("1H", NA, "2H")),
    c(1, NA_real_, 2)
  )
})

## Test 15: handles empty input ----
test_that("convert_xxtpt_to_hours Test 15: handles empty input", {
  expect_equal(
    convert_xxtpt_to_hours(character(0)),
    numeric(0)
  )
})

## Test 17: unrecognized patterns return NA ----
test_that("convert_xxtpt_to_hours Test 17: unrecognized patterns return NA", {
  result <- convert_xxtpt_to_hours(c(
    "Morning",
    "Evening",
    "UNKNOWN PATTERN"
  ))
  expect_true(all(is.na(result)))
  expect_equal(length(result), 3)
})

## Test 18: time ranges with direction use range_method ----
test_that("convert_xxtpt_to_hours Test 18: ranges with direction use range_method", {
  # With default midpoint
  expect_equal(
    convert_xxtpt_to_hours(c(
      "0-4H PRIOR START OF INFUSION",
      "8-16H POST START OF INFUSION",
      "0-4H BEFORE START OF TREATMENT",
      "8-16H AFTER START OF TREATMENT"
    )),
    c(-2, 12, -2, 12)
  )

  # With end method
  expect_equal(
    convert_xxtpt_to_hours(
      c(
        "0-4H PRIOR START OF INFUSION",
        "8-16H POST START OF INFUSION"
      ),
      range_method = "end"
    ),
    c(-4, 16)
  )

  # With start method
  expect_equal(
    convert_xxtpt_to_hours(
      c(
        "0-4H BEFORE START OF TREATMENT",
        "8-16H AFTER START OF TREATMENT"
      ),
      range_method = "start"
    ),
    c(0, 8)
  )
})

## Test 22: does not warn for time ranges ----
test_that("convert_xxtpt_to_hours Test 22: does not warn for time ranges", {
  expect_no_warning(
    convert_xxtpt_to_hours(c(
      "0-6h",
      "6-12h Post-dose",
      "0-4H PRIOR START OF INFUSION",
      "8-16H POST START OF INFUSION",
      "0-4H BEFORE START OF TREATMENT",
      "8-16H AFTER START OF TREATMENT"
    ))
  )
})

## Test 23: handles all patterns from original issue ----
test_that("convert_xxtpt_to_hours Test 23: handles all patterns with treatment_duration = 1", {
  expect_no_warning(
    result <- convert_xxtpt_to_hours(
      c(
        "60 MIN AFTER START INF",
        "PRE-INF",
        "5 MIN PREDOSE",
        "5 MIN BEFORE",
        "1 HOUR POST EOI",
        "1 HOUR POST",
        "1 HOUR AFTER EOT",
        "30 MIN POST EOI",
        "30 MIN POST",
        "30MIN AFTER END OF INFUSION",
        "30MIN AFTER END OF TREATMENT",
        "8H PRIOR START OF INFUSION",
        "8H POST START OF INFUSION",
        "8H BEFORE START OF TREATMENT",
        "30 DAYS AFTER LAST",
        "24 HR POST INF",
        "24 HR POST EOT",
        "24 HR POST EOI",
        "10MIN PRE EOI",
        "10MIN BEFORE EOT",
        "AFTER END OF INFUSION",
        "AFTER END OF TREATMENT",
        "Before"
      ),
      treatment_duration = 1
    )
  )

  expect_equal(
    result,
    c(
      1, # 60 MIN AFTER START INF
      0, # PRE-INF
      -5 / 60, # 5 MIN PREDOSE
      -5 / 60, # 5 MIN BEFORE
      2, # 1 HOUR POST EOI (1 + 1)
      1, # 1 HOUR POST
      2, # 1 HOUR AFTER EOT (1 + 1)
      1.5, # 30 MIN POST EOI (1 + 0.5)
      0.5, # 30 MIN POST
      1.5, # 30MIN AFTER END OF INFUSION (1 + 0.5)
      1.5, # 30MIN AFTER END OF TREATMENT (1 + 0.5)
      -8, # 8H PRIOR START OF INFUSION
      8, # 8H POST START OF INFUSION
      -8, # 8H BEFORE START OF TREATMENT
      720, # 30 DAYS AFTER LAST
      25, # 24 HR POST INF (1 + 24)
      25, # 24 HR POST EOT (1 + 24)
      25, # 24 HR POST EOI (1 + 24)
      1 - 10 / 60, # 10MIN PRE EOI (1 - 10/60)
      1 - 10 / 60, # 10MIN BEFORE EOT (1 - 10/60)
      1, # AFTER END OF INFUSION
      1, # AFTER END OF TREATMENT
      0 # Before
    )
  )
})

## Test 24: error if input is not character ----
test_that("convert_xxtpt_to_hours Test 24: error if input is not character", {
  expect_error(
    convert_xxtpt_to_hours(123),
    class = "assert_character_vector"
  )
})

## Test 25: handles whitespace variations ----
test_that("convert_xxtpt_to_hours Test 25: handles whitespace variations", {
  expect_equal(
    convert_xxtpt_to_hours(c(
      "1H",
      "1 H",
      "1  H",
      "1H ",
      " 1H"
    )),
    c(1, 1, 1, 1, 1)
  )
})

## Test 28: prioritizes specific patterns and handles ranges ----
test_that("convert_xxtpt_to_hours Test 28: prioritizes specific patterns correctly", {
  # Range with direction caught before simple range
  expect_equal(
    convert_xxtpt_to_hours("0-4H POST START OF INFUSION"),
    2
  )

  expect_equal(
    convert_xxtpt_to_hours("0-4H AFTER START OF TREATMENT"),
    2
  )

  # Simple range with midpoint
  expect_equal(
    convert_xxtpt_to_hours("0-6h"),
    3
  )

  # Ensure hours+minutes are caught before simple hours
  expect_equal(
    convert_xxtpt_to_hours("1H30M"),
    1.5
  )

  # Ensure predose/before is caught before simple minutes
  expect_equal(
    convert_xxtpt_to_hours("5 MIN PREDOSE"),
    -5 / 60
  )

  expect_equal(
    convert_xxtpt_to_hours("5 MIN BEFORE"),
    -5 / 60
  )
})

## Test 29: comprehensive integration test ----
test_that("convert_xxtpt_to_hours Test 29: comprehensive integration test", {
  input <- c(
    "Screening",
    "Pre-dose",
    "Pre-treatment",
    "PRE-INF",
    "Before",
    "30M",
    "1H",
    "1H30M",
    "2 hours Post-dose",
    "2 hours After-dose",
    "0-6h",
    "Day 1",
    "2D",
    "5 MIN PREDOSE",
    "5 MIN BEFORE",
    "1 HOUR POST EOI",
    "24 HR POST INF",
    "30 DAYS AFTER LAST",
    "AFTER END OF INFUSION",
    NA,
    "Morning"
  )

  # With default treatment_duration = 0
  expected_default <- c(
    0, # Screening
    0, # Pre-dose
    0, # Pre-treatment
    0, # PRE-INF
    0, # Before
    0.5, # 30M
    1, # 1H
    1.5, # 1H30M
    2, # 2 hours Post-dose
    2, # 2 hours After-dose
    3, # 0-6h (midpoint)
    24, # Day 1
    48, # 2D
    -5 / 60, # 5 MIN PREDOSE
    -5 / 60, # 5 MIN BEFORE
    1, # 1 HOUR POST EOI (0 + 1)
    24, # 24 HR POST INF (0 + 24)
    720, # 30 DAYS AFTER LAST
    0, # AFTER END OF INFUSION
    NA_real_, # NA
    NA_real_ # Morning
  )

  expect_no_warning(
    result <- convert_xxtpt_to_hours(input)
  )

  expect_equal(result, expected_default)

  # With treatment_duration = 1
  expected_treatment <- expected_default
  expected_treatment[16] <- 2 # 1 HOUR POST EOI (1 + 1)
  expected_treatment[17] <- 25 # 24 HR POST INF (1 + 24)
  expected_treatment[19] <- 1 # AFTER END OF INFUSION

  expect_no_warning(
    result_treatment <- convert_xxtpt_to_hours(input, treatment_duration = 1)
  )

  expect_equal(result_treatment, expected_treatment)
})

## Test 30: decimal ranges work correctly ----
test_that("convert_xxtpt_to_hours Test 30: handles decimal ranges", {
  expect_equal(
    convert_xxtpt_to_hours("0.5-6.5h", range_method = "midpoint"),
    3.5
  )

  expect_equal(
    convert_xxtpt_to_hours("0.5-6.5h", range_method = "start"),
    0.5
  )

  expect_equal(
    convert_xxtpt_to_hours("0.5-6.5h", range_method = "end"),
    6.5
  )
})

## Test 31: ranges with direction handle decimals ----
test_that("convert_xxtpt_to_hours Test 31: ranges with direction handle decimals", {
  expect_equal(
    convert_xxtpt_to_hours("0.5-2.5H PRIOR START OF INFUSION"),
    -1.5
  )

  expect_equal(
    convert_xxtpt_to_hours(
      "0.5-2.5H BEFORE START OF TREATMENT",
      range_method = "end"
    ),
    -2.5
  )
})

## Test 32: treatment_duration affects EOT and AFTER END patterns ----
test_that("convert_xxtpt_to_hours Test 32: treatment_duration in comprehensive context", {
  input <- c(
    "Pre-dose",
    "Pre-treatment",
    "PRE-INF",
    "Before",
    "EOI",
    "EOT",
    "End of Infusion",
    "End of Treatment",
    "AFTER END OF INFUSION",
    "AFTER END OF TREATMENT",
    "1 HOUR POST EOI",
    "1 HOUR AFTER EOT"
  )

  # With default treatment_duration = 0
  expect_equal(
    convert_xxtpt_to_hours(input),
    c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1)
  )

  # With treatment_duration = 2
  expect_equal(
    convert_xxtpt_to_hours(input, treatment_duration = 2),
    c(0, 0, 0, 0, 2, 2, 2, 2, 2, 2, 3, 3)
  )

  # With treatment_duration = 0.5
  expect_equal(
    convert_xxtpt_to_hours(input, treatment_duration = 0.5),
    c(0, 0, 0, 0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.5, 1.5)
  )
})

## Test 35: range_method affects all range patterns consistently ----
test_that("convert_xxtpt_to_hours Test 35: range_method affects all patterns consistently", {
  input <- c(
    "0-6h Post-dose",
    "0-4H PRIOR START OF INFUSION",
    "8-16H POST START OF INFUSION",
    "0-4H BEFORE START OF TREATMENT",
    "8-16H AFTER START OF TREATMENT"
  )

  # Midpoint (default)
  expect_equal(
    convert_xxtpt_to_hours(input, range_method = "midpoint"),
    c(3, -2, 12, -2, 12)
  )

  # Start
  expect_equal(
    convert_xxtpt_to_hours(input, range_method = "start"),
    c(0, 0, 8, 0, 8)
  )

  # End
  expect_equal(
    convert_xxtpt_to_hours(input, range_method = "end"),
    c(6, -4, 16, -4, 16)
  )
})

## Test 36: range_method with treatment_duration ----
test_that("convert_xxtpt_to_hours Test 36: range_method and treatment_duration work together", {
  input <- c(
    "0-6h Post-dose",
    "EOI",
    "1 HOUR POST EOI"
  )

  # With range_method = "end" and treatment_duration = 2
  expect_equal(
    convert_xxtpt_to_hours(
      input,
      treatment_duration = 2,
      range_method = "end"
    ),
    c(6, 2, 3)
  )

  # With range_method = "start" and treatment_duration = 2
  expect_equal(
    convert_xxtpt_to_hours(
      input,
      treatment_duration = 2,
      range_method = "start"
    ),
    c(0, 2, 3)
  )
})

## Test 37: EOI and EOT treated identically ----
test_that("convert_xxtpt_to_hours Test 37: EOI and EOT are equivalent", {
  eoi_patterns <- c(
    "EOI",
    "End of Infusion",
    "1 HOUR POST EOI",
    "30 MIN POST EOI",
    "10MIN PRE EOI",
    "AFTER END OF INFUSION"
  )

  eot_patterns <- c(
    "EOT",
    "End of Treatment",
    "1 HOUR POST EOT",
    "30 MIN POST EOT",
    "10MIN PRE EOT",
    "AFTER END OF TREATMENT"
  )

  # With treatment_duration = 0
  expect_equal(
    convert_xxtpt_to_hours(eoi_patterns, treatment_duration = 0),
    convert_xxtpt_to_hours(eot_patterns, treatment_duration = 0)
  )

  # With treatment_duration = 1
  expect_equal(
    convert_xxtpt_to_hours(eoi_patterns, treatment_duration = 1),
    convert_xxtpt_to_hours(eot_patterns, treatment_duration = 1)
  )

  # With treatment_duration = 2
  expect_equal(
    convert_xxtpt_to_hours(eoi_patterns, treatment_duration = 2),
    convert_xxtpt_to_hours(eot_patterns, treatment_duration = 2)
  )
})

## Test 38: MIN PRE EOI calculation with various treatment durations ----
test_that("convert_xxtpt_to_hours Test 38: MIN PRE EOI correctly uses treatment_duration", {
  # With 0 duration (oral): 10 minutes before end = -10 minutes
  expect_equal(
    convert_xxtpt_to_hours("10MIN PRE EOI", treatment_duration = 0),
    -10 / 60
  )

  # With 1 hour infusion: 10 minutes before end = 50 minutes from start
  expect_equal(
    convert_xxtpt_to_hours("10MIN PRE EOI", treatment_duration = 1),
    1 - 10 / 60
  )

  # With 2 hour infusion: 10 minutes before end = 110 minutes from start
  expect_equal(
    convert_xxtpt_to_hours("10MIN PRE EOI", treatment_duration = 2),
    2 - 10 / 60
  )

  # With 30 minute infusion: 10 minutes before end = 20 minutes from start
  expect_equal(
    convert_xxtpt_to_hours("10MIN PRE EOI", treatment_duration = 0.5),
    0.5 - 10 / 60
  )
})

## Test 39: POST INFUSION/INF patterns work ----
test_that("convert_xxtpt_to_hours Test 39: POST INFUSION/INF patterns work", {
  # With default treatment_duration = 0
  expect_equal(
    convert_xxtpt_to_hours(c(
      "30 MIN POST INFUSION",
      "30 MIN POST INF",
      "1 HOUR POST INFUSION",
      "2 HR POST INF"
    )),
    c(0.5, 0.5, 1, 2)
  )

  # With treatment_duration = 1
  expect_equal(
    convert_xxtpt_to_hours(
      c(
        "30 MIN POST INFUSION",
        "30 MIN POST INF"
      ),
      treatment_duration = 1
    ),
    c(1.5, 1.5)
  )

  # With treatment_duration = 2
  expect_equal(
    convert_xxtpt_to_hours(
      "1 HOUR POST INFUSION",
      treatment_duration = 2
    ),
    3
  )
})

## Test 40: vectorized treatment_duration works ----
test_that("convert_xxtpt_to_hours Test 40: vectorized treatment_duration works", {
  # Different durations for each timepoint
  expect_equal(
    convert_xxtpt_to_hours(
      c("EOI", "1 HOUR POST EOI", "EOI", "1 HOUR POST EOI"),
      treatment_duration = c(1, 1, 2, 2)
    ),
    c(1, 2, 2, 3)
  )

  # With PRE EOI patterns
  expect_equal(
    convert_xxtpt_to_hours(
      c("10MIN PRE EOI", "10MIN PRE EOI"),
      treatment_duration = c(1, 2)
    ),
    c(1 - 10 / 60, 2 - 10 / 60)
  )
})

## Test 41: treatment_duration vector length validation ----
test_that("convert_xxtpt_to_hours Test 41: treatment_duration vector length validated", {
  expect_error(
    convert_xxtpt_to_hours(
      c("EOI", "1 HOUR POST EOI", "EOI"),
      treatment_duration = c(1, 2) # Wrong length
    ),
    regexp = "must be either"
  )
})

## Test 42: treatment_duration with NA values ----
test_that("convert_xxtpt_to_hours Test 42: handles NA in treatment_duration", {
  result <- convert_xxtpt_to_hours(
    c("EOI", "1 HOUR POST EOI", "Pre-dose"),
    treatment_duration = c(1, NA, 1)
  )

  expect_equal(
    result,
    c(1, NA_real_, 0)
  )
})

## Test 43: ranges relative to EOI/EOT ----
test_that("convert_xxtpt_to_hours Test 43: ranges relative to EOI/EOT", {
  # With default treatment_duration = 0 and default midpoint
  expect_equal(
    convert_xxtpt_to_hours(c(
      "0-4H AFTER EOI",
      "0-4H EOI",
      "0-4H AFTER EOT",
      "0-4H EOT"
    )),
    c(2, 2, 2, 2) # midpoint of 0-4 is 2, plus 0 duration = 2
  )

  # With treatment_duration = 1
  expect_equal(
    convert_xxtpt_to_hours(
      c(
        "0-4H AFTER EOI",
        "0-4H EOI",
        "0-4H AFTER EOT",
        "0-4H EOT"
      ),
      treatment_duration = 1
    ),
    c(3, 3, 3, 3) # midpoint of 0-4 is 2, plus 1 duration = 3
  )

  # Long form: AFTER END OF INFUSION/TREATMENT
  expect_equal(
    convert_xxtpt_to_hours(
      c(
        "4-8H AFTER END OF INFUSION",
        "4-8H AFTER END OF TREATMENT",
        "4-8H AFTER EOI",
        "4-8H AFTER EOT"
      ),
      treatment_duration = 1
    ),
    c(7, 7, 7, 7) # midpoint of 4-8 is 6, plus 1 duration = 7
  )

  # With range_method = "start"
  expect_equal(
    convert_xxtpt_to_hours(
      c("0-4H AFTER EOI", "4-8H AFTER END OF INFUSION"),
      treatment_duration = 1,
      range_method = "start"
    ),
    c(1, 5) # start of 0-4 is 0, start of 4-8 is 4, plus 1 duration each
  )

  # With range_method = "end"
  expect_equal(
    convert_xxtpt_to_hours(
      c("0-4H AFTER EOI", "4-8H AFTER END OF TREATMENT"),
      treatment_duration = 1,
      range_method = "end"
    ),
    c(5, 9) # end of 0-4 is 4, end of 4-8 is 8, plus 1 duration each
  )

  # With vectorized treatment_duration
  expect_equal(
    convert_xxtpt_to_hours(
      c("0-4H EOI", "4-8H AFTER END OF INFUSION"),
      treatment_duration = c(1, 2)
    ),
    c(3, 8) # 1+2=3, 2+6=8
  )
})

## Test 44: POST vs POST EOI distinction ----
test_that("convert_xxtpt_to_hours Test 44: POST without EOI is relative to start", {
  # "POST" alone should NOT add treatment_duration (relative to start)
  expect_equal(
    convert_xxtpt_to_hours(
      c("1H POST", "30M AFTER"),
      treatment_duration = 2
    ),
    c(1, 0.5) # Just the time, no treatment_duration added
  )

  # "POST EOI" SHOULD add treatment_duration (relative to end)
  expect_equal(
    convert_xxtpt_to_hours(
      c("1H POST EOI", "30M AFTER EOT"),
      treatment_duration = 2
    ),
    c(3, 2.5)
  )

  # Same pattern, different durations to show the difference
  expect_equal(
    convert_xxtpt_to_hours(
      c("1H POST", "1H POST EOI"),
      treatment_duration = 2
    ),
    c(1, 3) # 1 vs (2+1)
  )

  # POST INFUSION should also add treatment_duration
  expect_equal(
    convert_xxtpt_to_hours(
      c("1H POST", "1H POST INFUSION"),
      treatment_duration = 2
    ),
    c(1, 3)
  )

  # Comprehensive comparison
  expect_equal(
    convert_xxtpt_to_hours(
      c(
        "Pre-dose",
        "1H POST",
        "2H POST",
        "EOI",
        "1H POST EOI",
        "2H POST EOI"
      ),
      treatment_duration = 2
    ),
    c(0, 1, 2, 2, 3, 4)
  )
})

## Test 45: bare numbers return NA (no unit ambiguity) ----
test_that("convert_xxtpt_to_hours: bare numbers without units return NA", {
  # Bare numbers without units should return NA (ambiguous)
  expect_equal(
    convert_xxtpt_to_hours(c("2", "1.5", "30")),
    c(NA_real_, NA_real_, NA_real_)
  )

  # With hour units, they work correctly
  expect_equal(
    convert_xxtpt_to_hours(c("2H", "1.5H", "30M")),
    c(2, 1.5, 0.5)
  )

  # Days patterns require at least one unit indicator
  expect_equal(
    convert_xxtpt_to_hours(c(
      "Day 1", # "Day" prefix
      "2D", # "D" suffix
      "2 days", # "days" suffix
      "30 DAYS AFTER LAST" # "DAYS AFTER LAST" context
    )),
    c(24, 48, 48, 720)
  )

  # More day variations
  expect_equal(
    convert_xxtpt_to_hours(c("Day 2", "3 day", "14D")),
    c(48, 72, 336)
  )
})
