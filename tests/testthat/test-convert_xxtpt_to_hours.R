# convert_xxtpt_to_hours ----

## Test 1: returns expected values for special cases ----
test_that("convert_xxtpt_to_hours Test 1: returns expected values for special cases", {
  expect_equal(
    convert_xxtpt_to_hours(c(
      "Screening",
      "SCREENING",
      "Pre-dose",
      "PREDOSE",
      "Pre-infusion",
      "Infusion",
      "0H",
      "0 H",
      "EOI",
      "End of Infusion",
      "Morning",
      "Evening"
    )),
    c(-1, -1, 0, 0, 0, 0, 0, 0, 1, 1, NA_real_, NA_real_)
  )
})

## Test 1b: infusion_duration parameter works ----
test_that("convert_xxtpt_to_hours Test 1b: infusion_duration parameter works", {
  expect_equal(
    convert_xxtpt_to_hours(c("EOI", "End of Infusion"), infusion_duration = 2),
    c(2, 2)
  )

  expect_equal(
    convert_xxtpt_to_hours(c("EOI", "End of Infusion"), infusion_duration = 0.5),
    c(0.5, 0.5)
  )

  expect_equal(
    convert_xxtpt_to_hours(c("EOI", "End of Infusion"), infusion_duration = 4),
    c(4, 4)
  )
})

## Test 1c: infusion_duration validation ----
test_that("convert_xxtpt_to_hours Test 1c: infusion_duration must be positive", {
  expect_error(
    convert_xxtpt_to_hours("EOI", infusion_duration = -1),
    regexp = "must be positive"
  )

  expect_error(
    convert_xxtpt_to_hours("EOI", infusion_duration = 0),
    regexp = "must be positive"
  )

  expect_error(
    convert_xxtpt_to_hours("EOI", infusion_duration = "1"),
    class = "assert_numeric_vector"
  )

  expect_error(
    convert_xxtpt_to_hours("EOI", infusion_duration = c(1, 2)),
    class = "assert_numeric_vector"
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
      "2HR15MIN"
    )),
    c(1.5, 1.5, 2.25)
  )
})

## Test 4: returns expected values for time ranges ----
test_that("convert_xxtpt_to_hours Test 4: returns expected values for time ranges", {
  expect_equal(
    convert_xxtpt_to_hours(c(
      "0-6h",
      "6-12h Post-dose",
      "0.5 - 6.5h"
    )),
    c(6, 12, 6.5)
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
      "1 HOUR POST"
    )),
    c(1, 2, 0.5, 4, 1)
  )
})

## Test 6: returns expected values for minutes only ----
test_that("convert_xxtpt_to_hours Test 6: returns expected values for minutes only", {
  expect_equal(
    convert_xxtpt_to_hours(c(
      "30M",
      "45 min",
      "2.5 Min",
      "30 MIN POST"
    )),
    c(0.5, 0.75, 2.5 / 60, 0.5)
  )
})

## Test 7: returns expected values for predose patterns (negative time) ----
test_that("convert_xxtpt_to_hours Test 7: returns expected values for predose patterns", {
  expect_equal(
    convert_xxtpt_to_hours(c(
      "5 MIN PREDOSE",
      "1 HOUR PREDOSE"
    )),
    c(-5 / 60, -1)
  )
})

## Test 8: returns expected values for post EOI/EOT patterns ----
test_that("convert_xxtpt_to_hours Test 8: returns expected values for post EOI/EOT", {
  expect_equal(
    convert_xxtpt_to_hours(c(
      "1 HOUR POST EOI",
      "30 MIN POST EOI",
      "24 HR POST EOI",
      "24 HR POST EOT",
      "1 HOUR POST EOT"
    )),
    c(1, 0.5, 24, 24, 1)
  )
})

## Test 9: returns expected values for after end patterns ----
test_that("convert_xxtpt_to_hours Test 9: returns expected values for after end patterns", {
  expect_equal(
    convert_xxtpt_to_hours(c(
      "30MIN AFTER END OF INFUSION",
      "30MIN AFTER END OF TREATMENT"
    )),
    c(0.5, 0.5)
  )
})

## Test 10: returns expected values for HR POST INF patterns ----
test_that("convert_xxtpt_to_hours Test 10: returns expected values for HR POST INF", {
  expect_equal(
    convert_xxtpt_to_hours(c(
      "24 HR POST INF",
      "24 HR POST EOI",
      "24 HR POST EOT"
    )),
    c(24, 24, 24)
  )
})

## Test 11: returns expected values for start of infusion patterns ----
test_that("convert_xxtpt_to_hours Test 11: returns expected values for start of infusion", {
  expect_equal(
    convert_xxtpt_to_hours(c(
      "8H PRIOR START OF INFUSION",
      "8H POST START OF INFUSION"
    )),
    c(-8, 8)
  )
})

## Test 12: returns expected values for MIN AFTER START INF ----
test_that("convert_xxtpt_to_hours Test 12: returns expected values for MIN AFTER START INF", {
  expect_equal(
    convert_xxtpt_to_hours("60 MIN AFTER START INF"),
    1
  )
})

## Test 13: returns expected values for MIN PRE EOI ----
test_that("convert_xxtpt_to_hours Test 13: returns expected values for MIN PRE EOI", {
  expect_equal(
    convert_xxtpt_to_hours("10MIN PRE EOI"),
    -10 / 60
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

## Test 16: is case-insensitive ----
test_that("convert_xxtpt_to_hours Test 16: is case-insensitive", {
  expect_equal(
    convert_xxtpt_to_hours(c(
      "1h post-dose",
      "1H POST-DOSE",
      "1H Post-Dose"
    )),
    c(1, 1, 1)
  )
})

## Test 17: warns for vague patterns and returns NA ----
test_that("convert_xxtpt_to_hours Test 17: warns for vague patterns and returns NA", {
  expect_warning(
    result <- convert_xxtpt_to_hours(c(
      "PRE-INF",
      "AFTER END OF INFUSION"
    )),
    regexp = "Vague timepoint descriptions cannot be converted to numeric hours"
  )
  expect_true(all(is.na(result)))
  expect_equal(length(result), 2)
})

## Test 18: warns for range patterns with direction and returns NA ----
test_that("convert_xxtpt_to_hours Test 18: warns for range patterns and returns NA", {
  expect_warning(
    result <- convert_xxtpt_to_hours(c(
      "0-4H PRIOR START OF INFUSION",
      "8-16H POST START OF INFUSION"
    )),
    regexp = "Time range patterns with direction cannot be converted to single numeric hours"
  )
  expect_true(all(is.na(result)))
  expect_equal(length(result), 2)
})

## Test 19: warns only once for multiple instances of same vague pattern ----
test_that("convert_xxtpt_to_hours Test 19: warns only once for duplicate vague patterns", {
  # unique() in the function should deduplicate warning values
  expect_warning(
    result <- convert_xxtpt_to_hours(c(
      "PRE-INF",
      "PRE-INF",
      "PRE-INF"
    )),
    regexp = "Vague timepoint descriptions cannot be converted to numeric hours"
  )
  expect_true(all(is.na(result)))
})

## Test 20: warns for both vague and range patterns when both present ----
test_that("convert_xxtpt_to_hours Test 20: warns for both vague and range patterns", {
  # Should get two warnings
  expect_warning(
    expect_warning(
      result <- convert_xxtpt_to_hours(c(
        "1 HOUR POST",
        "PRE-INF",
        "30 MIN POST EOI",
        "0-4H PRIOR START OF INFUSION"
      )),
      regexp = "Vague timepoint"
    ),
    regexp = "Time range"
  )
  expect_equal(result[c(1, 3)], c(1, 0.5))
  expect_true(all(is.na(result[c(2, 4)])))
})

## Test 21: does not warn for valid patterns ----
test_that("convert_xxtpt_to_hours Test 21: does not warn for valid patterns", {
  expect_no_warning(
    convert_xxtpt_to_hours(c(
      "1 HOUR POST",
      "30 MIN POST EOI",
      "5 MIN PREDOSE",
      "24 HR POST INF"
    ))
  )
})

## Test 22: does not warn for simple time ranges without direction ----
test_that("convert_xxtpt_to_hours Test 22: does not warn for simple time ranges", {
  expect_no_warning(
    convert_xxtpt_to_hours(c(
      "0-6h",
      "6-12h Post-dose"
    ))
  )
})

## Test 23: handles all patterns from original issue ----
test_that("convert_xxtpt_to_hours Test 23: handles all patterns from original issue", {
  # Valid conversions (no warning expected for these)
  expect_no_warning(
    result <- convert_xxtpt_to_hours(c(
      "60 MIN AFTER START INF",
      "5 MIN PREDOSE",
      "1 HOUR POST EOI",
      "1 HOUR POST",
      "30 MIN POST EOI",
      "30 MIN POST",
      "30MIN AFTER END OF INFUSION",
      "30MIN AFTER END OF TREATMENT",
      "8H PRIOR START OF INFUSION",
      "8H POST START OF INFUSION",
      "30 DAYS AFTER LAST",
      "24 HR POST INF",
      "24 HR POST EOT",
      "24 HR POST EOI",
      "10MIN PRE EOI"
    ))
  )

  expect_equal(
    result,
    c(
      1, # 60 MIN AFTER START INF
      -5 / 60, # 5 MIN PREDOSE
      1, # 1 HOUR POST EOI
      1, # 1 HOUR POST
      0.5, # 30 MIN POST EOI
      0.5, # 30 MIN POST
      0.5, # 30MIN AFTER END OF INFUSION
      0.5, # 30MIN AFTER END OF TREATMENT
      -8, # 8H PRIOR START OF INFUSION
      8, # 8H POST START OF INFUSION
      720, # 30 DAYS AFTER LAST
      24, # 24 HR POST INF
      24, # 24 HR POST EOT
      24, # 24 HR POST EOI
      -10 / 60 # 10MIN PRE EOI
    )
  )

  # Vague patterns (expect warning)
  expect_warning(
    result_vague <- convert_xxtpt_to_hours(c(
      "PRE-INF",
      "AFTER END OF INFUSION"
    )),
    regexp = "Vague timepoint"
  )
  expect_true(all(is.na(result_vague)))

  # Range patterns (expect warning)
  expect_warning(
    result_range <- convert_xxtpt_to_hours(c(
      "0-4H PRIOR START OF INFUSION",
      "8-16H POST START OF INFUSION"
    )),
    regexp = "Time range"
  )
  expect_true(all(is.na(result_range)))
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

## Test 26: handles decimal values ----
test_that("convert_xxtpt_to_hours Test 26: handles decimal values", {
  expect_equal(
    convert_xxtpt_to_hours(c(
      "0.5H",
      "1.5 HOURS",
      "2.25 hours Post-dose"
    )),
    c(0.5, 1.5, 2.25)
  )
})

## Test 27: handles post-dose variations ----
test_that("convert_xxtpt_to_hours Test 27: handles post-dose variations", {
  expect_equal(
    convert_xxtpt_to_hours(c(
      "1H Post-dose",
      "1H POST-DOSE",
      "1H POSTDOSE",
      "1H Post dose"
    )),
    c(1, 1, 1, 1)
  )
})

## Test 28: prioritizes specific patterns over general ones ----
test_that("convert_xxtpt_to_hours Test 28: prioritizes specific patterns correctly", {
  # Ensure time ranges are caught before simple hours
  expect_equal(
    convert_xxtpt_to_hours("0-6h"),
    6
  )

  # Ensure hours+minutes are caught before simple hours
  expect_equal(
    convert_xxtpt_to_hours("1H30M"),
    1.5
  )

  # Ensure predose is caught before simple minutes
  expect_equal(
    convert_xxtpt_to_hours("5 MIN PREDOSE"),
    -5 / 60
  )
})

## Test 29: comprehensive integration test ----
test_that("convert_xxtpt_to_hours Test 29: comprehensive integration test", {
  input <- c(
    "Screening",
    "Pre-dose",
    "30M",
    "1H",
    "1H30M",
    "2 hours Post-dose",
    "0-6h",
    "Day 1",
    "2D",
    "5 MIN PREDOSE",
    "1 HOUR POST EOI",
    "24 HR POST INF",
    "30 DAYS AFTER LAST",
    NA,
    "Morning"
  )

  expected <- c(
    -1,
    0,
    0.5,
    1,
    1.5,
    2,
    6,
    24,
    48,
    -5 / 60,
    1,
    24,
    720,
    NA_real_,
    NA_real_
  )

  expect_no_warning(
    result <- convert_xxtpt_to_hours(input)
  )

  expect_equal(result, expected)
})

## Test 30: warns with specific vague values in message ----
test_that("convert_xxtpt_to_hours Test 30: warning includes actual vague values", {
  expect_warning(
    convert_xxtpt_to_hours(c("PRE-INF", "AFTER END OF INFUSION")),
    regexp = "PRE-INF|AFTER END OF INFUSION"
  )
})

## Test 31: warns with specific range values in message ----
test_that("convert_xxtpt_to_hours Test 31: warning includes actual range values", {
  expect_warning(
    convert_xxtpt_to_hours(c("0-4H PRIOR START OF INFUSION")),
    regexp = "0-4H PRIOR START OF INFUSION"
  )
})

## Test 32: multiple vague patterns show unique values only ----
test_that("convert_xxtpt_to_hours Test 32: deduplicates vague values in warning", {
  # This tests that unique() is working in the warning collection
  expect_warning(
    result <- convert_xxtpt_to_hours(c(
      "PRE-INF",
      "PRE-INF",
      "AFTER END OF INFUSION",
      "AFTER END OF INFUSION"
    )),
    regexp = "Vague timepoint"
  )
  expect_equal(length(result), 4)
  expect_true(all(is.na(result)))
})

## Test 33: infusion_duration affects EOI in integration ----
test_that("convert_xxtpt_to_hours Test 33: infusion_duration in comprehensive context", {
  input <- c(
    "Pre-dose",
    "EOI",
    "End of Infusion",
    "1 HOUR POST EOI"
  )

  # With default infusion_duration = 1
  expect_equal(
    convert_xxtpt_to_hours(input),
    c(0, 1, 1, 1)
  )

  # With infusion_duration = 2
  expect_equal(
    convert_xxtpt_to_hours(input, infusion_duration = 2),
    c(0, 2, 2, 1)
  )

  # With infusion_duration = 0.5
  expect_equal(
    convert_xxtpt_to_hours(input, infusion_duration = 0.5),
    c(0, 0.5, 0.5, 1)
  )
})
