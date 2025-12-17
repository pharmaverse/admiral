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
      "PRE-INF",
      "Infusion",
      "0H",
      "0 H",
      "EOI",
      "End of Infusion",
      "AFTER END OF INFUSION",
      "AFTER END OF TREATMENT",
      "Morning",
      "Evening"
    )),
    c(-1, -1, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, NA_real_, NA_real_)
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
  # With default infusion_duration = 1
  expect_equal(
    convert_xxtpt_to_hours(c(
      "1 HOUR POST EOI",
      "30 MIN POST EOI",
      "24 HR POST EOI",
      "24 HR POST EOT",
      "1 HOUR POST EOT"
    )),
    c(2, 1.5, 25, 25, 2)
  )

  # With infusion_duration = 2
  expect_equal(
    convert_xxtpt_to_hours(
      c(
        "1 HOUR POST EOI",
        "30 MIN POST EOI"
      ),
      infusion_duration = 2
    ),
    c(3, 2.5)
  )
})

## Test 9: returns expected values for after end patterns ----
test_that("convert_xxtpt_to_hours Test 9: returns expected values for after end patterns", {
  # With default infusion_duration = 1
  expect_equal(
    convert_xxtpt_to_hours(c(
      "30MIN AFTER END OF INFUSION",
      "30MIN AFTER END OF TREATMENT"
    )),
    c(1.5, 1.5)
  )

  # With infusion_duration = 2
  expect_equal(
    convert_xxtpt_to_hours(
      c(
        "30MIN AFTER END OF INFUSION",
        "1 HOUR AFTER END OF TREATMENT"
      ),
      infusion_duration = 2
    ),
    c(2.5, 3)
  )
})

## Test 10: returns expected values for HR POST INF patterns ----
test_that("convert_xxtpt_to_hours Test 10: returns expected values for HR POST INF", {
  # With default infusion_duration = 1
  expect_equal(
    convert_xxtpt_to_hours(c(
      "24 HR POST INF",
      "24 HR POST EOI",
      "24 HR POST EOT"
    )),
    c(25, 25, 25)
  )

  # With infusion_duration = 0.5
  expect_equal(
    convert_xxtpt_to_hours(
      "24 HR POST INF",
      infusion_duration = 0.5
    ),
    24.5
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

## Test 19: warns only once for multiple instances of same range pattern ----
test_that("convert_xxtpt_to_hours Test 19: warns only once for duplicate range patterns", {
  expect_warning(
    result <- convert_xxtpt_to_hours(c(
      "0-4H PRIOR START OF INFUSION",
      "0-4H PRIOR START OF INFUSION",
      "0-4H PRIOR START OF INFUSION"
    )),
    regexp = "Time range patterns with direction cannot be converted to single numeric hours"
  )
  expect_true(all(is.na(result)))
})

## Test 20: warns for range patterns when present ----
test_that("convert_xxtpt_to_hours Test 20: warns for range patterns", {
  expect_warning(
    result <- convert_xxtpt_to_hours(c(
      "1 HOUR POST",
      "30 MIN POST EOI",
      "0-4H PRIOR START OF INFUSION"
    )),
    regexp = "Time range"
  )
  expect_equal(result[c(1, 2)], c(1, 1.5))
  expect_true(is.na(result[3]))
})

## Test 21: does not warn for valid patterns ----
test_that("convert_xxtpt_to_hours Test 21: does not warn for valid patterns", {
  expect_no_warning(
    convert_xxtpt_to_hours(c(
      "1 HOUR POST",
      "30 MIN POST EOI",
      "5 MIN PREDOSE",
      "24 HR POST INF",
      "PRE-INF",
      "AFTER END OF INFUSION"
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
  expect_no_warning(
    result <- convert_xxtpt_to_hours(c(
      "60 MIN AFTER START INF",
      "PRE-INF",
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
      "10MIN PRE EOI",
      "AFTER END OF INFUSION",
      "AFTER END OF TREATMENT"
    ))
  )

  expect_equal(
    result,
    c(
      1,
      0,
      -5 / 60,
      2,
      1,
      1.5,
      0.5,
      1.5,
      1.5,
      -8,
      8,
      720,
      25,
      25,
      25,
      -10 / 60,
      1,
      1
    )
  )

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
  expect_equal(
    convert_xxtpt_to_hours("0-6h"),
    6
  )

  expect_equal(
    convert_xxtpt_to_hours("1H30M"),
    1.5
  )

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
    "PRE-INF",
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
    "AFTER END OF INFUSION",
    NA,
    "Morning"
  )

  expected <- c(
    -1,
    0,
    0,
    0.5,
    1,
    1.5,
    2,
    6,
    24,
    48,
    -5 / 60,
    2,
    25,
    720,
    1,
    NA_real_,
    NA_real_
  )

  expect_no_warning(
    result <- convert_xxtpt_to_hours(input)
  )

  expect_equal(result, expected)
})

## Test 30: warns with specific range values in message ----
test_that("convert_xxtpt_to_hours Test 30: warning includes actual range values", {
  expect_warning(
    convert_xxtpt_to_hours(c("0-4H PRIOR START OF INFUSION")),
    regexp = "0-4H PRIOR START OF INFUSION"
  )
})

## Test 31: multiple range patterns show unique values only ----
test_that("convert_xxtpt_to_hours Test 31: deduplicates range values in warning", {
  expect_warning(
    result <- convert_xxtpt_to_hours(c(
      "0-4H PRIOR START OF INFUSION",
      "0-4H PRIOR START OF INFUSION",
      "8-16H POST START OF INFUSION",
      "8-16H POST START OF INFUSION"
    )),
    regexp = "Time range"
  )
  expect_equal(length(result), 4)
  expect_true(all(is.na(result)))
})

## Test 32: infusion_duration affects EOI and AFTER END patterns ----
test_that("convert_xxtpt_to_hours Test 32: infusion_duration in comprehensive context", {
  input <- c(
    "Pre-dose",
    "PRE-INF",
    "EOI",
    "End of Infusion",
    "AFTER END OF INFUSION",
    "AFTER END OF TREATMENT",
    "1 HOUR POST EOI"
  )

  # With default infusion_duration = 1
  expect_equal(
    convert_xxtpt_to_hours(input),
    c(0, 0, 1, 1, 1, 1, 2)
  )

  # With infusion_duration = 2
  expect_equal(
    convert_xxtpt_to_hours(input, infusion_duration = 2),
    c(0, 0, 2, 2, 2, 2, 3)
  )

  # With infusion_duration = 0.5
  expect_equal(
    convert_xxtpt_to_hours(input, infusion_duration = 0.5),
    c(0, 0, 0.5, 0.5, 0.5, 0.5, 1.5)
  )
})

## Test 33: PRE-INF treated as 0 like other pre-patterns ----
test_that("convert_xxtpt_to_hours Test 33: PRE-INF returns 0", {
  expect_equal(
    convert_xxtpt_to_hours(c(
      "PRE-INF",
      "Pre-inf",
      "pre-infusion",
      "PREDOSE"
    )),
    c(0, 0, 0, 0)
  )
})

## Test 34: AFTER END patterns equal infusion_duration ----
test_that("convert_xxtpt_to_hours Test 34: AFTER END patterns return infusion_duration", {
  # With default infusion_duration = 1
  expect_equal(
    convert_xxtpt_to_hours(c(
      "AFTER END OF INFUSION",
      "AFTER END OF TREATMENT",
      "EOI"
    )),
    c(1, 1, 1)
  )

  # With custom infusion_duration = 3
  expect_equal(
    convert_xxtpt_to_hours(
      c(
        "AFTER END OF INFUSION",
        "AFTER END OF TREATMENT"
      ),
      infusion_duration = 3
    ),
    c(3, 3)
  )
})
