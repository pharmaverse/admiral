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

## Test 4: returns expected values for time ranges with default midpoint ----
test_that("convert_xxtpt_to_hours Test 4: returns expected values for time ranges", {
  expect_equal(
    convert_xxtpt_to_hours(c(
      "0-6h",
      "6-12h Post-dose",
      "0.5 - 6.5h"
    )),
    c(3, 9, 3.5)
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
    regexp = "Argument `range_method` must be <character> with values"
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

## Test 18: time ranges with direction use range_method ----
test_that("convert_xxtpt_to_hours Test 18: ranges with direction use range_method", {
  # With default midpoint
  expect_equal(
    convert_xxtpt_to_hours(c(
      "0-4H PRIOR START OF INFUSION",
      "8-16H POST START OF INFUSION"
    )),
    c(-2, 12)
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
        "0-4H PRIOR START OF INFUSION",
        "8-16H POST START OF INFUSION"
      ),
      range_method = "start"
    ),
    c(0, 8)
  )
})

## Test 19: range_method parameter consistency ----
test_that("convert_xxtpt_to_hours Test 19: range_method works consistently", {
  input <- "2-8h Post-dose"

  expect_equal(
    convert_xxtpt_to_hours(input, range_method = "start"),
    2
  )

  expect_equal(
    convert_xxtpt_to_hours(input, range_method = "midpoint"),
    5
  )

  expect_equal(
    convert_xxtpt_to_hours(input, range_method = "end"),
    8
  )
})

## Test 20: does not warn for any patterns ----
test_that("convert_xxtpt_to_hours Test 20: does not warn for valid patterns", {
  expect_no_warning(
    convert_xxtpt_to_hours(c(
      "1 HOUR POST",
      "30 MIN POST EOI",
      "0-4H PRIOR START OF INFUSION"
    ))
  )
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

## Test 22: does not warn for time ranges ----
test_that("convert_xxtpt_to_hours Test 22: does not warn for time ranges", {
  expect_no_warning(
    convert_xxtpt_to_hours(c(
      "0-6h",
      "6-12h Post-dose",
      "0-4H PRIOR START OF INFUSION",
      "8-16H POST START OF INFUSION"
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

## Test 28: prioritizes specific patterns and handles ranges ----
test_that("convert_xxtpt_to_hours Test 28: prioritizes specific patterns correctly", {
  # Range with direction caught before simple range
  expect_equal(
    convert_xxtpt_to_hours("0-4H POST START OF INFUSION"),
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
    3,
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
      "0.5-2.5H PRIOR START OF INFUSION",
      range_method = "end"
    ),
    -2.5
  )
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

## Test 35: range_method affects all range patterns consistently ----
test_that("convert_xxtpt_to_hours Test 35: range_method affects all patterns consistently", {
  input <- c(
    "0-6h Post-dose",
    "0-4H PRIOR START OF INFUSION",
    "8-16H POST START OF INFUSION"
  )

  # Midpoint (default)
  expect_equal(
    convert_xxtpt_to_hours(input, range_method = "midpoint"),
    c(3, -2, 12)
  )

  # Start
  expect_equal(
    convert_xxtpt_to_hours(input, range_method = "start"),
    c(0, 0, 8)
  )

  # End
  expect_equal(
    convert_xxtpt_to_hours(input, range_method = "end"),
    c(6, -4, 16)
  )
})

## Test 36: range_method with infusion_duration ----
test_that("convert_xxtpt_to_hours Test 36: range_method and infusion_duration work together", {
  input <- c(
    "0-6h Post-dose",
    "EOI",
    "1 HOUR POST EOI"
  )

  # With range_method = "end" and infusion_duration = 2
  expect_equal(
    convert_xxtpt_to_hours(
      input,
      infusion_duration = 2,
      range_method = "end"
    ),
    c(6, 2, 3)
  )

  # With range_method = "start" and infusion_duration = 2
  expect_equal(
    convert_xxtpt_to_hours(
      input,
      infusion_duration = 2,
      range_method = "start"
    ),
    c(0, 2, 3)
  )
})
