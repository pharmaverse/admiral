# convert_xxtpt_to_hours ----
## Test 1: Special case values ----
test_that("convert_xxtpt_to_hours Test 1: Special case values", {
  expect_equal(
    convert_xxtpt_to_hours(c(
      "Screening", "SCREENING",
      "Pre-dose", "PREDOSE", "Pre-infusion", "Infusion", "0H", "0 h",
      "EOI", "End of Infusion",
      "Morning", "Evening"
    )),
    c(-1, -1, 0, 0, 0, 0, 0, 0, 1, 1, NA_real_, NA_real_)
  )
})
## Test 2: Days conversion ----
test_that("convert_xxtpt_to_hours Test 2: Days convert to hours (x24)", {
  expect_equal(
    convert_xxtpt_to_hours(c("Day 1", "2D", "1.5 days", "Day 7")),
    c(24, 48, 36, 168)
  )
})
## Test 3: Hours and minutes combinations ----
test_that("convert_xxtpt_to_hours Test 3: Hours+Minutes combinations", {
  expect_equal(
    convert_xxtpt_to_hours(c("1H30M", "1 hour 30 min", "2HR15MIN")),
    c(1.5, 1.5, 2.25)
  )
})
## Test 4: Time ranges return end value ----
test_that("convert_xxtpt_to_hours Test 4: Time ranges return end value", {
  expect_equal(
    convert_xxtpt_to_hours(c("0-6h", "6-12h Post-dose", "0.5 - 6.5h")),
    c(6, 12, 6.5)
  )
})
## Test 5: Hours only - various formats ----
test_that("convert_xxtpt_to_hours Test 5: Hours with format variations", {
  expect_equal(
    convert_xxtpt_to_hours(c(
      "1h", "2H", "4hr Post-dose", "8 hour", "12 HOURS", "0.5h", "24h"
    )),
    c(1, 2, 4, 8, 12, 0.5, 24)
  )
})
## Test 6: Minutes only - various formats ----
test_that("convert_xxtpt_to_hours Test 6: Minutes with format variations", {
  expect_equal(
    convert_xxtpt_to_hours(c(
      "5M", "15m Post-dose", "30 min", "45 MINUTES", "2.5 Min"
    )),
    c(5 / 60, 15 / 60, 30 / 60, 45 / 60, 2.5 / 60)
  )
})
## Test 7: Case insensitivity and spacing ----
test_that("convert_xxtpt_to_hours Test 7: Case and spacing variations", {
  expect_equal(
    convert_xxtpt_to_hours(c(
      "1h POST-DOSE", "1H post-dose", "1 h Post-Dose",
      "PRE-DOSE", "Predose", "PREDOSE"
    )),
    c(1, 1, 1, 0, 0, 0)
  )
})
## Test 8: NA and edge cases ----
test_that("convert_xxtpt_to_hours Test 8: NA and edge cases", {
  expect_equal(
    convert_xxtpt_to_hours(c(NA_character_, "Unknown", "EOS", "Pre-dose")),
    c(NA_real_, NA_real_, NA_real_, 0)
  )
  expect_equal(convert_xxtpt_to_hours(character(0)), numeric(0))
})
## Test 9: Comprehensive mixed input ----
test_that("convert_xxtpt_to_hours Test 9: Mixed comprehensive input", {
  input <- c(
    "Screening", "Pre-dose", "5 Min Post-dose", "30M", "1H30M",
    "0.5h", "1h Post-dose", "2h", "0-6h Post-dose", "Day 1",
    "24h", "EOI", "Morning", "Unknown"
  )
  expected <- c(
    -1, 0, 5 / 60, 30 / 60, 1.5,
    0.5, 1, 2, 6, 24,
    24, 1, NA_real_, NA_real_
  )
  expect_equal(convert_xxtpt_to_hours(input), expected)
})
