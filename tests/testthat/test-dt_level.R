# set up to create expected values:
expected_levels <- c("n", "D", "M", "Y")
add_class <- function(object, class = "dt_level") {
  class(object) <- c("dt_level", class(object))
  object
}

## Test 1: input is none; n ----
test_that("dt_level Test 1: input is none; n", {
  expected_n <- factor("n", levels = expected_levels, ordered = TRUE) %>%
    add_class("dt_level")

  expect_equal(dt_level("n"), expected_n)
})

## Test 2: input is day; D ----
test_that("dt_level Test 2: input is day; D", {
  expected_d <- factor("D", levels = expected_levels, ordered = TRUE) %>%
    add_class("dt_level")

  expect_equal(dt_level("D"), expected_d)
})

## Test 3: input is month; M ----
test_that("dt_level Test 3: input is month; M", {
  expected_m <- factor("M", levels = expected_levels, ordered = TRUE) %>%
    add_class("dt_level")

  expect_equal(dt_level("M"), expected_m)
})

## Test 4: input is year; Y ----
test_that("dt_level Test 4: input is year; Y", {
  expected_y <- factor("Y", levels = expected_levels, ordered = TRUE) %>%
    add_class("dt_level")

  expect_equal(dt_level("Y"), expected_y)
})

## Test 5: input is not scalar ----
test_that("dt_level Test 5: input is not scalar", {
  expect_snapshot(dt_level(c("D", "M", "Y")),
    error = TRUE
  )
})

## Test 6: input is scalar character but not in expected set ----
test_that("dt_level Test 6: input is scalar character but not in expected set", {
  expect_snapshot(dt_level("d"),
    error = TRUE
  )
})
