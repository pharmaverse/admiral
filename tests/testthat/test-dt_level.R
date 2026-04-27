## Test 1: input is none; n ----
test_that("dt_level Test 1: input is none; n", {
  expected_n <- factor("n", levels = c("n", "D", "M", "Y"), ordered = TRUE)
  class(expected_n) <- c("dt_level", class(expected_n))

  expect_equal(dt_level("n"), expected_n)
})

## Test 2: input is day; D ----
test_that("dt_level Test 2: input is day; D", {
  expected_d <- factor("D", levels = c("n", "D", "M", "Y"), ordered = TRUE)
  class(expected_d) <- c("dt_level", class(expected_d))

  expect_equal(dt_level("D"), expected_d)
})

## Test 3: input is month; M ----
test_that("dt_level Test 3: input is month; M", {
  expected_m <- factor("M", levels = c("n", "D", "M", "Y"), ordered = TRUE)
  class(expected_m) <- c("dt_level", class(expected_m))

  expect_equal(dt_level("M"), expected_m)
})

## Test 4: input is year; Y ----
test_that("dt_level Test 4: input is year; Y", {
  expected_y <- factor("Y", levels = c("n", "D", "M", "Y"), ordered = TRUE)
  class(expected_y) <- c("dt_level", class(expected_y))

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
