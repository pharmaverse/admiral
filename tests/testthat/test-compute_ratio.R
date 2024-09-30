test_that("compute_ratio Test 1: Positive numbers", {
  expect_equal(
    compute_ratio(
      x = c(182, 50.25),
      y = c(91, 201)
    ),
    c(2, 0.25)
  )
})

test_that("compute_ratio Test 2: Negative numbers", {
  expect_equal(
    compute_ratio(
      x = c(-63.9, 100, -100),
      y = c(21.3, -40, -10)
    ),
    c(-3, -2.5, 10)
  )
})

test_that("compute_ratio Test 3: Zeros", {
  expect_equal(
    compute_ratio(
      x = c(0, 95, 0),
      y = c(75, 0, 0)
    ),
    c(0, NA_real_, NA_real_)
  )
})

test_that("compute_ratio Test 4: NAs", {
  expect_equal(
    compute_ratio(
      x = c(44, NA, NA),
      y = c(NA, 140, NA)
    ),
    c(NA_real_, NA_real_, NA_real_)
  )
})

test_that("compute_ratio Test 5: Error is thrown when `x` and `y` are of different length", {
  x_input <- c(110, 120, 107)
  y_input <- c(125, 80)

  expected_output <- c(0.88, 1.5)

  expect_snapshot(
    error = TRUE,
    compute_ratio(x_input, y_input)
  )
})
