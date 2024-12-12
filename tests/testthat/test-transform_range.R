## Test 1: works as expected ----
test_that("transform_range Test 1: works as expected", {
  expect_equal(
    transform_range(
      c(5, 1, 6, 2, NA),
      source_range = c(1, 5),
      target_range = c(0, 100)
    ),
    c(100, 0, NA, 25, NA)
  )
})

## Test 2: range is flipped if flip_direction == TRUE ----
test_that("transform_range Test 2: range is flipped if flip_direction == TRUE", {
  expect_equal(
    transform_range(
      c(0, 4, 8, 11),
      c(0, 10),
      c(0, 100),
      flip_direction = TRUE
    ),
    c(100, 60, 20, NA)
  )
})

## Test 3: warning if outside range ----
test_that("transform_range Test 3: warning if outside range", {
  expect_snapshot(
    transform_range(
      c(5, 1, 6, 2, NA),
      source_range = c(1, 5),
      target_range = c(0, 100),
      outside_range = "warning"
    )
  )
})

## Test 4: error if outside range ----
test_that("transform_range Test 4: error if outside range", {
  expect_snapshot(
    transform_range(
      c(5, 1, 6, 2, 7),
      source_range = c(1, 5),
      target_range = c(0, 100),
      outside_range = "error"
    ),
    error = TRUE
  )
})
