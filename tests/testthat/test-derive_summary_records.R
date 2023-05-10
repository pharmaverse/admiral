## Test 1: creates new record per group and groups are retained ----
test_that("derive_summary_records Test 1: creates new record per group and groups are retained", {
  # group --> 4
  input <- tibble::tibble(x = rep(1:4, each = 4), y = rep(1:2, each = 8), z = runif(16))
  actual_output <- input %>%
    derive_summary_records(
      by_vars = exprs(x, y),
      analysis_var = z,
      summary_fun = mean
    )

  expect_equal(nrow(actual_output), nrow(input) + 4)
  expect_equal(dplyr::group_vars(actual_output), dplyr::group_vars(input))
})

## Test 2: `fns` as inlined ----
test_that("derive_summary_records Test 2: `fns` as inlined", {
  input <- tibble::tibble(x = rep(1:2, each = 2), y = 9:12, z = 101:104)
  actual_output <- derive_summary_records(
    input,
    by_vars = exprs(x),
    analysis_var = y,
    summary_fun = function(x) mean(x, na.rm = TRUE)
  )
  expected_output <- tibble::tibble(
    x = rep(1:2, each = 3),
    y = c(9:10, 9.5, 11:12, 11.5),
    z = c(101:102, NA, 103:104, NA)
  )

  expect_dfs_equal(actual_output, expected_output, keys = c("x", "y", "z"))
})

## Test 3: set new value to a derived record ----
test_that("derive_summary_records Test 3: set new value to a derived record", {
  input <- tibble::tibble(x = rep(1:2, each = 2), y = 9:12)
  actual_output <- derive_summary_records(
    input,
    by_vars = exprs(x),
    analysis_var = y,
    summary_fun = mean,
    set_values_to = exprs(z = "MEAN")
  )
  expected_output <- tibble::tibble(
    x = rep(1:2, each = 3),
    y = c(9:10, 9.5, 11:12, 11.5),
    z = c(NA, NA, "MEAN", NA, NA, "MEAN")
  )

  expect_dfs_equal(actual_output, expected_output, keys = c("x", "y", "z"))
})

## Test 4: check `set_values_to` mapping ----
test_that("derive_summary_records Test 4: check `set_values_to` mapping", {
  input <- tibble::tibble(x = rep(1:4, each = 4), y = rep(1:2, each = 8), z = runif(16))
  actual_output <- input %>%
    derive_summary_records(
      by_vars = exprs(x, y),
      analysis_var = z,
      summary_fun = mean,
      set_values_to = exprs(d = "MEAN")
    ) %>%
    derive_summary_records(
      by_vars = exprs(x, y),
      analysis_var = z,
      summary_fun = sum,
      set_values_to = exprs(d = "SUM")
    )
  tf <- rep(c(NA, "MEAN", "SUM"), c(16, 4, 4))

  expect_equal(actual_output$d, tf)

  actual_output <- input %>%
    derive_summary_records(
      by_vars = exprs(x, y),
      analysis_var = z,
      summary_fun = mean,
      set_values_to = exprs(d = "MEAN", p1 = "PARAM1", p2 = "PARAM2")
    )
  tf <- rep(c(NA, "MEAN"), c(16, 4))
  tp1 <- rep(c(NA, "PARAM1"), c(16, 4))
  tp2 <- rep(c(NA, "PARAM2"), c(16, 4))

  expect_equal(actual_output$d, tf)
  expect_equal(actual_output$p1, tp1)
  expect_equal(actual_output$p2, tp2)
})

## Test 5: Filter record within `by_vars` ----
test_that("derive_summary_records Test 5: Filter record within `by_vars`", {
  input <- tibble::tibble(x = c(rep(1:2, each = 2), 2), y = 9:13, z = c(1, 1, 2, 1, 1))

  actual_output <- derive_summary_records(
    input,
    by_vars = exprs(x),
    analysis_var = y,
    summary_fun = mean,
    filter = n() > 2,
    set_values_to = exprs(d = "MEAN")
  )
  expected_output <- tibble::tibble(
    x = c(rep(1, 2), rep(2, 4)),
    y = c(9:13, 12),
    z = c(1, 1, 2, 1, 1, NA),
    d = c(rep(NA, 5), "MEAN")
  )

  expect_dfs_equal(actual_output, expected_output, keys = c("x", "y", "z"))

  actual_output <- derive_summary_records(
    input,
    by_vars = exprs(x),
    analysis_var = y,
    summary_fun = mean,
    filter = z == 1,
    set_values_to = exprs(d = "MEAN")
  )
  expected_output <- tibble::tibble(
    x = c(rep(1, 3), rep(2, 4)),
    y = c(9:10, 9.5, 11:13, 12.5),
    z = c(1, 1, NA, 2, 1, 1, NA),
    d = c(rep(NA, 2), "MEAN", rep(NA, 3), "MEAN")
  )

  expect_dfs_equal(actual_output, expected_output, keys = c("x", "y", "z"))
})


## Test 6: Errors ----
test_that("derive_summary_records Test 6: Errors", {
  input <- tibble::tibble(x = rep(1:4, each = 4), y = rep(1:2, each = 8), z = runif(16))

  # Is by_vars quosures/`exprs()` object?
  expect_error(
    derive_summary_records(
      input,
      by_vars = "x",
      analysis_var = z,
      summary_fun = mean
    ),
    regexp = "`arg` must be an object of class 'list' but is `\"x\"`"
  )

  # Does by_vars exist in input dataset?
  expect_error(
    derive_summary_records(
      input,
      by_vars = exprs(a),
      analysis_var = z,
      summary_fun = mean
    ),
    regexp = "Required variable `a` is missing"
  )

  # summary_fun must be a single function
  expect_error(
    derive_summary_records(
      input,
      by_vars = exprs(x),
      analysis_var = y,
      summary_fun = list(mean, sum)
    ),
    regexp = "`summary_fun` must be an object of class 'function' but is a list"
  )

  # summary_fun must be a single function
  expect_error(
    derive_summary_records(
      input,
      by_vars = exprs(x),
      analysis_var = z,
      summary_fun = ~mean
    ),
    regexp = paste(
      "`summary_fun` must be an object of class 'function'",
      "but is an object of class 'formula'"
    )
  )
})
