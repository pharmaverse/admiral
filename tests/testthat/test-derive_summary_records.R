library(tibble)
library(dplyr)

test_that("creates a new record for each group and new data frame retains grouping", {
  # group --> 4
  input <- tibble(x = rep(1:4, each = 4), y = rep(1:2, each = 8), z = runif(16))
  actual_output <- input %>%
    derive_summary_records(
      by_vars = vars(x, y),
      analysis_var = z,
      summary_fun = mean
    )

  expect_equal(nrow(actual_output), nrow(input) + 4)
  expect_equal(group_vars(actual_output), group_vars(input))
})

test_that("`fns` as inlined", {
  input <- tibble(x = rep(1:2, each = 2), y = 9:12, z = 101:104)
  actual_output <- derive_summary_records(
    input,
    by_vars = vars(x),
    analysis_var = y,
    summary_fun = function(x) mean(x, na.rm = TRUE)
  )
  expected_output <- tibble(
    x = rep(1:2, each = 3),
    y = c(9:10, 9.5, 11:12, 11.5),
    z = c(101:102, NA, 103:104, NA)
  )

  expect_dfs_equal(actual_output, expected_output, keys = c("x", "y", "z"))
})

test_that("set new value to a derived record", {
  input <- tibble(x = rep(1:2, each = 2), y = 9:12)
  actual_output <- derive_summary_records(
    input,
    by_vars = vars(x),
    analysis_var = y,
    summary_fun = mean,
    set_values_to = vars(z = "MEAN")
  )
  expected_output <- tibble(
    x = rep(1:2, each = 3),
    y = c(9:10, 9.5, 11:12, 11.5),
    z = c(NA, NA, "MEAN", NA, NA, "MEAN")
  )

  expect_dfs_equal(actual_output, expected_output, keys = c("x", "y", "z"))
})

test_that("check `set_values_to` mapping", {
  input <- tibble(x = rep(1:4, each = 4), y = rep(1:2, each = 8), z = runif(16))
  actual_output <- input %>%
    derive_summary_records(
      by_vars = vars(x, y),
      analysis_var = z,
      summary_fun = mean,
      set_values_to = vars(d = "MEAN")
    ) %>%
    derive_summary_records(
      by_vars = vars(x, y),
      analysis_var = z,
      summary_fun = sum,
      set_values_to = vars(d = "SUM")
    )
  tf <- rep(c(NA, "MEAN", "SUM"), c(16, 4, 4))

  expect_equal(actual_output$d, tf)

  actual_output <- input %>%
    derive_summary_records(
      by_vars = vars(x, y),
      analysis_var = z,
      summary_fun = mean,
      set_values_to = vars(d = "MEAN", p1 = "PARAM1", p2 = "PARAM2")
    )
  tf <- rep(c(NA, "MEAN"), c(16, 4))
  tp1 <- rep(c(NA, "PARAM1"), c(16, 4))
  tp2 <- rep(c(NA, "PARAM2"), c(16, 4))

  expect_equal(actual_output$d, tf)
  expect_equal(actual_output$p1, tp1)
  expect_equal(actual_output$p2, tp2)
})

test_that("Filter record within `by_vars`", {
  input <- tibble(x = c(rep(1:2, each = 2), 2), y = 9:13, z = c(1, 1, 2, 1, 1))

  actual_output <- derive_summary_records(
    input,
    by_vars = vars(x),
    analysis_var = y,
    summary_fun = mean,
    filter = n() > 2,
    set_values_to = vars(d = "MEAN")
  )
  expected_output <- tibble(
    x = c(rep(1, 2), rep(2, 4)),
    y = c(9:13, 12),
    z = c(1, 1, 2, 1, 1, NA),
    d = c(rep(NA, 5), "MEAN")
  )

  expect_dfs_equal(actual_output, expected_output, keys = c("x", "y", "z"))

  actual_output <- derive_summary_records(
    input,
    by_vars = vars(x),
    analysis_var = y,
    summary_fun = mean,
    filter = z == 1,
    set_values_to = vars(d = "MEAN")
  )
  expected_output <- tibble(
    x = c(rep(1, 3), rep(2, 4)),
    y = c(9:10, 9.5, 11:13, 12.5),
    z = c(1, 1, NA, 2, 1, 1, NA),
    d = c(rep(NA, 2), "MEAN", rep(NA, 3), "MEAN")
  )

  expect_dfs_equal(actual_output, expected_output, keys = c("x", "y", "z"))
})

# Errors ---

test_that("Errors", {
  input <- tibble(x = rep(1:4, each = 4), y = rep(1:2, each = 8), z = runif(16))

  # Is by_vars quosures/`vars()` object?
  expect_error(
    derive_summary_records(
      input,
      by_vars = "x",
      analysis_var = z,
      summary_fun = mean
    ),
    regexp = "`by_vars` must be a list of unquoted variable names"
  )

  # Does by_vars exist in input dataset?
  expect_error(
    derive_summary_records(
      input,
      by_vars = vars(a),
      analysis_var = z,
      summary_fun = mean
    ),
    regexp = "Required variable `a` is missing"
  )

  # summary_fun must be a single function
  expect_error(
    derive_summary_records(
      input,
      by_vars = vars(x),
      analysis_var = y,
      summary_fun = list(mean, sum)
    ),
    regexp = "`summary_fun` must be an object of class 'function' but is a list"
  )

  # summary_fun must be a single function
  expect_error(
    derive_summary_records(
      input,
      by_vars = vars(x),
      analysis_var = z,
      summary_fun = ~mean
    ),
    regexp = paste(
      "`summary_fun` must be an object of class 'function'",
      "but is an object of class 'formula'"
    )
  )
})
