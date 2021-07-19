library(tibble)
library(dplyr)

test_that("creates a new record for each group and new data frame retains grouping", {
  # group --> 4
  input <- tibble(x = rep(1:4, each = 4), y = rep(1:2, each = 8), z = runif(16))
  input <- input %>% group_by(x)
  actual_output <- input %>%
    derive_summary_records(
      by_vars = vars(x, y),
      fns = list(z ~ mean)
    )

  expect_equal(nrow(actual_output), nrow(input) + 4)
  expect_equal(group_vars(actual_output), group_vars(input))
})

test_that("unique records within `by_vars` are retained, not set to <NA>", {
  input <- tibble(
    x = rep(1:2, each = 2),
    y = rep_len("a", 4),
    z = letters[1:4],
    a = sample(1:100, 4)
  )
  output <- derive_summary_records(
    input,
    by_vars = vars(x),
    fns = list(a ~ sum)
  )

  expect_equal(output$y, rep_len("a", 6))
  expect_equal(output$z, c("a", "b", NA, "c", "d", NA))
})


test_that("LHS of `fns` accepts multiple variable", {
  input <- tibble(x = rep(1:2, each = 2), y = 9:12, z = 101:104)
  actual_output <- derive_summary_records(
    input,
    by_vars = vars(x),
    fns = list(vars(y, z) ~ mean)
  )
  expected_output <- tibble(
    x = rep(1:2, each = 4),
    y = c(9:10, 9.5, NA, 11:12, 11.5, NA),
    z = c(101:102, NA, 101.5, 103:104, NA, 103.5)
  )

  expect_dfs_equal(actual_output, expected_output, keys = c("x", "y", "z"))
})

test_that("`fns` accepts single forumula without wrapping into list", {
  input <- tibble(x = rep(1:2, each = 2), y = 9:12, z = 101:104)
  actual_output <- derive_summary_records(
    input,
    by_vars = vars(x),
    fns = y ~ mean
  )
  expected_output <- tibble(
    x = rep(1:2, each = 3),
    y = c(9:10, 9.5, 11:12, 11.5),
    z = c(101:102, NA, 103:104, NA)
  )

  expect_dfs_equal(actual_output, expected_output, keys = c("x", "y", "z"))
})

test_that("`fns` as inlined", {
  input <- tibble(x = rep(1:2, each = 2), y = 9:12, z = 101:104)
  actual_output <- derive_summary_records(
    input,
    by_vars = vars(x),
    fns = list(y ~ mean(., na.rm = TRUE))
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
    fns = list(y ~ mean),
    set_values_to = vars(z = "MEAN")
  )
  expected_output <- tibble(
    x = rep(1:2, each = 3),
    y = c(9:10, 9.5, 11:12, 11.5),
    z = c(NA, NA, "MEAN", NA, NA, "MEAN")
  )

  expect_dfs_equal(actual_output, expected_output, keys = c("x", "y", "z"))
})

test_that("drop a value from derived record", {
  input <- tibble(x = rep(1:2, each = 2), y = 9:12, z = rep_len(1, 4))
  actual_output <- derive_summary_records(
    input,
    by_vars = vars(x),
    fns = list(y ~ mean),
    set_values_to = vars(d = "MEAN"),
    drop_values_from = vars(z)
  )
  expected_output <- tibble(
    x = rep(1:2, each = 3),
    y = c(9:10, 9.5, 11:12, 11.5),
    z = c(1, 1, NA, 1, 1, NA),
    d = c(NA, NA, "MEAN", NA, NA, "MEAN")
  )

  expect_dfs_equal(actual_output, expected_output, keys = c("x", "y", "z"))
})

test_that("check `set_values_to` mapping", {
  input <- tibble(x = rep(1:4, each = 4), y = rep(1:2, each = 8), z = runif(16))
  input <- input %>% group_by(x, y)
  actual_output <- derive_summary_records(
    input,
    by_vars = vars(x, y),
    fns = list(z ~ mean, z ~ sum),
    set_values_to = vars(d = c("MEAN", "SUM"))
  )
  tf <- rep(c(rep_len(NA, 4), "MEAN", "SUM"), 4)

  expect_equal(actual_output$d, tf)
})

test_that("Filter record within `by_vars`", {
  input <- tibble(x = c(rep(1:2, each = 2), 2), y = 9:13, z = c(1, 1, 2, 1, 1))

  actual_output <- derive_summary_records(
    input,
    by_vars = vars(x),
    fns = list(y ~ mean),
    filter = n() > 2,
    set_values_to = vars(d = "MEAN"),
    drop_values_from = vars(z)
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
    fns = list(y ~ mean),
    filter = z == 1,
    set_values_to = vars(d = "MEAN"),
    drop_values_from = vars(z)
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
  # Is by_vars and drop_values_from are quosures/`vars()` object?
  input <- tibble(x = rep(1:4, each = 4), y = rep(1:2, each = 8), z = runif(16))

  expect_error(
    derive_summary_records(
      input,
      by_vars = "x",
      fns = list(z ~ mean)
    )
  )

  expect_error(
    derive_summary_records(
      input,
      by_vars = vars(x),
      fns = list(z ~ mean),
      drop_values_from = "z"
    )
  )

  # Is by_vars and drop_values_from exits in input dataset?
  expect_error(
    derive_summary_records(
      input,
      by_vars = vars(a),
      fns = list(z ~ mean),
      drop_values_from = vars(b)
    )
  )

  # Can't have multiple function for a single analysis variable
  expect_error(
    derive_summary_records(
      input,
      by_vars = vars(x),
      fns = list(z ~ list(mean, sum))
    )
  )

  # check function must be a formula
  expect_error(
    derive_summary_records(
      input,
      by_vars = vars(x),
      fns = list(z = mean)
    )
  )

  # Is `set_values_to` is a quosures?
  expect_error(
    derive_summary_records(
      input,
      by_vars = vars(x),
      fns = list(z ~ mean),
      set_values_to = list(d = "a")
    )
  )

  # Is length of `set_values_to` equal to derived records within a `by_vars`?
  expect_error(
    derive_summary_records(
      input,
      by_vars = vars(x),
      fns = list(z ~ mean, z ~ sum),
      set_values_to = vars(d = "a")
    )
  )

  # Problem with handling RHS of `fns` argument
  expect_error(
    derive_summary_records(
      input,
      by_vars = vars(x),
      fns = list(z ~ vars(mean, sum)),
      set_values_to = vars(d = "a")
    )
  )
})
