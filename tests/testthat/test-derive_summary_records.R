library(tibble)
library(dplyr)

test_that("creates a new record for each group and new data frame retains grouping", {
  # group --> 4
  df <- tibble(x = rep(1:4, each = 4), y = rep(1:2, each = 8), z = runif(16))
  df <- df %>% group_by(x, y)
  cf <- df %>%
    derive_summary_records(
      by_vars = vars(x, y),
      fns = list(z ~ mean)
    )
  expect_equal(nrow(df) + 4, nrow(cf))
  expect_equal(group_vars(df), group_vars(cf))
})

test_that("unique records within `by_vars` are retained, not set to <NA>", {
  df <- tibble(x = rep(1:2, each = 2),
               y = rep_len("a", 4),
               z = letters[1:4],
               a = sample(1:100, 4))
  cf <- df %>%
    derive_summary_records(
      by_vars = vars(x),
      fns = list(a ~ sum)
    )
  expect_equal(cf$y, rep_len("a", 6))
  expect_equal(cf$z, c("a", "b", NA, "c", "d", NA))
})


test_that("LHS of `fns` accepts multiple variable", {
  df <- tibble(x = rep(1:2, each = 2), y = 9:12, z = 101:104)
  cf <- derive_summary_records(
    df,
    by_vars = vars(x),
    fns = list(vars(y, z) ~ mean)
  )
  tf <- tibble(x = rep(1:2, each = 4),
               y = c(9:10, 9.5, NA, 11:12, 11.5, NA),
               z = c(101:102, NA, 101.5, 103:104, NA, 103.5))
  expect_equal(tf, cf)
})

test_that("`fns` accepts single forumula without wrapping into list", {
  df <- tibble(x = rep(1:2, each = 2), y = 9:12, z = 101:104)
  cf <- derive_summary_records(
    df,
    by_vars = vars(x),
    fns = y ~ mean
  )
  tf <- tibble(x = rep(1:2, each = 3),
               y = c(9:10, 9.5, 11:12, 11.5),
               z = c(101:102, NA, 103:104, NA))
  expect_equal(tf, cf)
})

test_that("`fns` as inlined", {
  df <- tibble(x = rep(1:2, each = 2), y = 9:12, z = 101:104)
  cf <- derive_summary_records(
    df,
    by_vars = vars(x),
    fns = list(y ~ mean(., na.rm = TRUE))
  )
  tf <- tibble(x = rep(1:2, each = 3),
               y = c(9:10, 9.5, 11:12, 11.5),
               z = c(101:102, NA, 103:104, NA))
  expect_equal(tf, cf)
})

test_that("set new value to a derived record", {
  df <- tibble(x = rep(1:2, each = 2), y = 9:12)
  cf <- derive_summary_records(
    df,
    by_vars = vars(x),
    fns = list(y ~ mean),
    set_values_to = vars(z = "MEAN")
  )
  tf <- tibble(x = rep(1:2, each = 3),
               y = c(9:10, 9.5, 11:12, 11.5),
               z = c(NA, NA, "MEAN", NA, NA, "MEAN"))
  expect_equal(tf, cf)
})

test_that("drop a value from derived record", {
  df <- tibble(x = rep(1:2, each = 2), y = 9:12, z = rep_len(1, 4))
  cf <- derive_summary_records(
    df,
    by_vars = vars(x),
    fns = list(y ~ mean),
    set_values_to = vars(d = "MEAN"),
    drop_values_from = vars(z)
  )
  tf <- tibble(x = rep(1:2, each = 3),
               y = c(9:10, 9.5, 11:12, 11.5),
               z = c(1, 1, NA, 1, 1, NA),
               d = c(NA, NA, "MEAN", NA, NA, "MEAN"))
  expect_equal(tf, cf)
})

test_that("check `set_values_to` mapping", {
  df <- tibble(x = rep(1:4, each = 4), y = rep(1:2, each = 8), z = runif(16))
  df <- df %>% group_by(x, y)
  cf <- df %>%
    derive_summary_records(
      by_vars = vars(x, y),
      fns = list(z ~ mean, z ~ sum),
      set_values_to = vars(d = c("MEAN", "SUM"))
    )
  tf <- rep(c(rep_len(NA, 4), "MEAN", "SUM"), 4)
  expect_equal(cf$d, tf)
})

# Errors ---

test_that("Errors", {
  # Is by_vars and drop_values_from are quosures/`vars()` object?
  df <- tibble(x = rep(1:4, each = 4), y = rep(1:2, each = 8), z = runif(16))

  expect_error(
    derive_summary_records(
      df,
      by_vars = "x",
      fns = list(z ~ mean)))

  expect_error(
    derive_summary_records(
      df,
      by_vars = vars(x),
      fns = list(z ~ mean),
      drop_values_from = "z"))

  # Is by_vars and drop_values_from exits in input dataset?
  expect_error(
    derive_summary_records(
      df,
      by_vars = vars(a),
      fns = list(z ~ mean),
      drop_values_from = vars(b)))

  # Can't have multiple function for a single analysis variable
  expect_error(
    derive_summary_records(
      df,
      by_vars = vars(x),
      fns = list(z ~ list(mean, sum))))

  # check function must be a formula
  expect_error(
    derive_summary_records(
      df,
      by_vars = vars(x),
      fns = list(z = mean)))

  # Is `set_values_to` is a quosures?
  expect_error(
    derive_summary_records(
      df,
      by_vars = vars(x),
      fns = list(z ~ mean),
      set_values_to = list(d = "a")))

  # Is length of `set_values_to` equal to derived records within a `by_vars`?
  expect_error(
    derive_summary_records(
      df,
      by_vars = vars(x),
      fns = list(z ~ mean, z ~ sum),
      set_values_to = vars(d = "a")))

  # Problem with handling RHS of `fns` argument
  expect_error(
    derive_summary_records(
      df,
      by_vars = vars(x),
      fns = list(z ~ vars(mean, sum)),
      set_values_to = vars(d = "a")))

})

