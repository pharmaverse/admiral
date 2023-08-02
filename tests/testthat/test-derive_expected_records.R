## Test 1: missing values in `by_vars` ----
test_that("derive_expected_records Test 1: missing values in `by_vars`", {
  input <- tibble::tribble(
    ~USUBJID, ~PARAMCD, ~AVISITN, ~AVISIT, ~AVAL,
    "1", "a", 1, "WEEK 1", 10,
    "1", "b", 2, "WEEK 2", 11,
    NA_character_, "a", 1, "WEEK 1", 12,
    NA_character_, "a", 2, "WEEK 2", 13,
    NA_character_, "b", 2, "WEEK 2", 14
  )

  expected_obsv <- tibble::tribble(
    ~PARAMCD, ~AVISITN, ~AVISIT,
    "a", 1, "WEEK 1",
    "a", 2, "WEEK 2",
    "b", 1, "WEEK 1",
    "b", 2, "WEEK 2"
  )

  expected_output <- bind_rows(
    input,
    tibble::tribble(
      ~USUBJID, ~PARAMCD, ~AVISITN, ~AVISIT, ~AVAL,
      "1", "a", 2, "WEEK 2", NA_real_,
      "1", "b", 1, "WEEK 1", NA_real_,
      NA_character_, "b", 1, "WEEK 1", NA_real_
    ) %>%
      mutate(DTYPE = "DERIVED")
  )

  actual_output <- derive_expected_records(
    dataset = input,
    dataset_ref = expected_obsv,
    by_vars = exprs(USUBJID),
    set_values_to = exprs(DTYPE = "DERIVED")
  )

  expect_dfs_equal(
    base = expected_output,
    compare = actual_output,
    keys = c("USUBJID", "PARAMCD", "AVISITN", "AVISIT", "DTYPE")
  )
})


## Test 2: `by_vars` = NULL ----
test_that("derive_expected_records Test 2: `by_vars` = NULL", {
  input <- tibble::tribble(
    ~USUBJID, ~PARAMCD, ~AVISITN, ~AVISIT, ~AVAL,
    "1", "a", 1, "WEEK 1", 10,
    "1", "b", 2, "WEEK 2", 11
  )

  expected_obsv <- tibble::tribble(
    ~PARAMCD, ~AVISITN, ~AVISIT,
    "a", 1, "WEEK 1",
    "a", 2, "WEEK 2",
    "b", 1, "WEEK 1",
    "b", 2, "WEEK 2"
  )

  expected_output <- bind_rows(
    input,
    tibble::tribble(
      ~USUBJID, ~PARAMCD, ~AVISITN, ~AVISIT, ~AVAL,
      NA_character_, "a", 2, "WEEK 2", NA_real_,
      NA_character_, "b", 1, "WEEK 1", NA_real_
    ) %>%
      mutate(DTYPE = "DERIVED")
  )

  actual_output <- derive_expected_records(
    dataset = input,
    dataset_ref = expected_obsv,
    by_vars = NULL,
    set_values_to = exprs(DTYPE = "DERIVED")
  )

  expect_dfs_equal(
    base = expected_output,
    compare = actual_output,
    keys = c("USUBJID", "PARAMCD", "AVISITN", "AVISIT", "DTYPE")
  )
})


## Test 3: visit variables are parameter independent ----
test_that("derive_expected_records Test 3: visit variables are parameter independent", {
  input <- tibble::tribble(
    ~USUBJID, ~PARAMCD, ~AVISITN, ~AVISIT, ~AVAL,
    "1", "a", 1, "WEEK 1", 10,
    "1", "b", 1, "WEEK 1", 11,
    "2", "a", 2, "WEEK 2", 12,
    "2", "b", 2, "WEEK 2", 14
  )

  expected_obsv <- tibble::tribble(
    ~AVISITN, ~AVISIT,
    1, "WEEK 1",
    2, "WEEK 2"
  )

  expected_output <- bind_rows(
    input,
    tibble::tribble(
      ~USUBJID, ~PARAMCD, ~AVISITN, ~AVISIT, ~AVAL,
      "1", "a", 2, "WEEK 2", NA_real_,
      "1", "b", 2, "WEEK 2", NA_real_,
      "2", "a", 1, "WEEK 1", NA_real_,
      "2", "b", 1, "WEEK 1", NA_real_
    ) %>%
      mutate(DTYPE = "DERIVED")
  )

  actual_output <- derive_expected_records(
    dataset = input,
    dataset_ref = expected_obsv,
    by_vars = exprs(USUBJID, PARAMCD),
    set_values_to = exprs(DTYPE = "DERIVED")
  )

  expect_dfs_equal(
    base = expected_output,
    compare = actual_output,
    keys = c("USUBJID", "PARAMCD", "AVISITN", "AVISIT", "DTYPE")
  )
})


## Test 4: visit variables are parameter dependent ----
test_that("derive_expected_records Test 4: visit variables are parameter dependent", {
  input <- tibble::tribble(
    ~USUBJID, ~PARAMCD, ~AVISITN, ~AVISIT, ~AVAL,
    "1", "a", 1, "WEEK 1", 10,
    "1", "b", 1, "WEEK 1", 11,
    "2", "a", 2, "WEEK 2", 12,
    "2", "b", 2, "WEEK 2", 14
  )

  expected_obsv <- tibble::tribble(
    ~PARAMCD, ~AVISITN, ~AVISIT,
    "a", 1, "WEEK 1",
    "a", 2, "WEEK 2",
    "b", 1, "WEEK 1"
  )

  expected_output <- bind_rows(
    input,
    tibble::tribble(
      ~USUBJID, ~PARAMCD, ~AVISITN, ~AVISIT, ~AVAL,
      "1", "a", 2, "WEEK 2", NA_real_,
      "2", "a", 1, "WEEK 1", NA_real_,
      "2", "b", 1, "WEEK 1", NA_real_
    ) %>%
      mutate(DTYPE = "DERIVED")
  )

  actual_output <- derive_expected_records(
    dataset = input,
    dataset_ref = expected_obsv,
    by_vars = exprs(USUBJID),
    set_values_to = exprs(DTYPE = "DERIVED")
  )

  expect_dfs_equal(
    base = expected_output,
    compare = actual_output,
    keys = c("USUBJID", "PARAMCD", "AVISITN", "AVISIT", "DTYPE")
  )
})
