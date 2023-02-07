## Test 1: `by_vars` are missing ----
test_that("derive_expected_records Test 1: `by_vars` are missing", {
  input <- tibble::tribble(
    ~USUBJID, ~PARAMCD, ~AVISITN, ~AVISIT, ~AVAL,
    "1", "a", 1, "WEEK 1", 10,
    "1", "b", 2, "WEEK 2", 11,
    NA_character_, "a", 1, "WEEK 1", 12,
    NA_character_, "a", 2, "WEEK 2", 13,
    NA_character_, "b", 2, "WEEK 2", 14,
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
    dataset_expected_obs = expected_obsv,
    by_vars = vars(USUBJID),
    set_values_to = vars(DTYPE = "DERIVED")
  )

  expect_dfs_equal(
    base = expected_output,
    compare = actual_output,
    keys = c("USUBJID", "PARAMCD", "AVISITN", "AVISIT", "DTYPE")
  )
})
