# derive_extreme_records ----
## derive_extreme_records Test 1: add last observation for each group ----
test_that("derive_extreme_records Test 1: add last observation for each group", {
  input <- tibble::tribble(
    ~USUBJID, ~AVISITN, ~AVAL, ~LBSEQ,
    1, 1, 12, 1,
    1, 3, 9, 2,
    2, 2, 42, 1,
    3, 3, 14, 1,
    3, 3, 10, 2
  )

  expected_output <- bind_rows(
    input,
    tibble::tribble(
      ~USUBJID, ~AVISITN, ~AVAL, ~LBSEQ,
      1, 3, 9, 2,
      2, 2, 42, 1,
      3, 3, 10, 2
    ) %>%
      mutate(DTYPE = "LOV")
  )

  actual_output <- derive_extreme_records(
    input,
    order = exprs(AVISITN, LBSEQ),
    by_vars = exprs(USUBJID),
    mode = "last",
    set_values_to = exprs(DTYPE = "LOV")
  )

  expect_dfs_equal(
    base = expected_output,
    compare = actual_output,
    keys = c("USUBJID", "AVISITN", "LBSEQ", "DTYPE")
  )
})

# derive_extreme_records ----
## derive_extreme_records Test 2: keep variables specified in `keep_vars_source` in the new records ---- # nolint
test_that("derive_extreme_records Test 2: keep variables specified in `keep_vars_source` in the new records", { # nolint
  input <- tibble::tribble(
    ~USUBJID, ~AVISITN, ~AVAL, ~LBSEQ,
    1, 1, 12, 1,
    1, 3, 9, 2,
    2, 2, 42, 1,
    3, 3, 14, 1,
    3, 3, 10, 2
  )

  expected_output <- bind_rows(
    input,
    tibble::tribble(
      ~USUBJID, ~AVISITN, ~AVAL, ~LBSEQ,
      1, 3, 9, 2,
      2, 2, 42, 1,
      3, 3, 10, 2
    ) %>%
      select(USUBJID, AVISITN, AVAL) %>%
      mutate(DTYPE = "LOV")
  )

  actual_output <- derive_extreme_records(
    input,
    order = exprs(AVISITN, LBSEQ),
    by_vars = exprs(USUBJID),
    mode = "last",
    keep_vars_source = exprs(USUBJID, AVISITN, AVAL),
    set_values_to = exprs(DTYPE = "LOV")
  )

  expect_dfs_equal(
    base = expected_output,
    compare = actual_output,
    keys = c("USUBJID", "AVISITN", "LBSEQ", "DTYPE")
  )
})

# derive_extreme_records ----
## derive_extreme_records Test 3: keep all variables in the new records when `keep_vars_source` is NULL ---- # nolint
test_that("derive_extreme_records Test 3: keep all variables in the new records when `keep_vars_source` is NULL", { # nolint
  input <- tibble::tribble(
    ~USUBJID, ~AVISITN, ~AVAL, ~LBSEQ,
    1, 1, 12, 1,
    1, 3, 9, 2,
    2, 2, 42, 1,
    3, 3, 14, 1,
    3, 3, 10, 2
  )

  expected_output <- bind_rows(
    input,
    tibble::tribble(
      ~USUBJID, ~AVISITN, ~AVAL, ~LBSEQ,
      1, 3, 9, 2,
      2, 2, 42, 1,
      3, 3, 10, 2
    ) %>%
      mutate(DTYPE = "LOV")
  )

  actual_output <- derive_extreme_records(
    input,
    order = exprs(AVISITN, LBSEQ),
    by_vars = exprs(USUBJID),
    mode = "last",
    keep_vars_source = ,
    set_values_to = exprs(DTYPE = "LOV")
  )

  expect_dfs_equal(
    base = expected_output,
    compare = actual_output,
    keys = c("USUBJID", "AVISITN", "LBSEQ", "DTYPE")
  )
})
