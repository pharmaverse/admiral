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
    order = vars(AVISITN, LBSEQ),
    by_vars = vars(USUBJID),
    mode = "last",
    set_values_to = vars(DTYPE = "LOV")
  )

  expect_dfs_equal(
    base = expected_output,
    compare = actual_output,
    keys = c("USUBJID", "AVISITN", "LBSEQ", "DTYPE")
  )
})
