## Test 1: creates new record per group and groups are retained ----
test_that("derive_summary_records Test 1: creates new record per group and groups are retained", {
  input <- tibble::tribble(
    ~subj, ~visit,       ~val, ~seq,
    "1",        1,         10,    1,
    "1",        1,         14,    2,
    "1",        1,          9,    3,
    "1",        2,         11,    4,
    "2",        2,   NA_real_,    1
  )

  actual_output <- input %>%
    derive_summary_records(
      dataset_add = input,
      by_vars = exprs(subj, visit),
      set_values_to = exprs(
        val = mean(val),
        seq = max(seq),
        type = "AVERAGE"
      )
    )

  expected_output <- bind_rows(
    input,
    tibble::tribble(
      ~subj, ~visit,       ~val, ~seq,
      "1",        1,         11,    3,
      "1",        2,         11,    4,
      "2",        2,   NA_real_,    1
    ) %>%
      mutate(type = "AVERAGE")
  )

  expect_dfs_equal(
    base = expected_output,
    compare = actual_output,
    keys = c("subj", "visit", "seq", "type")
  )
})

## Test 2: Filter record within `by_vars` ----
test_that("derive_summary_records Test 2: Filter record within `by_vars`", {
  input <- tibble::tribble(
    ~subj, ~visit,       ~val, ~seq,
    "1",        1,         10,    1,
    "1",        1,         14,    2,
    "1",        1,          9,    3,
    "1",        2,         11,    4,
    "2",        2,   NA_real_,    1
  )

  actual_output <- input %>%
    derive_summary_records(
      dataset_add = input,
      by_vars = exprs(subj, visit),
      filter_add = n() > 2,
      set_values_to = exprs(
        val = mean(val),
        seq = max(seq),
        type = "AVERAGE"
      )
    )

  expected_output <- bind_rows(
    input,
    tibble::tribble(
      ~subj, ~visit,       ~val, ~seq,
      "1",        1,         11,    3,
    ) %>%
      mutate(type = "AVERAGE")
  )

  expect_dfs_equal(
    base = expected_output,
    compare = actual_output,
    keys = c("subj", "visit", "seq", "type")
  )
})

## Test 3: error if by_vars is not an exprs() object ----
test_that("derive_summary_records Test 3: error if by_vars is not an exprs() object", {
  input <- tibble::tribble(
    ~USUBJID, ~AVISIT,  ~AVAL,
    "1",      "WEEK 1",    10,
    "1",      "WEEK 1",    14,
    "1",      "WEEK 1",     9,
    "2",      "WEEK 1",    11,
    "2",      "WEEK 1",    12
  )

  expect_error(
    derive_summary_records(
      dataset_add = input,
      by_vars = "USUBJID",
      set_values_to = exprs(
        AVAL = mean(AVAL)
      )
    ),
    class = "assert_vars"
  )
})

## Test 4: error if by_vars does not exist in input dataset ----
test_that("derive_summary_records Test 4: error if by_vars does not exist in input dataset", {
  input <- tibble::tribble(
    ~USUBJID, ~AVISIT,  ~AVAL,
    "1",      "WEEK 1",    10,
    "1",      "WEEK 1",    14,
    "1",      "WEEK 1",     9,
    "2",      "WEEK 1",    11,
    "2",      "WEEK 1",    12
  )

  expect_error(
    derive_summary_records(
      dataset_add = input,
      by_vars = exprs(USUBJID, NONEXISTENT),
      set_values_to = exprs(
        AVAL = mean(AVAL)
      )
    ),
    class = "assert_data_frame"
  )
})

## Test 5: error if set_values_to returns more than one value per by group ----
test_that("derive_summary_records Test 5: error if set_values_to returns more than one value per by group", {
  input <- tibble::tribble(
    ~USUBJID, ~AVISIT,  ~AVAL, ~ASEQ,
    "1",      "WEEK 1",    10,     1,
    "1",      "WEEK 1",    14,     2,
    "1",      "WEEK 1",     9,     3,
    "2",      "WEEK 1",    11,     1,
    "2",      "WEEK 1",    12,     2
  )

  expect_error(
    derive_summary_records(
      dataset_add = input,
      by_vars = exprs(USUBJID, AVISIT),
      set_values_to = exprs(
        AVAL = AVAL
      )
    ),
    regexp = "Column\\(s\\) in `set_values_to` must return a single value"
  )
})

## Test 6: no error when set_values_to uses mean() properly ----
test_that("derive_summary_records Test 6: no error when set_values_to uses mean() properly", {
  input <- tibble::tribble(
    ~USUBJID, ~AVISIT,  ~AVAL, ~ASEQ,
    "1",      "WEEK 1",    10,     1,
    "1",      "WEEK 1",    14,     2,
    "1",      "WEEK 1",     9,     3,
    "2",      "WEEK 1",    11,     1,
    "2",      "WEEK 1",    12,     2
  )

  expected_output <- bind_rows(
    input,
    tibble::tribble(
      ~USUBJID, ~AVISIT,  ~AVAL,       ~DTYPE,
      "1",      "WEEK 1",    11,    "AVERAGE",
      "2",      "WEEK 1",  11.5,    "AVERAGE"
    )
  )

  actual_output <- input %>%
    derive_summary_records(
      dataset_add = input,
      by_vars = exprs(USUBJID, AVISIT),
      set_values_to = exprs(
        AVAL = mean(AVAL, na.rm = TRUE),
        DTYPE = "AVERAGE"
      )
    )

  expect_dfs_equal(
    base = expected_output,
    compare = actual_output,
    keys = c("USUBJID", "AVISIT", "ASEQ", "DTYPE")
  )
})

## Test 7: no error when set_values_to uses max() properly ----
test_that("derive_summary_records Test 7: no error when set_values_to uses max() properly", {
  input <- tibble::tribble(
    ~USUBJID, ~AVISIT,  ~AVAL, ~ASEQ,
    "1",      "WEEK 1",    10,     1,
    "1",      "WEEK 1",    14,     2,
    "1",      "WEEK 1",     9,     3,
    "2",      "WEEK 1",    11,     1,
    "2",      "WEEK 1",    12,     2
  )

  expected_output <- bind_rows(
    input,
    tibble::tribble(
      ~USUBJID, ~AVISIT,  ~AVAL,    ~DTYPE,
      "1",      "WEEK 1",    14,   "MAXIMUM",
      "2",      "WEEK 1",    12,   "MAXIMUM"
    )
  )

  actual_output <- input %>%
    derive_summary_records(
      dataset_add = input,
      by_vars = exprs(USUBJID, AVISIT),
      set_values_to = exprs(
        AVAL = max(AVAL, na.rm = TRUE),
        DTYPE = "MAXIMUM"
      )
    )

  expect_dfs_equal(
    base = expected_output,
    compare = actual_output,
    keys = c("USUBJID", "AVISIT", "ASEQ", "DTYPE")
  )
})

## Test 8: error when mean() is used but returns multiple values per group ----
test_that("derive_summary_records Test 8: error when mean() is used but returns multiple values per group", {
  input <- tibble::tribble(
    ~USUBJID, ~AVISIT,  ~AVAL, ~ASEQ,
    "1",      "WEEK 1",    10,     1,
    "1",      "WEEK 1",    14,     2,
    "1",      "WEEK 1",     9,     3,
    "2",      "WEEK 1",    11,     1,
    "2",      "WEEK 1",    12,     2
  )

  expect_error(
    derive_summary_records(
      dataset_add = input,
      by_vars = exprs(USUBJID, AVISIT),
      set_values_to = exprs(
        AVAL = c(mean(AVAL, na.rm = TRUE), max(AVAL, na.rm = TRUE))
      )
    ),
    regexp = "Column\\(s\\) in `set_values_to` must return a single value"
  )
})

## Test 9: error when max() is used but returns multiple values per group ----
test_that("derive_summary_records Test 9: error when max() is used but returns multiple values per group", {
  input <- tibble::tribble(
    ~USUBJID, ~AVISIT,  ~AVAL, ~ASEQ,
    "1",      "WEEK 1",    10,     1,
    "1",      "WEEK 1",    14,     2,
    "1",      "WEEK 1",     9,     3,
    "2",      "WEEK 1",    11,     1,
    "2",      "WEEK 1",    12,     2
  )

  expect_error(
    derive_summary_records(
      dataset_add = input,
      by_vars = exprs(USUBJID, AVISIT),
      set_values_to = exprs(
        AVAL = c(max(AVAL, na.rm = TRUE), min(AVAL, na.rm = TRUE))
      )
    ),
    regexp = "Column\\(s\\) in `set_values_to` must return a single value"
  )
})

## Test 10: make sure dataset_add works ----
test_that("derive_summary_records Test 10: make sure dataset_add works", {
  input <- tibble::tribble(
    ~subj, ~visit,       ~val, ~seq,
    "1",        1,         10,    1,
    "1",        1,         14,    2,
    "1",        1,          9,    3,
    "1",        2,         11,    4,
    "2",        2,   NA_real_,    1
  )
  input_add <- tibble::tribble(
    ~subj, ~visit,       ~add_val, ~seq,
    "1",        1,            100,    1,
    "1",        1,            140,    2,
    "1",        1,             90,    3
  )
  expected_output <- bind_rows(
    input,
    tibble::tribble(
      ~subj, ~visit, ~val, ~type,
      "1", 1, 110, "AVERAGE"
    )
  )
  actual_output <- input %>%
    derive_summary_records(
      dataset_add = input_add,
      by_vars = exprs(subj, visit),
      set_values_to = exprs(
        val = mean(add_val, na.rm = TRUE),
        type = "AVERAGE"
      )
    )
  expect_dfs_equal(
    base = expected_output,
    compare = actual_output,
    keys = c("subj", "visit", "seq", "type")
  )
})

## Test 11: test missing values ----
test_that("derive_summary_records Test 11: test missing values with dataset_ref", {
  input <- tibble::tribble(
    ~subj, ~visit,       ~val, ~seq,
    "1",        1,         10,    1,
    "1",        1,         14,    2,
    "1",        1,          9,    3,
    "1",        2,         11,    4,
    "2",        2,   NA_real_,    1
  )

  input_ref <- tibble::tribble(
    ~subj, ~visit,
    "1", 1,
    "1", 2,
    "2", 1,
    "2", 2,
  )

  expected_output <- bind_rows(
    input,
    tibble::tribble(
      ~subj, ~visit,       ~aval,         ~type,
      "1",        1,          11,     "AVERAGE",
      "1",        2,          11,     "AVERAGE",
      "2",        1,      999999,     "MISSING",
      "2",        2,    NA_real_,     "AVERAGE",
    )
  )

  actual_output <- input %>%
    derive_summary_records(
      dataset_add = input,
      dataset_ref = input_ref,
      by_vars = exprs(subj, visit),
      set_values_to = exprs(
        aval = mean(val, na.rm = TRUE),
        type = "AVERAGE"
      ),
      missing_values = exprs(aval = 999999, type = "MISSING")
    )

  expect_dfs_equal(
    base = expected_output,
    compare = actual_output,
    keys = c("subj", "visit", "seq", "type")
  )
})
