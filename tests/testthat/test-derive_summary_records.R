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

## Test 3: Errors ----
test_that("derive_summary_records Test 3: Errors", {
  input <- tibble::tibble(x = rep(1:4, each = 4), y = rep(1:2, each = 8), z = runif(16))

  # Is by_vars `exprs()` object?
  expect_error(
    derive_summary_records(
      input,
      by_vars = "x",
      set_values_to = exprs(
        z = mean(z)
      )
    ),
    class = "assert_vars"
  )

  # Does by_vars exist in input dataset?
  expect_error(
    derive_summary_records(
      dataset_add = input,
      by_vars = exprs(a),
      set_values_to = exprs(
        z = mean(z)
      )
    ),
    class = "assert_data_frame"
  )
})

## Test 4: make sure dataset_add works ----
test_that("derive_summary_records Test 4: make sure dataset_add works", {
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

## Test 5: test missing values with dataset_ref ----
test_that("derive_summary_records Test 5: test missing values with dataset_ref", {
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

## Test 6: test using constant_values ----
test_that("derive_summary_records Test 7: test using constant_values", {
  input <- tibble::tribble(
    ~subj, ~visit,       ~aval, ~seq,      ~param,
    "1",        1,          10,    1,     "VALUE",
    "1",        1,          14,    2,     "VALUE",
    "1",        1,           9,    3,     "VALUE",
    "1",        2,          11,    4,     "VALUE",
    "2",        2,          14,    1,     "VALUE"
  )

  expected_output <- bind_rows(
    input,
    tibble::tribble(
      ~subj, ~visit,       ~aval,         ~type,      ~param,
      "1",        1,          11,     "AVERAGE",    "RESULT",
      "1",        2,          11,     "AVERAGE",    "RESULT",
      "2",        2,          14,     "AVERAGE",    "RESULT"
    )
  )

  actual_output <- input %>%
    derive_summary_records(
      dataset_add = input,
      by_vars = exprs(subj, visit),
      constant_values = exprs(
        param = "RESULT"
      ),
      set_values_to = exprs(
        aval = mean(aval, na.rm = TRUE),
        type = "AVERAGE"
      )
    )

  expect_dfs_equal(
    base = expected_output,
    compare = actual_output,
    keys = c("subj", "visit", "seq", "type", "param")
  )
})

## Test 7: test missing values with constant values ----
test_that("derive_summary_records Test 8: test missing values using constant_values", {
  input <- tibble::tribble(
    ~subj, ~visit,       ~aval, ~seq,      ~param,
    "1",        1,          10,    1,     "VALUE",
    "1",        1,          14,    2,     "VALUE",
    "1",        1,           9,    3,     "VALUE",
    "1",        2,          11,    4,     "VALUE",
    "2",        2,    NA_real_,    1,     "VALUE"
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
      ~subj, ~visit,       ~aval,         ~type,      ~param,
      "1",        1,          11,     "AVERAGE",    "RESULT",
      "1",        2,          11,     "AVERAGE",    "RESULT",
      "2",        1,      999999,     "MISSING",    "RESULT",
      "2",        2,    NA_real_,     "AVERAGE",    "RESULT"
    )
  )

  actual_output <- input %>%
    derive_summary_records(
      dataset_add = input,
      dataset_ref = input_ref,
      by_vars = exprs(subj, visit),
      constant_values = exprs(
        param = "RESULT"
      ),
      set_values_to = exprs(
        aval = mean(aval, na.rm = TRUE),
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


## Test 6: error if no summary function ----
test_that("derive_summary_records Test 6: error if no summary function", {
  adbds <- tibble::tribble(
    ~AVISIT,  ~ASEQ, ~AVAL,
    "WEEK 1",     1,    10,
    "WEEK 1",     2,    NA,
    "WEEK 2",     3,    NA,
    "WEEK 3",     4,    42,
    "WEEK 4",     5,    12,
    "WEEK 4",     6,    12,
    "WEEK 4",     7,    15
  )

  expect_snapshot(
    derive_summary_records(
      dataset_add = adbds,
      by_vars = exprs(AVISIT),
      set_values_to = exprs(MEANVIS = AVAL / 2)
    ),
    error = TRUE
  )
})
