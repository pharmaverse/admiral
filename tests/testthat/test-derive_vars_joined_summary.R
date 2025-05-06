## Test 1: summing up values up to current record ----
test_that("derive_vars_joined_summary Test 1: summing up values up to current record", {
  adex <- tibble::tribble(
    ~USUBJID, ~ADY, ~AVAL,
    "1",         1,    10,
    "1",         8,    20,
    "1",        15,    10,
    "2",         8,     5
  )

  expected <- tibble::tribble(
    ~USUBJID, ~ADY, ~AEDECOD,     ~CUMDOSE,
    "1",         2, "Fatigue",          10,
    "1",         9, "Influenza",        30,
    "1",        15, "Theft",            40,
    "1",        15, "Fatigue",          40,
    "2",         4, "Parasomnia",       NA,
    "3",         2, "Truancy",          NA
  )

  adae <- select(expected, -CUMDOSE)

  actual <- derive_vars_joined_summary(
    dataset = adae,
    dataset_add = adex,
    by_vars = exprs(USUBJID),
    filter_join = ADY.join <= ADY,
    join_type = "all",
    join_vars = exprs(ADY),
    new_vars = exprs(CUMDOSE = sum(AVAL, na.rm = TRUE))
  )

  expect_dfs_equal(
    base = expected,
    compare = actual,
    keys = c("USUBJID", "ADY", "AEDECOD")
  )
})
## Test 2: missing_values works ----
test_that("derive_vars_joined_summary Test 2: missing_values works", {
  adex <- tibble::tribble(
    ~USUBJID, ~ADY, ~AVAL,
    "1",         1,    10,
    "1",         8,    20,
    "1",        15,    10,
    "2",         8,     5
  )

  expected <- tibble::tribble(
    ~USUBJID, ~ADY, ~AEDECOD,     ~CUMDOSE,
    "1",         2, "Fatigue",          10,
    "1",         9, "Influenza",        30,
    "1",        15, "Theft",            40,
    "2",         4, "Parasomnia",        0,
    "3",         2, "Truancy",           0
  )

  adae <- select(expected, -CUMDOSE)

  actual <- derive_vars_joined_summary(
    dataset = adae,
    dataset_add = adex,
    by_vars = exprs(USUBJID),
    filter_join = ADY.join <= ADY,
    join_type = "all",
    join_vars = exprs(ADY),
    new_vars = exprs(CUMDOSE = sum(AVAL, na.rm = TRUE)),
    missing_values = exprs(CUMDOSE = 0)
  )

  expect_dfs_equal(
    base = expected,
    compare = actual,
    keys = c("USUBJID", "ADY", "AEDECOD")
  )
})

## Test 3: error if new variable in input dataset ----
test_that("derive_vars_joined_summary Test 3: error if new variable in input dataset", {
  adex <- tibble::tribble(
    ~USUBJID, ~ADY, ~AVAL,
    "1",         1,    10,
    "1",         8,    20,
    "1",        15,    10,
    "2",         8,     5
  )

  adae <- tibble::tribble(
    ~USUBJID, ~ADY, ~AEDECOD,     ~CUMDOSE,
    "1",         2, "Fatigue",          10,
    "1",         9, "Influenza",        30,
    "1",        15, "Theft",            40,
    "1",        15, "Fatigue",          40,
    "2",         4, "Parasomnia",       NA,
    "3",         2, "Truancy",          NA
  )

  expect_snapshot(
    derive_vars_joined_summary(
      dataset = adae,
      dataset_add = adex,
      by_vars = exprs(USUBJID),
      filter_join = ADY.join <= ADY,
      join_type = "all",
      join_vars = exprs(ADY),
      new_vars = exprs(CUMDOSE = sum(AVAL, na.rm = TRUE))
    ),
    error = TRUE
  )
})
