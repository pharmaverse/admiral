## Test 1: works with defaults ----
test_that("derive_vars_crit_flag Test 1: works with defaults", {
  expected <- tibble::tribble(
    ~AVAL, ~CRIT1FL,      ~CRIT1,
    23,    NA_character_, NA_character_,
    42,    "Y",           "AVAL > 40",
    NA,    NA_character_, NA_character_
  )

  expect_dfs_equal(
    base = expected,
    compare = derive_vars_crit_flag(
      select(expected, AVAL),
      condition = AVAL > 40,
      description = "AVAL > 40"
    ),
    keys = "AVAL"
  )
})

## Test 2: create numeric flag ----
test_that("derive_vars_crit_flag Test 2: create numeric flag", {
  expected <- tibble::tribble(
    ~AVAL, ~CRIT1FL,      ~CRIT1FN,    ~CRIT1,
    23,    NA_character_, NA_integer_, NA_character_,
    42,    "Y",           1L,          "AVAL > 40",
    NA,    NA_character_, NA_integer_, NA_character_
  )

  expect_dfs_equal(
    base = expected,
    compare = derive_vars_crit_flag(
      select(expected, AVAL),
      condition = AVAL > 40,
      description = "AVAL > 40",
      create_numeric_flag = TRUE
    ),
    keys = "AVAL"
  )
})

## Test 3: using values Y and N ----
test_that("derive_vars_crit_flag Test 3: using values Y and N", {
  expected <- tibble::tribble(
    ~PARAMCD, ~AVAL, ~CRIT2FL,      ~CRIT2FN,    ~CRIT2,
    "AST",    23,    "N",           0L,          "AST > 40",
    "AST",    42,    "Y",           1L,          "AST > 40",
    "AST",    NA,    NA_character_, NA_integer_, "AST > 40",
    "ALT",    26,    "N",           0L,          "ALT > 40",
    "ALT",    56,    "Y",           1L,          "ALT > 40"
  )

  expect_dfs_equal(
    base = expected,
    compare = derive_vars_crit_flag(
      select(expected, PARAMCD, AVAL),
      crit_nr = 2,
      condition = AVAL > 40,
      description = paste(PARAMCD, "> 40"),
      values_yn = TRUE,
      create_numeric_flag = TRUE
    ),
    keys = "AVAL"
  )
})

## Test 4: error if invalid condition (var not in input) ----
test_that("derive_vars_crit_flag Test 4: error if invalid condition (var not in input)", {
  input <- tibble::tribble(
    ~AVAL,
    23,
    42,
    NA,
  )

  expect_snapshot(
    derive_vars_crit_flag(
      input,
      condition = AVAL > 3 * ANRHI,
      description = "> 3ULN"
    ),
    error = TRUE
  )
})

## Test 5: error if invalid description ----
test_that("derive_vars_crit_flag Test 5: error if invalid description (PARAMCD not in input)", {
  input <- tibble::tribble(
    ~PARAMCD, ~AVAL,
    "AST",    23,
    "AST",    42,
    "AST",    NA,
    "ALT",    26,
    "ALT",    56
  )

  expect_snapshot(
    derive_vars_crit_flag(
      select(input, AVAL),
      crit_nr = 2,
      condition = AVAL > 40,
      description = paste(PARAMCD, "> 40"),
      values_yn = TRUE,
      create_numeric_flag = TRUE
    ),
    error = TRUE
  )
})
