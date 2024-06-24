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
    keys = "AVAL")
})

## Test 2: create numeric flag ----
test_that("derive_vars_crit_flag Test 2: create numeric flag", {
  expected <- tibble::tribble(
    ~AVAL, ~CRIT1FL,      ~CRIT1FLN, ~CRIT1,
    23,    NA_character_, NA_real_,  NA_character_,
    42,    "Y",           1,         "AVAL > 40",
    NA,    NA_character_, NA_real_,  NA_character_
  )

  expect_dfs_equal(
    base = expected,
    compare = derive_vars_crit_flag(
      select(expected, AVAL),
      condition = AVAL > 40,
      description = "AVAL > 40",
      create_numeric_flag = TRUE
    ),
    keys = "AVAL")
})

## Test 3: using values Y and N ----
test_that("derive_vars_crit_flag Test 3: using values Y and N", {
  expected <- tibble::tribble(
    ~AVAL, ~CRIT2FL,      ~CRIT2FLN,
    23,    "N",           0,
    42,    "Y",           1,
    NA,    NA_character_, NA_real_
  ) %>%
    mutate(CRIT2 = "AVAL > 40")

  expect_dfs_equal(
    base = expected,
    compare = derive_vars_crit_flag(
      select(expected, AVAL),
      crit_nr = 2,
      condition = AVAL > 40,
      description = "AVAL > 40",
      values_yn = TRUE,
      create_numeric_flag = TRUE
    ),
    keys = "AVAL")
})
