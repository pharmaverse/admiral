expected_result <- tibble::tribble(
  ~USUBJID,       ~VSTEST,  ~AVAL, ~AVALCAT1, ~AVALCA1N,
  "01-701-1015", "Height", 147.32,    "<160",         2,
  "01-701-1023", "Height", 162.56,   ">=160",         1,
  "01-701-1028", "Height",  177.8,   ">=160",         1,
  "01-701-1033", "Height", 175.26,   ">=160",         1,
  "01-701-1034", "Height",     NA,        NA,        NA,
  "01-701-1047", "Height",     NA,        NA,        NA,
  "01-701-1097", "Height", 168.91,   ">=160",         1,
  "01-701-1111", "Height", 158.24,    "<160",         2,
  "01-701-1115", "Height", 181.61,   ">=160",         1,
  "01-701-1118", "Height", 180.34,   ">=160",         1,
  "01-701-1015", "Weight",  53.98,        NA,        NA,
  "01-701-1023", "Weight",  78.47,        NA,        NA,
  "01-701-1028", "Weight",  98.88,        NA,        NA,
  "01-701-1033", "Weight",  88.45,        NA,        NA,
  "01-701-1034", "Weight",     NA,        NA,        NA,
  "01-701-1047", "Weight",     NA,        NA,        NA,
  "01-701-1097", "Weight",  78.02,        NA,        NA,
  "01-701-1111", "Weight",  60.33,        NA,        NA,
  "01-701-1115", "Weight",   78.7,        NA,        NA,
  "01-701-1118", "Weight",  71.67,        NA,        NA
)

advs <- expected_result %>% select(
  USUBJID, VSTEST, AVAL
)
## Test 1: Basic functionality without by_vars ----
test_that("derive_vars_cat Test 1: Basic functionality without by_vars", {
  # Define the condition and categories
  definition <- exprs(
    ~condition,                     ~AVALCAT1, ~AVALCA1N,
    VSTEST == "Height" & AVAL >= 160, ">=160",         1,
    VSTEST == "Height" & AVAL < 160,   "<160",         2
  )

  expect_dfs_equal(
    base =
      derive_vars_cat(
        advs,
        definition
      ),
    compare = expected_result,
    keys = c("USUBJID", "VSTEST")
  )
})

## Test 2: Basic functionality with by_vars ----
test_that("derive_vars_cat Test 2: Basic functionality with by_vars", {
  # Define the condition and categories
  definition <- exprs(
    ~VSTEST,   ~condition, ~AVALCAT1, ~AVALCA1N,
    "Height", AVAL >= 160,   ">=160",         1,
    "Height",  AVAL < 160,    "<160",         2
  )



  expect_dfs_equal(
    base =
      derive_vars_cat(
        advs,
        definition,
        by_vars = exprs(VSTEST)
      ),
    compare = expected_result,
    keys = c("USUBJID", "VSTEST")
  )
})

## Test 3: Forgot to specify by_vars ----
test_that("derive_vars_cat Test 3: Forgot to specify by_vars", {
  definition <- exprs(
    ~VSTEST,   ~condition, ~AVALCAT1, ~AVALCA1N,
    "Height", AVAL >= 160,   ">=160",         1,
    "Height",  AVAL < 160,    "<160",         2
  )

  expect_snapshot_warning(derive_vars_cat(advs, definition))
})

## Test 4: Error when dataset is not a dataframe ----
test_that("derive_vars_cat Test 4: Error when dataset is not a dataframe", {
  # Define the condition and categories
  definition <- exprs(
    ~condition,                     ~AVALCAT1, ~AVALCA1N,
    VSTEST == "Height" & AVAL >= 160, ">=160",         1,
    VSTEST == "Height" & AVAL < 160,   "<160",         2
  )

  # Snapshot the error message
  expect_error(
    derive_vars_cat(list(1, 2, 3), definition),
    class = "assert_data_frame"
  )
})

## Test 5: Error when definition is not an exprs object ----
test_that("derive_vars_cat Test 5: Error when definition is not an exprs object", {
  definition <- tribble(
    ~condition,  ~AVALCAT1, ~AVALCA1N,
    "AVAL >= 160", ">=160",         1,
    "AVAL < 160",   "<160",         2
  )
  # Snapshot the error message
  expect_error(
    derive_vars_cat(advs, definition),
    class = "assert_expr_list"
  )
})

## Test 6: Error when required columns are missing from dataset ----
test_that("derive_vars_cat Test 6: Error when required columns are missing from dataset", {
  # Define the condition and categories (without VSTEST in the dataset)
  definition <- exprs(
    ~VSTEST,   ~condition, ~AVALCAT1, ~AVALCA1N,
    "Height", AVAL >= 160,   ">=160",         1,
    "Height",  AVAL < 160,    "<160",         2
  )

  # Remove VSTEST column from dataset
  advs_missing_col <- advs %>% select(-VSTEST)

  # Snapshot the error message
  expect_error(
    derive_vars_cat(advs_missing_col, definition, by_vars = exprs(VSTEST)),
    class = "assert_data_frame"
  )
})

## Test 7: Correct behavior when no conditions are met ----
test_that("derive_vars_cat Test 7: Correct behavior when no conditions are met", {
  # Define conditions that do not match any rows
  definition <- exprs(
    ~condition,                    ~AVALCAT1, ~AVALCA1N,
    VSTEST == "Height" & AVAL < 0,      "<0",         1
  )

  expected_result <- tibble::tribble(
    ~USUBJID,       ~VSTEST,  ~AVAL,     ~AVALCAT1, ~AVALCA1N,
    "01-701-1015", "Height", 147.32, NA_character_,  NA_real_,
    "01-701-1023", "Height", 162.56, NA_character_,  NA_real_,
    "01-701-1028", "Height",  177.8, NA_character_,  NA_real_,
    "01-701-1033", "Height", 175.26, NA_character_,  NA_real_,
    "01-701-1034", "Height",     NA, NA_character_,  NA_real_,
    "01-701-1047", "Height",     NA, NA_character_,  NA_real_,
    "01-701-1097", "Height", 168.91, NA_character_,  NA_real_,
    "01-701-1111", "Height", 158.24, NA_character_,  NA_real_,
    "01-701-1115", "Height", 181.61, NA_character_,  NA_real_,
    "01-701-1118", "Height", 180.34, NA_character_,  NA_real_,
    "01-701-1015", "Weight",  53.98, NA_character_,  NA_real_,
    "01-701-1023", "Weight",  78.47, NA_character_,  NA_real_,
    "01-701-1028", "Weight",  98.88, NA_character_,  NA_real_,
    "01-701-1033", "Weight",  88.45, NA_character_,  NA_real_,
    "01-701-1034", "Weight",     NA, NA_character_,  NA_real_,
    "01-701-1047", "Weight",     NA, NA_character_,  NA_real_,
    "01-701-1097", "Weight",  78.02, NA_character_,  NA_real_,
    "01-701-1111", "Weight",  60.33, NA_character_,  NA_real_,
    "01-701-1115", "Weight",   78.7, NA_character_,  NA_real_,
    "01-701-1118", "Weight",  71.67, NA_character_,  NA_real_
  )

  expect_dfs_equal(
    base = derive_vars_cat(advs, definition),
    compare = expected_result,
    keys = c("USUBJID", "VSTEST")
  )
})

## Test 8: Overlapping conditions handled correctly ----
test_that("derive_vars_cat Test 8: Overlapping conditions handled correctly", {
  # Define overlapping conditions
  definition <- exprs(
    ~VSTEST,   ~condition, ~AVALCAT1, ~AVALCA1N,
    "Height",  AVAL < 160,    "<160",         3,
    "Height",  AVAL < 170,    "<170",         2,
    "Height", AVAL >= 170,   ">=170",         1
  )

  expected_result <- tibble::tribble(
    ~USUBJID,       ~VSTEST,  ~AVAL, ~AVALCAT1, ~AVALCA1N,
    "01-701-1015", "Height", 147.32,    "<160",         3,
    "01-701-1023", "Height", 162.56,    "<170",         2,
    "01-701-1028", "Height",  177.8,   ">=170",         1,
    "01-701-1033", "Height", 175.26,   ">=170",         1,
    "01-701-1034", "Height",     NA,        NA,        NA,
    "01-701-1047", "Height",     NA,        NA,        NA,
    "01-701-1097", "Height", 168.91,    "<170",         2,
    "01-701-1111", "Height", 158.24,    "<160",         3,
    "01-701-1115", "Height", 181.61,   ">=170",         1,
    "01-701-1118", "Height", 180.34,   ">=170",         1,
    "01-701-1015", "Weight",  53.98,        NA,        NA,
    "01-701-1023", "Weight",  78.47,        NA,        NA,
    "01-701-1028", "Weight",  98.88,        NA,        NA,
    "01-701-1033", "Weight",  88.45,        NA,        NA,
    "01-701-1034", "Weight",     NA,        NA,        NA,
    "01-701-1047", "Weight",     NA,        NA,        NA,
    "01-701-1097", "Weight",  78.02,        NA,        NA,
    "01-701-1111", "Weight",  60.33,        NA,        NA,
    "01-701-1115", "Weight",   78.7,        NA,        NA,
    "01-701-1118", "Weight",  71.67,        NA,        NA
  )

  expect_dfs_equal(
    base = derive_vars_cat(advs, definition, by_vars = exprs(VSTEST)),
    compare = expected_result,
    keys = c("USUBJID", "VSTEST")
  )
})


## Test 9: Error when condition is missing from `definition` ----
test_that("derive_vars_cat Test 9: Error when condition is missing from `definition`", {
  # Define the condition but omit the 'condition' column from the definition
  definition <- exprs(
    ~AVALCAT1, ~AVALCA1N,
    ">=160",           1,
    "<160",            2
  )

  # Snapshot the error message
  expect_error(
    derive_vars_cat(advs, definition),
    class = "assert_data_frame"
  )
})

## Test 10: Conditions for multiple VSTESTs (Height and Weight) ----
test_that("derive_vars_cat Test 10: Conditions for multiple VSTESTs (Height and Weight)", {
  # Define conditions for two different VSTEST values: Height and BILI
  definition <- exprs(
    ~VSTEST,     ~condition,         ~AVALCAT1, ~AVALCA1N,
    "Height",   AVAL >= 160,   "Height >= 160",         1,
    "Height",    AVAL < 160,    "Height < 160",         2,
    "Weight", AVAL >= 66.68, "Weight >= 66.68",         1,
    "Weight",  AVAL < 66.68,  "Weight < 66.68",         2
  )

  expected_result <- tibble::tribble(
    ~USUBJID,       ~VSTEST,  ~AVAL,         ~AVALCAT1, ~AVALCA1N,
    "01-701-1015", "Height", 147.32,    "Height < 160",         2,
    "01-701-1023", "Height", 162.56,   "Height >= 160",         1,
    "01-701-1028", "Height",  177.8,   "Height >= 160",         1,
    "01-701-1033", "Height", 175.26,   "Height >= 160",         1,
    "01-701-1034", "Height",     NA,                NA,        NA,
    "01-701-1047", "Height",     NA,                NA,        NA,
    "01-701-1097", "Height", 168.91,   "Height >= 160",         1,
    "01-701-1111", "Height", 158.24,    "Height < 160",         2,
    "01-701-1115", "Height", 181.61,   "Height >= 160",         1,
    "01-701-1118", "Height", 180.34,   "Height >= 160",         1,
    "01-701-1015", "Weight",  53.98,  "Weight < 66.68",         2,
    "01-701-1023", "Weight",  78.47, "Weight >= 66.68",         1,
    "01-701-1028", "Weight",  98.88, "Weight >= 66.68",         1,
    "01-701-1033", "Weight",  88.45, "Weight >= 66.68",         1,
    "01-701-1034", "Weight",     NA,                NA,        NA,
    "01-701-1047", "Weight",     NA,                NA,        NA,
    "01-701-1097", "Weight",  78.02, "Weight >= 66.68",         1,
    "01-701-1111", "Weight",  60.33,  "Weight < 66.68",         2,
    "01-701-1115", "Weight",   78.7, "Weight >= 66.68",         1,
    "01-701-1118", "Weight",  71.67, "Weight >= 66.68",         1
  )
  expect_dfs_equal(
    base = derive_vars_cat(advs, definition, by_vars = exprs(VSTEST)),
    compare = expected_result,
    keys = c("USUBJID", "VSTEST")
  )
})

## Test 11: Adding an extra variable (flag) to the dataset ----
test_that("derive_vars_cat Test 11: Adding an extra variable (flag) to the dataset", {
  # Define conditions and add a third variable (flag) that is TRUE or FALSE
  definition <- exprs(
    ~VSTEST,   ~condition, ~AVALCAT1, ~AVALCA1N, ~extra_var,
    "Height", AVAL >= 160,   ">=160",         1,       TRUE,
    "Height",  AVAL < 160,    "<160",         2,      FALSE
  )

  expected_result <- tibble::tribble(
    ~USUBJID,       ~VSTEST,  ~AVAL, ~AVALCAT1, ~AVALCA1N, ~extra_var,
    "01-701-1015", "Height", 147.32,    "<160",         2,      FALSE,
    "01-701-1023", "Height", 162.56,   ">=160",         1,       TRUE,
    "01-701-1028", "Height",  177.8,   ">=160",         1,       TRUE,
    "01-701-1033", "Height", 175.26,   ">=160",         1,       TRUE,
    "01-701-1034", "Height",     NA,        NA,        NA,         NA,
    "01-701-1047", "Height",     NA,        NA,        NA,         NA,
    "01-701-1097", "Height", 168.91,   ">=160",         1,       TRUE,
    "01-701-1111", "Height", 158.24,    "<160",         2,      FALSE,
    "01-701-1115", "Height", 181.61,   ">=160",         1,       TRUE,
    "01-701-1118", "Height", 180.34,   ">=160",         1,       TRUE,
    "01-701-1015", "Weight",  53.98,        NA,        NA,         NA,
    "01-701-1023", "Weight",  78.47,        NA,        NA,         NA,
    "01-701-1028", "Weight",  98.88,        NA,        NA,         NA,
    "01-701-1033", "Weight",  88.45,        NA,        NA,         NA,
    "01-701-1034", "Weight",     NA,        NA,        NA,         NA,
    "01-701-1047", "Weight",     NA,        NA,        NA,         NA,
    "01-701-1097", "Weight",  78.02,        NA,        NA,         NA,
    "01-701-1111", "Weight",  60.33,        NA,        NA,         NA,
    "01-701-1115", "Weight",   78.7,        NA,        NA,         NA,
    "01-701-1118", "Weight",  71.67,        NA,        NA,         NA
  )
  expect_dfs_equal(
    base = derive_vars_cat(advs, definition, by_vars = exprs(VSTEST)),
    compare = expected_result,
    keys = c("USUBJID", "VSTEST")
  )
})

## Test 12: Wrong input for by_vars ----
test_that("derive_vars_cat Test 12: Wrong input for by_vars", {
  # Define conditions
  definition <- exprs(
    ~VSTEST,   ~condition, ~AVALCAT1, ~AVALCA1N,
    "Height", AVAL >= 160,   ">=160",         1,
    "Height",  AVAL < 160,    "<160",         2
  )

  expect_error(derive_vars_cat(advs, definition, by_vars = exprs(VSTEST == "Height")),
    class = "assert_vars"
  )
})

## Test 13: definition has wrong shape ----
test_that("derive_vars_cat Test 13: definition has wrong shape", {
  # Define conditions
  definition_wrong_shape <- exprs(
    ~VSTEST,   ~condition, ~AVALCAT1, ~AVALCA1N,
    "Height", AVAL >= 160,   ">=160",         1,
    "Height",  AVAL < 160,    "<160"
  )

  expect_snapshot_error(derive_vars_cat(advs, definition_wrong_shape, by_vars = exprs(VSTEST)))
})

## Test 14: two by_vars variables ----
test_that("derive_vars_cat Test 14: two by_vars variables", {
  # Define conditions
  definition <- exprs(
    ~VISIT,     ~VSTEST,  ~condition, ~AVALCAT1, ~AVALCA1N,
    "Week 24", "Height", AVAL >= 160,   ">=160",         1,
    "Week 24", "Height",  AVAL < 160,    "<160",         2,
  )

  advs_visit <- advs %>% mutate(
    VISIT = "Week 24"
  )

  expect_snapshot_error(derive_vars_cat(advs_visit, definition, by_vars = exprs(VSTEST, VISIT)))
})
