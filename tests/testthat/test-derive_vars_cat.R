# Load the advs dataset
advs <- tibble::tribble(
  ~USUBJID, ~VSTEST, ~AVAL,
  "01-701-1015", "Height", 147.32,
  "01-701-1015", "Weight", 53.98,
  "01-701-1023", "Height", 162.56,
  "01-701-1023", "Weight", 78.47,
  "01-701-1028", "Height", 177.8,
  "01-701-1028", "Weight", 98.88,
  "01-701-1033", "Height", 175.26,
  "01-701-1033", "Weight", 88.45,
  "01-701-1034", "Height", NA,
  "01-701-1034", "Weight", NA,
  "01-701-1047", "Height", NA,
  "01-701-1047", "Weight", NA,
  "01-701-1097", "Height", 168.91,
  "01-701-1097", "Weight", 78.02,
  "01-701-1111", "Height", 158.24,
  "01-701-1111", "Weight", 60.33,
  "01-701-1115", "Height", 181.61,
  "01-701-1115", "Weight", 78.7,
  "01-701-1118", "Height", 180.34,
  "01-701-1118", "Weight", 71.67
) %>% arrange(VSTEST)

## Test 1: Basic functionality with advs dataset ----
test_that("derive_vars_cat Test 1: Basic functionality with advs dataset", {
  # Define the condition and categories
  definition <- exprs(
    ~condition, ~AVALCAT1, ~AVALCA1N,
    VSTEST == "Height" & AVAL >= 160, ">=160", 1,
    VSTEST == "Height" & AVAL < 160, "<160", 2
  )

  result <- derive_vars_cat(advs, definition)

  # using by_vars
  definition2 <- exprs(
    ~VSTEST, ~condition, ~AVALCAT1, ~AVALCA1N,
    "Height", AVAL >= 160, ">=160", 1,
    "Height", AVAL < 160, "<160", 2
  )

  result2 <- derive_vars_cat(advs, definition2, by_vars = exprs(VSTEST))

  expect_snapshot(result)
  expect_snapshot(result2)
  expect_dfs_equal(base = result, compare = result2, keys = c("USUBJID", "VSTEST"))
})

## Test 2: Error when dataset is not a dataframe ----
test_that("derive_vars_cat Test 2: Error when dataset is not a dataframe", {
  # Define the condition and categories
  definition <- exprs(
    ~condition, ~AVALCAT1, ~AVALCA1N,
    VSTEST == "Height" & AVAL >= 160, ">=160", 1,
    VSTEST == "Height" & AVAL < 160, "<160", 2
  )

  # Snapshot the error message
  expect_error(
    derive_vars_cat(list(1, 2, 3), definition),
    class = "assert_data_frame"
  )
})

## Test 3: Error when definition is not an exprs object ----
test_that("derive_vars_cat Test 3: Error when definition is not an exprs object", {
  definition <- tribble(
    ~condition, ~AVALCAT1, ~AVALCA1N,
    "AVAL >= 160", ">=160", 1,
    "AVAL < 160", "<160", 2
  )
  # Snapshot the error message
  expect_snapshot_error(
    derive_vars_cat(advs, definition)
  )
})

## Test 4: Error when required columns are missing from dataset ----
test_that("derive_vars_cat Test 4: Error when required columns are missing from dataset", {
  # Define the condition and categories (without VSTEST in the dataset)
  definition <- exprs(
    ~VSTEST, ~condition, ~AVALCAT1, ~AVALCA1N,
    "Height", AVAL >= 160, ">=160", 1,
    "Height", AVAL < 160, "<160", 2
  )

  # Remove VSTEST column from dataset
  advs_missing_col <- advs %>% select(-VSTEST)

  # Snapshot the error message
  expect_error(
    derive_vars_cat(advs_missing_col, definition, by_vars = exprs(VSTEST)),
    class = "assert_data_frame"
  )
})

## Test 5: Correct behavior when no conditions are met ----
test_that("derive_vars_cat Test 5: Correct behavior when no conditions are met", {
  # Define conditions that do not match any rows
  definition <- exprs(
    ~condition, ~AVALCAT1, ~AVALCA1N,
    VSTEST == "Height" & AVAL < 0, "<0", 1
  )

  result <- derive_vars_cat(advs, definition)

  expect_snapshot(result)
})

## Test 6: Overlapping conditions handled correctly ----
test_that("derive_vars_cat Test 6: Overlapping conditions handled correctly", {
  # Define overlapping conditions
  definition <- exprs(
    ~VSTEST, ~condition, ~AVALCAT1, ~AVALCA1N,
    "Height", AVAL < 160, "<160", 3,
    "Height", AVAL < 170, "<170", 2,
    "Height", AVAL >= 170, ">=170", 1
  )

  result <- derive_vars_cat(advs, definition, by_vars = exprs(VSTEST))

  expect_snapshot(result)
})


## Test 7: Error when condition is missing from `definition` ----
test_that("derive_vars_cat Test 7: Error when condition is missing from `definition`", {
  # Define the condition but omit the 'condition' column from the definition
  definition <- exprs(
    ~AVALCAT1, ~AVALCA1N,
    ">=160", 1,
    "<160", 2
  )

  # Snapshot the error message
  expect_error(
    derive_vars_cat(advs, definition),
    class = "assert_data_frame"
  )
})

## Test 8: Conditions for multiple VSTESTs (Height and Weight) ----
test_that("derive_vars_cat Test 8: Conditions for multiple VSTESTs (Height and Weight)", {
  # Define conditions for two different VSTEST values: Height and BILI
  definition <- exprs(
    ~VSTEST, ~condition, ~AVALCAT1, ~AVALCA1N,
    "Height", AVAL >= 160, "Height >= 160", 1,
    "Height", AVAL < 160, "Height < 160", 2,
    "Weight", AVAL >= 66.68, "Weight >= 66.68", 1,
    "Weight", AVAL < 66.68, "Weight < 66.68", 2
  )

  result <- derive_vars_cat(advs, definition, by_vars = exprs(VSTEST))

  expect_snapshot(result)
})

## Test 9: Adding an extra variable (flag) to the dataset ----
test_that("derive_vars_cat Test 9: Adding an extra variable (flag) to the dataset", {
  # Define conditions and add a third variable (flag) that is TRUE or FALSE
  definition <- exprs(
    ~VSTEST, ~condition, ~AVALCAT1, ~AVALCA1N, ~extra_var,
    "Height", AVAL >= 160, ">=160", 1, TRUE,
    "Height", AVAL < 160, "<160", 2, FALSE
  )

  result <- derive_vars_cat(advs, definition, by_vars = exprs(VSTEST))

  expect_snapshot(result)
})

## Test 10: Wrong input for by_vars ----
test_that("derive_vars_cat Test 10: Wrong input for by_vars", {
  # Define conditions
  definition <- exprs(
    ~VSTEST, ~condition, ~AVALCAT1, ~AVALCA1N,
    "Height", AVAL >= 160, ">=160", 1,
    "Height", AVAL < 160, "<160", 2
  )

  expect_error(derive_vars_cat(advs, definition, by_vars = exprs(VSTEST == "Height")),
               class = "assert_vars")
})

## Test 11: definition has wrong shape ----
test_that("derive_vars_cat Test 11: definition has wrong shape", {
  # Define conditions
  definition_wrong_shape <- exprs(
    ~VSTEST, ~condition, ~AVALCAT1, ~AVALCA1N,
    "Height", AVAL >= 160, ">=160", 1,
    "Height", AVAL < 160, "<160"
  )

  expect_snapshot_error(derive_vars_cat(advs, definition_wrong_shape, by_vars = exprs(VSTEST)))
})
