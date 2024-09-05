library(testthat)
library(dplyr)
library(pharmaverseadam)

# Load the advs dataset
advs <- pharmaverseadam::advs

## Test 1: Basic functionality with advs dataset ----
test_that("derive_vars_cat Test 1: Basic functionality with advs dataset", {
  # Define the condition and categories
  definition <- exprs(
    ~condition, ~AVALCAT1, ~AVALCA1N,
    VSTEST == "Height" & AVAL >= 160, ">=160", 1,
    VSTEST == "Height" & AVAL < 160, "<160", 2
  )

  result <- derive_vars_cat(advs, definition) %>%
    select(USUBJID, VSTEST, AVAL, AVALCAT1, AVALCA1N) %>%
    filter(VSTEST == "Height")

  expect_snapshot(result)
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
  expect_snapshot_error(
    derive_vars_cat(list(1, 2, 3), definition)
  )
})

## Test 3: Error when definition is not an exprs object ----
test_that("derive_vars_cat Test 3: Error when definition is not an exprs object", {
  # Snapshot the error message
  expect_snapshot_error(
    derive_vars_cat(advs, list(condition = "Height", AVALCAT1 = 1))
  )
})

## Test 4: Error when required columns are missing from dataset ----
test_that("derive_vars_cat Test 4: Error when required columns are missing from dataset", {
  # Define the condition and categories (without VSTEST in the dataset)
  definition <- exprs(
    ~condition, ~AVALCAT1, ~AVALCA1N,
    VSTEST == "Height" & AVAL >= 160, ">=160", 1,
    VSTEST == "Height" & AVAL < 160, "<160", 2
  )

  # Remove VSTEST column from dataset
  advs_missing_col <- advs %>% select(-VSTEST)

  # Snapshot the error message
  expect_snapshot_error(
    derive_vars_cat(advs_missing_col, definition)
  )
})

## Test 5: Correct behavior when no conditions are met ----
test_that("derive_vars_cat Test 5: Correct behavior when no conditions are met", {
  # Define conditions that do not match any rows
  definition <- exprs(
    ~condition, ~AVALCAT1, ~AVALCA1N,
    VSTEST == "Height" & AVAL < 0, "<0", 1
  )

  result <- derive_vars_cat(advs, definition) %>%
    filter(VSTEST == "Height") %>%
    select(USUBJID, VSTEST, AVAL, AVALCAT1, AVALCA1N)

  expect_snapshot(result)
})

## Test 6: Overlapping conditions handled correctly ----
test_that("derive_vars_cat Test 6: Overlapping conditions handled correctly", {
  # Define overlapping conditions
  definition <- exprs(
    ~condition, ~AVALCAT1, ~AVALCA1N,
    VSTEST == "Height" & AVAL >= 160, ">=160", 1,
    VSTEST == "Height" & AVAL > 155, ">155", 2
  )

  result <- derive_vars_cat(advs, definition) %>%
    select(USUBJID, VSTEST, AVAL, AVALCAT1, AVALCA1N) %>%
    filter(VSTEST == "Height")

  expect_snapshot(result)
})

## Test 7: Handles missing values in dataset correctly ----
test_that("derive_vars_cat Test 7: Handles missing values in dataset correctly", {
  # Introduce missing values in AVAL
  advs_missing <- advs %>% filter(VSTEST == "Height")
  advs_missing$AVAL[1:5] <- NA

  # Define the condition and categories
  definition <- exprs(
    ~condition, ~AVALCAT1, ~AVALCA1N,
    VSTEST == "Height" & AVAL >= 160, ">=160", 1,
    VSTEST == "Height" & AVAL < 160, "<160", 2
  )

  result <- derive_vars_cat(advs_missing, definition) %>%
    select(USUBJID, VSTEST, AVAL, AVALCAT1, AVALCA1N)

  expect_snapshot(result)
})

## too slow, not sure if necessary.
# test_that("Large dataset is handled efficiently", {
#   # Create a large dataset by replicating advs
#   large_advs <- bind_rows(replicate(1000, advs, simplify = FALSE))
#
#   # Define the condition and categories
#   definition <- exprs(
#     ~condition, ~AVALCAT1, ~AVALCA1N,
#     VSTEST == "Height" & AVAL >= 160, ">=160", 1,
#     VSTEST == "Height" & AVAL < 160, "<160", 2
#   )
#
#   expect_snapshot(derive_vars_cat(large_advs, definition))
# })

## Test 8: Error when condition is missing from exprs() definition object ----
test_that("derive_vars_cat Test 8: Error when condition is missing from exprs() definition object", {
  # Define the condition but omit the 'condition' column from the definition
  definition <- exprs(
    ~AVALCAT1, ~AVALCA1N,
    ">=160", 1,
    "<160", 2
  )

  # Snapshot the error message
  expect_snapshot_error(
    derive_vars_cat(advs, definition)
  )
})

## Test 9: Conditions for multiple VSTESTs (Height and BILI) ----
test_that("derive_vars_cat Test 9: Conditions for multiple VSTESTs (Height and BILI)", {
  # Define conditions for two different VSTEST values: Height and BILI
  definition <- exprs(
    ~condition, ~AVALCAT1, ~AVALCA1N,
    VSTEST == "Height" & AVAL >= 160, "Height >= 160", 1,
    VSTEST == "Height" & AVAL < 160, "Height < 160", 2,
    VSTEST == "Weight" & AVAL >= 66.68, "Weight >= 66.68", 3,
    VSTEST == "Weight" & AVAL < 66.68, "Weight < 66.68", 4
  )

  result <- derive_vars_cat(advs, definition) %>%
    select(USUBJID, VSTEST, AVAL, AVALCAT1, AVALCA1N) %>%
    filter(VSTEST %in% c("Height", "Weight"))

  expect_snapshot(result)
})

## Test 10: Adding an extra variable (flag) to the dataset ----
test_that("derive_vars_cat Test 10: Adding an extra variable (flag) to the dataset", {
  # Define conditions and add a third variable (flag) that is TRUE or FALSE
  definition <- exprs(
    ~condition, ~AVALCAT1, ~AVALCA1N, ~extra_var,
    VSTEST == "Height" & AVAL >= 160, ">=160", 1, TRUE,
    VSTEST == "Height" & AVAL < 160, "<160", 2, FALSE
  )

  result <- derive_vars_cat(advs, definition) %>%
    select(USUBJID, VSTEST, AVAL, AVALCAT1, AVALCA1N, extra_var) %>%
    filter(VSTEST == "Height")

  expect_snapshot(result)
})
