## Test 1: add variables ----
test_that("process_set_values_to Test 1: add variables", {
  bds <- tibble::tribble(
    ~USUBJID, ~AVAL,
    "1",         20,
    "2",         35
  )
  expected <- bds %>%
    mutate(
      PARAMCD = "BMI",
      PARAM = "Body-Mass-Index",
      PARAMN = 1
    )
  expect_dfs_equal(
    base = expected,
    compare = process_set_values_to(
      bds,
      set_values_to = exprs(
        PARAMCD = "BMI",
        PARAM = "Body-Mass-Index",
        PARAMN = 1
      ),
      expected_types = c(
        PARAMCD = "character",
        PARAM = "character",
        PARAMN = "numeric"
      )
    ),
    keys = c("USUBJID")
  )
})

## Test 2: catch error ----
test_that("process_set_values_to Test 2: catch error", {
  bds <- tibble::tribble(
    ~USUBJID, ~AVAL,
    "1",         20,
    "2",         35
  )

  expect_error(
    process_set_values_to(
      bds,
      set_values_to = exprs(
        PARAMCD = BMI,
        PARAM = "Body-Mass-Index",
        PARAMN = 1
      )
    ),
    "Assigning variables failed!\nset_values_to = (",
    fixed = TRUE
  )
})

## Test 3: check types ----
test_that("process_set_values_to Test 3: check types", {
  bds <- tibble::tribble(
    ~USUBJID, ~AVAL,
    "1",         20,
    "2",         35
  )

  expect_error(
    process_set_values_to(
      bds,
      set_values_to = exprs(
        PARAMCD = 1,
        PARAM = "Body-Mass-Index",
        PARAMN = "BMI"
      ),
      expected_types = c(
        PARAMCD = "character",
        PARAM = "character",
        PARAMN = "numeric"
      )
    ),
    paste(
      "The following variables have an unexpected type:",
      "PARAMCD: expected: character, actual: numeric",
      "PARAMN: expected: numeric, actual: character",
      sep = "\n"
    ),
    fixed = TRUE
  )
})
