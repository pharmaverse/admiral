data <- tibble::tribble(
  ~USUBJID, ~AVISITN, ~AVALC,
  "1",      1,        "PR",
  "1",      2,        "CR",
  "1",      3,        "CR",
  "1",      4,        "SD",
  "1",      5,        "NE",
  "2",      1,        "SD",
  "2",      2,        "PD",
  "2",      3,        "PD",
  "3",      1,        "SD",
  "4",      1,        "PR",
  "4",      2,        "PD",
  "4",      3,        "SD",
  "4",      4,        "SD",
  "4",      5,        "PR"
)

# filter_relative ----
## filter_relative Test 1: mode = first, selection = before, inclusive = TRUE ----
test_that("filter_relative Test 1: mode = first, selection = before, inclusive = TRUE", {
  actual <- filter_relative(
    data,
    by_vars = vars(USUBJID),
    order = vars(AVISITN),
    condition = AVALC == "PD",
    mode = "first",
    selection = "before",
    inclusive = TRUE
  )

  expected <- tibble::tribble(
    ~USUBJID, ~AVISITN, ~AVALC,
    "1",      1,        "PR",
    "1",      2,        "CR",
    "1",      3,        "CR",
    "1",      4,        "SD",
    "1",      5,        "NE",
    "2",      1,        "SD",
    "2",      2,        "PD",
    "3",      1,        "SD",
    "4",      1,        "PR",
    "4",      2,        "PD"
  )

  expect_dfs_equal(
    base = expected,
    comp = actual,
    keys = c("USUBJID", "AVISITN")
  )
})

## filter_relative Test 2: mode = first, selection = before, inclusive = FALSE ----
test_that("filter_relative Test 2: mode = first, selection = before, inclusive = FALSE", {
  actual <- filter_relative(
    data,
    by_vars = vars(USUBJID),
    order = vars(AVISITN),
    condition = AVALC == "PD",
    mode = "first",
    selection = "before",
    inclusive = FALSE,
    keep_no_ref_group = FALSE
  )

  expected <- tibble::tribble(
    ~USUBJID, ~AVISITN, ~AVALC,
    "2",      1,        "SD",
    "4",      1,        "PR"
  )

  expect_dfs_equal(
    base = expected,
    comp = actual,
    keys = c("USUBJID", "AVISITN")
  )
})

## filter_relative Test 3: mode = first, selection = after, inclusive = TRUE ----
test_that("filter_relative Test 3: mode = first, selection = after, inclusive = TRUE", {
  actual <- filter_relative(
    data,
    by_vars = vars(USUBJID),
    order = vars(AVISITN),
    condition = AVALC == "PD",
    mode = "first",
    selection = "after",
    inclusive = TRUE,
    keep_no_ref_group = FALSE
  )

  expected <- tibble::tribble(
    ~USUBJID, ~AVISITN, ~AVALC,
    "2",      2,        "PD",
    "2",      3,        "PD",
    "4",      2,        "PD",
    "4",      3,        "SD",
    "4",      4,        "SD",
    "4",      5,        "PR"
  )

  expect_dfs_equal(
    base = expected,
    comp = actual,
    keys = c("USUBJID", "AVISITN")
  )
})

## filter_relative Test 4: mode = first, selection = after, inclusive = FALSE ----
test_that("filter_relative Test 4: mode = first, selection = after, inclusive = FALSE", {
  actual <- filter_relative(
    data,
    by_vars = vars(USUBJID),
    order = vars(AVISITN),
    condition = AVALC == "PD",
    mode = "first",
    selection = "after",
    inclusive = FALSE
  )

  expected <- tibble::tribble(
    ~USUBJID, ~AVISITN, ~AVALC,
    "1",      1,        "PR",
    "1",      2,        "CR",
    "1",      3,        "CR",
    "1",      4,        "SD",
    "1",      5,        "NE",
    "2",      3,        "PD",
    "3",      1,        "SD",
    "4",      3,        "SD",
    "4",      4,        "SD",
    "4",      5,        "PR"
  )

  expect_dfs_equal(
    base = expected,
    comp = actual,
    keys = c("USUBJID", "AVISITN")
  )
})

## filter_relative Test 5: mode = last, selection = before, inclusive = TRUE ----
test_that("filter_relative Test 1: mode = last, selection = before, inclusive = TRUE", {
  actual <- filter_relative(
    data,
    by_vars = vars(USUBJID),
    order = vars(AVISITN),
    condition = AVALC == "SD",
    mode = "last",
    selection = "before",
    inclusive = TRUE,
    keep_no_ref_group = FALSE
  )

  expected <- tibble::tribble(
    ~USUBJID, ~AVISITN, ~AVALC,
    "1",      1,        "PR",
    "1",      2,        "CR",
    "1",      3,        "CR",
    "1",      4,        "SD",
    "2",      1,        "SD",
    "3",      1,        "SD",
    "4",      1,        "PR",
    "4",      2,        "PD",
    "4",      3,        "SD",
    "4",      4,        "SD",
  )

  expect_dfs_equal(
    base = expected,
    comp = actual,
    keys = c("USUBJID", "AVISITN")
  )
})

## filter_relative Test 6: mode = last, selection = before, inclusive = FALSE ----
test_that("filter_relative Test 2: mode = last, selection = before, inclusive = FALSE", {
  actual <- filter_relative(
    data,
    by_vars = vars(USUBJID),
    order = vars(AVISITN),
    condition = AVALC == "SD",
    mode = "last",
    selection = "before",
    inclusive = FALSE,
    keep_no_ref_group = FALSE
  )

  expected <- tibble::tribble(
    ~USUBJID, ~AVISITN, ~AVALC,
    "1",      1,        "PR",
    "1",      2,        "CR",
    "1",      3,        "CR",
    "4",      1,        "PR",
    "4",      2,        "PD",
    "4",      3,        "SD",
  )

  expect_dfs_equal(
    base = expected,
    comp = actual,
    keys = c("USUBJID", "AVISITN")
  )
})

## filter_relative Test 7: mode = last, selection = after, inclusive = TRUE ----
test_that("filter_relative Test 7: mode = last, selection = after, inclusive = TRUE", {
  actual <- filter_relative(
    data,
    by_vars = vars(USUBJID),
    order = vars(AVISITN),
    condition = AVALC == "SD",
    mode = "last",
    selection = "after",
    inclusive = TRUE,
    keep_no_ref_group = FALSE
  )

  expected <- tibble::tribble(
    ~USUBJID, ~AVISITN, ~AVALC,
    "1",      4,        "SD",
    "1",      5,        "NE",
    "2",      1,        "SD",
    "2",      2,        "PD",
    "2",      3,        "PD",
    "3",      1,        "SD",
    "4",      4,        "SD",
    "4",      5,        "PR"
  )

  expect_dfs_equal(
    base = expected,
    comp = actual,
    keys = c("USUBJID", "AVISITN")
  )
})

## filter_relative Test 8: mode = last, selection = after, inclusive = FALSE ----
test_that("filter_relative Test 8: mode = last, selection = after, inclusive = FALSE", {
  actual <- filter_relative(
    data,
    by_vars = vars(USUBJID),
    order = vars(AVISITN),
    condition = AVALC == "SD",
    mode = "last",
    selection = "after",
    inclusive = FALSE,
    keep_no_ref_group = FALSE
  )

  expected <- tibble::tribble(
    ~USUBJID, ~AVISITN, ~AVALC,
    "1",      5,        "NE",
    "2",      2,        "PD",
    "2",      3,        "PD",
    "4",      5,        "PR"
  )

  expect_dfs_equal(
    base = expected,
    comp = actual,
    keys = c("USUBJID", "AVISITN")
  )
})
