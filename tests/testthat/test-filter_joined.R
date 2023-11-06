data <- tibble::tribble(
  ~USUBJID, ~AVISITN, ~AVALC,
  "1",      1,        "PR",
  "1",      2,        "CR",
  "1",      3,        "CR",
  "1",      4,        "SD",
  "1",      5,        "NE",
  "2",      1,        "SD",
  "2",      2,        "PR",
  "2",      3,        "PD",
  "3",      1,        "SD",
  "4",      1,        "PR",
  "4",      2,        "PD",
  "4",      3,        "SD",
  "4",      4,        "SD",
  "4",      5,        "PR"
)

# filter_joined ----
## Test 1: filter without first_cond_upper ----
test_that("filter_joined Test 1: filter without first_cond_upper", {
  actual <-
    filter_joined(
      data,
      dataset_add = data,
      by_vars = exprs(USUBJID),
      join_vars = exprs(AVISITN, AVALC),
      join_type = "after",
      order = exprs(AVISITN),
      filter_join = AVALC == "PR" & AVALC.join %in% c("CR", "PR") &
        AVISITN < AVISITN.join
    )

  expected <- tibble::tribble(
    ~USUBJID, ~AVISITN, ~AVALC,
    "1",      1,        "PR",
    "4",      1,        "PR"
  )

  expect_dfs_equal(
    base = expected,
    compare = actual,
    keys = c("USUBJID", "AVISITN")
  )
})

## Test 2: filter with first_cond ----
test_that("filter_joined Test 2: filter with first_cond", {
  actual <-
    filter_joined(
      data,
      dataset_add = data,
      by_vars = exprs(USUBJID),
      join_vars = exprs(AVALC),
      join_type = "after",
      first_cond_upper = AVALC == "CR" &
        AVALC.join == "CR",
      order = exprs(AVISITN),
      filter_join = TRUE
    )

  expected <- tibble::tribble(
    ~USUBJID, ~AVISITN, ~AVALC,
    "1",      2,        "CR"
  )

  expect_dfs_equal(
    base = expected,
    compare = actual,
    keys = c("USUBJID", "AVISITN")
  )
})

## Test 3: filter with first_cond and summary function ----
test_that("filter_joined Test 3: filter with first_cond and summary function", {
  actual <-
    filter_joined(
      data,
      dataset_add = data,
      by_vars = exprs(USUBJID),
      join_vars = exprs(AVALC),
      join_type = "after",
      first_cond_upper = AVALC == "PR" &
        AVALC.join %in% c("CR", "PR"),
      order = exprs(AVISITN),
      filter_join = count_vals(AVALC.join, "SD") <= 1
    )

  expected <- tibble::tribble(
    ~USUBJID, ~AVISITN, ~AVALC,
    "1",      1,        "PR"
  )

  expect_dfs_equal(
    base = expected,
    compare = actual,
    keys = c("USUBJID", "AVISITN")
  )
})

## Test 4: join_type = 'all' ----
test_that("filter_joined Test 4: join_type = 'all'", {
  adae <- tibble::tribble(
    ~USUBJID, ~ADY, ~ACOVFL, ~ADURN,
    "1",        10, "N",          1,
    "1",        21, "N",         50,
    "1",        23, "Y",         14,
    "1",        32, "N",         31,
    "1",        42, "N",         20,
    "2",        11, "Y",         13,
    "2",        23, "N",          2,
    "3",        13, "Y",         12,
    "4",        14, "N",         32,
    "4",        21, "N",         41
  )

  actual <- filter_joined(
    adae,
    dataset_add = adae,
    by_vars = exprs(USUBJID),
    join_vars = exprs(ACOVFL, ADY),
    join_type = "all",
    order = exprs(ADY),
    filter_join = ADURN > 30 & ACOVFL.join == "Y" & ADY >= ADY.join - 7
  )

  expected <- tibble::tribble(
    ~USUBJID, ~ADY, ~ACOVFL, ~ADURN,
    "1",        21, "N",         50,
    "1",        32, "N",         31,
  )

  expect_dfs_equal(
    base = expected,
    compare = actual,
    keys = c("USUBJID", "ADY")
  )
})

## Test 5: deprecation of `filter` ----
test_that("filter_joined Test 5: deprecation of `filter`", {
  expect_warning(
    actual <-
      filter_joined(
        data,
        dataset_add = data,
        by_vars = exprs(USUBJID),
        join_vars = exprs(AVISITN, AVALC),
        join_type = "after",
        order = exprs(AVISITN),
        filter = AVALC == "PR" & AVALC.join %in% c("CR", "PR") &
          AVISITN < AVISITN.join
      ),
    class = "lifecycle_warning_deprecated"
  )

  expected <- tibble::tribble(
    ~USUBJID, ~AVISITN, ~AVALC,
    "1",      1,        "PR",
    "4",      1,        "PR"
  )

  expect_dfs_equal(
    base = expected,
    compare = actual,
    keys = c("USUBJID", "AVISITN")
  )
})

## Test 6: deprecation of `first_cond` ----
test_that("filter_joined Test 6: deprecation of `first_cond`", {
  expect_warning(
    actual <-
      filter_joined(
        data,
        dataset_add = data,
        by_vars = exprs(USUBJID),
        join_vars = exprs(AVALC),
        join_type = "after",
        first_cond = AVALC == "CR" &
          AVALC.join == "CR",
        order = exprs(AVISITN),
        filter_join = TRUE
      ),
    class = "lifecycle_warning_deprecated"
  )

  expected <- tibble::tribble(
    ~USUBJID, ~AVISITN, ~AVALC,
    "1",      2,        "CR"
  )

  expect_dfs_equal(
    base = expected,
    compare = actual,
    keys = c("USUBJID", "AVISITN")
  )
})

# min_cond ----
## Test 7: minimum is derived correctly ----
test_that("min_cond Test 7: minimum is derived correctly", {
  data <- tibble::tribble(
    ~USUBJID, ~AVISITN, ~AVALC,
    "1",      1,        "PR",
    "1",      2,        "CR",
    "1",      3,        "NE",
    "1",      4,        "CR",
    "1",      5,        "NE",
    "2",      1,        "CR",
    "2",      2,        "PR",
    "2",      3,        "CR",
  )

  actual <- group_by(data, USUBJID) %>% mutate(
    first_cr_vis = min_cond(var = AVISITN, cond = AVALC == "CR")
  )

  expected <- tibble::tribble(
    ~USUBJID, ~AVISITN, ~AVALC, ~first_cr_vis,
    "1",      1,        "PR",               2,
    "1",      2,        "CR",               2,
    "1",      3,        "NE",               2,
    "1",      4,        "CR",               2,
    "1",      5,        "NE",               2,
    "2",      1,        "CR",               1,
    "2",      2,        "PR",               1,
    "2",      3,        "CR",               1,
  )

  expect_dfs_equal(
    base = expected,
    compare = actual,
    keys = c("USUBJID", "AVISITN")
  )
})

# max_cond ----
## Test 8: maximum is derived correctly ----
test_that("max_cond Test 8: maximum is derived correctly", {
  data <- tibble::tribble(
    ~USUBJID, ~AVISITN, ~AVALC,
    "1",      1,        "PR",
    "1",      2,        "CR",
    "1",      3,        "NE",
    "1",      4,        "CR",
    "1",      5,        "NE",
    "2",      1,        "CR",
    "2",      2,        "PR",
    "2",      3,        "CR",
  )

  actual <- group_by(data, USUBJID) %>% mutate(
    last_pr_vis = max_cond(var = AVISITN, cond = AVALC == "PR")
  )

  expected <- tibble::tribble(
    ~USUBJID, ~AVISITN, ~AVALC, ~last_pr_vis,
    "1",      1,        "PR",              1,
    "1",      2,        "CR",              1,
    "1",      3,        "NE",              1,
    "1",      4,        "CR",              1,
    "1",      5,        "NE",              1,
    "2",      1,        "CR",              2,
    "2",      2,        "PR",              2,
    "2",      3,        "CR",              2,
  )

  expect_dfs_equal(
    base = expected,
    compare = actual,
    keys = c("USUBJID", "AVISITN")
  )
})
