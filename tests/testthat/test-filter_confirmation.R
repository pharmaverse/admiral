library(tibble)
data <- tribble(
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

# filter_confirmation ----
## Test 1: filter without first_cond ----
test_that("filter_confirmation Test 1: filter without first_cond", {
  actual <-
    filter_confirmation(
      data,
      by_vars = vars(USUBJID),
      join_vars = vars(AVISITN, AVALC),
      join_type = "after",
      order = vars(AVISITN),
      filter = AVALC == "PR" & AVALC.join %in% c("CR", "PR") &
        AVISITN < AVISITN.join
    )

  expected <- tribble(
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
test_that("filter_confirmation Test 2: filter with first_cond", {
  actual <-
    filter_confirmation(
      data,
      by_vars = vars(USUBJID),
      join_vars = vars(AVALC),
      join_type = "after",
      first_cond = AVALC == "CR" &
        AVALC.join == "CR",
      order = vars(AVISITN),
      filter = TRUE
    )

  expected <- tribble(
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
test_that("filter_confirmation Test 3: filter with first_cond and summary function", {
  actual <-
    filter_confirmation(
      data,
      by_vars = vars(USUBJID),
      join_vars = vars(AVALC),
      join_type = "after",
      first_cond = AVALC == "PR" &
        AVALC.join %in% c("CR", "PR"),
      order = vars(AVISITN),
      filter = count_vals(AVALC.join, "SD") <= 1
    )

  expected <- tribble(
    ~USUBJID, ~AVISITN, ~AVALC,
    "1",      1,        "PR"
  )

  expect_dfs_equal(
    base = expected,
    compare = actual,
    keys = c("USUBJID", "AVISITN")
  )
})

## Test 4: join_type = "all" ----
test_that("filter_confirmation Test 4: join_type = 'all'", {
  adae <- tribble(
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

  actual <- filter_confirmation(
    adae,
    by_vars = vars(USUBJID),
    join_vars = vars(ACOVFL, ADY),
    join_type = "all",
    order = vars(ADY),
    filter = ADURN > 30 & ACOVFL.join == "Y" & ADY >= ADY.join - 7
  )

  expected <- tribble(
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

# min_cond ----
## Test 1: test it ----
test_that("min_cond, Test 1: test it", {
  data <- tribble(
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

  expected <- tribble(
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
## Test 1: test it ----
test_that("max_cond, Test 1: test it", {
  data <- tribble(
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

  expected <- tribble(
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
