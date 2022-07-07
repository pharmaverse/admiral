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
## filter_confirmation Test 1: filter without first_cond ----
test_that("filter_confirmation Test 1: filter without first_cond", {
  actual <-
    filter_confirmation(
      data,
      by_vars = vars(USUBJID),
      join_vars = vars(AVISITN, AVALC),
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

## filter_confirmation Test 2: filter with first_cond ----
test_that("filter_confirmation Test 2: filter with first_cond", {
  actual <-
    filter_confirmation(
      data,
      by_vars = vars(USUBJID),
      join_vars = vars(AVALC),
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

## filter_confirmation Test 3: filter with first_cond and summary function ----
test_that("filter_confirmation Test 3: filter with first_cond and summary function", {
  actual <-
    filter_confirmation(
      data,
      by_vars = vars(USUBJID),
      join_vars = vars(AVALC),
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
