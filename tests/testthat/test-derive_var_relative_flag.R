
## Test 1: flag observations up to first PD ----
test_that("derive_var_relative_flag Test 1: flag observations up to first PD", {
  expected <- tibble::tribble(
    ~USUBJID, ~AVISITN, ~AVALC, ~ANL02FL,
    "1",      0,        "PR",   "Y",
    "1",      1,        "CR",   "Y",
    "1",      2,        "CR",   "Y",
    "1",      3,        "SD",   "Y",
    "1",      4,        "NE",   "Y",
    "2",      0,        "SD",   "Y",
    "2",      1,        "PD",   "Y",
    "2",      2,        "PD",   NA,
    "3",      0,        "SD",   "Y",
    "4",      0,        "SD",   "Y",
    "4",      1,        "PR",   "Y",
    "4",      2,        "PD",   "Y",
    "4",      3,        "SD",   NA,
    "4",      4,        "PR",   NA
  )

  response <- select(expected, -ANL02FL)

  expect_dfs_equal(
    base = expected,
    compare = derive_var_relative_flag(
      response,
      by_vars = vars(USUBJID),
      order = vars(AVISITN),
      new_var = ANL02FL,
      condition = AVALC == "PD",
      mode = "first",
      selection = "before",
      inclusive = TRUE
    ),
    keys = c("USUBJID", "AVISITN")
  )
})
