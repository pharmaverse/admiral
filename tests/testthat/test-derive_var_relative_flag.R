
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
      by_vars = exprs(USUBJID),
      order = exprs(AVISITN),
      new_var = ANL02FL,
      condition = AVALC == "PD",
      mode = "first",
      selection = "before",
      inclusive = TRUE
    ),
    keys = c("USUBJID", "AVISITN")
  )
})

## Test 2: Flag AEs after COVID AE ----
test_that("derive_var_relative_flag Test 2: Flag AEs after COVID AE", {
  expected <- tibble::tribble(
    ~USUBJID, ~ASTDY, ~ACOVFL, ~AESEQ, ~PSTCOVFL,
    "1",           2, NA,           1, NA,
    "1",           5, "Y",          2, NA,
    "1",           5, NA,           3, "Y",
    "1",          17, NA,           4, "Y",
    "1",          27, "Y",          5, "Y",
    "1",          32, NA,           6, "Y",
    "2",           8, NA,           1, NA,
    "2",          11, NA,           2, NA
  )

  adae <- select(expected, -PSTCOVFL)

  expect_dfs_equal(
    base = expected,
    compare = derive_var_relative_flag(
      adae,
      by_vars = exprs(USUBJID),
      order = exprs(ASTDY, AESEQ),
      new_var = PSTCOVFL,
      condition = ACOVFL == "Y",
      mode = "first",
      selection = "after",
      inclusive = FALSE,
      flag_no_ref_groups = FALSE
    ),
    keys = c("USUBJID", "AESEQ")
  )
})
