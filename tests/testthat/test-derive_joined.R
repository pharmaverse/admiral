library(tibble)

# derive_vars_joined ----
test_that("no by_vars, no order", {
  expected <- tribble(
    ~USUBJID, ~ADY, ~AVISIT,    ~AWLO, ~AWHI,
    "1",        -2, "BASELINE",   -30,     1,
    "1",         3, "WEEK 1",       2,     7,
    "1",        24, "WEEK 4",      23,    30,
    "2",        NA, NA,            NA,    NA
  )

  windows <- tribble(
    ~AWLO, ~AWHI, ~AVISIT,
      -30,     1, "BASELINE",
        2,     7, "WEEK 1",
        8,    15, "WEEK 2",
       16,    22, "WEEK 3",
       23,    30, "WEEK 4"
  )

  expect_dfs_equal(
    base = expected,
    comp = derive_vars_joined(
      select(expected, USUBJID, ADY),
      dataset_add = windows,
      join_vars = vars(AWHI, AWLO),
      filter_join = AWLO <= ADY & ADY <= AWHI
    ),
    keys = c("USUBJID", "ADY")
  )
})
