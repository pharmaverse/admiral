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

# derive_var_confirmation_flag ----
## Test 1: filter without first_cond ----
## Flagging any patient PR value that is followed subsequently by a CR or PR
test_that("derive_var_confirmation_flag Test 1: filter without first_cond", {
  actual <-
    derive_var_confirmation_flag(
      data,
      new_var = CONFFL,
      by_vars = vars(USUBJID),
      join_vars = vars(AVISITN, AVALC),
      join_type = "after",
      order = vars(AVISITN),
      filter = AVALC == "PR" & AVALC.join %in% c("CR", "PR")
    )

  expected <- tribble(
    ~USUBJID, ~AVISITN, ~AVALC, ~CONFFL,
    "1",      1,        "PR",   "Y",
    "1",      2,        "CR",   NA_character_,
    "1",      3,        "CR",   NA_character_,
    "1",      4,        "SD",   NA_character_,
    "1",      5,        "NE",   NA_character_,
    "2",      1,        "SD",   NA_character_,
    "2",      2,        "PR",   NA_character_,
    "2",      3,        "PD",   NA_character_,
    "3",      1,        "SD",   NA_character_,
    "4",      1,        "PR",   "Y",
    "4",      2,        "PD",   NA_character_,
    "4",      3,        "SD",   NA_character_,
    "4",      4,        "SD",   NA_character_,
    "4",      5,        "PR",   NA_character_
  )

  expect_dfs_equal(
    base = expected,
    compare = actual,
    keys = c("USUBJID", "AVISITN")
  )
})

## Test 2: filter with first_cond ----
## Flagging any patient CR value that is followed subsequently by a CR
test_that("derive_var_confirmation_flag Test 2: filter with first_cond", {
  actual <-
    derive_var_confirmation_flag(
      data,
      new_var = CONFFL,
      by_vars = vars(USUBJID),
      join_vars = vars(AVALC),
      first_cond = AVALC == "CR" &
        AVALC.join == "CR",
      order = vars(AVISITN),
      filter = TRUE
    )

  expected <- tribble(
    ~USUBJID, ~AVISITN, ~AVALC, ~CONFFL,
    "1",      1,        "PR",   NA_character_,
    "1",      2,        "CR",   "Y",
    "1",      3,        "CR",   NA_character_,
    "1",      4,        "SD",   NA_character_,
    "1",      5,        "NE",   NA_character_,
    "2",      1,        "SD",   NA_character_,
    "2",      2,        "PR",   NA_character_,
    "2",      3,        "PD",   NA_character_,
    "3",      1,        "SD",   NA_character_,
    "4",      1,        "PR",   NA_character_,
    "4",      2,        "PD",   NA_character_,
    "4",      3,        "SD",   NA_character_,
    "4",      4,        "SD",   NA_character_,
    "4",      5,        "PR",   NA_character_
  )

  expect_dfs_equal(
    base = expected,
    compare = actual,
    keys = c("USUBJID", "AVISITN")
  )
})

## Test 3: filter with first_cond and summary function ----
## Flagging any patient PR value that is followed subsequently by a CR or PR
## and at most one SD in between
test_that("derive_var_confirmation_flag Test 3: filter with first_cond and summary function", {
  actual <-
    derive_var_confirmation_flag(
      data,
      new_var = CONFFL,
      by_vars = vars(USUBJID),
      join_vars = vars(AVALC),
      first_cond = AVALC == "PR" &
        AVALC.join %in% c("CR", "PR"),
      order = vars(AVISITN),
      filter = count_vals(AVALC.join, "SD") <= 1,
      false_value = "N"
    )

  expected <- tribble(
    ~USUBJID, ~AVISITN, ~AVALC, ~CONFFL,
    "1",      1,        "PR",   "Y",
    "1",      2,        "CR",   "N",
    "1",      3,        "CR",   "N",
    "1",      4,        "SD",   "N",
    "1",      5,        "NE",   "N",
    "2",      1,        "SD",   "N",
    "2",      2,        "PR",   "N",
    "2",      3,        "PD",   "N",
    "3",      1,        "SD",   "N",
    "4",      1,        "PR",   "N",
    "4",      2,        "PD",   "N",
    "4",      3,        "SD",   "N",
    "4",      4,        "SD",   "N",
    "4",      5,        "PR",   "N"
  )

  expect_dfs_equal(
    base = expected,
    compare = actual,
    keys = c("USUBJID", "AVISITN")
  )
})

## Test 4: join_type = "all" ----
## Flagging observations with a duration longer than 30 and
## on or before 7 days before a COVID AE (ACOVFL == "Y")
test_that("derive_var_confirmation_flag, Test 4: join_type = 'all'", {
  adae <- tribble(
    ~USUBJID, ~ADY, ~ACOVFL, ~ADURN,
    "1",        10, "N",          1,
    "1",        21, "N",         50,
    "1",        23, "Y",         14,
    "1",        32, "N",         31,
    "1",        42, "N",         20,
    "2",         2, "N",         33,
    "2",        11, "Y",         13,
    "2",        23, "N",          2,
    "3",        13, "Y",         12,
    "4",        14, "N",         32,
    "4",        21, "N",         41
  )

  actual <- derive_var_confirmation_flag(
    adae,
    by_vars = vars(USUBJID),
    new_var = ALCOVFL,
    join_vars = vars(ACOVFL, ADY),
    join_type = "all",
    order = vars(ADY),
    filter = ADURN > 30 & ACOVFL.join == "Y" & ADY >= ADY.join - 7
  )

  expected <- tribble(
    ~USUBJID, ~ADY, ~ACOVFL, ~ADURN, ~ALCOVFL,
    "1",        10, "N",          1, NA_character_,
    "1",        21, "N",         50, "Y",
    "1",        23, "Y",         14, NA_character_,
    "1",        32, "N",         31, "Y",
    "1",        42, "N",         20, NA_character_,
    "2",         2, "N",         33, NA_character_,
    "2",        11, "Y",         13, NA_character_,
    "2",        23, "N",          2, NA_character_,
    "3",        13, "Y",         12, NA_character_,
    "4",        14, "N",         32, NA_character_,
    "4",        21, "N",         41, NA_character_
  )

  expect_dfs_equal(
    base = expected,
    compare = actual,
    keys = c("USUBJID", "ADY")
  )
})

## Test 5: join_type = "before" ----
## Flagging observations with AVALC = Y and an observation with CRIT1FL = Y before
test_that("derive_var_confirmation_flag, Test 5: join_type = 'before'", {
  data <- tribble(
    ~USUBJID, ~ASEQ, ~AVALC, ~CRIT1FL,
    "1",          1, "Y",    "Y",
    "1",          2, "N",    "N",
    "1",          3, "Y",    "N",
    "2",          1, "Y",    "Y",
    "3",          1, "N",    "Y"
  )

  actual <- derive_var_confirmation_flag(
    data,
    by_vars = vars(USUBJID),
    order = vars(ASEQ),
    new_var = CONFFL,
    join_vars = vars(CRIT1FL),
    join_type = "before",
    filter = AVALC == "Y" & CRIT1FL.join == "Y",
    false_value = "N"
  )

  expected <- tribble(
    ~USUBJID, ~ASEQ, ~AVALC, ~CRIT1FL, ~CONFFL,
    "1",          1, "Y",    "Y",      "N",
    "1",          2, "N",    "N",      "N",
    "1",          3, "Y",    "N",      "Y",
    "2",          1, "Y",    "Y",      "N",
    "3",          1, "N",    "Y",      "N"
  )

  expect_dfs_equal(
    base = expected,
    compare = actual,
    keys = c("USUBJID", "ASEQ")
  )
})

## Test 6: join_type = "at_before" ----
## Flagging observations with AVALC = Y and an observation with CRIT1FL = Y at or before
test_that("derive_var_confirmation_flag, Test 6: join_type = 'at_before'", {
  data <- tribble(
    ~USUBJID, ~ASEQ, ~AVALC, ~CRIT1FL,
    "1",          1, "Y",    "Y",
    "1",          2, "N",    "N",
    "1",          3, "Y",    "N",
    "2",          1, "Y",    "Y",
    "3",          1, "N",    "Y"
  )

  actual <- derive_var_confirmation_flag(
    data,
    by_vars = vars(USUBJID),
    order = vars(ASEQ),
    new_var = CONFFL,
    join_vars = vars(CRIT1FL),
    join_type = "at_before",
    filter = AVALC == "Y" & CRIT1FL.join == "Y",
    false_value = "N"
  )

  expected <- tribble(
    ~USUBJID, ~ASEQ, ~AVALC, ~CRIT1FL, ~CONFFL,
    "1",          1, "Y",    "Y",      "Y",
    "1",          2, "N",    "N",      "N",
    "1",          3, "Y",    "N",      "Y",
    "2",          1, "Y",    "Y",      "Y",
    "3",          1, "N",    "Y",      "N"
  )

  expect_dfs_equal(
    base = expected,
    compare = actual,
    keys = c("USUBJID", "ASEQ")
  )
})
