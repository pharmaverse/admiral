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


# derive_var_joined_exist_flag ----
## Flagging any patient PR value that is followed by a CR or PR
## Test 1: filter without first_cond ----
test_that("derive_var_joined_exist_flag Test 1: filter without first_cond", {
  actual <-
    derive_var_joined_exist_flag(
      data,
      dataset_add = data,
      new_var = CONFFL,
      by_vars = exprs(USUBJID),
      join_vars = exprs(AVALC),
      join_type = "after",
      order = exprs(AVISITN),
      filter_join = AVALC == "PR" & AVALC.join %in% c("CR", "PR")
    )

  expected <- tibble::tribble(
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

## Flagging any patient CR value that is followed by a CR
## Test 2: filter with first_cond ----
test_that("derive_var_joined_exist_flag Test 2: filter with first_cond", {
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
    "3",      1,        "CR",
    "4",      1,        "CR",
    "4",      2,        "SD",
    "4",      3,        "CR",
    "4",      4,        "CR"
  )
  actual <-
    derive_var_joined_exist_flag(
      data,
      dataset_add = data,
      new_var = CONFFL,
      by_vars = exprs(USUBJID),
      join_vars = exprs(AVALC),
      join_type = "after",
      first_cond_upper = AVALC == "CR" &
        AVALC.join == "CR",
      order = exprs(AVISITN),
      filter_join = TRUE
    )

  expected <- tibble::tribble(
    ~USUBJID, ~AVISITN, ~AVALC, ~CONFFL,
    "1",      1,        "PR",   NA_character_,
    "1",      2,        "CR",   "Y",
    "1",      3,        "CR",   NA_character_,
    "1",      4,        "SD",   NA_character_,
    "1",      5,        "NE",   NA_character_,
    "2",      1,        "SD",   NA_character_,
    "2",      2,        "PR",   NA_character_,
    "2",      3,        "PD",   NA_character_,
    "3",      1,        "CR",   NA_character_,
    "4",      1,        "CR",   "Y",
    "4",      2,        "SD",   NA_character_,
    "4",      3,        "CR",   "Y",
    "4",      4,        "CR",   NA_character_
  )

  expect_dfs_equal(
    base = expected,
    compare = actual,
    keys = c("USUBJID", "AVISITN")
  )
})

## Flagging any patient PR value that is followed by a CR or PR
## and at most one SD in between
## Test 3: filter with first_cond and summary function ----
test_that("derive_var_joined_exist_flag Test 3: filter with first_cond and summary function", {
  actual <-
    derive_var_joined_exist_flag(
      data,
      dataset_add = data,
      new_var = CONFFL,
      by_vars = exprs(USUBJID),
      join_vars = exprs(AVALC),
      join_type = "after",
      first_cond_upper = AVALC == "PR" &
        AVALC.join %in% c("CR", "PR"),
      order = exprs(AVISITN),
      filter_join = count_vals(AVALC.join, "SD") <= 1,
      false_value = "N"
    )

  expected <- tibble::tribble(
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

## Flagging observations with a duration longer than 30 and
## on or after 7 days of a COVID AE (ACOVFL == "Y")
## Test 4: join_type = 'all' ----
test_that("derive_var_joined_exist_flag Test 4: join_type = 'all'", {
  adae <- tibble::tribble(
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

  actual <- derive_var_joined_exist_flag(
    adae,
    dataset_add = adae,
    by_vars = exprs(USUBJID),
    new_var = ALCOVFL,
    join_vars = exprs(ACOVFL, ADY),
    join_type = "all",
    order = exprs(ADY),
    filter_join = ADURN > 30 & ACOVFL.join == "Y" & ADY >= ADY.join - 7
  )

  expected <- tibble::tribble(
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

## Flagging observations with AVALC = Y and an observation with CRIT1FL = Y before
## Test 5: join_type = 'before' ----
test_that("derive_var_joined_exist_flag Test 5: join_type = 'before'", {
  data <- tibble::tribble(
    ~USUBJID, ~ASEQ, ~AVALC, ~CRIT1FL,
    "1",          1, "Y",    "Y",
    "1",          2, "N",    "N",
    "1",          3, "Y",    "N",
    "2",          1, "Y",    "Y",
    "3",          1, "N",    "Y"
  )

  actual <- derive_var_joined_exist_flag(
    data,
    dataset_add = data,
    by_vars = exprs(USUBJID),
    order = exprs(ASEQ),
    new_var = CONFFL,
    join_vars = exprs(CRIT1FL),
    join_type = "before",
    filter_join = AVALC == "Y" & CRIT1FL.join == "Y",
    false_value = "N"
  )

  expected <- tibble::tribble(
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


## Test 6: tmp_obs_nr_var argument works ----
test_that("derive_var_joined_exist_flag Test 6: tmp_obs_nr_var argument works", {
  expected <- tibble::tribble(
    ~USUBJID, ~AVISITN, ~CRIT1FL, ~CONFFL,
    "1",      1,        "Y",      "N",
    "1",      2,        "N",      "N",
    "1",      3,        "Y",      "N",
    "1",      5,        "N",      "N",
    "2",      1,        "Y",      "Y",
    "2",      3,        "Y",      "N",
    "2",      5,        "N",      "N",
    "3",      1,        "Y",      "Y",
    "4",      1,        "Y",      "N",
    "4",      2,        "N",      "N"
  )

  input <- select(expected, -CONFFL)
  expect_dfs_equal(
    base = expected,
    compare = derive_var_joined_exist_flag(
      input,
      dataset_add = input,
      by_vars = exprs(USUBJID),
      new_var = CONFFL,
      tmp_obs_nr_var = tmp_obs_nr,
      join_vars = exprs(CRIT1FL),
      join_type = "all",
      order = exprs(AVISITN),
      filter_join = CRIT1FL == "Y" & CRIT1FL.join == "Y" &
        (tmp_obs_nr + 1 == tmp_obs_nr.join | tmp_obs_nr == max(tmp_obs_nr.join)),
      false_value = "N"
    ),
    keys = c("USUBJID", "AVISITN")
  )
})

## Test 7: deprecation of `filter` ----
test_that("derive_var_joined_exist_flag Test 7: deprecation of `filter`", {
  expect_error(
    actual <-
      derive_var_joined_exist_flag(
        data,
        dataset_add = data,
        new_var = CONFFL,
        by_vars = exprs(USUBJID),
        join_vars = exprs(AVALC),
        join_type = "after",
        order = exprs(AVISITN),
        filter = AVALC == "PR" & AVALC.join %in% c("CR", "PR")
      ),
    class = "lifecycle_error_deprecated"
  )
})

## Test 8: deprecation of `first_cond` ----
test_that("derive_var_joined_exist_flag Test 8: deprecation of `first_cond`", {
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
    "3",      1,        "CR",
    "4",      1,        "CR",
    "4",      2,        "SD",
    "4",      3,        "CR",
    "4",      4,        "CR"
  )

  expect_error(
    actual <-
      derive_var_joined_exist_flag(
        data,
        dataset_add = data,
        new_var = CONFFL,
        by_vars = exprs(USUBJID),
        join_vars = exprs(AVALC),
        join_type = "after",
        first_cond = AVALC == "CR" &
          AVALC.join == "CR",
        order = exprs(AVISITN),
        filter_join = TRUE
      ),
    class = "lifecycle_error_deprecated"
  )
})
