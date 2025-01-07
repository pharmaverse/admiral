# derive_vars_joined ----
## Test 1: no by_vars, no order, no new_vars ----
test_that("derive_vars_joined Test 1: no by_vars, no order, no new_vars", {
  expected <- tibble::tribble(
    ~USUBJID, ~ADY, ~AVISIT,    ~AWLO, ~AWHI,
    "1",        -2, "BASELINE",   -30,     1,
    "1",         3, "WEEK 1",       2,     7,
    "1",        24, "WEEK 4",      23,    30,
    "2",        NA, NA,            NA,    NA
  )

  windows <- tibble::tribble(
    ~AVISIT,    ~AWLO, ~AWHI,
    "BASELINE",   -30,     1,
    "WEEK 1",       2,     7,
    "WEEK 2",       8,    15,
    "WEEK 3",      16,    22,
    "WEEK 4",      23,    30
  )

  expect_dfs_equal(
    base = expected,
    comp = derive_vars_joined(
      select(expected, USUBJID, ADY),
      dataset_add = windows,
      join_vars = exprs(AWHI, AWLO),
      join_type = "all",
      filter_join = AWLO <= ADY & ADY <= AWHI
    ),
    keys = c("USUBJID", "ADY")
  )
})

## Test 2: new_vars with rename ----
test_that("derive_vars_joined Test 2: new_vars with rename", {
  expected <- tibble::tribble(
    ~USUBJID, ~ADY, ~AVAL, ~NADIR,
    "1",        -7,    10,     NA,
    "1",         1,    12,     NA,
    "1",         8,    11,     12,
    "1",        15,     9,     11,
    "1",        20,    14,      9,
    "1",        24,    12,      9,
    "2",        13,     8,     NA
  )

  adbds <- select(expected, -NADIR)

  expect_dfs_equal(
    base = expected,
    comp = derive_vars_joined(
      adbds,
      dataset_add = adbds,
      by_vars = exprs(USUBJID),
      order = exprs(AVAL),
      new_vars = exprs(NADIR = AVAL),
      join_vars = exprs(ADY),
      join_type = "all",
      filter_add = ADY > 0,
      filter_join = ADY.join < ADY,
      mode = "first",
      check_type = "none"
    ),
    keys = c("USUBJID", "ADY")
  )
})

## Test 3: by_vars with rename ----
test_that("derive_vars_joined Test 3: by_vars with rename", {
  adae <- tibble::tribble(
    ~AEGRPID,
    "1",
    "2"
  ) %>%
    mutate(
      TRTSDTM = ymd_hms("2020-01-06T12:00:00")
    )

  faae <- tibble::tribble(
    ~FAGRPID, ~FADT,        ~FAORRES,
    "1",      "2020-01-01", "1",
    "1",      "2020-01-03", "2",
    "1",      "2020-01-05", "3",
    "1",      "2020-01-08", "4"
  ) %>%
    mutate(FADT = ymd(FADT))
  expect_dfs_equal(
    base = mutate(adae, ATOXGR_pre = c("3", NA)),
    comp = derive_vars_joined(
      adae,
      dataset_add = faae,
      by_vars = exprs(AEGRPID = FAGRPID),
      order = exprs(FADT),
      new_vars = exprs(ATOXGR_pre = FAORRES),
      join_vars = exprs(FADT),
      join_type = "all",
      filter_join = FADT < TRTSDTM,
      mode = "last"
    ),
    keys = c("AEGRPID")
  )
})

## Test 4: order with expression ----
test_that("derive_vars_joined Test 4: order with expression", {
  adae <- tibble::tribble(
    ~AEGRPID,
    "1",
    "2"
  ) %>%
    mutate(
      TRTSDTM = ymd_hms("2020-01-06T12:00:00")
    )

  faae <- tibble::tribble(
    ~FAGRPID, ~FADTC,       ~FAORRES,
    "1",      "2020-01-01", "1",
    "1",      "2020-01-03", "2",
    "1",      "2020-01-05", "3",
    "1",      "2020-01-08", "4"
  )
  expect_dfs_equal(
    base = mutate(adae, ATOXGR_pre = c("3", NA)),
    comp = derive_vars_joined(
      adae,
      dataset_add = faae,
      by_vars = exprs(AEGRPID = FAGRPID),
      order = exprs(FADT = convert_dtc_to_dt(FADTC)),
      new_vars = exprs(ATOXGR_pre = FAORRES),
      join_vars = exprs(FADT),
      join_type = "all",
      filter_join = FADT < TRTSDTM,
      mode = "last"
    ),
    keys = c("AEGRPID")
  )
})

## Test 5: join_vars with expression ----
test_that("derive_vars_joined Test 5: join_vars with expression", {
  add <- tibble::tribble(
    ~USUBJID, ~TRDTC,       ~TRSTRESN,
    "1",      "2020-02-01",        10,
    "1",      "2020-02-04",        12,
    "1",      "2020-02-08",        11,
    "1",      "2020-02-13",         9,
    "1",      "2020-02-24",        14,
    "1",      "2020-03-01",        12,
    "2",      "2021-01-13",         8
  )

  expected <- tibble::tribble(
    ~USUBJID, ~ADT,         ~AVAL,
    "1",      "2020-02-09",    10,
    "1",      "2020-02-13",     9,
    "1",      "2020-02-24",     9,
    "1",      "2020-03-01",     9,
    "2",      "2021-01-13",     8
  ) %>%
    mutate(
      ADT = ymd(ADT)
    )

  expect_dfs_equal(
    base = expected,
    comp = derive_vars_joined(
      select(expected, -AVAL),
      dataset_add = add,
      by_vars = exprs(USUBJID),
      order = exprs(TRSTRESN),
      new_vars = exprs(AVAL = TRSTRESN),
      join_vars = exprs(TRDT = convert_dtc_to_dt(TRDTC)),
      join_type = "all",
      filter_join = TRDT <= ADT,
      mode = "first",
      check_type = "none"
    ),
    keys = c("USUBJID", "ADT")
  )
})


## Test 6: no join_vars, no filter_join ----
test_that("derive_vars_joined Test 6: no join_vars, no filter_join", {
  adae <- tibble::tribble(
    ~AEGRPID,
    "1",
    "2"
  ) %>%
    mutate(
      TRTSDTM = ymd_hms("2020-01-06T12:00:00")
    )

  faae <- tibble::tribble(
    ~FAGRPID, ~FADT,        ~FAORRES,
    "1",      "2020-01-01", "1",
    "1",      "2020-01-03", "2",
    "1",      "2020-01-05", "3",
    "1",      "2020-01-08", "4"
  ) %>%
    mutate(FADT = ymd(FADT))
  expect_dfs_equal(
    base = mutate(adae, ATOXGR_pre = c("1", NA)),
    comp = derive_vars_joined(
      adae,
      dataset_add = faae,
      by_vars = exprs(AEGRPID = FAGRPID),
      order = exprs(FAORRES),
      join_type = "all",
      new_vars = exprs(ATOXGR_pre = FAORRES),
      mode = "first"
    ),
    keys = c("AEGRPID")
  )
})

## Test 7: new_vars expressions using variables from both datasets ----
test_that("derive_vars_joined Test 7: new_vars expressions using variables from both datasets", {
  expected <- tibble::tribble(
    ~USUBJID, ~ASTDT,       ~AESEQ, ~LSTDSDUR,
    "1",      "2020-02-02",      1,        14,
    "1",      "2020-02-04",      2,         2
  ) %>%
    mutate(ASTDT = ymd(ASTDT))

  ex <- tibble::tribble(
    ~USUBJID, ~EXSDTC,
    "1",      "2020-01-10",
    "1",      "2020-01",
    "1",      "2020-01-20",
    "1",      "2020-02-03"
  )

  expect_dfs_equal(
    base = expected,
    compare = derive_vars_joined(
      select(expected, -LSTDSDUR),
      dataset_add = ex,
      by_vars = exprs(USUBJID),
      order = exprs(EXSDT = convert_dtc_to_dt(EXSDTC)),
      join_type = "all",
      new_vars = exprs(LSTDSDUR = compute_duration(
        start_date = EXSDT, end_date = ASTDT
      )),
      filter_add = !is.na(EXSDT),
      filter_join = EXSDT <= ASTDT,
      mode = "last"
    ),
    keys = c("USUBJID", "AESEQ")
  )
})

## Test 8: error if new_vars are already in dataset ----
test_that("derive_vars_joined Test 8: error if new_vars are already in dataset", {
  myd <- data.frame(day = c(1, 2, 3), val = c(0, 17, 21))
  expect_snapshot(
    derive_vars_joined(
      myd,
      dataset_add = myd,
      order = exprs(day),
      join_type = "all",
      mode = "last",
      filter_join = day < day.join
    ),
    error = TRUE
  )
})

## Test 9: fixing a bug from issue 1966 ----
test_that("derive_vars_joined Test 9: fixing a bug from issue 1966", { # nolint
  adlb_ast <- tribble(
    ~ADT,         ~ASEQ,
    "2002-01-01", 1,
    "2002-02-02", 2,
    "2002-02-02", 3
  ) %>%
    mutate(
      STUDYID = "ABC",
      USUBJID = "1",
      ADT = ymd(ADT),
      ADTM = as_datetime(ADT)
    )

  adlb_tbili_pbl <- tribble(
    ~ADT,         ~ASEQ,
    "2002-01-01", 4,
    "2002-02-02", 5,
    "2002-02-02", 6
  ) %>%
    mutate(
      STUDYID = "ABC",
      USUBJID = "1",
      ADT = ymd(ADT),
      ADTM = as_datetime(ADT)
    )

  adlb_joined <- derive_vars_joined(
    adlb_ast,
    dataset_add = adlb_tbili_pbl,
    by_vars = exprs(STUDYID, USUBJID),
    order = exprs(ADTM, ASEQ),
    join_type = "all",
    new_vars = exprs(TBILI_ADT = ADT),
    filter_join = ADT <= ADT.join,
    mode = "first"
  )

  expected <- adlb_ast %>%
    mutate(TBILI_ADT = as.Date(c("2002-01-01", "2002-02-02", "2002-02-02"), "%Y-%m-%d"))

  expect_dfs_equal(
    base = expected,
    compare = adlb_joined,
    keys = c("ADT", "ASEQ", "STUDYID", "USUBJID", "ADTM", "TBILI_ADT")
  )
})

## Test 10: order vars are selected properly in function body ----
test_that("derive_vars_joined Test 10: order vars are selected properly in function body", {
  myd <- data.frame(day = c(1, 2, 3), val = c(0, 17, 21))
  actual <- derive_vars_joined(
    myd,
    dataset_add = myd,
    new_vars = exprs(first_val = val),
    join_vars = exprs(day),
    join_type = "all",
    order = exprs(-day),
    mode = "last",
    filter_join = day < day.join
  )
  expected <- tribble(
    ~day, ~val, ~first_val,
    1,       0,         17,
    2,      17,         21,
    3,      21,         NA
  )

  expect_dfs_equal(
    base = expected,
    compare = actual,
    keys = c("day", "val", "first_val")
  )
})


## Test 11: Ensure exist_flag, true/false value arguments work ----
test_that("derive_vars_joined Test 11: Ensure exist_flag, true/false value arguments work", {
  expected <- tibble::tribble(
    ~USUBJID, ~ADY, ~AVISIT,    ~AWLO, ~AWHI,  ~flag,
    "1",        -2, "BASELINE",   -30,     1,  "Yes",
    "1",         3, "WEEK 1",       2,     7,  "Yes",
    "1",        24, "WEEK 4",      23,    30,  "Yes",
    "2",        NA, NA,            NA,    NA,   "No"
  )

  windows <- tibble::tribble(
    ~AVISIT,    ~AWLO, ~AWHI,
    "BASELINE",   -30,     1,
    "WEEK 1",       2,     7,
    "WEEK 2",       8,    15,
    "WEEK 3",      16,    22,
    "WEEK 4",      23,    30
  )

  expect_dfs_equal(
    base = expected,
    comp = derive_vars_joined(
      select(expected, USUBJID, ADY),
      dataset_add = windows,
      join_vars = exprs(AWHI, AWLO),
      join_type = "all",
      filter_join = AWLO <= ADY & ADY <= AWHI,
      exist_flag = flag,
      true_value = "Yes",
      false_value = "No"
    ),
    keys = c("USUBJID", "ADY")
  )
})

## Test 12: by_vars with rename, no new_vars ----
test_that("derive_vars_joined Test 12: by_vars with rename, no new_vars", {
  expected <- tibble::tribble(
    ~USUBJID, ~EXLINKID, ~FA,
    "1",      "001",       1,
    "1",      "002",      NA,
    "2",      "001",       1
  )

  fa <- tibble::tribble(
    ~USUBJID, ~FALINKID, ~FA,
    "1",      "001",     1,
    "1",      "002",     0,
    "2",      "001",     1,
  )

  ex <- select(expected, -FA)

  expect_dfs_equal(
    base = expected,
    compare = derive_vars_joined(
      ex,
      dataset_add = fa,
      by_vars = exprs(USUBJID, EXLINKID = FALINKID),
      join_type = "all",
      filter_add = FA == 1
    ),
    keys = c("USUBJID", "EXLINKID")
  )
})

## Test 13: with save_memory and no by_vars ----
test_that("derive_vars_joined Test 13: with save_memory and no by_vars", {
  expected <- tibble::tribble(
    ~USUBJID, ~ADY, ~AVISIT,    ~AWLO, ~AWHI,
    "1",        -2, "BASELINE",   -30,     1,
    "1",         3, "WEEK 1",       2,     7,
    "1",        24, "WEEK 4",      23,    30,
    "2",        NA, NA,            NA,    NA
  )

  windows <- tibble::tribble(
    ~AVISIT,    ~AWLO, ~AWHI,
    "BASELINE",   -30,     1,
    "WEEK 1",       2,     7,
    "WEEK 2",       8,    15,
    "WEEK 3",      16,    22,
    "WEEK 4",      23,    30
  )

  save_memory <- get_admiral_option("save_memory")
  set_admiral_options(save_memory = TRUE)

  expect_dfs_equal(
    base = expected,
    comp = derive_vars_joined(
      select(expected, USUBJID, ADY),
      dataset_add = windows,
      join_vars = exprs(AWHI, AWLO),
      join_type = "all",
      filter_join = AWLO <= ADY & ADY <= AWHI
    ),
    keys = c("USUBJID", "ADY")
  )
  set_admiral_options(save_memory = save_memory)
})

## Test 14: save_memory with by_vars ----
test_that("derive_vars_joined Test 14: save_memory with by_vars", {
  expected <- tibble::tribble(
    ~USUBJID, ~ADY, ~AVAL, ~NADIR,
    "1",        -7,    10,     NA,
    "1",         1,    12,     NA,
    "1",         8,    11,     12,
    "1",        15,     9,     11,
    "1",        20,    14,      9,
    "1",        24,    12,      9,
    "2",        13,     8,     NA
  )

  adbds <- select(expected, -NADIR)

  save_memory <- get_admiral_option("save_memory")
  set_admiral_options(save_memory = TRUE)

  expect_dfs_equal(
    base = expected,
    comp = derive_vars_joined(
      adbds,
      dataset_add = adbds,
      by_vars = exprs(USUBJID),
      order = exprs(AVAL),
      new_vars = exprs(NADIR = AVAL),
      join_vars = exprs(ADY),
      join_type = "all",
      filter_add = ADY > 0,
      filter_join = ADY.join < ADY,
      mode = "first",
      check_type = "none"
    ),
    keys = c("USUBJID", "ADY")
  )
  set_admiral_options(save_memory = save_memory)
})

# get_joined_data ----
## Test 15: `first_cond_lower` works ----
test_that("get_joined_data Test 15: `first_cond_lower` works", {
  data <- tribble(
    ~subj, ~day, ~val,
    "1",      1, "++",
    "1",      2, "-",
    "1",      3, "0",
    "1",      4, "+",
    "1",      5, "++",
    "1",      6, "-",
    "2",      1, "-",
    "2",      2, "++",
    "2",      3, "+",
    "2",      4, "0",
    "2",      5, "-",
    "2",      6, "++"
  )

  expected <- tibble::tribble(
    ~day.join, ~val.join,
    2,         "++",
    3,         "+"
  ) %>%
    mutate(
      subj = "2",
      day = 4,
      val = "0"
    )

  expect_dfs_equal(
    base = expected,
    compare = get_joined_data(
      data,
      dataset_add = data,
      by_vars = exprs(subj),
      order = exprs(day),
      join_vars = exprs(val),
      join_type = "before",
      first_cond_lower = val.join == "++",
      filter_join = val == "0" & all(val.join %in% c("+", "++"))
    ),
    keys = c("subj", "day.join")
  )
})
