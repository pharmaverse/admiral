adsl <- tibble::tribble(
  ~USUBJID, ~SEX, ~COUNTRY,
  "ST42-1", "F",  "AUT",
  "ST42-2", "M",  "MWI",
  "ST42-3", "M",  "NOR",
  "ST42-4", "F",  "UGA"
) %>% mutate(STUDYID = "ST42")

adsl1 <- tibble::tribble(
  ~ID, ~SEX, ~COUNTRY,
  "ST42-1", "F", "AUT",
  "ST42-2", "M", "MWI",
  "ST42-3", "M", "NOR",
  "ST42-4", "F", "UGA"
) %>% mutate(STUDYID = "ST42")

adsl2 <- tibble::tribble(
  ~ID, ~SEX, ~COUNTRY,
  "ST42-1", "F", "AUT",
  "ST42-1", "F", "NOR",
  "ST42-2", "M", "MWI",
  "ST42-3", "M", "NOR",
  "ST42-4", "F", "UGA"
) %>% mutate(STUDYID = "ST42")


advs <- tibble::tribble(
  ~USUBJID, ~PARAMCD, ~AVISIT,    ~AVAL,
  "ST42-1", "WEIGHT", "BASELINE", 66,
  "ST42-1", "WEIGHT", "Week 2",   68,
  "ST42-2", "WEIGHT", "BASELINE", 88,
  "ST42-3", "WEIGHT", "Week 2",   55,
  "ST42-3", "WEIGHT", "Week 4",   50
) %>% mutate(STUDYID = "ST42")

advs1 <- tibble::tribble(
  ~ID, ~PARAMCD, ~AVISIT, ~AVAL,
  "ST42-1", "WEIGHT", "BASELINE", 66,
  "ST42-1", "WEIGHT", "Week 2", 68,
  "ST42-2", "WEIGHT", "BASELINE", 88,
  "ST42-3", "WEIGHT", "Week 2", 55,
  "ST42-3", "WEIGHT", "Week 4", 50
) %>% mutate(STUDYID = "ST42")


ex <- tibble::tribble(
  ~USUBJID, ~EXSTDTC,
  "ST42-1", "2020-12-07",
  "ST42-1", "2020-12-14",
  "ST42-2", "2021-01-12T12:00:00",
  "ST42-2", "2021-01-26T13:21",
  "ST42-3", "2021-03-02"
) %>% mutate(STUDYID = "ST42")

vs <- tibble::tribble(
  ~USUBJID, ~VSTESTCD, ~VSTEST, ~VSORRES, ~VSSEQ,
  "ST42-1", "DIABP", "Diastolic Blood Pressure", 64, 1,
  "ST42-1", "DIABP", "Diastolic Blood Pressure", 83, 2,
  "ST42-1", "WEIGHT", "Weight", 120, 3,
  "ST42-2", "WEIGHT", "Weight", 110, 1,
  "ST42-2", "HEIGHT", "Height", 58, 2
) %>% mutate(STUDYID = "ST42")

# derive_vars_merged ----
## Test 1: merge all variables ----
test_that("derive_vars_merged Test 1: merge all variables", {
  actual <- derive_vars_merged(advs,
    dataset_add = adsl,
    by_vars = exprs(STUDYID, USUBJID)
  )

  expected <- left_join(advs, adsl, by = c("STUDYID", "USUBJID"))

  expect_dfs_equal(
    base = expected,
    compare = actual,
    keys = c("USUBJID", "AVISIT")
  )
})

## Test 2: merge selected variables ----
test_that("derive_vars_merged Test 2: merge selected variables", {
  actual <- derive_vars_merged(advs,
    dataset_add = adsl,
    by_vars = exprs(USUBJID),
    new_vars = exprs(SEX)
  )

  expected <- left_join(advs, select(adsl, USUBJID, SEX), by = "USUBJID")

  expect_dfs_equal(
    base = expected,
    compare = actual,
    keys = c("USUBJID", "AVISIT")
  )
})

## Test 3: merge last value and flag matched by groups ----
test_that("derive_vars_merged Test 3: merge last value and flag matched by groups", {
  actual <- derive_vars_merged(adsl,
    dataset_add = advs,
    order = exprs(AVAL),
    by_vars = exprs(STUDYID, USUBJID),
    new_vars = exprs(WEIGHTBL = AVAL),
    mode = "last",
    exist_flag = matched
  )
  expected <- adsl %>% mutate(
    WEIGHTBL = c(68, 88, 55, NA),
    matched = c("Y", "Y", "Y", NA_character_)
  )

  expect_dfs_equal(
    base = expected,
    compare = actual,
    keys = c("USUBJID")
  )
})

## Test 4: merge last value and flag matched by groups ----
test_that("derive_vars_merged Test 4: merge last value and flag matched by groups", {
  actual <- derive_vars_merged(adsl,
    dataset_add = advs,
    order = exprs(AVAL),
    by_vars = exprs(STUDYID, USUBJID),
    new_vars = exprs(WEIGHTBL = AVAL),
    mode = "last",
    exist_flag = matched,
    true_value = "Y",
    false_value = "N"
  )
  expected <- adsl %>% mutate(
    WEIGHTBL = c(68, 88, 55, NA),
    matched = c("Y", "Y", "Y", "N")
  )

  expect_dfs_equal(
    base = expected,
    compare = actual,
    keys = c("USUBJID")
  )
})

## Test 5: error if variable in both datasets ----
test_that("derive_vars_merged Test 5: error if variable in both datasets", {
  expect_snapshot(
    derive_vars_merged(advs,
      dataset_add = adsl,
      by_vars = exprs(USUBJID)
    ),
    error = TRUE
  )
})

## Test 6: by_vars with rename ----
test_that("derive_vars_merged Test 6: by_vars with rename", {
  actual <- derive_vars_merged(advs,
    dataset_add = adsl1,
    by_vars = exprs(STUDYID, USUBJID = ID),
    filter_add = SEX == "F"
  )

  adsl_1 <- adsl1 %>% filter(SEX == "F")
  expected <- left_join(advs, adsl_1, by = c("STUDYID", "USUBJID" = "ID"))

  expect_dfs_equal(
    base = expected,
    compare = actual,
    keys = c("USUBJID", "AVISIT")
  )
})

## Test 7: expressions for new_vars and missing_values ----
test_that("derive_vars_merged Test 7: expressions for new_vars and missing_values", {
  actual <- derive_vars_merged(
    adsl,
    dataset_add = advs,
    by_vars = exprs(USUBJID),
    order = exprs(AVISIT),
    new_vars = exprs(LASTVIS = str_to_upper(AVISIT)),
    mode = "last",
    missing_values = exprs(LASTVIS = "UNKNOWN")
  )

  expected <-
    mutate(adsl, LASTVIS = c("WEEK 2", "BASELINE", "WEEK 4", "UNKNOWN"))


  expect_dfs_equal(
    base = expected,
    compare = actual,
    keys = "USUBJID"
  )
})


## Test 8: Use of missing_values and exist_flags ----
test_that("derive_vars_merged Test 8: Use of missing_values and exist_flags", {
  actual <- derive_vars_merged(
    adsl,
    dataset_add = advs,
    by_vars = exprs(USUBJID),
    order = exprs(AVISIT),
    new_vars = exprs(LASTVIS = str_to_upper(AVISIT)),
    mode = "last",
    missing_values = exprs(LASTVIS = "UNKNOWN"),
    exist_flag = matched,
    true_value = NA,
    false_value = "No"
  )

  expected <- adsl %>%
    mutate(
      LASTVIS = c("WEEK 2", "BASELINE", "WEEK 4", "UNKNOWN"),
      matched = c(NA, NA, NA, "No")
    )


  expect_dfs_equal(
    base = expected,
    compare = actual,
    keys = "USUBJID"
  )
})

## Test 9: use new variables in filter_add and order ----
test_that("derive_vars_merged Test 9: use new variables in filter_add and order", {
  expected <- tibble::tribble(
    ~USUBJID, ~TRTSDT,      ~TRTSSEQ,
    "ST42-1", "2020-12-14",        2,
    "ST42-2", "2021-01-26",        2,
    "ST42-3", NA,                 NA,
    "ST42-4", NA,                 NA
  ) %>% mutate(
    STUDYID = "ST42",
    TRTSDT = ymd(TRTSDT)
  )

  ex <- tibble::tribble(
    ~USUBJID, ~EXSTDTC,              ~EXSEQ,
    "ST42-1", "2020-12-07",               1,
    "ST42-1", "2020-12-14",               2,
    "ST42-2", "2021-01-12T12:00:00",      1,
    "ST42-2", "2021-01-26T13:21",         2,
    "ST42-3", "2021-03",                  1
  ) %>% mutate(STUDYID = "ST42")

  actual <- derive_vars_merged(
    select(adsl, STUDYID, USUBJID),
    dataset_add = ex,
    by_vars = exprs(USUBJID),
    order = exprs(TRTSDT),
    new_vars = exprs(TRTSDT = convert_dtc_to_dt(EXSTDTC), TRTSSEQ = EXSEQ),
    filter_add = !is.na(TRTSDT),
    mode = "last"
  )

  expect_dfs_equal(
    base = expected,
    compare = actual,
    keys = "USUBJID"
  )
})

## Test 10: warning if not unique w.r.t the by variables and the order ----
test_that("derive_vars_merged Test 10: warning if not unique w.r.t the by variables and the order", { # nolint
  expect_warning(
    actual <- derive_vars_merged(advs,
      dataset_add = adsl2,
      by_vars = exprs(STUDYID, USUBJID = ID),
      order = exprs(ID),
      mode = "last",
      check_type = "warning"
    ),
    regexp = ""
  )
})

## Test 11: error if not unique w.r.t the by variables and the order ----
test_that("derive_vars_merged Test 11: error if not unique w.r.t the by variables and the order", {
  expect_snapshot(
    actual <- derive_vars_merged(advs,
      dataset_add = adsl2,
      by_vars = exprs(STUDYID, USUBJID = ID),
      order = exprs(ID),
      mode = "last",
      check_type = "error",
      duplicate_msg = "Duplicate records present!"
    ),
    error = TRUE
  )
})

## Test 12: error if variables in missing_values but not in new_vars ----
test_that("derive_vars_merged Test 12: error if variables in missing_values but not in new_vars", {
  expect_snapshot(
    derive_vars_merged(
      adsl,
      dataset_add = advs,
      by_vars = exprs(USUBJID),
      order = exprs(AVISIT),
      new_vars = exprs(LASTVIS = str_to_upper(AVISIT)),
      mode = "last",
      missing_values = exprs(LASTVIS = "UNKNOWN", LASTVISN = -1)
    ),
    error = TRUE
  )
})

## Test 13: error if not unique, no order, check_type = error ----
test_that("derive_vars_merged Test 13: error if not unique, no order, check_type = error", {
  expect_snapshot(
    actual <- derive_vars_merged(advs,
      dataset_add = adsl2,
      by_vars = exprs(STUDYID, USUBJID = ID),
      order = NULL,
      check_type = "error"
    ),
    error = TRUE
  )
})

## Test 14: error if not unique, no order, check_type = warning ----
test_that("derive_vars_merged Test 14: error if not unique, no order, check_type = warning", {
  expect_snapshot(
    actual <- derive_vars_merged(
      advs,
      dataset_add = adsl2,
      by_vars = exprs(STUDYID, USUBJID = ID),
      order = NULL,
      check_type = "warning"
    ),
    error = TRUE
  )
})

## Test 15: error if not unique, no order, check_type = NULL ----
test_that("derive_vars_merged Test 15: error if not unique, no order, check_type = NULL", {
  expect_snapshot(
    actual <- derive_vars_merged(
      advs,
      dataset_add = adsl2,
      by_vars = exprs(STUDYID, USUBJID = ID),
      order = NULL,
      check_type = NULL
    ),
    error = TRUE
  )
})


# derive_var_merged_exist_flag ----
## Test 16: merge existence flag ----
test_that("derive_var_merged_exist_flag Test 16: merge existence flag", {
  actual <- derive_var_merged_exist_flag(
    adsl,
    advs,
    by_vars = exprs(USUBJID),
    new_var = VSEVALFL,
    condition = AVISIT == "BASELINE"
  )
  expected <-
    mutate(adsl, VSEVALFL = c("Y", "Y", NA_character_, NA_character_))
  expect_dfs_equal(
    base = expected,
    compare = actual,
    keys = "USUBJID"
  )
})

## Test 17: by_vars with rename ----
test_that("derive_var_merged_exist_flag Test 17: by_vars with rename", {
  actual <- derive_var_merged_exist_flag(
    adsl,
    dataset_add = advs1,
    by_vars = exprs(USUBJID = ID),
    new_var = VSEVALFL,
    condition = AVISIT == "BASELINE"
  )

  expected <-
    mutate(adsl, VSEVALFL = c("Y", "Y", NA_character_, NA_character_))


  expect_dfs_equal(
    base = expected,
    compare = actual,
    keys = "USUBJID"
  )
})

# derive_vars_merged_lookup ----
## Test 18: merge lookup table ----
test_that("derive_vars_merged_lookup Test 18: merge lookup table", {
  param_lookup <- tibble::tribble(
    ~VSTESTCD, ~VSTEST, ~PARAMCD, ~DESCRIPTION,
    "WEIGHT", "Weight", "WEIGHT", "Weight (kg)",
    "HEIGHT", "Height", "HEIGHT", "Height (cm)",
    "BMI", "Body Mass Index", "BMI", "Body Mass Index(kg/m^2)"
  )

  attr(param_lookup$VSTESTCD, "label") <- "Vital Signs Test Short Name"
  attr(param_lookup$VSTEST, "label") <- "Vital Signs Test Name"

  expect_snapshot(
    actual <- derive_vars_merged_lookup(
      vs,
      dataset_add = param_lookup,
      by_vars = exprs(VSTESTCD, VSTEST),
      new_var = exprs(PARAMCD, PARAM = DESCRIPTION),
      print_not_mapped = TRUE
    )
  )

  expected <-
    left_join(vs, param_lookup, by = c("VSTESTCD", "VSTEST")) %>%
    rename(PARAM = DESCRIPTION)


  expect_dfs_equal(
    base = expected,
    compare = actual,
    keys = c("USUBJID", "VSSEQ", "VSTESTCD")
  )
})



## the lookup table
## Test 19:  all by_vars have records in the lookup table ----
test_that("derive_vars_merged_lookup Test 19:  all by_vars have records in the lookup table", {
  param_lookup <- tibble::tribble(
    ~VSTESTCD, ~VSTEST, ~PARAMCD, ~DESCRIPTION,
    "WEIGHT", "Weight", "WEIGHT", "Weight (kg)",
    "HEIGHT", "Height", "HEIGHT", "Height (cm)",
    "BMI", "Body Mass Index", "BMI", "Body Mass Index(kg/m^2)",
    "DIABP", "Diastolic Blood Pressure", "DIABP", "Diastolic Blood Pressure (mmHg)"
  )

  attr(param_lookup$VSTESTCD, "label") <- "Vital Signs Test Short Name"
  attr(param_lookup$VSTEST, "label") <- "Vital Signs Test Name"

  expect_message(
    actual <- derive_vars_merged_lookup(
      vs,
      dataset_add = param_lookup,
      by_vars = exprs(VSTESTCD, VSTEST),
      new_var = exprs(PARAMCD, PARAM = DESCRIPTION),
      print_not_mapped = TRUE
    ),
    regex = "All `VSTESTCD` and `VSTEST` are mapped."
  )

  expected <-
    left_join(vs, param_lookup, by = c("VSTESTCD", "VSTEST")) %>%
    rename(PARAM = DESCRIPTION)


  expect_dfs_equal(
    base = expected,
    compare = actual,
    keys = c("USUBJID", "VSSEQ", "VSTESTCD")
  )
})

## Test 20: by_vars with rename ----
test_that("derive_vars_merged_lookup Test 20: by_vars with rename", {
  param_lookup <- tibble::tribble(
    ~TESTCD, ~VSTEST, ~PARAMCD, ~DESCRIPTION,
    "WEIGHT", "Weight", "WEIGHT", "Weight (kg)",
    "HEIGHT", "Height", "HEIGHT", "Height (cm)",
    "BMI", "Body Mass Index", "BMI", "Body Mass Index(kg/m^2)"
  )

  attr(param_lookup$TESTCD, "label") <- "Vital Signs Test Short Name"
  attr(param_lookup$VSTEST, "label") <- "Vital Signs Test Name"

  expect_snapshot(
    actual <- derive_vars_merged_lookup(
      vs,
      dataset_add = param_lookup,
      by_vars = exprs(VSTESTCD = TESTCD, VSTEST),
      new_var = exprs(PARAMCD, PARAM = DESCRIPTION),
      print_not_mapped = TRUE
    )
  )

  expected <-
    left_join(vs, param_lookup, by = c("VSTESTCD" = "TESTCD", "VSTEST")) %>%
    rename(PARAM = DESCRIPTION)

  expect_dfs_equal(
    base = expected,
    compare = actual,
    keys = c("USUBJID", "VSSEQ", "VSTESTCD")
  )
})


# get_not_mapped ----
## Test 21: not all by_vars have records in the lookup table ----
test_that("get_not_mapped Test 21: not all by_vars have records in the lookup table", {
  param_lookup <- tibble::tribble(
    ~VSTESTCD, ~VSTEST, ~PARAMCD, ~DESCRIPTION,
    "WEIGHT", "Weight", "WEIGHT", "Weight (kg)",
    "HEIGHT", "Height", "HEIGHT", "Height (cm)",
    "BMI", "Body Mass Index", "BMI", "Body Mass Index(kg/m^2)"
  )

  attr(param_lookup$VSTESTCD, "label") <- "Vital Signs Test Short Name"
  attr(param_lookup$VSTEST, "label") <- "Vital Signs Test Name"

  expect_snapshot(
    act_vs_param <- derive_vars_merged_lookup(
      vs,
      dataset_add = param_lookup,
      by_vars = exprs(VSTESTCD, VSTEST),
      new_var = exprs(PARAMCD, PARAM = DESCRIPTION),
      print_not_mapped = TRUE
    )
  )

  actual <- get_not_mapped()

  expected <- left_join(vs, param_lookup, by = c("VSTESTCD", "VSTEST")) %>%
    rename(PARAM = DESCRIPTION) %>%
    filter(is.na(PARAMCD)) %>%
    select(VSTESTCD, VSTEST) %>%
    distinct()

  expect_dfs_equal(
    base = expected,
    compare = actual,
    keys = c("VSTESTCD", "VSTEST")
  )
})

# derive_var_merged_summary ----
## Test 22: dataset == dataset_add, no filter ----
test_that("derive_var_merged_summary Test 22: dataset == dataset_add, no filter", {
  expected <- tibble::tribble(
    ~AVISIT,  ~ASEQ, ~AVAL, ~MEANVIS,
    "WEEK 1",     1,    10,       10,
    "WEEK 1",     2,    NA,       10,
    "WEEK 2",     3,    NA,       NA,
    "WEEK 3",     4,    42,       42,
    "WEEK 4",     5,    12,       13,
    "WEEK 4",     6,    12,       13,
    "WEEK 4",     7,    15,       13
  )

  adbds <- select(expected, -MEANVIS)

  expect_dfs_equal(
    base = expected,
    compare = derive_var_merged_summary(
      adbds,
      dataset_add = adbds,
      by_vars = exprs(AVISIT),
      new_vars = exprs(MEANVIS = mean(AVAL, na.rm = TRUE))
    ),
    keys = c("AVISIT", "ASEQ")
  )
})

## Test 23: dataset != dataset_add, filter ----
test_that("derive_var_merged_summary Test 23: dataset != dataset_add, filter", {
  expected <- tibble::tribble(
    ~USUBJID, ~MEANPBL,
    "1",          13.5,
    "2",            NA,
    "3",          42.0
  )

  adbds <- tibble::tribble(
    ~USUBJID, ~ADY, ~AVAL,
    "1",        -3,    10,
    "1",         2,    12,
    "1",         8,    15,
    "3",         4,    42
  )

  adsl <- select(expected, -MEANPBL)

  expect_dfs_equal(
    base = expected,
    compare = derive_var_merged_summary(
      adsl,
      dataset_add = adbds,
      by_vars = exprs(USUBJID),
      new_vars = exprs(MEANPBL = mean(AVAL, na.rm = TRUE)),
      filter_add = ADY > 0
    ),
    keys = c("USUBJID")
  )
})

## Test 24: by_vars with rename ----
test_that("derive_var_merged_summary Test 24: by_vars with rename", {
  expected <- tibble::tribble(
    ~AVISIT,  ~ASEQ, ~AVAL, ~MEANVIS,
    "WEEK 1",     1,    10,       10,
    "WEEK 1",     2,    NA,       10,
    "WEEK 2",     3,    NA,       NA,
    "WEEK 3",     4,    42,       42,
    "WEEK 4",     5,    12,       13,
    "WEEK 4",     6,    12,       13,
    "WEEK 4",     7,    15,       13
  )
  adbds <- select(expected, -MEANVIS)
  adbds1 <- select(expected, -MEANVIS) %>%
    rename(VISIT = AVISIT)

  expect_dfs_equal(
    base = expected,
    compare = derive_var_merged_summary(
      adbds,
      dataset_add = adbds1,
      by_vars = exprs(AVISIT = VISIT),
      new_vars = exprs(MEANVIS = mean(AVAL, na.rm = TRUE))
    ),
    keys = c("AVISIT", "ASEQ")
  )
})

## Test 25: merge relationship as 'many-to-one' ----
test_that("derive_var_merged_summary Test 25: merge relationship as 'many-to-one'", {
  actual <- derive_vars_merged(advs,
    dataset_add = adsl,
    by_vars = exprs(USUBJID),
    new_vars = exprs(SEX),
    relationship = "many-to-one"
  )

  expected <- left_join(advs, select(adsl, USUBJID, SEX), by = "USUBJID")

  expect_dfs_equal(
    base = expected,
    compare = actual,
    keys = c("USUBJID", "AVISIT")
  )
})

## Test 26: error incorrect 'one-to-one' ----
test_that("derive_var_merged_summary Test 26: error incorrect 'one-to-one'", {
  expect_snapshot(
    derive_vars_merged(advs,
      dataset_add = adsl,
      by_vars = exprs(USUBJID),
      new_vars = exprs(SEX),
      relationship = "one-to-one"
    ),
    error = TRUE
  )
})

## Test 27: merge sel vars 'one-to-one' ----
test_that("derive_var_merged_summary Test 27: merge sel vars 'one-to-one'", {
  expect_snapshot(
    derive_vars_merged(adsl,
      dataset_add = advs,
      by_vars = exprs(USUBJID),
      new_vars = exprs(WEIGHTBL = AVAL),
      filter_add = AVISIT == "BASELINE",
      relationship = "one-to-one"
    )
  )
})

## Test 28: error if no summary function ----
test_that("derive_var_merged_summary Test 28: error if no summary function", {
  adbds <- tibble::tribble(
    ~AVISIT,  ~ASEQ, ~AVAL,
    "WEEK 1",     1,    10,
    "WEEK 1",     2,    NA,
    "WEEK 2",     3,    NA,
    "WEEK 3",     4,    42,
    "WEEK 4",     5,    12,
    "WEEK 4",     6,    12,
    "WEEK 4",     7,    15
  )

  expect_snapshot(
    derive_var_merged_summary(
      adbds,
      dataset_add = adbds,
      by_vars = exprs(AVISIT),
      new_vars = exprs(MEANVIS = AVAL / 2)
    ),
    error = TRUE
  )
})

## Test 29: test get_not_mapped with unmapped records ----
test_that("derive_vars_merged_lookup Test 29: test get_not_mapped with unmapped records", {
  # Create a lookup table that doesn't include BMI
  param_lookup <- tibble::tribble(
    ~VSTESTCD, ~VSTEST, ~PARAMCD,  ~DESCRIPTION,
    "HEIGHT", "Height", "HEIGHT", "Height (cm)",
    "WEIGHT", "Weight", "WEIGHT", "Weight (kg)",
  )

  # Run derive_vars_merged_lookup with print_not_mapped = TRUE
  actual <- derive_vars_merged_lookup(
    vs,
    dataset_add = param_lookup,
    by_vars = exprs(VSTESTCD, VSTEST),
    new_vars = exprs(PARAMCD, PARAM = DESCRIPTION),
    print_not_mapped = TRUE
  )

  # Get the not mapped records
  not_mapped <- get_not_mapped()

  # Verify the not mapped records
  expected_not_mapped <- tibble::tribble(
    ~VSTESTCD,                    ~VSTEST,
    "DIABP",   "Diastolic Blood Pressure"
  )

  expect_dfs_equal(
    base = expected_not_mapped,
    compare = not_mapped,
    keys = c("VSTESTCD", "VSTEST")
  )
})
