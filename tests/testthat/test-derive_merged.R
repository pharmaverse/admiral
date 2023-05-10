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
    match_flag = matched
  )
  expected <- adsl %>% mutate(
    WEIGHTBL = c(68, 88, 55, NA),
    matched = c(TRUE, TRUE, TRUE, NA)
  )

  expect_dfs_equal(
    base = expected,
    compare = actual,
    keys = c("USUBJID")
  )
})

## Test 4: error if variable in both datasets ----
test_that("derive_vars_merged Test 4: error if variable in both datasets", {
  expect_error(
    derive_vars_merged(advs,
      dataset_add = adsl,
      by_vars = exprs(USUBJID)
    ),
    regexp = ""
  )
})

## Test 5: by_vars with rename ----
test_that("derive_vars_merged Test 5: by_vars with rename", {
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

## Test 6: expressions for new_vars and missing_values ----
test_that("derive_vars_merged Test 6: expressions for new_vars and missing_values", {
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

## Test 7: use new variables in filter_add and order ----
test_that("derive_vars_merged Test 7: use new variables in filter_add and order", {
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

## Test 8: warning if not unique w.r.t the by variables and the order ----
test_that("derive_vars_merged Test 8: warning if not unique w.r.t the by variables and the order", {
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

## Test 9: error if not unique w.r.t the by variables and the order ----
test_that("derive_vars_merged Test 9: error if not unique w.r.t the by variables and the order", {
  expect_error(
    actual <- derive_vars_merged(advs,
      dataset_add = adsl2,
      by_vars = exprs(STUDYID, USUBJID = ID),
      order = exprs(ID),
      mode = "last",
      check_type = "error",
      duplicate_msg = "Duplicate records present!"
    ),
    regexp = ""
  )
})

## Test 10: error if variables in missing_values but not in new_vars ----
test_that("derive_vars_merged Test 10: error if variables in missing_values but not in new_vars", {
  expect_error(
    derive_vars_merged(
      adsl,
      dataset_add = advs,
      by_vars = exprs(USUBJID),
      order = exprs(AVISIT),
      new_vars = exprs(LASTVIS = str_to_upper(AVISIT)),
      mode = "last",
      missing_values = exprs(LASTVIS = "UNKNOWN", LASTVISN = -1)
    ),
    regexp = "The variables `LASTVISN` were specified for `missing_values` but not for `new_vars`.",
    fixed = TRUE
  )
})

# derive_var_merged_cat ----

## Test 11: deprecation warning ----
test_that("derive_var_merged_cat Test 11: deprecation warning", {
  get_vscat <- function(x) {
    if_else(x == "BASELINE", "BASELINE", "POST-BASELINE")
  }

  expect_warning(
    derive_var_merged_cat(
      adsl,
      dataset_add = advs,
      by_vars = exprs(USUBJID),
      new_var = LSTVSCAT,
      source_var = AVISIT,
      cat_fun = get_vscat,
      order = exprs(AVISIT),
      mode = "last",
      missing_value = "MISSING"
    ),
    class = "lifecycle_warning_deprecated"
  )
})

## Test 12: define value for non-matched by groups ----
test_that("derive_var_merged_cat Test 12: define value for non-matched by groups", {
  get_vscat <- function(x) {
    if_else(x == "BASELINE", "BASELINE", "POST-BASELINE")
  }

  actual <- suppress_warning(
    derive_var_merged_cat(
      adsl,
      dataset_add = advs,
      by_vars = exprs(USUBJID),
      new_var = LSTVSCAT,
      source_var = AVISIT,
      cat_fun = get_vscat,
      order = exprs(AVISIT),
      mode = "last",
      missing_value = "MISSING"
    ),
    regexpr = "was deprecated"
  )

  expected <-
    mutate(adsl,
      LSTVSCAT = c("POST-BASELINE", "BASELINE", "POST-BASELINE", "MISSING")
    )

  expect_dfs_equal(
    base = expected,
    compare = actual,
    keys = c("USUBJID")
  )
})

## Test 13: by_vars with rename ----
test_that("derive_var_merged_cat Test 13: by_vars with rename", {
  get_region <- function(x) {
    if_else(x %in% c("AUT", "NOR"), "EUROPE", "AFRICA")
  }

  actual <- suppress_warning(
    derive_var_merged_cat(
      advs,
      dataset_add = adsl1,
      by_vars = exprs(USUBJID = ID),
      new_var = REGION,
      source_var = COUNTRY,
      cat_fun = get_region,
      filter_add = SEX == "M"
    ),
    regexpr = "was deprecated"
  )

  adsl_1 <- adsl1 %>% filter(SEX == "M")
  expected <- left_join(advs, select(adsl_1, ID, COUNTRY), by = c("USUBJID" = "ID")) %>%
    mutate(REGION = get_region(COUNTRY)) %>%
    mutate(REGION = case_when(
      !is.na(COUNTRY) ~ REGION,
      TRUE ~ NA_character_
    )) %>%
    select(-COUNTRY)


  expect_dfs_equal(
    base = expected,
    compare = actual,
    keys = c("USUBJID", "AVISIT")
  )
})

# derive_var_merged_exist_flag ----
## Test 14: merge existence flag ----
test_that("derive_var_merged_exist_flag Test 14: merge existence flag", {
  actual <- derive_var_merged_exist_flag(
    adsl,
    dataset_add = advs,
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

## Test 15: by_vars with rename ----
test_that("derive_var_merged_exist_flag Test 15: by_vars with rename", {
  actual <- derive_var_merged_exist_flag(
    adsl,
    dataset_add = advs1,
    by_vars = exprs(USUBJID = ID),
    new_var = VSEVALFL,
    condition = AVISIT == "BASELINE",
  )

  expected <-
    mutate(adsl, VSEVALFL = c("Y", "Y", NA_character_, NA_character_))


  expect_dfs_equal(
    base = expected,
    compare = actual,
    keys = "USUBJID"
  )
})

# derive_var_merged_character ----

## Test 16: deprecation warning ----
test_that("derive_var_merged_character Test 16: deprecation warning", {
  expect_warning(
    derive_var_merged_character(
      adsl,
      dataset_add = advs,
      by_vars = exprs(USUBJID),
      order = exprs(AVISIT),
      new_var = LASTVIS,
      source_var = AVISIT,
      mode = "last"
    ),
    class = "lifecycle_warning_deprecated"
  )
})

## Test 17: no transformation ----
test_that("derive_var_merged_character Test 17: no transformation", {
  actual <- suppress_warning(
    derive_var_merged_character(
      adsl,
      dataset_add = advs,
      by_vars = exprs(USUBJID),
      order = exprs(AVISIT),
      new_var = LASTVIS,
      source_var = AVISIT,
      mode = "last"
    ),
    regexpr = "was deprecated"
  )

  expected <-
    mutate(adsl, LASTVIS = c("Week 2", "BASELINE", "Week 4", NA_character_))


  expect_dfs_equal(
    base = expected,
    compare = actual,
    keys = "USUBJID"
  )
})

## Test 18: upper case ----
test_that("derive_var_merged_character Test 18: upper case", {
  actual <- suppress_warning(
    derive_var_merged_character(
      adsl,
      dataset_add = advs,
      by_vars = exprs(USUBJID),
      order = exprs(AVISIT),
      new_var = LASTVIS,
      source_var = AVISIT,
      mode = "last",
      case = "upper",
      missing_value = "UNKNOWN"
    ),
    regexpr = "was deprecated"
  )

  expected <-
    mutate(adsl, LASTVIS = c("WEEK 2", "BASELINE", "WEEK 4", "UNKNOWN"))


  expect_dfs_equal(
    base = expected,
    compare = actual,
    keys = "USUBJID"
  )
})

## Test 19: lower case ----
test_that("derive_var_merged_character Test 19: lower case", {
  actual <- suppress_warning(
    derive_var_merged_character(
      adsl,
      dataset_add = advs,
      by_vars = exprs(USUBJID),
      order = exprs(AVISIT),
      new_var = LASTVIS,
      source_var = AVISIT,
      mode = "last",
      case = "lower"
    ),
    regexpr = "was deprecated"
  )

  expected <-
    mutate(adsl, LASTVIS = c("week 2", "baseline", "week 4", NA_character_))


  expect_dfs_equal(
    base = expected,
    compare = actual,
    keys = "USUBJID"
  )
})

## Test 20: title case ----
test_that("derive_var_merged_character Test 20: title case", {
  actual <- suppress_warning(
    derive_var_merged_character(
      adsl,
      dataset_add = advs,
      by_vars = exprs(USUBJID),
      order = exprs(AVISIT),
      new_var = LASTVIS,
      source_var = AVISIT,
      mode = "last",
      case = "title"
    ),
    regexpr = "was deprecated"
  )

  expected <-
    mutate(adsl, LASTVIS = c("Week 2", "Baseline", "Week 4", NA_character_))


  expect_dfs_equal(
    base = expected,
    compare = actual,
    keys = "USUBJID"
  )
})

## Test 21: by_vars with rename ----
test_that("derive_var_merged_character Test 21: by_vars with rename", {
  actual <- suppress_warning(
    derive_var_merged_character(
      adsl,
      dataset_add = advs1,
      by_vars = exprs(USUBJID = ID),
      order = exprs(AVISIT),
      new_var = LASTVIS,
      source_var = AVISIT,
      mode = "last"
    ),
    regexpr = "was deprecated"
  )

  expected <-
    mutate(adsl, LASTVIS = c("Week 2", "BASELINE", "Week 4", NA_character_))


  expect_dfs_equal(
    base = expected,
    compare = actual,
    keys = "USUBJID"
  )
})


# derive_vars_merged_lookup ----
## Test 22: merge lookup table ----
test_that("derive_vars_merged_lookup Test 22: merge lookup table", {
  param_lookup <- tibble::tribble(
    ~VSTESTCD, ~VSTEST, ~PARAMCD, ~DESCRIPTION,
    "WEIGHT", "Weight", "WEIGHT", "Weight (kg)",
    "HEIGHT", "Height", "HEIGHT", "Height (cm)",
    "BMI", "Body Mass Index", "BMI", "Body Mass Index(kg/m^2)"
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
    regex = "^List of `VSTESTCD` and `VSTEST` not mapped: .*$"
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
## Test 23:  all by_vars have records in the lookup table ----
test_that("derive_vars_merged_lookup Test 23:  all by_vars have records in the lookup table", {
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

## Test 24: by_vars with rename ----
test_that("derive_vars_merged_lookup Test 24: by_vars with rename", {
  param_lookup <- tibble::tribble(
    ~TESTCD, ~VSTEST, ~PARAMCD, ~DESCRIPTION,
    "WEIGHT", "Weight", "WEIGHT", "Weight (kg)",
    "HEIGHT", "Height", "HEIGHT", "Height (cm)",
    "BMI", "Body Mass Index", "BMI", "Body Mass Index(kg/m^2)"
  )

  attr(param_lookup$TESTCD, "label") <- "Vital Signs Test Short Name"
  attr(param_lookup$VSTEST, "label") <- "Vital Signs Test Name"

  expect_message(
    actual <- derive_vars_merged_lookup(
      vs,
      dataset_add = param_lookup,
      by_vars = exprs(VSTESTCD = TESTCD, VSTEST),
      new_var = exprs(PARAMCD, PARAM = DESCRIPTION),
      print_not_mapped = TRUE
    ),
    regex = "^List of `VSTESTCD` and `VSTEST` not mapped: .*$"
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
## Test 25: not all by_vars have records in the lookup table ----
test_that("get_not_mapped Test 25: not all by_vars have records in the lookup table", {
  param_lookup <- tibble::tribble(
    ~VSTESTCD, ~VSTEST, ~PARAMCD, ~DESCRIPTION,
    "WEIGHT", "Weight", "WEIGHT", "Weight (kg)",
    "HEIGHT", "Height", "HEIGHT", "Height (cm)",
    "BMI", "Body Mass Index", "BMI", "Body Mass Index(kg/m^2)"
  )

  attr(param_lookup$VSTESTCD, "label") <- "Vital Signs Test Short Name"
  attr(param_lookup$VSTEST, "label") <- "Vital Signs Test Name"

  expect_message(
    act_vs_param <- derive_vars_merged_lookup(
      vs,
      dataset_add = param_lookup,
      by_vars = exprs(VSTESTCD, VSTEST),
      new_var = exprs(PARAMCD, PARAM = DESCRIPTION),
      print_not_mapped = TRUE
    ),
    regex = "^List of `VSTESTCD` and `VSTEST` not mapped: .*$"
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
## Test 26: dataset == dataset_add, no filter ----
test_that("derive_var_merged_summary Test 26: dataset == dataset_add, no filter", {
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
      new_var = MEANVIS,
      analysis_var = AVAL,
      summary_fun = function(x) mean(x, na.rm = TRUE)
    ),
    keys = c("AVISIT", "ASEQ")
  )
})

## Test 27: dataset != dataset_add, filter ----
test_that("derive_var_merged_summary Test 27: dataset != dataset_add, filter", {
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
      new_var = MEANPBL,
      filter_add = ADY > 0,
      analysis_var = AVAL,
      summary_fun = function(x) mean(x, na.rm = TRUE)
    ),
    keys = c("USUBJID")
  )
})

## Test 28: by_vars with rename ----
test_that("derive_var_merged_summary Test 28: by_vars with rename", {
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
      new_var = MEANVIS,
      analysis_var = AVAL,
      summary_fun = function(x) mean(x, na.rm = TRUE)
    ),
    keys = c("AVISIT", "ASEQ")
  )
})
