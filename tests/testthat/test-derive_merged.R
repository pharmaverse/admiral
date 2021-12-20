
adsl <- tibble::tribble(
  ~USUBJID, ~SEX, ~COUNTRY,
  "ST42-1", "F",  "AUT",
  "ST42-2", "M",  "MWI",
  "ST42-3", "M",  "NOR",
  "ST42-4", "F",  "UGA"
) %>% mutate(STUDYID = "ST42")

advs <- tibble::tribble(
  ~USUBJID, ~PARAMCD, ~AVISIT,    ~AVAL,
  "ST42-1", "WEIGHT", "BASELINE", 66,
  "ST42-1", "WEIGHT", "WEEK 2",   68,
  "ST42-2", "WEIGHT", "BASELINE", 88,
  "ST42-2", "WEIGHT", "WEEK 2",   85,
  "ST42-3", "WEIGHT", "WEEK 2",   55
) %>% mutate(STUDYID = "ST42")

ex <- tibble::tribble(
  ~USUBJID, ~EXSTDTC,
  "ST42-1", "2020-12-07",
  "ST42-1", "2020-12-14",
  "ST42-2", "2021-01-12T12:00:00",
  "ST42-2", "2021-01-26T13:21",
  "ST42-3", "2021-03-02"
) %>% mutate(STUDYID = "ST42")
# derive_vars_merged ----
## derive_vars_merged: merge all variables ----
test_that("derive_vars_merged: merge all variables", {
  actual <- derive_vars_merged(advs,
                               dataset_add = adsl,
                               by_vars = vars(STUDYID, USUBJID))

  expected <- left_join(advs, adsl, by = c("STUDYID", "USUBJID"))

  expect_dfs_equal(base = expected,
                   compare = actual,
                   keys = c("USUBJID", "AVISIT"))
})

## derive_vars_merged: merge selected variables ----
test_that("derive_vars_merged: merge selected variables", {
  actual <- derive_vars_merged(advs,
                               dataset_add = adsl,
                               by_vars = vars(USUBJID),
                               new_vars = vars(SEX))

  expected <- left_join(advs, select(adsl, USUBJID, SEX), by = "USUBJID")

  expect_dfs_equal(base = expected,
                   compare = actual,
                   keys = c("USUBJID", "AVISIT"))
})

## derive_vars_merged: merge last value ----
test_that("derive_vars_merged: merge last value", {
  actual <- derive_vars_merged(adsl,
                               dataset_add = advs,
                               order = vars(AVAL),
                               by_vars = vars(STUDYID, USUBJID),
                               new_vars = vars(WEIGHTBL = AVAL),
                               mode = "last")
  expected <- adsl %>% mutate(WEIGHTBL = c(68, 88, 55, NA))

  expect_dfs_equal(base = expected,
                   compare = actual,
                   keys = c("USUBJID"))
})

## derive_vars_merged: error if variable in both datasets ----
test_that("derive_vars_merged: error if variable in both datasets", {
  expect_error(derive_vars_merged(advs,
                                  dataset_add = adsl,
                                  by_vars = vars(USUBJID)),
               regexp = "")
})

# derive_vars_merged_dt ----
## derive_vars_merged_dt: merge first date ----
test_that("derive_vars_merged_dt: merge first date", {
  actual <- derive_vars_merged_dt(adsl,
                               dataset_add = ex,
                               order = vars(TRTSDT),
                               by_vars = vars(STUDYID, USUBJID),
                               dtc = EXSTDTC,
                               new_vars_prefix = "TRTS",
                               mode = "first")
  expected <- adsl %>% mutate(TRTSDT = ymd(c("2020-12-07", "2021-01-12", "2021-03-02", NA)),
                              TRTSDTF = NA_character_)

  expect_dfs_equal(base = expected,
                   compare = actual,
                   keys = c("USUBJID"))
})

# derive_vars_merged_dtm ----
## derive_vars_merged_dtm: merge first date ----
test_that("derive_vars_merged_dt: merge first date", {
  actual <- derive_vars_merged_dtm(
    adsl,
    dataset_add = ex,
    order = vars(TRTSDTM),
    by_vars = vars(STUDYID, USUBJID),
    dtc = EXSTDTC,
    new_vars_prefix = "TRTS",
    time_imputation = "first",
    mode = "first"
  )
  expected <- adsl %>% mutate(TRTSDTM = as_iso_dtm(ymd_hms(c("2020-12-07T00:00:00", "2021-01-12T12:00:00", "2021-03-02T00:00:00", NA))),
                              TRTSTMF = c("H", NA, "H", NA))

  expect_dfs_equal(base = expected,
                   compare = actual,
                   keys = c("USUBJID"))
})

# derive_var_merged_cat ----
## derive_var_merged_cat: merge categorized variable ----
test_that("derive_vars_merged_cat: merge categorized variable", {
  get_region <- function(x) {
    if_else(x %in% c("AUT", "NOR"), "EUROPE", "AFRICA")
  }

  actual <- derive_var_merged_cat(
    advs,
    dataset_add = adsl,
    by_vars = vars(USUBJID),
    new_var = REGION,
    source_var = COUNTRY,
    cat_fun = get_region
  )

  expected <- left_join(advs, select(adsl, USUBJID, COUNTRY), by = "USUBJID") %>%
    mutate(REGION = get_region(COUNTRY)) %>%
    select(-COUNTRY)


  expect_dfs_equal(base = expected,
                   compare = actual,
                   keys = c("USUBJID", "AVISIT"))
})

# derive_var_merged_exist_flag ----
## derive_var_merged_exist_flag: merge existence flag ----
test_that("derive_vars_merged_cat: merge existence flag", {
  actual <- derive_var_merged_exist_flag(
    adsl,
    dataset_add = advs,
    by_vars = vars(USUBJID),
    new_var = VSEVALFL,
    condition = AVISIT == 'BASELINE'
  )

  expected <-
    mutate(adsl, VSEVALFL = c("Y", "Y", NA_character_, NA_character_))


  expect_dfs_equal(base = expected,
                   compare = actual,
                   keys = "USUBJID")
})
