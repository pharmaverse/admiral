# derive_vars_dt Test 8: min_dates length mismatch provides error

    Code
      impute_dtc_dt(input, min_dates = list(c(ymd("2019-07-06")), c(ymd("2019-06-06"))),
      highest_imputation = "Y", date_imputation = "first")
    Condition
      Error in `restrict_imputed_dtc_dt()`:
      ! Length of `min_dates` do not match length of dates to be imputed.

# derive_vars_dt Test 9: max_dates length mismatch provides error

    Code
      impute_dtc_dt(input, max_dates = list(c(ymd("2019-07-06")), c(ymd("2019-06-06"))),
      highest_imputation = "Y", date_imputation = "last")
    Condition
      Error in `restrict_imputed_dtc_dt()`:
      ! Length of `max_dates` do not match length of dates to be imputed.

# derive_vars_dt Test 19: NA imputation for highest_imputation = Y & max_dates but date_imputation = first

    Code
      data.frame(AESTDTC = c(NA_character_, NA_character_), TRTSDT = c(ymd(
        "2022-01-01"), NA)) %>% mutate(AESTDTC = as.character(AESTDTC)) %>%
        derive_vars_dt(dtc = AESTDTC, new_vars_prefix = "AST", highest_imputation = "Y",
          date_imputation = "first", flag_imputation = "auto", max_dates = exprs(
            TRTSDT))
    Condition
      Warning:
      If `highest_impuation = "Y"` and `date_imputation = "first"` is specified, `min_dates` should be specified.
    Output
        AESTDTC     TRTSDT ASTDT ASTDTF
      1    <NA> 2022-01-01  <NA>   <NA>
      2    <NA>       <NA>  <NA>   <NA>

# derive_vars_dt Test 21: NA imputation for highest_imputation = Y & min_dates but date_imputation = last

    Code
      data.frame(AESTDTC = c(NA_character_, NA_character_), TRTSDT = c(ymd(
        "2022-01-01"), NA)) %>% mutate(AESTDTC = as.character(AESTDTC)) %>%
        derive_vars_dt(dtc = AESTDTC, new_vars_prefix = "AST", highest_imputation = "Y",
          date_imputation = "last", flag_imputation = "auto", min_dates = exprs(
            TRTSDT))
    Condition
      Warning:
      If `highest_impuation = "Y"` and `date_imputation = "last"` is specified, `max_dates` should be specified.
    Output
        AESTDTC     TRTSDT ASTDT ASTDTF
      1    <NA> 2022-01-01  <NA>   <NA>
      2    <NA>       <NA>  <NA>   <NA>

# derive_vars_dt Test 22: NA imputation for highest_imputation = Y but null min/max dates fails

    Code
      data.frame(AESTDTC = c(NA_character_, NA_character_), TRTSDT = c(ymd(
        "2022-01-01"), NA)) %>% mutate(AESTDTC = as.character(AESTDTC)) %>%
        derive_vars_dt(dtc = AESTDTC, new_vars_prefix = "AST", highest_imputation = "Y",
          date_imputation = "first", flag_imputation = "auto")
    Condition
      Error in `derive_vars_dt()`:
      ! If `highest_impuation = "Y"` is specified, `min_dates` or `max_dates` must be specified respectively.

