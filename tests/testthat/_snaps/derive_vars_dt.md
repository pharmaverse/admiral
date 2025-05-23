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

# derive_vars_dt Test 10: Error if null min/max_dates when highest_imputation = Y

    Code
      impute_dtc_dt(input, highest_imputation = "Y")
    Condition
      Error in `assert_highest_imputation()`:
      ! If `highest_imputation = "Y"` is specified, `min_dates` or `max_dates` must be specified respectively.

# derive_vars_dt Test 12: wrong input to `date_imputation`

    Code
      impute_dtc_dt(dtc = input, highest_imputation = "M", date_imputation = "13-01")
    Condition
      Error in `assert_date_imputation()`:
      ! If `highest_imputation = "M"` is specified, `date_imputation` must be one of "first", "mid", "last" or a format with month and day specified as "mm-dd": e.g. "06-15"

---

    Code
      impute_dtc_dt(dtc = input, highest_imputation = "M", date_imputation = "12:01")
    Condition
      Error in `assert_date_imputation()`:
      ! If `highest_imputation = "M"` is specified, `date_imputation` must be one of "first", "mid", "last" or a format with month and day specified as "mm-dd": e.g. "06-15"

---

    Code
      impute_dtc_dt(dtc = input, highest_imputation = "M", date_imputation = "01")
    Condition
      Error in `assert_date_imputation()`:
      ! If `highest_imputation = "M"` is specified, `date_imputation` must be one of "first", "mid", "last" or a format with month and day specified as "mm-dd": e.g. "06-15"

---

    Code
      impute_dtc_dt(dtc = input, highest_imputation = "D", date_imputation = "01")
    Condition
      Error in `assert_date_imputation()`:
      ! Argument `date_imputation` must be equal to one of "first", "mid", or "last".

---

    Code
      impute_dtc_dt(dtc = input, highest_imputation = "Y", date_imputation = "2006-01-01")
    Condition
      Error in `assert_highest_imputation()`:
      ! Argument `date_imputation` must be equal to one of "first" or "last".

# derive_vars_dt Test 21: NA imputation for highest_imputation = Y & max_dates but date_imputation = first

    Code
      data.frame(AESTDTC = c(NA_character_, NA_character_), TRTSDT = c(ymd(
        "2022-01-01"), NA)) %>% mutate(AESTDTC = as.character(AESTDTC)) %>%
        derive_vars_dt(dtc = AESTDTC, new_vars_prefix = "AST", highest_imputation = "Y",
          date_imputation = "first", flag_imputation = "auto", max_dates = exprs(
            TRTSDT))
    Condition
      Error in `assert_highest_imputation()`:
      ! If `highest_imputation = "Y"` and `date_imputation = "first"` is specified, `min_dates` must be specified.

# derive_vars_dt Test 23: Error for highest_imputation = Y & min_dates but date_imputation = last

    Code
      data.frame(AESTDTC = c(NA_character_, NA_character_), TRTSDT = c(ymd(
        "2022-01-01"), NA)) %>% mutate(AESTDTC = as.character(AESTDTC)) %>%
        derive_vars_dt(dtc = AESTDTC, new_vars_prefix = "AST", highest_imputation = "Y",
          date_imputation = "last", flag_imputation = "auto", min_dates = exprs(
            TRTSDT))
    Condition
      Error in `assert_highest_imputation()`:
      ! If `highest_imputation = "Y"` and `date_imputation = "last"` is specified, `max_dates` must be specified.

# derive_vars_dt Test 24: Error for highest_imputation = Y but null min/max dates fails

    Code
      data.frame(AESTDTC = c(NA_character_, NA_character_), TRTSDT = c(ymd(
        "2022-01-01"), NA)) %>% mutate(AESTDTC = as.character(AESTDTC)) %>%
        derive_vars_dt(dtc = AESTDTC, new_vars_prefix = "AST", highest_imputation = "Y",
          date_imputation = "first", flag_imputation = "auto")
    Condition
      Error in `assert_highest_imputation()`:
      ! If `highest_imputation = "Y"` is specified, `min_dates` or `max_dates` must be specified respectively.

