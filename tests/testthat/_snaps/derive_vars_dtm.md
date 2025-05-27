# impute_dtc_dtm Test 10: min_dates length mismatch provides error

    Code
      impute_dtc_dtm(c("2020-12", NA_character_), min_dates = list(c(ymd_hms(
        "2020-12-06T12:12:12")), c(ymd_hms("2020-11-11T11:11:11"))),
      highest_imputation = "Y")
    Condition
      Error in `restrict_imputed_dtc_dtm()`:
      ! Length of `min_dates` do not match length of dates to be imputed.

# impute_dtc_dtm Test 11: max_dates length mismatch provides error

    Code
      impute_dtc_dtm(c("2020-12", NA_character_), max_dates = list(c(ymd_hms(
        "2020-12-06T12:12:12")), c(ymd_hms("2020-11-11T11:11:11"))),
      highest_imputation = "Y", date_imputation = "last")
    Condition
      Error in `restrict_imputed_dtc_dtm()`:
      ! Length of `max_dates` do not match length of dates to be imputed.

# impute_dtc_dtm Test 12: Error if null min/max_dates when highest_imputation = Y

    Code
      impute_dtc_dtm(c("2020-12", NA_character_), highest_imputation = "Y")
    Condition
      Error in `assert_highest_imputation()`:
      ! If `highest_imputation = "Y"` is specified, `min_dates` or `max_dates` must be specified respectively.

# impute_dtc_dtm Test 13: wrong input to `date_imputation`

    Code
      impute_dtc_dtm(dtc = c("2020-12", "2020-11", NA_character_, "2020-02-02"),
      highest_imputation = "D", date_imputation = "15", time_imputation = "last")
    Condition
      Error in `assert_date_imputation()`:
      ! Argument `date_imputation` must be equal to one of "first", "mid", or "last".

---

    Code
      impute_dtc_dtm(dtc = c("2020-12", "2020-11", NA_character_, "2020-02-02",
        "2020"), highest_imputation = "M", date_imputation = "10:12",
      time_imputation = "last")
    Condition
      Error in `assert_date_imputation()`:
      ! If `highest_imputation = "M"` is specified, `date_imputation` must be one of "first", "mid", "last" or a format with month and day specified as "mm-dd": e.g. "06-15"

---

    Code
      impute_dtc_dtm(dtc = c("2020-12", "2020-11", NA_character_, "2020-02-02",
        "2020"), highest_imputation = "M", date_imputation = "10", time_imputation = "last")
    Condition
      Error in `assert_date_imputation()`:
      ! If `highest_imputation = "M"` is specified, `date_imputation` must be one of "first", "mid", "last" or a format with month and day specified as "mm-dd": e.g. "06-15"

---

    Code
      impute_dtc_dtm(dtc = c("2020-12", "2020-11", NA_character_, "2020-02-02",
        "2020"), highest_imputation = "D", date_imputation = "15", time_imputation = "last")
    Condition
      Error in `assert_date_imputation()`:
      ! Argument `date_imputation` must be equal to one of "first", "mid", or "last".

---

    Code
      impute_dtc_dtm(dtc = c("2020-12", "2020-11", NA_character_, "2020-02-02",
        "2020"), highest_imputation = "Y", date_imputation = "2020-02-02",
      time_imputation = "last")
    Condition
      Error in `assert_highest_imputation()`:
      ! Argument `date_imputation` must be equal to one of "first" or "last".

---

    Code
      impute_dtc_dtm(dtc = c("2020-12", "2020-11", NA_character_, "2020-02-02",
        "2020"), highest_imputation = "M", date_imputation = "first",
      time_imputation = "WRONG")
    Condition
      Error in `assert_time_imputation()`:
      ! `time_imputation` must be one of "first", "last" or time specified as "hh:mm:ss": e.g. "12:00:00"

---

    Code
      impute_dtc_dtm(dtc = c("2020-12", "2020-11", NA_character_, "2020-02-02",
        "2020"), highest_imputation = "M", date_imputation = "first",
      time_imputation = "12:12:61")
    Condition
      Error in `assert_time_imputation()`:
      ! `time_imputation` must be one of "first", "last" or time specified as "hh:mm:ss": e.g. "12:00:00"

# compute_tmf Test 16: throws ERROR when ignore_seconds_flag  = T and seconds are present

    Code
      compute_tmf(dtc = c("2020-11-11T11:11:11", "2020-11-11T11:11"), dtm = ymd_hms(c(
        "2020-11-11T11:11:11", "2020-11-11T11:11:00")), ignore_seconds_flag = TRUE)
    Condition
      Error in `compute_tmf()`:
      ! Seconds detected in data while `ignore_seconds_flag` is invoked

# derive_vars_dtm Test 23: No re-derivation is done if --DTF variable already exists

    Code
      actual_output <- derive_vars_dtm(mutate(input, ASTDTF = c(NA, NA, NA, NA, "D",
        "MD", "M")), new_vars_prefix = "AST", dtc = XXSTDTC, highest_imputation = "M",
      date_imputation = "first")
    Message
      The `ASTDTF` variable is already present in the input dataset and will not be re-derived.

# derive_vars_dtm Test 26: Error for highest_imputation = Y & max_dates but date_imputation = first

    Code
      (data.frame(AESTDTC = c(NA_character_, NA_character_), TRTSDTM = c(ymd_hms(
        "2022-01-01 23:59:59"), NA)) %>% mutate(AESTDTC = as.character(AESTDTC)) %>%
        derive_vars_dtm(dtc = AESTDTC, new_vars_prefix = "AST", highest_imputation = "Y",
          date_imputation = "first", time_imputation = "first", flag_imputation = "both",
          max_dates = exprs(TRTSDTM)))
    Condition
      Error in `assert_highest_imputation()`:
      ! If `highest_imputation = "Y"` and `date_imputation = "first"` is specified, `min_dates` must be specified.

# derive_vars_dtm Test 28: Error for highest_imputation = Y & min_dates but date_imputation = last

    Code
      data.frame(AESTDTC = c(NA_character_, NA_character_), TRTSDTM = c(ymd_hms(
        "2022-01-01 23:59:59"), NA)) %>% mutate(AESTDTC = as.character(AESTDTC)) %>%
        derive_vars_dtm(dtc = AESTDTC, new_vars_prefix = "AST", highest_imputation = "Y",
          date_imputation = "last", time_imputation = "last", flag_imputation = "both",
          min_dates = exprs(TRTSDTM))
    Condition
      Error in `assert_highest_imputation()`:
      ! If `highest_imputation = "Y"` and `date_imputation = "last"` is specified, `max_dates` must be specified.

# derive_vars_dtm Test 29: NA imputation for highest_imputation = Y but null min/max dates fails

    Code
      data.frame(AESTDTC = c(NA_character_, NA_character_), TRTSDTM = c(ymd_hms(
        "2022-01-01 23:59:59"), NA)) %>% mutate(AESTDTC = as.character(AESTDTC)) %>%
        derive_vars_dtm(dtc = AESTDTC, new_vars_prefix = "AST", highest_imputation = "Y",
          date_imputation = "first", time_imputation = "first", flag_imputation = "both")
    Condition
      Error in `assert_highest_imputation()`:
      ! If `highest_imputation = "Y"` is specified, `min_dates` or `max_dates` must be specified respectively.

# derive_vars_dtm Test 32: catch ignore_seconds_flag error

    Code
      derive_vars_dtm(input, new_vars_prefix = "AST", dtc = XXSTDTC,
        highest_imputation = "M", ignore_seconds_flag = TRUE)
    Condition
      Error in `derive_vars_dtm()`:
      ! Seconds detected in data while `ignore_seconds_flag` is invoked

