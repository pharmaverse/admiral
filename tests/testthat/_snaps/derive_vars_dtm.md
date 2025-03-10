# compute_tmf Test 15: throws ERROR when ignore_seconds_flag  = T and seconds are present

    Code
      compute_tmf(dtc = c("2020-11-11T11:11:11", "2020-11-11T11:11"), dtm = ymd_hms(c(
        "2020-11-11T11:11:11", "2020-11-11T11:11:00")), ignore_seconds_flag = TRUE)
    Condition
      Error in `compute_tmf()`:
      ! Seconds detected in data while `ignore_seconds_flag` is invoked

# derive_vars_dtm Test 22: No re-derivation is done if --DTF variable already exists

    Code
      actual_output <- derive_vars_dtm(mutate(input, ASTDTF = c(NA, NA, NA, "D", "MD",
        "M")), new_vars_prefix = "AST", dtc = XXSTDTC, highest_imputation = "M",
      date_imputation = "first")
    Message
      The `ASTDTF` variable is already present in the input dataset and will not be re-derived.

# derive_vars_dtm Test 25: NA imputation for highest_imputation = Y & max_dates but date_imputation = first

    Code
      (data.frame(AESTDTC = c(NA_character_, NA_character_), TRTSDTM = c(ymd_hms(
        "2022-01-01 23:59:59"), NA)) %>% mutate(AESTDTC = as.character(AESTDTC)) %>%
        derive_vars_dtm(dtc = AESTDTC, new_vars_prefix = "AST", highest_imputation = "Y",
          date_imputation = "first", time_imputation = "first", flag_imputation = "both",
          max_dates = exprs(TRTSDTM)))
    Condition
      Warning:
      If `highest_impuation = "Y"` and `date_imputation = "first"` is specified, `min_dates` should be specified.
    Output
        AESTDTC             TRTSDTM ASTDTM ASTDTF ASTTMF
      1    <NA> 2022-01-01 23:59:59   <NA>   <NA>   <NA>
      2    <NA>                <NA>   <NA>   <NA>   <NA>

# derive_vars_dtm Test 27: NA imputation for highest_imputation = Y & min_dates but date_imputation = last

    Code
      data.frame(AESTDTC = c(NA_character_, NA_character_), TRTSDTM = c(ymd_hms(
        "2022-01-01 23:59:59"), NA)) %>% mutate(AESTDTC = as.character(AESTDTC)) %>%
        derive_vars_dtm(dtc = AESTDTC, new_vars_prefix = "AST", highest_imputation = "Y",
          date_imputation = "last", time_imputation = "last", flag_imputation = "both",
          min_dates = exprs(TRTSDTM))
    Condition
      Warning:
      If `highest_impuation = "Y"` and `date_imputation = "last"` is specified, `max_dates` should be specified.
    Output
        AESTDTC             TRTSDTM ASTDTM ASTDTF ASTTMF
      1    <NA> 2022-01-01 23:59:59   <NA>   <NA>   <NA>
      2    <NA>                <NA>   <NA>   <NA>   <NA>

# derive_vars_dtm Test 28: NA imputation for highest_imputation = Y but null min/max dates fails

    Code
      data.frame(AESTDTC = c(NA_character_, NA_character_), TRTSDTM = c(ymd_hms(
        "2022-01-01 23:59:59"), NA)) %>% mutate(AESTDTC = as.character(AESTDTC)) %>%
        derive_vars_dtm(dtc = AESTDTC, new_vars_prefix = "AST", highest_imputation = "Y",
          date_imputation = "first", time_imputation = "first", flag_imputation = "both")
    Condition
      Error in `derive_vars_dtm()`:
      ! If `highest_impuation = "Y"` is specified, `min_dates` or `max_dates` must be specified respectively.

# derive_vars_dtm Test 31: catch ignore_seconds_flag error

    Code
      derive_vars_dtm(input, new_vars_prefix = "AST", dtc = XXSTDTC,
        ignore_seconds_flag = TRUE)
    Output
      # A tibble: 6 x 3
        XXSTDTC          ASTDTM              ASTTMF
        <chr>            <dttm>              <chr> 
      1 2019-07-18T15:25 2019-07-18 15:25:00 <NA>  
      2 2019-07-18T15    2019-07-18 15:00:00 M     
      3 2019-07-18       2019-07-18 00:00:00 H     
      4 2019-02          NA                  <NA>  
      5 2019             NA                  <NA>  
      6 2019---07        NA                  <NA>  

