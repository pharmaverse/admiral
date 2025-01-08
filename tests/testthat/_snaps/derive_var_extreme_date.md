# derive_var_extreme_dt Test 1: Message sent to users

    Code
      derive_var_extreme_dt(adsl, new_var = LSTALVDT, source_datasets = list(ae = ae,
        adsl = adsl), ae_start, ae_end, adsl_trtdate, adsl_dthdate, mode = "last")
    Message
      `derive_var_extreme_dt()` was deprecated in admiral 1.2.0.
      i Please use `derive_vars_extreme_event()` instead.
      x This message will turn into a warning in the next release.
      i https://pharmaverse.github.io/admiral/reference/derive_vars_extreme_event.html
      `derive_var_extreme_dtm()` was deprecated in admiral 1.2.0.
      i Please use `derive_vars_extreme_event()` instead.
      x This message will turn into a warning in the next release.
      i https://pharmaverse.github.io/admiral/reference/derive_vars_extreme_event.html
    Output
      # A tibble: 3 x 6
        STUDYID USUBJID TRTEDTM             DTHDTC  DTHDT      LSTALVDT  
        <chr>   <chr>   <dttm>              <chr>   <date>     <date>    
      1 STUDY01 1       2020-01-01 12:00:00 <NA>    NA         2020-02-01
      2 STUDY01 2       NA                  2020-06 2020-06-01 NA        
      3 STUDY01 3       2020-04-12 13:15:00 <NA>    NA         2020-04-12

# derive_var_extreme_dtm Test 5: Message sent to users

    Code
      derive_var_extreme_dtm(adsl, new_var = LSTALVDTM, source_datasets = list(ae = ae,
        adsl = adsl), ae_start, ae_end, adsl_trtdate, adsl_dthdate, mode = "last")
    Message
      `derive_var_extreme_dtm()` was deprecated in admiral 1.2.0.
      i Please use `derive_vars_extreme_event()` instead.
      x This message will turn into a warning in the next release.
      i https://pharmaverse.github.io/admiral/reference/derive_vars_extreme_event.html
    Output
      # A tibble: 3 x 9
        STUDYID USUBJID TRTEDTM             DTHDTC  DTHDT      LALVDOM LALVSEQ LALVVAR
        <chr>   <chr>   <dttm>              <chr>   <date>     <chr>     <dbl> <chr>  
      1 STUDY01 1       2020-01-01 12:00:00 <NA>    NA         AE            2 AEENDTC
      2 STUDY01 2       NA                  2020-06 2020-06-01 <NA>         NA <NA>   
      3 STUDY01 3       2020-04-12 13:15:00 <NA>    NA         ADSL         NA TRTEDTM
      # i 1 more variable: LSTALVDTM <dttm>

# derive_var_extreme_dtm Test 8: error if source dataset is not available

    Code
      derive_var_extreme_dtm(adsl, new_var = LSTALVDTM, source_datasets = list(ea = ae),
      ae_start, mode = "last")
    Condition
      Error in `derive_var_extreme_dtm()`:
      ! The dataset names must be included in the list specified for the `source_datasets` argument.
      i Following names were provided by `source_datasets`: ea
      i But, `sources[[1]]$dataset_name = ae`

