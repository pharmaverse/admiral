# derive_var_extreme_dt Test 1: deprecation message if function is called

    Code
      ae_start <- date_source(dataset_name = "ae", date = AESTDTM)
    Message
      `date_source()` was deprecated in admiral 1.2.0.
      i Please use `event()` instead.
      x This message will turn into a warning at the beginning of 2026.
      i See admiral's deprecation guidance: https://pharmaverse.github.io/admiraldev/dev/articles/programming_strategy.html#deprecation
    Code
      ae_end <- date_source(dataset_name = "ae", date = AEENDTM)
    Message
      `date_source()` was deprecated in admiral 1.2.0.
      i Please use `event()` instead.
      x This message will turn into a warning at the beginning of 2026.
      i See admiral's deprecation guidance: https://pharmaverse.github.io/admiraldev/dev/articles/programming_strategy.html#deprecation
    Code
      adsl_trtdate <- date_source(dataset_name = "adsl", date = TRTEDTM)
    Message
      `date_source()` was deprecated in admiral 1.2.0.
      i Please use `event()` instead.
      x This message will turn into a warning at the beginning of 2026.
      i See admiral's deprecation guidance: https://pharmaverse.github.io/admiraldev/dev/articles/programming_strategy.html#deprecation
    Code
      adsl_dthdate <- date_source(dataset_name = "adsl", date = DTHDT, filter = nchar(
        DTHDTC) >= 10)
    Message
      `date_source()` was deprecated in admiral 1.2.0.
      i Please use `event()` instead.
      x This message will turn into a warning at the beginning of 2026.
      i See admiral's deprecation guidance: https://pharmaverse.github.io/admiraldev/dev/articles/programming_strategy.html#deprecation
    Code
      derive_var_extreme_dt(adsl, new_var = LSTALVDT, source_datasets = list(ae = ae,
        adsl = adsl), ae_start, ae_end, adsl_trtdate, adsl_dthdate, mode = "last")
    Message
      `derive_var_extreme_dt()` was deprecated in admiral 1.2.0.
      i Please use `derive_vars_extreme_event()` instead.
      x This message will turn into a warning at the beginning of 2026.
      i See admiral's deprecation guidance: https://pharmaverse.github.io/admiraldev/dev/articles/programming_strategy.html#deprecation
      `derive_var_extreme_dtm()` was deprecated in admiral 1.2.0.
      i Please use `derive_vars_extreme_event()` instead.
      x This message will turn into a warning at the beginning of 2026.
      i See admiral's deprecation guidance: https://pharmaverse.github.io/admiraldev/dev/articles/programming_strategy.html#deprecation
    Output
      # A tibble: 3 x 6
        STUDYID USUBJID TRTEDTM             DTHDTC  DTHDT      LSTALVDT  
        <chr>   <chr>   <dttm>              <chr>   <date>     <date>    
      1 STUDY01 1       2020-01-01 12:00:00 <NA>    NA         2020-02-01
      2 STUDY01 2       NA                  2020-06 2020-06-01 NA        
      3 STUDY01 3       2020-04-12 13:15:00 <NA>    NA         2020-04-12

# derive_var_extreme_dtm Test 5: Message sent to users

    Code
      ae_start <- date_source(dataset_name = "ae", date = convert_dtc_to_dtm(AESTDTC),
      set_values_to = exprs(LALVDOM = "AE", LALVSEQ = AESEQ, LALVVAR = "AESTDTC"))
    Message
      `date_source()` was deprecated in admiral 1.2.0.
      i Please use `event()` instead.
      x This message will turn into a warning at the beginning of 2026.
      i See admiral's deprecation guidance: https://pharmaverse.github.io/admiraldev/dev/articles/programming_strategy.html#deprecation
    Code
      ae_end <- date_source(dataset_name = "ae", date = AEENDTM, set_values_to = exprs(
        LALVDOM = "AE", LALVSEQ = AESEQ, LALVVAR = "AEENDTC"))
    Message
      `date_source()` was deprecated in admiral 1.2.0.
      i Please use `event()` instead.
      x This message will turn into a warning at the beginning of 2026.
      i See admiral's deprecation guidance: https://pharmaverse.github.io/admiraldev/dev/articles/programming_strategy.html#deprecation
    Code
      adsl_trtdate <- date_source(dataset_name = "adsl", date = TRTEDTM,
        set_values_to = exprs(LALVDOM = "ADSL", LALVSEQ = NA_integer_, LALVVAR = "TRTEDTM"))
    Message
      `date_source()` was deprecated in admiral 1.2.0.
      i Please use `event()` instead.
      x This message will turn into a warning at the beginning of 2026.
      i See admiral's deprecation guidance: https://pharmaverse.github.io/admiraldev/dev/articles/programming_strategy.html#deprecation
    Code
      adsl_dthdate <- date_source(dataset_name = "adsl", date = DTHDT, filter = nchar(
        DTHDTC) >= 10, set_values_to = exprs(LALVDOM = "ADSL", LALVSEQ = NA_integer_,
        LALVVAR = "DTHDTC"))
    Message
      `date_source()` was deprecated in admiral 1.2.0.
      i Please use `event()` instead.
      x This message will turn into a warning at the beginning of 2026.
      i See admiral's deprecation guidance: https://pharmaverse.github.io/admiraldev/dev/articles/programming_strategy.html#deprecation
    Code
      derive_var_extreme_dtm(adsl, new_var = LSTALVDTM, source_datasets = list(ae = ae,
        adsl = adsl), ae_start, ae_end, adsl_trtdate, adsl_dthdate, mode = "last")
    Message
      `derive_var_extreme_dtm()` was deprecated in admiral 1.2.0.
      i Please use `derive_vars_extreme_event()` instead.
      x This message will turn into a warning at the beginning of 2026.
      i See admiral's deprecation guidance: https://pharmaverse.github.io/admiraldev/dev/articles/programming_strategy.html#deprecation
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
      ae_start <- date_source(dataset_name = "ae", date = AESTDT, set_values_to = exprs(
        LALVDOM = "AE", LALVSEQ = AESEQ, LALVVAR = "AESTDTC"))
      derive_var_extreme_dtm(adsl, new_var = LSTALVDTM, source_datasets = list(ea = ae),
      ae_start, mode = "last")
    Condition
      Error in `derive_var_extreme_dtm()`:
      ! The dataset names must be included in the list specified for the `source_datasets` argument.
      i Following names were provided by `source_datasets`: ea
      i But, `sources[[1]]$dataset_name = ae`

