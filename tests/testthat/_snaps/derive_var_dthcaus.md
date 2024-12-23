# derive_var_dthcaus Test 1: Message sent to users

    Code
      src_ae <- dthcaus_source(dataset_name = "ae", filter = AEOUT == "FATAL", date = AEDTHDT,
      mode = "first", dthcaus = AEDECOD)
    Message
      `date_source()` was deprecated in admiral 1.2.0.
      i Please use `derive_vars_extreme_event()` instead.
      x This message will turn into a warning with release of 1.3.0
      i https://pharmaverse.github.io/admiral/reference/derive_vars_extreme_event.html
    Code
      src_ds <- dthcaus_source(dataset_name = "ds", filter = DSDECOD == "DEATH" &
        grepl("DEATH DUE TO", DSTERM), date = convert_dtc_to_dt(DSSTDTC), mode = "first",
      dthcaus = str_to_upper(DSTERM))
    Message
      `date_source()` was deprecated in admiral 1.2.0.
      i Please use `derive_vars_extreme_event()` instead.
      x This message will turn into a warning with release of 1.3.0
      i https://pharmaverse.github.io/admiral/reference/derive_vars_extreme_event.html
    Code
      derive_var_dthcaus(adsl, source_datasets = list(ae = ae, ds = ds), src_ae,
      src_ds)
    Message
      `derive_var_dthcaus()` was deprecated in admiral 1.2.0.
      i Please use `derive_vars_extreme_event()` instead.
      x This message will turn into a warning with release of 1.3.0
      i https://pharmaverse.github.io/admiral/reference/derive_vars_extreme_event.html
    Output
      # A tibble: 3 x 3
        STUDYID USUBJID DTHCAUS                            
        <chr>   <chr>   <chr>                              
      1 TEST01  PAT01   DEATH DUE TO PROGRESSION OF DISEASE
      2 TEST01  PAT02   <NA>                               
      3 TEST01  PAT03   SUDDEN DEATH                       

# derive_var_dthcaus Test 12: error if source dataset is not available

    Code
      derive_var_dthcaus(adsl, source_datasets = list(ae = ae, dd = ds), src_ae,
      src_ds)
    Condition
      Error in `derive_var_dthcaus()`:
      ! The dataset names must be included in the list specified for the `source_datasets` argument.
      i Following names were provided by `source_datasets`: ae and dd
      i But, `sources[[2]]$dataset_name = ds`

