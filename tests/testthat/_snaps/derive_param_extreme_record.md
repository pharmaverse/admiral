# derive_param_extreme_record Test 1: Message Sent to users

    Code
      derive_param_extreme_record(dataset = aevent, sources = list(records_source(
        dataset_name = "cm", filter = CMDECOD == "ACT", new_vars = exprs(ADT = convert_dtc_to_dt(
          CMSTDTC), AVALC = CMDECOD)), records_source(dataset_name = "pr", filter = PRDECOD ==
        "ACS", new_vars = exprs(ADT = convert_dtc_to_dt(PRSTDTC), AVALC = PRDECOD))),
      source_datasets = list(cm = cm, pr = pr), by_vars = exprs(STUDYID, USUBJID),
      order = exprs(ADT), mode = "first", set_values_to = exprs(PARAMCD = "FIRSTACT",
        PARAM = "First Anti-Cancer Therapy"))
    Message
      `derive_param_extreme_record()` was deprecated in admiral 1.2.0.
      i Please use `derive_extreme_event()` instead.
      x This message will turn into a warning with release of 1.3.0
      i https://pharmaverse.github.io/admiral/reference/derive_extreme_event.html
    Output
      # A tibble: 6 x 7
        STUDYID USUBJID LBSTDTC    PARAMCD  PARAM                     ADT        AVALC
        <chr>   <chr>   <chr>      <chr>    <chr>                     <date>     <chr>
      1 1001    1       2023-01-01 TST      TEST                      NA         <NA> 
      2 1001    2       2023-01-01 TST      TEST                      NA         <NA> 
      3 1001    3       2023-01-01 TST      TEST                      NA         <NA> 
      4 1001    1       <NA>       FIRSTACT First Anti-Cancer Therapy 2020-12-25 ACT  
      5 1001    2       <NA>       FIRSTACT First Anti-Cancer Therapy 2021-12-25 ACS  
      6 1001    3       <NA>       FIRSTACT First Anti-Cancer Therapy 2022-12-25 ACS  

