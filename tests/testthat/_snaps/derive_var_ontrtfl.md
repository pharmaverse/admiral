# derive_var_ontrtfl Test 15: if trt end date is missing, the obs may still be flagged

    Code
      derive_var_ontrtfl(adcm, start_date = ASTDT, end_date = AENDT, ref_start_date = TRTSDT,
        ref_end_date = TRTEDT, span_period = TRUE)
    Output
        USUBJID      ASTDT     TRTSDT TRTEDT      AENDT ONTRTFL
      1     P01 2018-03-15 2019-01-01     NA 2022-12-01       Y
      2     P02 2020-04-30 2019-01-01     NA 2022-03-15       Y
      3     P03 2020-04-30 2019-01-01     NA       <NA>       Y
      4     P04 2020-04-30       <NA>     NA       <NA>    <NA>

---

    Code
      derive_var_ontrtfl(adcm, start_date = ASTDT, end_date = AENDT, ref_start_date = TRTSDT,
        ref_end_date = TRTEDT)
    Output
        USUBJID      ASTDT     TRTSDT TRTEDT      AENDT ONTRTFL
      1     P01 2018-03-15 2019-01-01     NA 2022-12-01    <NA>
      2     P02 2020-04-30 2019-01-01     NA 2022-03-15       Y
      3     P03 2020-04-30 2019-01-01     NA       <NA>       Y
      4     P04 2020-04-30       <NA>     NA       <NA>    <NA>

