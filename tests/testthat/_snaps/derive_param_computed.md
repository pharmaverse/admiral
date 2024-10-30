# derive_param_computed Test 11: error if keep_nas is invalid

    Code
      derive_param_computed(advs, by_vars = exprs(USUBJID, AVISIT), parameters = c(
        "SYSBP", "DIABP"), set_values_to = exprs(AVAL = (AVAL.SYSBP + 2 *
      AVAL.DIABP) / 3, PARAMCD = "MAP", PARAM = "Mean Arterial Pressure (mmHg)"),
      keep_nas = 3)
    Condition
      Error in `derive_param_computed()`:
      ! Argument `keep_nas` must be either TRUE, FALSE, or a list of <symbol>, but is a number.

# derive_param_computed Test 12: inform if no new records due to NAs

    Code
      derive_param_computed(advs, by_vars = exprs(USUBJID, AVISIT), parameters = c(
        "SYSBP", "DIABP"), set_values_to = exprs(AVAL = (AVAL.SYSBP + 2 *
      AVAL.DIABP) / 3, PARAMCD = "MAP", PARAM = "Mean Arterial Pressure (mmHg)", ADT = ADT.SYSBP,
      ADTF = ADTF.SYSBP))
    Message
      No computed records were added because for all potential computed records at least one of the contributing values was NA.
      If this is not expected, please check the input data and the value of the `keep_nas` argument.
    Output
      # A tibble: 8 x 7
        USUBJID     PARAMCD PARAM                         AVAL AVISIT ADT        ADTF 
        <chr>       <chr>   <chr>                        <dbl> <chr>  <date>     <chr>
      1 01-701-1015 DIABP   Diastolic Blood Pressure (m~    51 BASEL~ 2024-01-10 <NA> 
      2 01-701-1015 DIABP   Diastolic Blood Pressure (m~    50 WEEK 2 2024-01-24 <NA> 
      3 01-701-1015 SYSBP   Systolic Blood Pressure (mm~   121 BASEL~ 2024-01-10 <NA> 
      4 01-701-1015 SYSBP   Systolic Blood Pressure (mm~   121 WEEK 2 2024-01-24 <NA> 
      5 01-701-1028 DIABP   Diastolic Blood Pressure (m~    79 BASEL~ 2024-01-10 <NA> 
      6 01-701-1028 DIABP   Diastolic Blood Pressure (m~    80 WEEK 2 2024-01-24 <NA> 
      7 01-701-1028 SYSBP   Systolic Blood Pressure (mm~   130 BASEL~ 2024-01-10 <NA> 
      8 01-701-1028 SYSBP   Systolic Blood Pressure (mm~    NA WEEK 2 2024-01-24 <NA> 

# assert_parameters_argument Test 13: error if argument is of wrong type

    Code
      assert_parameters_argument(myparameters <- c(1, 2, 3))
    Condition
      Error in `assert_parameters_argument()`:
      ! `myparameters <- c(1, 2, 3)` must be a character vector or a list of expressions but it is a double vector.

# get_hori_data Test 14: error if variables with more than one dot

    Code
      get_hori_data(input, parameters = exprs(SYSBP, DIABP), by_vars = exprs(USUBJID,
        VISIT), set_values_to = exprs(AVAL = (AVAL.SYSBP + 2 * AVAL.DIA.BP) / 3),
      filter = NULL)
    Condition
      Error in `get_hori_data()`:
      ! The `set_values_to` argument contains variable names with more than one dot:
      `AVAL.DIA.BP`

