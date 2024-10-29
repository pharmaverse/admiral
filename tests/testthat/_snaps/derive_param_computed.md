# derive_param_computed Test 11: error if keep_nas is invalid

    Code
      derive_param_computed(advs, by_vars = exprs(USUBJID, AVISIT), parameters = c(
        "SYSBP", "DIABP"), set_values_to = exprs(AVAL = (AVAL.SYSBP + 2 *
      AVAL.DIABP) / 3, PARAMCD = "MAP", PARAM = "Mean Arterial Pressure (mmHg)"),
      keep_nas = 3)
    Condition
      Error in `derive_param_computed()`:
      ! Argument `keep_nas` must be either TRUE, FALSE, or a list of <symbol>, but is a number.

# assert_parameters_argument Test 12: error if argument is of wrong type

    Code
      assert_parameters_argument(myparameters <- c(1, 2, 3))
    Condition
      Error in `assert_parameters_argument()`:
      ! `myparameters <- c(1, 2, 3)` must be a character vector or a list of expressions but it is a double vector.

# get_hori_data Test 13: error if variables with more than one dot

    Code
      get_hori_data(input, parameters = exprs(SYSBP, DIABP), by_vars = exprs(USUBJID,
        VISIT), set_values_to = exprs(AVAL = (AVAL.SYSBP + 2 * AVAL.DIA.BP) / 3),
      filter = NULL)
    Condition
      Error in `get_hori_data()`:
      ! The `set_values_to` argument contains variable names with more than one dot:
      `AVAL.DIA.BP`

