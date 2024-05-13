# assert_parameters_argument Test 10: error if argument is of wrong type

    Code
      assert_parameters_argument(myparameters <- c(1, 2, 3))
    Condition
      Error in `assert_parameters_argument()`:
      ! `myparameters <- c(1, 2, 3)` must be a character vector or a list of expressions but it is a double vector.

# get_hori_data Test 11: error if variables with more than one dot

    Code
      get_hori_data(input, parameters = exprs(SYSBP, DIABP), by_vars = exprs(USUBJID,
        VISIT), set_values_to = exprs(AVAL = (AVAL.SYSBP + 2 * AVAL.DIA.BP) / 3),
      filter = NULL)
    Condition
      Error in `get_hori_data()`:
      ! The `set_values_to` argument contains variable names with more than one dot:
      `AVAL.DIA.BP`

