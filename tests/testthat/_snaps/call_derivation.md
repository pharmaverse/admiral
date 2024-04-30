# call_derivation Test 3: Test that Error is thrown if ... has no arguments

    Code
      call_derivation(dataset = input, derivation = derive_vars_dt, variable_params = list(
        params(dtc = AESTDTC, date_imputation = "first", new_vars_prefix = "AST"),
        params(dtc = AEENDTC, date_imputation = "last", new_vars_prefix = "AEN")))
    Condition
      Error in `call_derivation()`:
      ! At least one argument must be set inside `...`.

# call_derivation Test 4: Error is thrown if ... arguments are not properly named

    Code
      call_derivation(dataset = input, derivation = derive_vars_dt, variable_params = list(
        params(dtc = AESTDTC, date_imputation = "first", new_vars_prefix = "AST"),
        params(dtc = AEENDTC, date_imputation = "last", new_vars_prefix = "AEN")),
      XYZSDT, XYZEDT)
    Condition
      Error in `call_derivation()`:
      ! All arguments inside `...` must be named.

# call_derivation Test 5: Error is thrown if params is empty

    Code
      call_derivation(dataset = input, derivation = derive_vars_dt, variable_params = list(
        params(), params()), min_dates = exprs(TRTSDT), max_dates = exprs(TRTEDT))
    Condition
      Error in `params()`:
      ! At least one argument must be provided.

# call_derivation Test 6: Error is thrown if passed params are not properly named

    Code
      call_derivation(dataset = input, derivation = derive_vars_dt, variable_params = list(
        params(XYZ), params(XYZ)), min_dates = exprs(TRTSDT), max_dates = exprs(
        TRTEDT))
    Condition
      Error in `params()`:
      ! All arguments passed to `params()` must be named.

# call_derivation Test 7: Error is thrown if `...` arguments are not properly named

    Code
      call_derivation(dataset = input, derivation = derive_vars_dt, variable_params = list(
        params(dtc = AESTDTC, date_imputation = "first", new_vars_prefix = "AST"),
        params(dtc = AEENDTC, date_imputation = "last", new_vars_prefix = "AEN")),
      XYZSDT, XYZEDT)
    Condition
      Error in `call_derivation()`:
      ! All arguments inside `...` must be named.

# call_derivation Test 8: Error is thrown if duplicate parameters

    Code
      params(dtc = VSDTC, dtc = VSDTC, new_vars_prefix = "A")
    Condition
      Error in `params()`:
      ! The following argument has been specified more than once: "dtc".

