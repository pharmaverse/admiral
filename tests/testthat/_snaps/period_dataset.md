# create_period_dataset Test 4: error if no period/phase variable on RHS

    Code
      create_period_dataset(adsl, new_vars = exprs(USUBJ = USUBJID))
    Condition
      Error in `create_period_dataset()`:
      ! The right hand side values of `new_vars` have to be CDISC style subperiod, period, or phase variables.
      I.e., they must contain the xx or w fragment, e.g., APxxSDT, PxxSwSDT, or PHwSDT.

# create_period_dataset Test 5: error if different type of RHSs

    Code
      create_period_dataset(adsl, new_vars = exprs(APERSDT = APxxSDT, ASPRSDT = PxxSwSDT))
    Condition
      Error in `create_period_dataset()`:
      ! More than one type of subperiod, period, or phase variables
      is specified for `new_vars`:
      subperiod: `PxxSwSDT`
      period: `APxxSDT`

# create_period_dataset Test 6: error if RHS variable not in input dataset

    Code
      create_period_dataset(adsl, new_vars = exprs(PHSDT = PHwSDT))
    Condition
      Error in `create_period_dataset()`:
      ! No variables of the form `PHwSDT` were found in the input dataset.

# derive_vars_period Test 10: error if no period/phase variable on LHS

    Code
      derive_vars_period(adsl, dataset_ref = period_ref, new_vars = exprs(USUBJ = USUBJID))
    Condition
      Error in `derive_vars_period()`:
      ! The left hand side values of `new_vars` have to be CDISC style subperiod, period, or phase variables.
      I.e., they must contain the xx or w fragment, e.g., APxxSDT, PxxSwSDT, or PHwSDT.

# derive_vars_period Test 11: error if different type of LHSs

    Code
      derive_vars_period(adsl, dataset_ref = period_ref, new_vars = exprs(APxxSDT = APERSDT,
        PxxSwSDT = ASPRSDT))
    Condition
      Error in `derive_vars_period()`:
      ! More than one type of subperiod, period, or phase variables
      is specified for `new_vars`:
      subperiod: `PxxSwSDT`
      period: `APxxSDT`

