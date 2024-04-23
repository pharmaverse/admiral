# call_derivation Test 8: Error is thrown if duplicate parameters

    Code
      params(dtc = VSDTC, dtc = VSDTC, new_vars_prefix = "A")
    Condition
      Error in `params()`:
      ! The following parameters have been specified more than once: "dtc".

