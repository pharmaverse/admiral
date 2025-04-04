# slice_derivation Test 7: Error if a mandatory argument is not in arg or all slices

    Code
      actual <- slice_derivation(advs, derivation = derive_vars_dtm, args = params(
        dtc = VSDTC), derivation_slice(filter = str_detect(VSTPT, "PRE|BEFORE"),
      args = params(time_imputation = "first", new_vars_prefix = "A")),
      derivation_slice(filter = TRUE, args = params(time_imputation = "last")))
    Condition
      Error in `slice_derivation()`:
      ! Issue with the mandatory argument `new_vars_prefix` of derivation function `derive_vars_dtm`. It must: (1) be passed to the `args` argument, or (2) be passed to all derivation slices.

