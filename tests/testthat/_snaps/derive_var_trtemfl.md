# derive_var_trtemfl Test 6: error if `end_window` without `trt_end_date`

    Code
      derive_var_trtemfl(adae, end_window = 10)
    Condition
      Error in `derive_var_trtemfl()`:
      ! `end_window` argument was specified but not `trt_end_date`
      Either both or none of them must be specified.

# derive_var_trtemfl Test 7: error if `initial_intensity` without `intensity`

    Code
      derive_var_trtemfl(adae, initial_intensity = AEITOXGR)
    Condition
      Error in `derive_var_trtemfl()`:
      ! `initial_intensity` argument was specified but not `intensity`
      Either both or none of them must be specified.

# derive_var_trtemfl Test 8: error if `intensity` without `initial_intensity`

    Code
      derive_var_trtemfl(adae, intensity = AETOXGR)
    Condition
      Error in `derive_var_trtemfl()`:
      ! `intensity` argument was specified but not `initial_intensity`
      Either both or none of them must be specified.

# derive_var_trtemfl Test 9: error if `intensity` without `initial_intensity`

    Code
      derive_var_trtemfl(adae2, intensity = AETOXGR)
    Condition
      Error in `derive_var_trtemfl()`:
      ! `intensity` argument was specified but not `initial_intensity`
      Either both or none of them must be specified.

