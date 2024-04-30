# derive_vars_dy Test 8: error if source variables do not end in DT or DTM

    Code
      derive_vars_dy(datain, reference_date = TRTSDTW, source_vars = exprs(TRTSDTW,
        ASTDTW, AENDTW))
    Condition
      Error in `derive_vars_dy()`:
      ! `source_vars` must end in DT or DTM or be explicitly and uniquely named.
      i Please name or rename the following source_vars: `TRTSDTW`, `ASTDTW`, and `AENDTW`

