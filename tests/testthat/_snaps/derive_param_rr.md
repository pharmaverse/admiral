# derive_param_rr Test 2: Message if no new records

    Code
      actual <- derive_param_rr(input, by_vars = exprs(USUBJID, VISIT),
      set_values_to = exprs(PARAMCD = "RRR", PARAM = "RR Duration Rederived (ms)",
        AVALU = "ms"), get_unit_expr = AVALU)
    Message
      No computed records were added because for all potential computed records at least one of the contributing values was NA.
      If this is not expected, please check the input data.

