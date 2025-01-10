# derive_param_map Test 11: MAP parameter NOT added

    Code
      result <- derive_param_map(input, by_vars = exprs(USUBJID, VISIT), hr_code = "PULSE",
      get_unit_expr = extract_unit(PARAM))
    Message
      No computed records were added because for all potential computed records at least one of the contributing values was NA.
      If this is not expected, please check the input data.

