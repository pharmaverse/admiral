# derive_param_bmi Test 36: BMI parameter NOT added

    Code
      result <- derive_param_bmi(input, by_vars = exprs(USUBJID, VISIT),
      get_unit_expr = VSSTRESU)
    Message
      No computed records were added because for all potential computed records at least one of the contributing values was NA.
      If this is not expected, please check the input data.

# derive_param_bsa Test 43: BSA parameter NOT added

    Code
      result <- derive_param_bsa(input, by_vars = exprs(USUBJID, VISIT), method = "Mosteller",
      get_unit_expr = VSSTRESU)
    Message
      No computed records were added because for all potential computed records at least one of the contributing values was NA.
      If this is not expected, please check the input data.

# derive_param_map Test 56: MAP parameter NOT added

    Code
      result <- derive_param_map(input, by_vars = exprs(USUBJID, VISIT), hr_code = "PULSE",
      get_unit_expr = extract_unit(PARAM))
    Message
      No computed records were added because for all potential computed records at least one of the contributing values was NA.
      If this is not expected, please check the input data.

