# derive_param_bmi Test 8: BMI parameter NOT added

    Code
      result <- derive_param_bmi(input, by_vars = exprs(USUBJID, VISIT),
      get_unit_expr = VSSTRESU)
    Message
      No computed records were added because for all potential computed records at least one of the contributing values was NA.
      If this is not expected, please check the input data.

