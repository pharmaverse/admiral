# derive_param_qtc Test 4: Message if no new records

    Code
      actual <- derive_param_qtc(input, by_vars = exprs(USUBJID, VISIT), method = "Bazett",
      get_unit_expr = AVALU)
    Message
      No computed records were added because for all potential computed records at least one of the contributing values was NA.
      If this is not expected, please check the input data.

