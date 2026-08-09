# derive_var_atoxgr_dir Test 130: CTCAEv6  Blood bilirubin increased

    Code
      actual_bili_ctcv6 <- derive_var_atoxgr_dir(input_bili_ctcv6, new_var = ATOXGRH,
        meta_criteria = atoxgr_criteria_ctcv6, tox_description_var = ATOXDSCH,
        criteria_direction = "H", abnormal_indicator = "HIGH", get_unit_expr = AVALU)
    Condition
      Warning:
      `derive_var_atoxgr_dir()` was deprecated in admiral 1.4.0.
      i The `abnormal_indicator` argument is deprecated. Please use `low_indicator` or `high_indicator` instead.
      x This message will turn into an error at the beginning of 2028.
      i See admiral's deprecation guidance: https://pharmaverse.github.io/admiraldev/dev/articles/programming_strategy.html#deprecation
      x The argument is mapped to `high_indicator`, i.e., `high_indicator = "HIGH"`

