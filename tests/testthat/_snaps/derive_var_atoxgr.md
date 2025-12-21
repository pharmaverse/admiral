# derive_var_atoxgr_dir Test 130: CTCAEv6  Blood bilirubin increased

    Code
      actual_bili_ctcv6 <- derive_var_atoxgr_dir(input_bili_ctcv6, new_var = ATOXGRH,
        meta_criteria = atoxgr_criteria_ctcv6, tox_description_var = ATOXDSCH,
        criteria_direction = "H", abnormal_indicator = "HIGH", get_unit_expr = AVALU)
    Message
      The `abnormal_indicator` argument of `derive_var_atoxgr_dir()` is deprecated as of admiral 1.4.0.
      x This message will turn into a warning at the beginning of 2027.
      i See admiral's deprecation guidance: https://pharmaverse.github.io/admiraldev/dev/articles/programming_strategy.html#deprecation
      x The argument is mapped to `high_indicator`, i.e., `high_indicator = "HIGH"`

