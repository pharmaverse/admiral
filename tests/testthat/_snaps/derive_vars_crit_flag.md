# derive_vars_crit_flag Test 4: error if invalid condition (var not in input)

    Code
      derive_vars_crit_flag(input, condition = AVAL > 3 * ANRHI, description = "> 3ULN")
    Condition
      Error in `derive_vars_crit_flag()`:
      ! Evaluating `condition` (`AVAL > 3 * ANRHI`) in `dataset` failed:
        object 'ANRHI' not found

# derive_vars_crit_flag Test 5: error if invalid description (PARAMCD not in input)

    Code
      derive_vars_crit_flag(select(input, AVAL), crit_nr = 2, condition = AVAL > 40,
      description = paste(PARAMCD, "> 40"), values_yn = TRUE, create_numeric_flag = TRUE)
    Condition
      Error in `derive_vars_crit_flag()`:
      ! Evaluating `description` (`paste(PARAMCD, "> 40")`) in `dataset` failed:
        object 'PARAMCD' not found

