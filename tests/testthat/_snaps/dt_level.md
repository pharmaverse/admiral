# dt_level Test 5: input is not scalar

    Code
      dt_level(c("D", "M", "Y"))
    Condition
      Error in `dt_level()`:
      ! Argument `level` must be a scalar of class <character>, but is length 3

# dt_level Test 6: input is scalar character but not in expected set

    Code
      dt_level("d")
    Condition
      Error in `dt_level()`:
      ! Argument `level` must be equal to one of "n", "D", "M", or "Y".

