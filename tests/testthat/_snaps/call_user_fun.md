# Test 1: Error is issued to function call if function errors

    Code
      call_user_fun(compute_bmi(height = 172, weight = "hallo"))
    Message
      `call_user_fun()` was deprecated in admiral 1.3.0.
      i 
    Condition
      Error in `call_user_fun()`:
      ! Calling `compute_bmi(height = 172, weight = "hallo")` caused the following error:
      Argument `weight` must be a numeric vector, but it is a string.

