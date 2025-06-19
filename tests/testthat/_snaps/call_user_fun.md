# Test 1: Error is issued to function call if function errors

    Code
      call_user_fun(compute_bmi(height = 172, weight = "hallo"))
    Message
      `call_user_fun()` was deprecated in admiral 1.3.0.
      i `call_user_fun()` is no longer supported and no replacement is provided;
      i The original code for this function is here:
      i https://github.com/pharmaverse/admiral/blob/v1.2.0/R/call_user_fun.R#L26-L39
    Condition
      Error in `call_user_fun()`:
      ! Calling `compute_bmi(height = 172, weight = "hallo")` caused the following error:
      Argument `weight` must be a numeric vector, but it is a string.

