# Test 1: deprecation error

    Code
      call_user_fun(compute_bmi(height = 172, weight = "hallo"))
    Condition
      Error:
      ! `call_user_fun()` was deprecated in admiral 1.3.0 and is now defunct.
      i `call_user_fun()` is no longer supported and no replacement is provided;
      i The original code for this function is here:
      i https://github.com/pharmaverse/admiral/blob/v1.2.0/R/call_user_fun.R#L26-L39

