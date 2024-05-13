# derive_var_age_years Test 6: Error is issued if age_unit is missing

    Code
      derive_var_age_years(input, AGE, new_var = AAGE)
    Condition
      Error in `derive_var_age_years()`:
      ! There is no unit variable (`AGEU`) associated with `AGE`. Please specify a value for `age_unit`.

# derive_var_age_years Test 7: warn if `age_unit` doesn't match units in data

    Code
      derive_var_age_years(input, AGE, age_unit = "months", new_var = AAGE)
    Condition
      Warning:
      The unit variable `AGEU` is associated with `AGE` and contains multiple values but the argument `age_unit` has been specified with a single different value. The `age_unit` argument is ignored and the conversion will be based on `AGEU`.
    Output
      # A tibble: 5 x 3
          AGE AGEU    AAGE
        <dbl> <chr>  <dbl>
      1    25 years   25  
      2   312 months  26  
      3    51 years   51  
      4   402 months  33.5
      5   432 months  36  

# derive_var_age_years Test 10: warn if unit in data differs from `age_unit`

    Code
      derive_var_age_years(input, AGE, age_unit = "years", new_var = AAGE)
    Condition
      Warning:
      The unit variable `AGEU` is associated with `AGE` but the argument `age_unit` has been specified with a different value. The `age_unit` argument is ignored and the conversion will be based on `AGEU`.
    Output
      # A tibble: 5 x 3
          AGE AGEU    AAGE
        <dbl> <chr>  <dbl>
      1   459 months  38.2
      2   312 months  26  
      3   510 months  42.5
      4   402 months  33.5
      5   432 months  36  

