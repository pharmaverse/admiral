# compute_qual_imputation Test 33: inform about deprecated arguments

    Code
      compute_qual_imputation("<40", imputation_type = 1)
    Message
      The `imputation_type` argument of `compute_qual_imputation()` must be a string as of admiral 1.6.0.
      i Use `imputation_type = "strip"` instead of 1 and `imputation_type = "shift"` instead of 2.
    Output
      [1] 40
    Code
      compute_qual_imputation("<40", imputation_type = 2)
    Message
      The `imputation_type` argument of `compute_qual_imputation()` must be a string as of admiral 1.6.0.
      i Use `imputation_type = "strip"` instead of 1 and `imputation_type = "shift"` instead of 2.
    Output
      [1] 39
    Code
      compute_qual_imputation("<40", factor = 0.5)
    Message
      The `factor` argument of `compute_qual_imputation()` is deprecated as of admiral 1.6.0.
      i Please use the `offset` argument instead.
    Output
      [1] 39.5

