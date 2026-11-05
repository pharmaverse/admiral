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

# compute_qual_imputation Test 30: validate argument types and method

    Code
      compute_qual_imputation("<40", imputation_type = 0)
    Condition <compute_qual_imputation_method>
      Error in `compute_qual_imputation()`:
      ! Numeric `imputation_type` must be 1 or 2.

# compute_qual_imputation Test 31: validate offsets

    Code
      compute_qual_imputation("<40", offset = -1)
    Condition <compute_qual_imputation_offset>
      Error in `compute_qual_imputation()`:
      ! Argument `offset` must be finite and nonnegative.

---

    Code
      compute_qual_imputation("<40", factor = 1, offset = 1)
    Condition <compute_qual_imputation_arguments>
      Error in `compute_qual_imputation()`:
      ! Supply only `offset`, not both `factor` and `offset`.

