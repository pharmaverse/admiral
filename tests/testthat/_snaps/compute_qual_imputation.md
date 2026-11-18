# compute_qual_imputation Test 30: validate argument types and method

    Code
      compute_qual_imputation("<40", imputation_type = 0)
    Condition <compute_qual_imputation_method>
      Error in `compute_qual_imputation()`:
      ! Argument `imputation_type` must be "strip" or "shift".

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

# compute_qual_imputation Test 32: preserve deprecated calls and messages

    Code
      expect_equal(compute_qual_imputation(values, imputation_type = 1), c(40, 40.1,
        40, 40.1, 3))
    Message
      The `imputation_type` argument of `compute_qual_imputation()` must be a string as of admiral 1.6.0.
      i Use `imputation_type = "strip"` instead of 1 and `imputation_type = "shift"` instead of 2.
    Code
      expect_equal(compute_qual_imputation(values, imputation_type = 2), c(39, 40.2,
        40, 40.1, 3))
    Message
      The `imputation_type` argument of `compute_qual_imputation()` must be a string as of admiral 1.6.0.
      i Use `imputation_type = "strip"` instead of 1 and `imputation_type = "shift"` instead of 2.
    Code
      expect_equal(compute_qual_imputation(values, 1, 0.5), c(39.5, 40.6, 40, 40.1, 3))
    Message
      The `imputation_type` argument of `compute_qual_imputation()` must be a string as of admiral 1.6.0.
      i Use `imputation_type = "strip"` instead of 1 and `imputation_type = "shift"` instead of 2.
      The `factor` argument of `compute_qual_imputation()` is deprecated as of admiral 1.6.0.
      i Please use the `offset` argument instead.
    Code
      expect_equal(compute_qual_imputation(values, 2, 100), c(39, 40.2, 40, 40.1, 3))
    Message
      The `imputation_type` argument of `compute_qual_imputation()` must be a string as of admiral 1.6.0.
      i Use `imputation_type = "strip"` instead of 1 and `imputation_type = "shift"` instead of 2.
      The `factor` argument of `compute_qual_imputation()` is deprecated as of admiral 1.6.0.
      i Please use the `offset` argument instead.
    Code
      expect_equal(compute_qual_imputation(values, factor = 0.5), c(39.5, 40.6, 40,
        40.1, 3))
    Message
      The `factor` argument of `compute_qual_imputation()` is deprecated as of admiral 1.6.0.
      i Please use the `offset` argument instead.
    Code
      expect_equal(compute_qual_imputation(values, 2L), c(39, 40.2, 40, 40.1, 3))
    Message
      The `imputation_type` argument of `compute_qual_imputation()` must be a string as of admiral 1.6.0.
      i Use `imputation_type = "strip"` instead of 1 and `imputation_type = "shift"` instead of 2.

