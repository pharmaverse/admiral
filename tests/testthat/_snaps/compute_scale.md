# compute_scale Test 6: error if source_range is supplied, but not target_range

    Code
      compute_scale(input, source_range = c(1, 5), min_n = 2)
    Condition
      Error in `compute_scale()`:
      ! Argument `target_range` is missing with no default and `source_range` is not missing.
      i Either both or neither arguments should be specified.

# compute_scale Test 7: error if target_range is supplied, but not source_range

    Code
      compute_scale(input, target_range = c(0, 100), min_n = 2)
    Condition
      Error in `compute_scale()`:
      ! Argument `source_range` is missing with no default and `target_range` is not missing.
      i Either both or neither arguments should be specified.

