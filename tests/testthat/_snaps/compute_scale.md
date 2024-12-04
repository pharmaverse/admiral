# compute_scale Test 6: error if source_range is supplied, but not target_range

    Code
      compute_scale(input, source_range = c(1, 5), min_n = 2)
    Condition
      Error in `compute_scale()`:
      ! could not find function "compute_scale"

---

    Code
      compute_scale(input, target_range = c(0, 100), min_n = 2)
    Condition
      Error in `compute_scale()`:
      ! could not find function "compute_scale"

