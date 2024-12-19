# transform_range Test 3: warning if outside range

    Code
      transform_range(c(5, 1, 6, 2, NA), source_range = c(1, 5), target_range = c(0,
        100), outside_range = "warning")
    Condition
      Warning:
      `source` contains values outside the range of 1 to 5:
      source[[3]] = 6
    Output
      [1] 100   0  NA  25  NA

# transform_range Test 4: error if outside range

    Code
      transform_range(c(5, 1, 6, 2, 7), source_range = c(1, 5), target_range = c(0,
        100), outside_range = "error")
    Condition
      Error in `transform_range()`:
      ! `source` contains values outside the range of 1 to 5:
      source[[3]] = 6
      source[[5]] = 7

