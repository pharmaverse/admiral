# process_set_values_to Test 2: catch error

    Code
      process_set_values_to(bds, set_values_to = exprs(PARAMCD = BMI, PARAM = "Body-Mass-Index",
        PARAMN = 1))
    Condition
      Error in `process_set_values_to()`:
      ! Assigning variables failed!
      * `set_values_to = exprs(PARAMCD = BMI, PARAM = Body-Mass-Index, PARAMN = 1)`
      See error message below:
      i In argument: `PARAMCD = BMI`. Caused by error: ! object 'BMI' not found

# process_set_values_to Test 3: check types

    Code
      process_set_values_to(bds, set_values_to = exprs(PARAMCD = 1, PARAM = "Body-Mass-Index",
        PARAMN = "BMI"), expected_types = c(PARAMCD = "character", PARAM = "character",
        PARAMN = "numeric"))
    Condition
      Error in `process_set_values_to()`:
      ! The following variables have an unexpected type:
      * PARAMCD: expected is <character>, but it is <numeric>.
      * PARAMN: expected is <numeric>, but it is <character>.

