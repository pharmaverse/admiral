# set_admiral_keys Test 10: an error is issued for invalid arguments

    Code
      set_admiral_keys(input, keys = c("USUBJID", "PARAMCD"))
    Condition
      Error in `set_admiral_keys()`:
      ! Argument `keys` must be a list of variables created by `exprs()`, but is a character vector.
      i Use `chr2vars()` to convert a character vector, e.g. `chr2vars(metacore::get_keys(spec, "ADVS")$variable)`.

