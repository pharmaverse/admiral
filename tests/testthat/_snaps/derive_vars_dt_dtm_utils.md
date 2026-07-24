# get_dt_dtm_range Test 13: get_dt_dtm_range correctly imputes date ranges

    Code
      get_dt_dtm_range(invalid_dtc, create_datetime = FALSE)
    Condition
      Warning:
      Dataset contains incorrect datetime format: --DTC may be incorrectly imputed on row(s)
      * Row 1 : --DTC = invalid-date
      * ISO representations of the form YYYY-MM-DDThh:mm:ss.ddd are expected, e.g., 2003-12-15T13:15:17.123. Missing parts at the end can be omitted. Missing parts in the middle must be represented by a dash, e.g., 2003---15.
      Warning:
       1 failed to parse.
    Output
      $lower
      [1] "0000-01-01" "2021-13-40"
      
      $upper
      [1] "9999-12-31" "2021-13-40"
      

# is_partial_datetime Test 33: correctly identifies datetime and date partials

    Code
      is_partial_datetime(partial_invalid)
    Condition
      Error in `is_partial_datetime()`:
      ! `partial` must be a named list containing either all date components or all datetime components

---

    Code
      is_partial_datetime(list())
    Condition
      Error in `is_partial_datetime()`:
      ! `partial` must be a named list containing either all date components or all datetime components

---

    Code
      is_partial_datetime(partial_time_only)
    Condition
      Error in `is_partial_datetime()`:
      ! `partial` must be a named list containing either all date components or all datetime components

