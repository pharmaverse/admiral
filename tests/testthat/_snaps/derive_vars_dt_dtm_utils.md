# is_partial_datetime correctly identifies datetime and date partials

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

