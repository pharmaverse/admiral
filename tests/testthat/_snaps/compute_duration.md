# compute_duration Test 18: error if time zones differ

    Code
      compute_duration(ymd_hms("2024-01-01T10:00:00", tz = "Asia/Tokyo"), ymd(
        "2024-01-03"))
    Condition
      Error in `compute_duration()`:
      ! `start_date` and `end_date` must have the same time zone but have "Asia/Tokyo" and "UTC".
      i Please use `lubridate::force_tz()` to set the time zone.

# compute_duration Test 20: error if a date-time has no time zone

    Code
      compute_duration(as.POSIXct("2024-01-01 10:00:00", tz = ""), ymd_hms(
        "2024-01-03T10:00:00"))
    Condition
      Error in `compute_duration()`:
      ! `start_date` and `end_date` must have a time zone.
      i Please use `lubridate::force_tz()` to set the time zone.

