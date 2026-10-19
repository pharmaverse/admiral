# compute_duration Test 18: warning if time zones differ

    Code
      compute_duration(ymd_hms("2024-01-01T10:00:00", tz = "Asia/Tokyo"), ymd(
        "2024-01-03"))
    Condition
      Warning:
      `start_date` and `end_date` have different time zones: "Asia/Tokyo" and "UTC".
      Date-times are expected in UTC. Please use `lubridate::force_tz()` to set the time zone to UTC.
    Output
      [1] 3.375

