# compute_duration ----

## Test 1: Default duration, i.e., relative day ----
test_that("compute_duration Test 1: Default duration, i.e., relative day", {
  expect_equal(
    compute_duration(
      ymd_hms("2020-12-06T15:00:00"),
      ymd_hms("2020-12-24T08:15:00")
    ),
    19
  )
})

## Test 2: Fractional duration ----
test_that("compute_duration Test 2: Fractional duration", {
  expect_equal(
    compute_duration(
      ymd_hms("2020-12-06T15:00:00"),
      ymd_hms("2020-12-24T09:00:00"),
      floor_in = FALSE,
      add_one = FALSE
    ),
    17.75
  )
})

## Test 3: Age in years ----
test_that("compute_duration Test 3: Age in years", {
  expect_equal(
    compute_duration(
      ymd("1984-09-06"),
      ymd("2020-02-24"),
      trunc_out = TRUE,
      out_unit = "years",
      add_one = FALSE
    ),
    35
  )
})

## Test 4: Age in months ----
test_that("compute_duration Test 4: Age in months", {
  expect_equal(
    compute_duration(
      ymd("1984-09-06"),
      ymd("2020-02-24"),
      trunc_out = TRUE,
      out_unit = "months",
      add_one = FALSE
    ),
    425
  )
})

## Test 5: Age in weeks ----
test_that("compute_duration Test 5: Age in weeks", {
  expect_equal(
    compute_duration(
      ymd("2020-02-03"),
      ymd("2020-02-24"),
      trunc_out = TRUE,
      out_unit = "weeks",
      add_one = FALSE
    ),
    3
  )
})

## Test 6: Duration in hours ----
test_that("compute_duration Test 6: Duration in hours", {
  expect_equal(
    compute_duration(
      ymd_hms("2020-12-06T9:00:00"),
      ymd_hms("2020-12-06T13:30:00"),
      out_unit = "hours",
      floor_in = FALSE,
      add_one = FALSE
    ),
    4.5
  )
})


## Test 7: Duration in minutes (minutes option) ----
test_that("compute_duration Test 7: Duration in minutes (minutes option)", {
  expect_equal(
    compute_duration(
      ymd_hms("2020-12-06T13:00:00"),
      ymd_hms("2020-12-06T13:30:00"),
      out_unit = "minutes",
      floor_in = FALSE,
      add_one = FALSE
    ),
    30
  )
})

## Test 8: Duration in minutes (min option) ----
test_that("compute_duration Test 8: Duration in minutes (min option)", {
  expect_equal(
    compute_duration(
      ymd_hms("2020-12-06T13:00:00"),
      ymd_hms("2020-12-06T13:30:00"),
      out_unit = "min",
      floor_in = FALSE,
      add_one = FALSE
    ),
    30
  )
})


## Test 9: Duration in seconds (seconds option) ----
test_that("compute_duration Test 9: Duration in seconds (seconds option)", {
  expect_equal(
    compute_duration(
      ymd_hms("2020-12-06T13:30:00"),
      ymd_hms("2020-12-06T13:30:29"),
      out_unit = "seconds",
      floor_in = FALSE,
      add_one = FALSE
    ),
    29
  )
})



## Test 10: Duration in seconds (sec option) ----
test_that("compute_duration Test 10: Duration in seconds (sec option)", {
  expect_equal(
    compute_duration(
      ymd_hms("2020-12-06T13:30:00"),
      ymd_hms("2020-12-06T13:30:29"),
      out_unit = "sec",
      floor_in = FALSE,
      add_one = FALSE
    ),
    29
  )
})


## Test 11: Duration (instead of interval) ----
test_that("compute_duration Test 11: Duration (instead of interval)", {
  expect_equal(
    compute_duration(
      ymd("2000-02-01"),
      ymd("2000-03-01"),
      out_unit = "months",
      add_one = FALSE,
      type = "duration"
    ),
    29 / (365.25 / 12) # 29 days divided by the average month length
  )

  expect_equal(
    compute_duration(
      ymd("2000-02-01"),
      ymd("2001-02-01"),
      out_unit = "years",
      add_one = FALSE,
      type = "duration"
    ),
    366 / 365.25 # 366 days in this leap year divided by the average year length
  )
})

## Test 12: Interval (instead of duration) ----
test_that("compute_duration Test 12: Interval (instead of duration)", {
  expect_equal(
    compute_duration(
      ymd("2000-02-01"),
      ymd("2000-03-01"),
      out_unit = "months",
      add_one = FALSE,
      type = "interval"
    ),
    1
  )

  expect_equal(
    compute_duration(
      ymd("2000-02-01"),
      ymd("2001-02-01"),
      out_unit = "years",
      add_one = FALSE,
      type = "interval"
    ),
    1
  )
})

## Test 13: Interval with duration/interval invariant units ----
test_that("compute_duration Test 13: Interval with duration/interval invariant units", {
  expect_equal(
    compute_duration(
      ymd("2000-02-01"),
      ymd("2000-03-01"),
      out_unit = "days",
      add_one = FALSE,
      type = "interval"
    ),
    compute_duration(
      ymd("2000-02-01"),
      ymd("2000-03-01"),
      out_unit = "days",
      add_one = FALSE,
      type = "duration"
    )
  )

  expect_equal(
    compute_duration(
      ymd("2000-02-01"),
      ymd("2001-02-01"),
      out_unit = "weeks",
      add_one = FALSE,
      type = "interval"
    ),
    compute_duration(
      ymd("2000-02-01"),
      ymd("2001-02-01"),
      out_unit = "weeks",
      add_one = FALSE,
      type = "duration"
    )
  )
})
