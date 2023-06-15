test_that("compute_duration Test 1: Default duration, i.e., relative day", {
  expect_equal(
    compute_duration(
      ymd_hms("2020-12-06T15:00:00"),
      ymd_hms("2020-12-24T08:15:00")
    ),
    19
  )
})

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


test_that("compute_duration Test 7a: Duration in minutes (minutes option)", {
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

test_that("compute_duration Test 7b: Duration in minutes (min option)", {
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


test_that("compute_duration Test 8a: Duration in seconds (seconds option)", {
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



test_that("compute_duration Test 8b: Duration in seconds (sec option)", {
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
