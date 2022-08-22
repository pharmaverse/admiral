test_that("default duration, i.e., relative day", {
  expect_equal(
    compute_duration(
      ymd_hms("2020-12-06T15:00:00"),
      ymd_hms("2020-12-24T08:15:00")
    ),
    19
  )
})

test_that("fractional duration", {
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

test_that("age in years", {
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

test_that("age in months", {
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

test_that("age in weeks", {
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

test_that("duration in hours", {
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

test_that("duration in days after imputation", {
  actual <- compute_duration(
    convert_dtc_to_dt("2020-12", date_imputation = "first"),
    ymd("2020-12-10"),
    out_unit = "days"
  )
  expected <- 10
  expect_equal(actual, expected)
})
