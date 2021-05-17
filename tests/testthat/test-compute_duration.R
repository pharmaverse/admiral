context("test-compute_duration")

test_that("default duration, i.e., relative day", {
  expect_equal(compute_duration(
    ymd_hms("2020-12-06T15:00:00"),
    ymd_hms("2020-12-24T08:15:00")
  ),
  19)
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
