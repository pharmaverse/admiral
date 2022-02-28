context("test-create_single_dose_dataset")


test_that("derive_vars_single_dose works as expected", {
  input <- tibble::tribble(
    ~USUBJID, ~EXDOSFRQ, ~ASTDT, ~AENDT,
    "P01", "Q2D", lubridate::ymd("2021-01-01"), lubridate::ymd("2021-01-03"),
    "P01", "Q3D", lubridate::ymd("2021-01-08"), lubridate::ymd("2021-01-12"),
    "P01", "EVERY 2 WEEKS", lubridate::ymd("2021-01-15"), lubridate::ymd("2021-01-29")
  )
  expected_output <- tibble::tribble(
    ~USUBJID, ~EXDOSFRQ, ~ASTDT, ~AENDT,
    "P01", "ONCE", lubridate::ymd("2021-01-01"), lubridate::ymd("2021-01-01"),
    "P01", "ONCE", lubridate::ymd("2021-01-03"), lubridate::ymd("2021-01-03"),
    "P01", "ONCE", lubridate::ymd("2021-01-08"), lubridate::ymd("2021-01-08"),
    "P01", "ONCE", lubridate::ymd("2021-01-11"), lubridate::ymd("2021-01-11"),
    "P01", "ONCE", lubridate::ymd("2021-01-15"), lubridate::ymd("2021-01-15"),
    "P01", "ONCE", lubridate::ymd("2021-01-29"), lubridate::ymd("2021-01-29")
  )

  expect_equal(create_single_dose_dataset(input, by_vars = vars(USUBJID)), expected_output)
})

test_that("derive_vars_single_dose works for different treatments", {
  input <- tibble::tribble(
    ~USUBJID, ~EXDOSFRQ, ~ASTDT, ~AENDT, ~EXTRT,
    "P01", "Q2D", lubridate::ymd("2021-01-01"), lubridate::ymd("2021-01-03"), "XANOMELINE",
    "P01", "EVERY 2 DAYS", lubridate::ymd("2021-01-01"), lubridate::ymd("2021-01-05"), "PLACEBO"
  )
  expected_output <- tibble::tribble(
    ~USUBJID, ~EXDOSFRQ, ~ASTDT, ~AENDT, ~EXTRT,
    "P01", "ONCE", lubridate::ymd("2021-01-01"), lubridate::ymd("2021-01-01"), "XANOMELINE",
    "P01", "ONCE", lubridate::ymd("2021-01-03"), lubridate::ymd("2021-01-03"), "XANOMELINE",
    "P01", "ONCE", lubridate::ymd("2021-01-01"), lubridate::ymd("2021-01-01"), "PLACEBO",
    "P01", "ONCE", lubridate::ymd("2021-01-03"), lubridate::ymd("2021-01-03"), "PLACEBO",
    "P01", "ONCE", lubridate::ymd("2021-01-05"), lubridate::ymd("2021-01-05"), "PLACEBO"
  )

  expect_equal(create_single_dose_dataset(input, by_vars = vars(USUBJID, EXTRT)), expected_output)
})
