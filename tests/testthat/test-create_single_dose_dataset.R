test_that("create_single_dose_dataset works as expected for Q#/EVERY # cases", {
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

  expect_equal(create_single_dose_dataset(input), expected_output)
})


test_that("create_single_dose_dataset works as expected for # TIMES PER cases", {
  input <- tibble::tribble(
    ~USUBJID, ~EXDOSFRQ, ~ASTDT, ~AENDT,
    "P01", "2 TIMES PER YEAR", lubridate::ymd("2021-01-01"), lubridate::ymd("2021-07-01"),
    "P02", "2 TIMES PER YEAR", lubridate::ymd("2021-01-01"), lubridate::ymd("2021-12-31"),
    "P03", "4 TIMES PER MONTH", lubridate::ymd("2021-02-01"), lubridate::ymd("2021-03-01"),
    "P04", "4 TIMES PER MONTH", lubridate::ymd("2021-01-01"), lubridate::ymd("2021-01-20"),
    "P05", "5 TIMES PER WEEK", lubridate::ymd("2021-01-15"), lubridate::ymd("2021-01-17"),
    "P06", "5 TIMES PER WEEK", lubridate::ymd("2021-01-15"), lubridate::ymd("2021-01-21")
  )
  expected_output <- tibble::tribble(
    ~USUBJID, ~EXDOSFRQ, ~ASTDT, ~AENDT,
    "P01", "ONCE", lubridate::ymd("2021-01-01"), lubridate::ymd("2021-01-01"),
    "P02", "ONCE", lubridate::ymd("2021-01-01"), lubridate::ymd("2021-01-01"),
    "P02", "ONCE", lubridate::ymd("2021-07-02"), lubridate::ymd("2021-07-02"),
    "P03", "ONCE", lubridate::ymd("2021-02-01"), lubridate::ymd("2021-02-01"),
    "P03", "ONCE", lubridate::ymd("2021-02-08"), lubridate::ymd("2021-02-08"),
    "P03", "ONCE", lubridate::ymd("2021-02-16"), lubridate::ymd("2021-02-16"),
    "P03", "ONCE", lubridate::ymd("2021-02-23"), lubridate::ymd("2021-02-23"),
    "P04", "ONCE", lubridate::ymd("2021-01-01"), lubridate::ymd("2021-01-01"),
    "P04", "ONCE", lubridate::ymd("2021-01-08"), lubridate::ymd("2021-01-08"),
    "P04", "ONCE", lubridate::ymd("2021-01-16"), lubridate::ymd("2021-01-16"),
    "P05", "ONCE", lubridate::ymd("2021-01-15"), lubridate::ymd("2021-01-15"),
    "P05", "ONCE", lubridate::ymd("2021-01-16"), lubridate::ymd("2021-01-16"),
    "P05", "ONCE", lubridate::ymd("2021-01-17"), lubridate::ymd("2021-01-17"),
    "P06", "ONCE", lubridate::ymd("2021-01-15"), lubridate::ymd("2021-01-15"),
    "P06", "ONCE", lubridate::ymd("2021-01-16"), lubridate::ymd("2021-01-16"),
    "P06", "ONCE", lubridate::ymd("2021-01-17"), lubridate::ymd("2021-01-17"),
    "P06", "ONCE", lubridate::ymd("2021-01-19"), lubridate::ymd("2021-01-19"),
    "P06", "ONCE", lubridate::ymd("2021-01-20"), lubridate::ymd("2021-01-20"),
  )

  expect_equal(create_single_dose_dataset(input), expected_output)
})

test_that("create_single_dose_dataset works for different treatments", {
  input <- tibble::tribble(
    ~USUBJID, ~EXDOSFRQ, ~ASTDT, ~AENDT, ~EXTRT,
    "P01", "Q2D", lubridate::ymd("2021-01-01"), lubridate::ymd("2021-01-03"), "XANOMELINE",
    "P01", "QOD", lubridate::ymd("2021-01-01"), lubridate::ymd("2021-01-05"), "PLACEBO"
  )
  expected_output <- tibble::tribble(
    ~USUBJID, ~EXDOSFRQ, ~ASTDT, ~AENDT, ~EXTRT,
    "P01", "ONCE", lubridate::ymd("2021-01-01"),
    lubridate::ymd("2021-01-01"), "XANOMELINE",
    "P01", "ONCE", lubridate::ymd("2021-01-03"),
    lubridate::ymd("2021-01-03"), "XANOMELINE",
    "P01", "ONCE", lubridate::ymd("2021-01-01"),
    lubridate::ymd("2021-01-01"), "PLACEBO",
    "P01", "ONCE", lubridate::ymd("2021-01-03"),
    lubridate::ymd("2021-01-03"), "PLACEBO",
    "P01", "ONCE", lubridate::ymd("2021-01-05"),
    lubridate::ymd("2021-01-05"), "PLACEBO"
  )

  expect_equal(create_single_dose_dataset(input), expected_output)
})

test_that("custom lookup works", {
  custom_lookup <- tibble::tribble(
    ~VALUE, ~DOSE_COUNT, ~DOSE_WINDOW, ~CONVERSION_FACTOR,
    "Q30MIN", (1 / 30), "MINUTE", 1,
    "Q90MIN", (1 / 90), "MINUTE", 1
  )

  input <- tibble::tribble(
    ~USUBJID, ~EXDOSFRQ, ~ASTDTM, ~AENDTM,
    "P01", "Q30MIN", lubridate::ymd_hms("2021-01-01T06:00:00"),
    lubridate::ymd_hms("2021-01-01T07:00:00"),
    "P02", "Q90MIN", lubridate::ymd_hms("2021-01-01T06:00:00"),
    lubridate::ymd_hms("2021-01-01T09:00:00")
  )

  expected_output <- tibble::tribble(
    ~USUBJID, ~EXDOSFRQ, ~ASTDTM, ~AENDTM,
    "P01", "ONCE", lubridate::ymd_hms("2021-01-01T06:00:00"),
    lubridate::ymd_hms("2021-01-01T06:00:00"),
    "P01", "ONCE", lubridate::ymd_hms("2021-01-01T06:30:00"),
    lubridate::ymd_hms("2021-01-01T06:30:00"),
    "P01", "ONCE", lubridate::ymd_hms("2021-01-01T07:00:00"),
    lubridate::ymd_hms("2021-01-01T07:00:00"),
    "P02", "ONCE", lubridate::ymd_hms("2021-01-01T06:00:00"),
    lubridate::ymd_hms("2021-01-01T06:00:00"),
    "P02", "ONCE", lubridate::ymd_hms("2021-01-01T07:30:00"),
    lubridate::ymd_hms("2021-01-01T07:30:00"),
    "P02", "ONCE", lubridate::ymd_hms("2021-01-01T09:00:00"),
    lubridate::ymd_hms("2021-01-01T09:00:00")
  )

  expect_equal(
    create_single_dose_dataset(input,
      lookup_table = custom_lookup,
      lookup_column = VALUE,
      start_date = ASTDTM,
      end_date = AENDTM
    ),
    expected_output
  )
})

test_that("Warning is returned when values in EXDOSFRQ does not appear in lookup table", {
  input <- tibble::tribble(
    ~USUBJID, ~EXDOSFRQ, ~ASTDT, ~AENDT,
    "P01", "1", lubridate::ymd("2021-01-01"), lubridate::ymd("2021-01-03"),
    "P01", "1", lubridate::ymd("2021-01-08"), lubridate::ymd("2021-01-12"),
    "P01", "1", lubridate::ymd("2021-01-15"), lubridate::ymd("2021-01-29")
  )
  expect_error(
    create_single_dose_dataset(input)
  )
})

test_that("Error is returned when a date variable is supplied rather than
          a datetime variable", {
  custom_lookup <- tibble::tribble(
    ~VALUE, ~DOSE_COUNT, ~DOSE_WINDOW, ~CONVERSION_FACTOR,
    "Q30MIN", (1 / 30), "MINUTE", 1,
    "Q90MIN", (1 / 90), "MINUTE", 1
  )

  input <- tibble::tribble(
    ~USUBJID, ~EXDOSFRQ, ~ASTDTM, ~AENDTM,
    "P01", "Q30MIN", lubridate::ymd("2021-01-01"),
    lubridate::ymd_hms("2021-01-01T07:00:00"),
    "P02", "Q90MIN", lubridate::ymd("2021-01-01"),
    lubridate::ymd_hms("2021-01-01T09:00:00")
  )

  expect_error(create_single_dose_dataset(input,
    lookup_table = custom_lookup,
    lookup_column = VALUE,
    start_date = ASTDTM,
    end_date = AENDTM
  ))
})

test_that("Error is returned when a date variable contains NA values", {

input <- tibble::tribble(
              ~USUBJID, ~EXDOSFRQ, ~ASTDT, ~AENDT,
              "P01", "Q2D", ymd("2021-01-01"), NA,
              "P01", "Q3D", ymd("2021-01-08"), ymd("2021-01-15"),
              "P01", "EVERY 2 WEEKS", ymd("2021-01-15"), ymd("2021-01-29")
            )
expect_error(
  create_single_dose_dataset(input), regexp = "cannot contain `NA`")

})
