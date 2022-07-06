test_that("create_single_dose_dataset works as expected for Q#/EVERY # cases", {
  input <- tibble::tribble(
    ~USUBJID, ~EXDOSFRQ, ~ASTDT, ~ASTDTM, ~AENDT, ~AENDTM,
    "P01", "Q2D", ymd("2021-01-01"), ymd_hms("2021-01-01 10:30:00"),
    ymd("2021-01-07"), ymd_hms("2021-01-07 11:30:00"),
    "P01", "Q3D", ymd("2021-01-08"), ymd_hms("2021-01-08 12:00:00"),
    ymd("2021-01-14"), ymd_hms("2021-01-14 14:00:00"),
    "P01", "EVERY 2 WEEKS", ymd("2021-01-15"), ymd_hms("2021-01-15 09:57:00"),
    ymd("2021-01-29"), ymd_hms("2021-01-29 10:57:00")
  )
  expected_output <- tibble::tribble(
    ~USUBJID, ~EXDOSFRQ, ~ASTDT, ~ASTDTM, ~AENDT, ~AENDTM,
    "P01", "ONCE", lubridate::ymd("2021-01-01"), lubridate::ymd_hms("2021-01-01 10:30:00"),
    lubridate::ymd("2021-01-01"), lubridate::ymd_hms("2021-01-01 11:30:00"),
    "P01", "ONCE", lubridate::ymd("2021-01-03"), lubridate::ymd_hms("2021-01-03 10:30:00"),
    lubridate::ymd("2021-01-03"), lubridate::ymd_hms("2021-01-03 11:30:00"),
    "P01", "ONCE", lubridate::ymd("2021-01-05"), lubridate::ymd_hms("2021-01-05 10:30:00"),
    lubridate::ymd("2021-01-05"), lubridate::ymd_hms("2021-01-05 11:30:00"),
    "P01", "ONCE", lubridate::ymd("2021-01-07"), lubridate::ymd_hms("2021-01-07 10:30:00"),
    lubridate::ymd("2021-01-07"), lubridate::ymd_hms("2021-01-07 11:30:00"),
    "P01", "ONCE", lubridate::ymd("2021-01-08"), lubridate::ymd_hms("2021-01-08 12:00:00"),
    lubridate::ymd("2021-01-08"), lubridate::ymd_hms("2021-01-08 14:00:00"),
    "P01", "ONCE", lubridate::ymd("2021-01-11"), lubridate::ymd_hms("2021-01-11 12:00:00"),
    lubridate::ymd("2021-01-11"), lubridate::ymd_hms("2021-01-11 14:00:00"),
    "P01", "ONCE", lubridate::ymd("2021-01-14"), lubridate::ymd_hms("2021-01-14 12:00:00"),
    lubridate::ymd("2021-01-14"), lubridate::ymd_hms("2021-01-14 14:00:00"),
    "P01", "ONCE", lubridate::ymd("2021-01-15"), lubridate::ymd_hms("2021-01-15 09:57:00"),
    lubridate::ymd("2021-01-15"), lubridate::ymd_hms("2021-01-15 10:57:00"),
    "P01", "ONCE", lubridate::ymd("2021-01-29"), lubridate::ymd_hms("2021-01-29 09:57:00"),
    lubridate::ymd("2021-01-29"), lubridate::ymd_hms("2021-01-29 10:57:00")
  )

  expect_equal(create_single_dose_dataset(input), expected_output)
})


test_that("create_single_dose_dataset works as expected for # TIMES PER cases", {
  input <- tibble::tribble(
    ~USUBJID, ~EXDOSFRQ, ~ASTDT, ~ASTDTM, ~AENDT, ~AENDTM,
    "P01", "2 TIMES PER YEAR",
    lubridate::ymd("2021-01-01"), lubridate::ymd_hms("2021-01-01 10:00:00"),
    lubridate::ymd("2021-07-01"), lubridate::ymd_hms("2021-07-01 10:00:00"),
    "P02", "2 TIMES PER YEAR",
    lubridate::ymd("2021-01-01"), lubridate::ymd_hms("2021-01-01 10:30:00"),
    lubridate::ymd("2021-12-31"), lubridate::ymd_hms("2021-12-31 10:30:00"),
    "P03", "4 TIMES PER MONTH",
    lubridate::ymd("2021-02-01"), lubridate::ymd_hms("2021-02-01 11:00:00"),
    lubridate::ymd("2021-03-01"), lubridate::ymd_hms("2021-03-01 11:00:00"),
    "P04", "4 TIMES PER MONTH",
    lubridate::ymd("2021-01-01"), lubridate::ymd_hms("2021-01-01 11:30:00"),
    lubridate::ymd("2021-01-20"), lubridate::ymd_hms("2021-01-20 11:30:00"),
    "P05", "5 TIMES PER WEEK",
    lubridate::ymd("2021-01-15"), lubridate::ymd_hms("2021-01-15 12:00:00"),
    lubridate::ymd("2021-01-17"), lubridate::ymd_hms("2021-01-17 12:00:00"),
    "P06", "5 TIMES PER WEEK",
    lubridate::ymd("2021-01-15"), lubridate::ymd_hms("2021-01-15 12:30:00"),
    lubridate::ymd("2021-01-21"), lubridate::ymd_hms("2021-01-21 12:30:00"),
  )
  expected_output <- tibble::tribble(
    ~USUBJID, ~EXDOSFRQ, ~ASTDT, ~ASTDTM, ~AENDT, ~AENDTM,
    "P01", "ONCE", lubridate::ymd("2021-01-01"), lubridate::ymd_hms("2021-01-01 10:00:00"),
    lubridate::ymd("2021-01-01"), lubridate::ymd_hms("2021-01-01 10:00:00"),
    "P02", "ONCE", lubridate::ymd("2021-01-01"), lubridate::ymd_hms("2021-01-01 10:30:00"),
    lubridate::ymd("2021-01-01"), lubridate::ymd_hms("2021-01-01 10:30:00"),
    "P02", "ONCE", lubridate::ymd("2021-07-02"), lubridate::ymd_hms("2021-07-02 10:30:00"),
    lubridate::ymd("2021-07-02"), lubridate::ymd_hms("2021-07-02 10:30:00"),
    "P03", "ONCE", lubridate::ymd("2021-02-01"), lubridate::ymd_hms("2021-02-01 11:00:00"),
    lubridate::ymd("2021-02-01"), lubridate::ymd_hms("2021-02-01 11:00:00"),
    "P03", "ONCE", lubridate::ymd("2021-02-08"), lubridate::ymd_hms("2021-02-08 11:00:00"),
    lubridate::ymd("2021-02-08"), lubridate::ymd_hms("2021-02-08 11:00:00"),
    "P03", "ONCE", lubridate::ymd("2021-02-16"), lubridate::ymd_hms("2021-02-16 11:00:00"),
    lubridate::ymd("2021-02-16"), lubridate::ymd_hms("2021-02-16 11:00:00"),
    "P03", "ONCE", lubridate::ymd("2021-02-23"), lubridate::ymd_hms("2021-02-23 11:00:00"),
    lubridate::ymd("2021-02-23"), lubridate::ymd_hms("2021-02-23 11:00:00"),
    "P04", "ONCE", lubridate::ymd("2021-01-01"), lubridate::ymd_hms("2021-01-01 11:30:00"),
    lubridate::ymd("2021-01-01"), lubridate::ymd_hms("2021-01-01 11:30:00"),
    "P04", "ONCE", lubridate::ymd("2021-01-08"), lubridate::ymd_hms("2021-01-08 11:30:00"),
    lubridate::ymd("2021-01-08"), lubridate::ymd_hms("2021-01-08 11:30:00"),
    "P04", "ONCE", lubridate::ymd("2021-01-16"), lubridate::ymd_hms("2021-01-16 11:30:00"),
    lubridate::ymd("2021-01-16"), lubridate::ymd_hms("2021-01-16 11:30:00"),
    "P05", "ONCE", lubridate::ymd("2021-01-15"), lubridate::ymd_hms("2021-01-15 12:00:00"),
    lubridate::ymd("2021-01-15"), lubridate::ymd_hms("2021-01-15 12:00:00"),
    "P05", "ONCE", lubridate::ymd("2021-01-16"), lubridate::ymd_hms("2021-01-16 12:00:00"),
    lubridate::ymd("2021-01-16"), lubridate::ymd_hms("2021-01-16 12:00:00"),
    "P05", "ONCE", lubridate::ymd("2021-01-17"), lubridate::ymd_hms("2021-01-17 12:00:00"),
    lubridate::ymd("2021-01-17"), lubridate::ymd_hms("2021-01-17 12:00:00"),
    "P06", "ONCE", lubridate::ymd("2021-01-15"), lubridate::ymd_hms("2021-01-15 12:30:00"),
    lubridate::ymd("2021-01-15"), lubridate::ymd_hms("2021-01-15 12:30:00"),
    "P06", "ONCE", lubridate::ymd("2021-01-16"), lubridate::ymd_hms("2021-01-16 12:30:00"),
    lubridate::ymd("2021-01-16"), lubridate::ymd_hms("2021-01-16 12:30:00"),
    "P06", "ONCE", lubridate::ymd("2021-01-17"), lubridate::ymd_hms("2021-01-17 12:30:00"),
    lubridate::ymd("2021-01-17"), lubridate::ymd_hms("2021-01-17 12:30:00"),
    "P06", "ONCE", lubridate::ymd("2021-01-19"), lubridate::ymd_hms("2021-01-19 12:30:00"),
    lubridate::ymd("2021-01-19"), lubridate::ymd_hms("2021-01-19 12:30:00"),
    "P06", "ONCE", lubridate::ymd("2021-01-20"), lubridate::ymd_hms("2021-01-20 12:30:00"),
    lubridate::ymd("2021-01-20"), lubridate::ymd_hms("2021-01-20 12:30:00")
  )

  expect_equal(create_single_dose_dataset(input), expected_output)
})

test_that("create_single_dose_dataset works for different treatments", {
  input <- tibble::tribble(
    ~USUBJID, ~EXDOSFRQ, ~ASTDT, ~ASTDTM, ~AENDT, ~AENDTM, ~EXTRT,
    "P01", "Q2D", lubridate::ymd("2021-01-01"), lubridate::ymd_hms("2021-01-01 09:00:00"),
    lubridate::ymd("2021-01-03"), lubridate::ymd_hms("2021-01-03 09:00:00"), "XANOMELINE",
    "P01", "QOD", lubridate::ymd("2021-01-01"), lubridate::ymd_hms("2021-01-01 09:15:00"),
    lubridate::ymd("2021-01-05"), lubridate::ymd_hms("2021-01-05 09:15:00"), "PLACEBO"
  )
  expected_output <- tibble::tribble(
    ~USUBJID, ~EXDOSFRQ, ~ASTDT, ~ASTDTM, ~AENDT, ~AENDTM, ~EXTRT,
    "P01", "ONCE", lubridate::ymd("2021-01-01"), lubridate::ymd_hms("2021-01-01 09:00:00"),
    lubridate::ymd("2021-01-01"), lubridate::ymd_hms("2021-01-01 09:00:00"), "XANOMELINE",
    "P01", "ONCE", lubridate::ymd("2021-01-03"), lubridate::ymd_hms("2021-01-03 09:00:00"),
    lubridate::ymd("2021-01-03"), lubridate::ymd_hms("2021-01-03 09:00:00"), "XANOMELINE",
    "P01", "ONCE", lubridate::ymd("2021-01-01"), lubridate::ymd_hms("2021-01-01 09:15:00"),
    lubridate::ymd("2021-01-01"), lubridate::ymd_hms("2021-01-01 09:15:00"), "PLACEBO",
    "P01", "ONCE", lubridate::ymd("2021-01-03"), lubridate::ymd_hms("2021-01-03 09:15:00"),
    lubridate::ymd("2021-01-03"), lubridate::ymd_hms("2021-01-03 09:15:00"), "PLACEBO",
    "P01", "ONCE", lubridate::ymd("2021-01-05"), lubridate::ymd_hms("2021-01-05 09:15:00"),
    lubridate::ymd("2021-01-05"), lubridate::ymd_hms("2021-01-05 09:15:00"), "PLACEBO"
  )

  expect_equal(
    create_single_dose_dataset(input,
      keep_source_vars = vars(USUBJID, EXDOSFRQ, ASTDT, ASTDTM, AENDT, AENDTM, EXTRT)
    ),
    expected_output
  )
})

test_that("custom lookup works", {
  custom_lookup <- tibble::tribble(
    ~VALUE, ~DOSE_COUNT, ~DOSE_WINDOW, ~CONVERSION_FACTOR,
    "Q30MIN", (1 / 30), "MINUTE", 1,
    "Q90MIN", (1 / 90), "MINUTE", 1
  )

  input <- tibble::tribble(
    ~USUBJID, ~EXDOSFRQ, ~ASTDT, ~ASTDTM, ~AENDT, ~AENDTM,
    "P01", "Q30MIN", lubridate::ymd("2021-01-01"), lubridate::ymd_hms("2021-01-01T06:00:00"),
    lubridate::ymd("2021-01-01"), lubridate::ymd_hms("2021-01-01T07:00:00"),
    "P02", "Q90MIN", lubridate::ymd("2021-01-01"), lubridate::ymd_hms("2021-01-01T06:00:00"),
    lubridate::ymd("2021-01-01"), lubridate::ymd_hms("2021-01-01T09:00:00")
  )

  expected_output <- tibble::tribble(
    ~USUBJID, ~EXDOSFRQ, ~ASTDT, ~ASTDTM, ~AENDT, ~AENDTM,
    "P01", "ONCE", lubridate::ymd("2021-01-01"), lubridate::ymd_hms("2021-01-01T06:00:00"),
    lubridate::ymd("2021-01-01"), lubridate::ymd_hms("2021-01-01T06:00:00"),
    "P01", "ONCE", lubridate::ymd("2021-01-01"), lubridate::ymd_hms("2021-01-01T06:30:00"),
    lubridate::ymd("2021-01-01"), lubridate::ymd_hms("2021-01-01T06:30:00"),
    "P01", "ONCE", lubridate::ymd("2021-01-01"), lubridate::ymd_hms("2021-01-01T07:00:00"),
    lubridate::ymd("2021-01-01"), lubridate::ymd_hms("2021-01-01T07:00:00"),
    "P02", "ONCE", lubridate::ymd("2021-01-01"), lubridate::ymd_hms("2021-01-01T06:00:00"),
    lubridate::ymd("2021-01-01"), lubridate::ymd_hms("2021-01-01T06:00:00"),
    "P02", "ONCE", lubridate::ymd("2021-01-01"), lubridate::ymd_hms("2021-01-01T07:30:00"),
    lubridate::ymd("2021-01-01"), lubridate::ymd_hms("2021-01-01T07:30:00"),
    "P02", "ONCE", lubridate::ymd("2021-01-01"), lubridate::ymd_hms("2021-01-01T09:00:00"),
    lubridate::ymd("2021-01-01"), lubridate::ymd_hms("2021-01-01T09:00:00")
  )

  expect_equal(
    create_single_dose_dataset(input,
      lookup_table = custom_lookup,
      lookup_column = VALUE
    ),
    expected_output
  )
})

test_that("Warning is returned when values in EXDOSFRQ does not appear in lookup table", {
  input <- tibble::tribble(
    ~USUBJID, ~EXDOSFRQ, ~ASTDT, ~ASTDTM, ~AENDT, ~AENDTM,
    "P01", "1", lubridate::ymd("2021-01-01"), lubridate::ymd_hms("2021-01-01T09:00:00"),
    lubridate::ymd("2021-01-03"), lubridate::ymd_hms("2021-01-03T09:00:00"),
    "P01", "1", lubridate::ymd("2021-01-08"), lubridate::ymd_hms("2021-01-08T09:00:00"),
    lubridate::ymd("2021-01-12"), lubridate::ymd_hms("2021-01-12T09:00:00"),
    "P01", "1", lubridate::ymd("2021-01-15"), lubridate::ymd_hms("2021-01-15T09:00:00"),
    lubridate::ymd("2021-01-29"), lubridate::ymd_hms("2021-01-29T09:00:00")
  )
  expect_error(
    create_single_dose_dataset(input)
  )
})

test_that("Error is returned when a date variable contains NA values", {
  input <- tibble::tribble(
    ~USUBJID, ~EXDOSFRQ, ~ASTDT, ~ASTDTM, ~AENDT, ~AENDTM,
    "P01", "Q2D", ymd("2021-01-01"), lubridate::ymd_hms("2021-01-01T09:00:00"), NA, NA,
    "P01", "Q3D", ymd("2021-01-08"), lubridate::ymd_hms("2021-01-08T09:00:00"),
    ymd("2021-01-15"), lubridate::ymd_hms("2021-01-15T09:00:00"),
    "P01", "EVERY 2 WEEKS", ymd("2021-01-15"), lubridate::ymd_hms("2021-01-15T09:00:00"),
    ymd("2021-01-29"), lubridate::ymd_hms("2021-01-29T09:00:00")
  )
  expect_error(
    create_single_dose_dataset(input),
    regexp = "cannot contain `NA`"
  )
})
