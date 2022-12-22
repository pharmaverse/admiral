# create_single_dose_dataset ----
## Test 1: Works as expected for Q#/EVERY # cases ----
test_that("cases Test 1: Works as expected for Q#/EVERY # cases", {
  input <- tibble::tribble(
    ~USUBJID, ~EXDOSFRQ,       ~ASTDT,            ~AENDT,
    "P01",    "Q2D",           ymd("2021-01-01"), ymd("2021-01-07"),
    "P01",    "Q3D",           ymd("2021-01-08"), ymd("2021-01-14"),
    "P01",    "EVERY 2 WEEKS", ymd("2021-01-15"), ymd("2021-01-29"),
    "P02",    "ONCE",          ymd("2021-02-02"), ymd("2021-02-02")
  )
  expected_output <- tibble::tribble(
    ~USUBJID, ~EXDOSFRQ, ~ASTDT,            ~AENDT,
    "P01",    "ONCE",    ymd("2021-01-01"), ymd("2021-01-01"),
    "P01",    "ONCE",    ymd("2021-01-03"), ymd("2021-01-03"),
    "P01",    "ONCE",    ymd("2021-01-05"), ymd("2021-01-05"),
    "P01",    "ONCE",    ymd("2021-01-07"), ymd("2021-01-07"),
    "P01",    "ONCE",    ymd("2021-01-08"), ymd("2021-01-08"),
    "P01",    "ONCE",    ymd("2021-01-11"), ymd("2021-01-11"),
    "P01",    "ONCE",    ymd("2021-01-14"), ymd("2021-01-14"),
    "P01",    "ONCE",    ymd("2021-01-15"), ymd("2021-01-15"),
    "P01",    "ONCE",    ymd("2021-01-29"), ymd("2021-01-29"),
    "P02",    "ONCE",    ymd("2021-02-02"), ymd("2021-02-02")
  )

  expect_dfs_equal(
    create_single_dose_dataset(input),
    expected_output,
    keys = "ASTDT"
  )
})



## Test 2: Works as expected for # TIMES PER cases ----
test_that("cases Test 2: Works as expected for # TIMES PER cases", {
  input <- tibble::tribble(
    ~USUBJID, ~DOSFREQ, ~EXSTDT, ~EXSTDTM, ~EXENDT, ~EXENDTM,
    "P01", "2 TIMES PER YEAR",
    ymd("2021-01-01"), ymd_hms("2021-01-01 10:00:00"),
    ymd("2021-07-01"), ymd_hms("2021-07-01 10:00:00"),
    "P02", "2 TIMES PER YEAR",
    ymd("2021-01-01"), ymd_hms("2021-01-01 10:30:00"),
    ymd("2021-12-31"), ymd_hms("2021-12-31 10:30:00"),
    "P03", "4 TIMES PER MONTH",
    ymd("2021-02-01"), ymd_hms("2021-02-01 11:00:00"),
    ymd("2021-03-01"), ymd_hms("2021-03-01 11:00:00"),
    "P04", "4 TIMES PER MONTH",
    ymd("2021-01-01"), ymd_hms("2021-01-01 11:30:00"),
    ymd("2021-01-20"), ymd_hms("2021-01-20 11:30:00"),
    "P05", "5 TIMES PER WEEK",
    ymd("2021-01-15"), ymd_hms("2021-01-15 12:00:00"),
    ymd("2021-01-17"), ymd_hms("2021-01-17 12:00:00"),
    "P06", "5 TIMES PER WEEK",
    ymd("2021-01-15"), ymd_hms("2021-01-15 12:30:00"),
    ymd("2021-01-21"), ymd_hms("2021-01-21 12:30:00"),
  )
  expected_output <- tibble::tribble(
    ~USUBJID, ~DOSFREQ, ~EXSTDT, ~EXSTDTM, ~EXENDT, ~EXENDTM,
    "P01", "ONCE", ymd("2021-01-01"), ymd_hms("2021-01-01 10:00:00"),
    ymd("2021-01-01"), ymd_hms("2021-01-01 10:00:00"),
    "P02", "ONCE", ymd("2021-01-01"), ymd_hms("2021-01-01 10:30:00"),
    ymd("2021-01-01"), ymd_hms("2021-01-01 10:30:00"),
    "P02", "ONCE", ymd("2021-07-02"), ymd_hms("2021-07-02 10:30:00"),
    ymd("2021-07-02"), ymd_hms("2021-07-02 10:30:00"),
    "P03", "ONCE", ymd("2021-02-01"), ymd_hms("2021-02-01 11:00:00"),
    ymd("2021-02-01"), ymd_hms("2021-02-01 11:00:00"),
    "P03", "ONCE", ymd("2021-02-08"), ymd_hms("2021-02-08 11:00:00"),
    ymd("2021-02-08"), ymd_hms("2021-02-08 11:00:00"),
    "P03", "ONCE", ymd("2021-02-16"), ymd_hms("2021-02-16 11:00:00"),
    ymd("2021-02-16"), ymd_hms("2021-02-16 11:00:00"),
    "P03", "ONCE", ymd("2021-02-23"), ymd_hms("2021-02-23 11:00:00"),
    ymd("2021-02-23"), ymd_hms("2021-02-23 11:00:00"),
    "P04", "ONCE", ymd("2021-01-01"), ymd_hms("2021-01-01 11:30:00"),
    ymd("2021-01-01"), ymd_hms("2021-01-01 11:30:00"),
    "P04", "ONCE", ymd("2021-01-08"), ymd_hms("2021-01-08 11:30:00"),
    ymd("2021-01-08"), ymd_hms("2021-01-08 11:30:00"),
    "P04", "ONCE", ymd("2021-01-16"), ymd_hms("2021-01-16 11:30:00"),
    ymd("2021-01-16"), ymd_hms("2021-01-16 11:30:00"),
    "P05", "ONCE", ymd("2021-01-15"), ymd_hms("2021-01-15 12:00:00"),
    ymd("2021-01-15"), ymd_hms("2021-01-15 12:00:00"),
    "P05", "ONCE", ymd("2021-01-16"), ymd_hms("2021-01-16 12:00:00"),
    ymd("2021-01-16"), ymd_hms("2021-01-16 12:00:00"),
    "P05", "ONCE", ymd("2021-01-17"), ymd_hms("2021-01-17 12:00:00"),
    ymd("2021-01-17"), ymd_hms("2021-01-17 12:00:00"),
    "P06", "ONCE", ymd("2021-01-15"), ymd_hms("2021-01-15 12:30:00"),
    ymd("2021-01-15"), ymd_hms("2021-01-15 12:30:00"),
    "P06", "ONCE", ymd("2021-01-16"), ymd_hms("2021-01-16 12:30:00"),
    ymd("2021-01-16"), ymd_hms("2021-01-16 12:30:00"),
    "P06", "ONCE", ymd("2021-01-17"), ymd_hms("2021-01-17 12:30:00"),
    ymd("2021-01-17"), ymd_hms("2021-01-17 12:30:00"),
    "P06", "ONCE", ymd("2021-01-19"), ymd_hms("2021-01-19 12:30:00"),
    ymd("2021-01-19"), ymd_hms("2021-01-19 12:30:00"),
    "P06", "ONCE", ymd("2021-01-20"), ymd_hms("2021-01-20 12:30:00"),
    ymd("2021-01-20"), ymd_hms("2021-01-20 12:30:00")
  )

  expect_dfs_equal(
    create_single_dose_dataset(
      input,
      dose_freq = DOSFREQ,
      start_date = EXSTDT,
      start_datetime = EXSTDTM,
      end_date = EXENDT,
      end_datetime = EXENDTM
    ),
    expected_output,
    keys = c("USUBJID", "EXSTDT")
  )
})

## Test 3: Works for different treatments ----
test_that("cases Test 3: Works for different treatments", {
  input <- tibble::tribble(
    ~USUBJID, ~EXDOSFRQ, ~ASTDT, ~ASTDTM, ~AENDT, ~AENDTM, ~EXTRT,
    "P01", "Q2D", ymd("2021-01-01"), ymd_hms("2021-01-01 09:00:00"),
    ymd("2021-01-03"), ymd_hms("2021-01-03 09:00:00"), "XANOMELINE",
    "P01", "QOD", ymd("2021-01-01"), ymd_hms("2021-01-01 09:15:00"),
    ymd("2021-01-05"), ymd_hms("2021-01-05 09:15:00"), "PLACEBO"
  )
  expected_output <- tibble::tribble(
    ~USUBJID, ~EXDOSFRQ, ~ASTDT, ~ASTDTM, ~AENDT, ~AENDTM, ~EXTRT,
    "P01", "ONCE", ymd("2021-01-01"), ymd_hms("2021-01-01 09:00:00"),
    ymd("2021-01-01"), ymd_hms("2021-01-01 09:00:00"), "XANOMELINE",
    "P01", "ONCE", ymd("2021-01-03"), ymd_hms("2021-01-03 09:00:00"),
    ymd("2021-01-03"), ymd_hms("2021-01-03 09:00:00"), "XANOMELINE",
    "P01", "ONCE", ymd("2021-01-01"), ymd_hms("2021-01-01 09:15:00"),
    ymd("2021-01-01"), ymd_hms("2021-01-01 09:15:00"), "PLACEBO",
    "P01", "ONCE", ymd("2021-01-03"), ymd_hms("2021-01-03 09:15:00"),
    ymd("2021-01-03"), ymd_hms("2021-01-03 09:15:00"), "PLACEBO",
    "P01", "ONCE", ymd("2021-01-05"), ymd_hms("2021-01-05 09:15:00"),
    ymd("2021-01-05"), ymd_hms("2021-01-05 09:15:00"), "PLACEBO"
  )

  expect_dfs_equal(
    create_single_dose_dataset(
      input,
      start_datetime = ASTDTM,
      end_datetime = AENDTM,
      keep_source_vars = vars(USUBJID, EXDOSFRQ, ASTDT, ASTDTM, AENDT, AENDTM, EXTRT)
    ),
    expected_output,
    keys = c("EXTRT", "ASTDT")
  )
})

## Test 4: Custom lookup works ----
test_that("cases Test 4: Custom lookup works", {
  custom_lookup <- tibble::tribble(
    ~VALUE, ~DOSE_COUNT, ~DOSE_WINDOW, ~CONVERSION_FACTOR,
    "Q30MIN", (1 / 30), "MINUTE", 1,
    "Q90MIN", (1 / 90), "MINUTE", 1
  )

  input <- tibble::tribble(
    ~USUBJID, ~EXDOSFRQ, ~ASTDT, ~ASTDTM, ~AENDT, ~AENDTM,
    "P01", "Q30MIN", ymd("2021-01-01"), ymd_hms("2021-01-01T06:00:00"),
    ymd("2021-01-01"), ymd_hms("2021-01-01T07:00:00"),
    "P02", "Q90MIN", ymd("2021-01-01"), ymd_hms("2021-01-01T06:00:00"),
    ymd("2021-01-01"), ymd_hms("2021-01-01T09:00:00")
  )

  expected_output <- tibble::tribble(
    ~USUBJID, ~EXDOSFRQ, ~ASTDT, ~ASTDTM, ~AENDT, ~AENDTM,
    "P01", "ONCE", ymd("2021-01-01"), ymd_hms("2021-01-01T06:00:00"),
    ymd("2021-01-01"), ymd_hms("2021-01-01T06:00:00"),
    "P01", "ONCE", ymd("2021-01-01"), ymd_hms("2021-01-01T06:30:00"),
    ymd("2021-01-01"), ymd_hms("2021-01-01T06:30:00"),
    "P01", "ONCE", ymd("2021-01-01"), ymd_hms("2021-01-01T07:00:00"),
    ymd("2021-01-01"), ymd_hms("2021-01-01T07:00:00"),
    "P02", "ONCE", ymd("2021-01-01"), ymd_hms("2021-01-01T06:00:00"),
    ymd("2021-01-01"), ymd_hms("2021-01-01T06:00:00"),
    "P02", "ONCE", ymd("2021-01-01"), ymd_hms("2021-01-01T07:30:00"),
    ymd("2021-01-01"), ymd_hms("2021-01-01T07:30:00"),
    "P02", "ONCE", ymd("2021-01-01"), ymd_hms("2021-01-01T09:00:00"),
    ymd("2021-01-01"), ymd_hms("2021-01-01T09:00:00")
  )

  expect_dfs_equal(
    create_single_dose_dataset(
      input,
      start_datetime = ASTDTM,
      end_datetime = AENDTM,
      lookup_table = custom_lookup,
      lookup_column = VALUE
    ),
    expected_output,
    keys = c("USUBJID", "ASTDTM")
  )
})

## Test 5: Warning is returned when values in EXDOSFRQ does not appear in lookup table ----
test_that("cases Test 5: Warning is returned when values in EXDOSFRQ does not appear in lookup table", { # nolint
  input <- tibble::tribble(
    ~USUBJID, ~EXDOSFRQ, ~ASTDT, ~ASTDTM, ~AENDT, ~AENDTM,
    "P01", "1", ymd("2021-01-01"), ymd_hms("2021-01-01T09:00:00"),
    ymd("2021-01-03"), ymd_hms("2021-01-03T09:00:00"),
    "P01", "1", ymd("2021-01-08"), ymd_hms("2021-01-08T09:00:00"),
    ymd("2021-01-12"), ymd_hms("2021-01-12T09:00:00"),
    "P01", "1", ymd("2021-01-15"), ymd_hms("2021-01-15T09:00:00"),
    ymd("2021-01-29"), ymd_hms("2021-01-29T09:00:00")
  )
  expect_error(
    create_single_dose_dataset(input)
  )
})

## Test 6: Error is returned when a date variable contains NA values ----
test_that("cases Test 6: Error is returned when a date variable contains NA values", { # nolint
  input <- tibble::tribble(
    ~USUBJID, ~EXDOSFRQ, ~ASTDT, ~ASTDTM, ~AENDT, ~AENDTM,
    "P01", "Q2D", ymd("2021-01-01"), ymd_hms("2021-01-01T09:00:00"), NA, NA,
    "P01", "Q3D", ymd("2021-01-08"), ymd_hms("2021-01-08T09:00:00"),
    ymd("2021-01-15"), ymd_hms("2021-01-15T09:00:00"),
    "P01", "EVERY 2 WEEKS", ymd("2021-01-15"), ymd_hms("2021-01-15T09:00:00"),
    ymd("2021-01-29"), ymd_hms("2021-01-29T09:00:00")
  )
  expect_error(
    create_single_dose_dataset(input),
    regexp = "cannot contain `NA`"
  )
})

## Test 7: Message for improper DT column names, ASTDT ----
test_that("cases Test 7: Message for improper DT column names, ASTDT", {
  input <- tibble::tribble(
    ~USUBJID, ~EXDOSFRQ, ~ADTSTD, ~ASTDTM, ~AENDT, ~AENDTM,
    "P01", "Q2D", ymd("2021-01-01"), ymd_hms("2021-01-01 10:30:00"),
    ymd("2021-01-07"), ymd_hms("2021-01-07 11:30:00"),
    "P01", "Q3D", ymd("2021-01-01"), ymd_hms("2021-01-08 12:00:00"),
    ymd("2021-01-14"), ymd_hms("2021-01-14 14:00:00"),
    "P01", "EVERY 2 WEEKS", ymd("2021-01-15"), ymd_hms("2021-01-15 09:57:00"),
    ymd("2021-01-29"), ymd_hms("2021-01-29 10:57:00")
  )
  expect_error(
    create_single_dose_dataset(input,
      start_date = ADTSTD,
      keep_source_vars = vars(
        USUBJID, EXDOSFRQ,
        ADTSTD, ASTDTM,
        AENDT, AENDTM
      )
    ),
    regexp = paste0(
      "The argument start_date is expected to have a name like xxxDT.\n",
      "Please check as it does not follow the expected naming convention"
    )
  )
})

## Test 8: Message for improper DT column names, AENDT ----
test_that("cases Test 8: Message for improper DT column names, AENDT", {
  input <- tibble::tribble(
    ~USUBJID, ~EXDOSFRQ, ~ASTDT, ~ASTDTM, ~ADTEND, ~AENDTM,
    "P01", "Q2D", ymd("2021-01-01"), ymd_hms("2021-01-01 10:30:00"),
    ymd("2021-01-07"), ymd_hms("2021-01-07 11:30:00"),
    "P01", "Q3D", ymd("2021-01-01"), ymd_hms("2021-01-08 12:00:00"),
    ymd("2021-01-14"), ymd_hms("2021-01-14 14:00:00"),
    "P01", "EVERY 2 WEEKS", ymd("2021-01-15"), ymd_hms("2021-01-15 09:57:00"),
    ymd("2021-01-29"), ymd_hms("2021-01-29 10:57:00")
  )
  expect_error(
    create_single_dose_dataset(input,
      end_date = ADTEND,
    ),
    regexp = paste0(
      "The argument end_date is expected to have a name like xxxDT.\n",
      "Please check as it does not follow the expected naming convention"
    )
  )
})

## Test 9: error if no datetime specified and freq more than QD ----
test_that("cases Test 9: error if no datetime specified and freq more than QD", {
  input <- tibble::tribble(
    ~USUBJID, ~EXDOSFRQ, ~ASTDT,            ~AENDT,
    "P01",    "Q12H",    ymd("2021-01-01"), ymd("2021-01-01"),
    "P02",    "Q12H",    ymd("2021-01-01"), ymd("2021-01-01")
  )

  expect_error(
    create_single_dose_dataset(input),
    regexp = paste(
      "There are dose frequencies more frequent than once a day.",
      "Thus `start_datetime` and `end_datetime` must be specified.",
      sep = "\n"
    ),
    fixed = TRUE
  )
})
