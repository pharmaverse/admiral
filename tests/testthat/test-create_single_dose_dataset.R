# create_single_dose_dataset ----
## Test 1: Works as expected for Q*/EVERY * cases ----
test_that("create_single_dose_dataset Test 1: Works as expected for Q*/EVERY * cases", {
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
test_that("create_single_dose_dataset Test 2: Works as expected for # TIMES PER cases", {
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
test_that("create_single_dose_dataset Test 3: Works for different treatments", {
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
      keep_source_vars = exprs(USUBJID, EXDOSFRQ, ASTDT, ASTDTM, AENDT, AENDTM, EXTRT)
    ),
    expected_output,
    keys = c("EXTRT", "ASTDT")
  )
})

## Test 4: Custom lookup works ----
test_that("create_single_dose_dataset Test 4: Custom lookup works", {
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
test_that("create_single_dose_dataset Test 5: Warning is returned when values in EXDOSFRQ does not appear in lookup table", { # nolint
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

## Test 6: Error when a date variable contains NA values ----
test_that("create_single_dose_dataset Test 6: Error when a date variable contains NA values", {
  input <- tibble::tribble(
    ~USUBJID, ~EXDOSFRQ, ~ASTDT, ~ASTDTM, ~AENDT, ~AENDTM,
    "P01", "Q2D", ymd("2021-01-01"), ymd_hms("2021-01-01T09:00:00"), NA, NA,
    "P01", "Q3D", ymd("2021-01-08"), ymd_hms("2021-01-08T09:00:00"),
    ymd("2021-01-15"), ymd_hms("2021-01-15T09:00:00"),
    "P01", "EVERY 2 WEEKS", ymd("2021-01-15"), ymd_hms("2021-01-15T09:00:00"),
    ymd("2021-01-29"), ymd_hms("2021-01-29T09:00:00")
  )
  expect_snapshot(
    error = TRUE,
    create_single_dose_dataset(input)
  )
})

## Test 7: Message for improper DT column names, ASTDT ----
test_that("create_single_dose_dataset Test 7: Message for improper DT column names, ASTDT", {
  input <- tibble::tribble(
    ~USUBJID, ~EXDOSFRQ, ~ADTSTD, ~ASTDTM, ~AENDT, ~AENDTM,
    "P01", "Q2D", ymd("2021-01-01"), ymd_hms("2021-01-01 10:30:00"),
    ymd("2021-01-07"), ymd_hms("2021-01-07 11:30:00"),
    "P01", "Q3D", ymd("2021-01-01"), ymd_hms("2021-01-08 12:00:00"),
    ymd("2021-01-14"), ymd_hms("2021-01-14 14:00:00"),
    "P01", "EVERY 2 WEEKS", ymd("2021-01-15"), ymd_hms("2021-01-15 09:57:00"),
    ymd("2021-01-29"), ymd_hms("2021-01-29 10:57:00")
  )
  expect_snapshot(
    error = TRUE,
    create_single_dose_dataset(input,
      start_date = ADTSTD,
      keep_source_vars = exprs(
        USUBJID, EXDOSFRQ,
        ADTSTD, ASTDTM,
        AENDT, AENDTM
      )
    )
  )
})

## Test 8: Message for improper DT column names, AENDT ----
test_that("create_single_dose_dataset Test 8: Message for improper DT column names, AENDT", {
  input <- tibble::tribble(
    ~USUBJID, ~EXDOSFRQ, ~ASTDT, ~ASTDTM, ~ADTEND, ~AENDTM,
    "P01", "Q2D", ymd("2021-01-01"), ymd_hms("2021-01-01 10:30:00"),
    ymd("2021-01-07"), ymd_hms("2021-01-07 11:30:00"),
    "P01", "Q3D", ymd("2021-01-01"), ymd_hms("2021-01-08 12:00:00"),
    ymd("2021-01-14"), ymd_hms("2021-01-14 14:00:00"),
    "P01", "EVERY 2 WEEKS", ymd("2021-01-15"), ymd_hms("2021-01-15 09:57:00"),
    ymd("2021-01-29"), ymd_hms("2021-01-29 10:57:00")
  )
  expect_snapshot(
    error = TRUE,
    create_single_dose_dataset(input,
      end_date = ADTEND,
    )
  )
})

## Test 9: error if no datetime and freq more than QD ----
test_that("create_single_dose_dataset Test 9: error if no datetime and freq more than QD", {
  input <- tibble::tribble(
    ~USUBJID, ~EXDOSFRQ, ~ASTDT,            ~AENDT,
    "P01",    "Q12H",    ymd("2021-01-01"), ymd("2021-01-01"),
    "P02",    "Q12H",    ymd("2021-01-01"), ymd("2021-01-01")
  )

  expect_snapshot(
    error = TRUE,
    create_single_dose_dataset(input)
  )
})

## Test 10: Works as expected for BID cases ----
test_that("create_single_dose_dataset Test 10: Works as expected for BID cases", {
  input <- tibble::tribble(
    ~USUBJID, ~EXDOSFRQ, ~ASTDT, ~ASTDTM, ~AENDT, ~AENDTM,
    "P01", "BID", ymd("2021-01-01"), ymd_hms("2021-01-01 08:00:00"),
    ymd("2021-01-03"), ymd_hms("2021-01-03 20:00:00"),
    "P01", "BID", ymd("2021-01-04"), ymd_hms("2021-01-04 08:00:00"),
    ymd("2021-01-06"), ymd_hms("2021-01-06 20:00:00")
  )
  expected_output <- tibble::tribble(
    ~USUBJID, ~EXDOSFRQ, ~ASTDT, ~ASTDTM, ~AENDT, ~AENDTM,
    "P01", "ONCE", ymd("2021-01-01"), ymd_hms("2021-01-01 08:00:00"),
    ymd("2021-01-01"), ymd_hms("2021-01-01 08:00:00"),
    "P01", "ONCE", ymd("2021-01-01"), ymd_hms("2021-01-01 20:00:00"),
    ymd("2021-01-01"), ymd_hms("2021-01-01 20:00:00"),
    "P01", "ONCE", ymd("2021-01-02"), ymd_hms("2021-01-02 08:00:00"),
    ymd("2021-01-02"), ymd_hms("2021-01-02 08:00:00"),
    "P01", "ONCE", ymd("2021-01-02"), ymd_hms("2021-01-02 20:00:00"),
    ymd("2021-01-02"), ymd_hms("2021-01-02 20:00:00"),
    "P01", "ONCE", ymd("2021-01-03"), ymd_hms("2021-01-03 08:00:00"),
    ymd("2021-01-03"), ymd_hms("2021-01-03 08:00:00"),
    "P01", "ONCE", ymd("2021-01-03"), ymd_hms("2021-01-03 20:00:00"),
    ymd("2021-01-03"), ymd_hms("2021-01-03 20:00:00"),
    "P01", "ONCE", ymd("2021-01-04"), ymd_hms("2021-01-04 08:00:00"),
    ymd("2021-01-04"), ymd_hms("2021-01-04 08:00:00"),
    "P01", "ONCE", ymd("2021-01-04"), ymd_hms("2021-01-04 20:00:00"),
    ymd("2021-01-04"), ymd_hms("2021-01-04 20:00:00"),
    "P01", "ONCE", ymd("2021-01-05"), ymd_hms("2021-01-05 08:00:00"),
    ymd("2021-01-05"), ymd_hms("2021-01-05 08:00:00"),
    "P01", "ONCE", ymd("2021-01-05"), ymd_hms("2021-01-05 20:00:00"),
    ymd("2021-01-05"), ymd_hms("2021-01-05 20:00:00"),
    "P01", "ONCE", ymd("2021-01-06"), ymd_hms("2021-01-06 08:00:00"),
    ymd("2021-01-06"), ymd_hms("2021-01-06 08:00:00"),
    "P01", "ONCE", ymd("2021-01-06"), ymd_hms("2021-01-06 20:00:00"),
    ymd("2021-01-06"), ymd_hms("2021-01-06 20:00:00"),
  )

  expect_dfs_equal(
    create_single_dose_dataset(
      dataset = input,
      dose_freq = EXDOSFRQ,
      start_date = ASTDT,
      start_datetime = ASTDTM,
      end_date = AENDT,
      end_datetime = AENDTM,
      lookup_table = dose_freq_lookup,
      lookup_column = CDISC_VALUE,
      keep_source_vars = exprs(
        USUBJID, EXDOSFRQ, ASTDT, ASTDTM, AENDT, AENDTM
      )
    ),
    expected_output,
    keys = "ASTDTM"
  )
})

## Test 11: Works as expected for cases with nominal time ----
test_that("create_single_dose_dataset Test 11: Works as expected for cases with nominal time", {
  input <- tibble::tribble(
    ~USUBJID, ~EXDOSFRQ, ~ASTDT, ~ASTDTM, ~AENDT, ~AENDTM, ~NFRLT,
    "P01", "QD", ymd("2021-01-01"), ymd_hms("2021-01-01 08:00:00"),
    ymd("2021-01-07"), ymd_hms("2021-01-07 08:00:00"), 0,
    "P01", "QD", ymd("2021-01-08"), ymd_hms("2021-01-08 08:00:00"),
    ymd("2021-01-14"), ymd_hms("2021-01-14 08:00:00"), 168
  )
  expected_output <- tibble::tribble(
    ~USUBJID, ~EXDOSFRQ, ~ASTDT, ~ASTDTM, ~AENDT, ~AENDTM, ~NFRLT,
    "P01", "ONCE", ymd("2021-01-01"), ymd_hms("2021-01-01 08:00:00"),
    ymd("2021-01-01"), ymd_hms("2021-01-01 08:00:00"), 0,
    "P01", "ONCE", ymd("2021-01-02"), ymd_hms("2021-01-02 08:00:00"),
    ymd("2021-01-02"), ymd_hms("2021-01-02 08:00:00"), 24,
    "P01", "ONCE", ymd("2021-01-03"), ymd_hms("2021-01-03 08:00:00"),
    ymd("2021-01-03"), ymd_hms("2021-01-03 08:00:00"), 48,
    "P01", "ONCE", ymd("2021-01-04"), ymd_hms("2021-01-04 08:00:00"),
    ymd("2021-01-04"), ymd_hms("2021-01-04 08:00:00"), 72,
    "P01", "ONCE", ymd("2021-01-05"), ymd_hms("2021-01-05 08:00:00"),
    ymd("2021-01-05"), ymd_hms("2021-01-05 08:00:00"), 96,
    "P01", "ONCE", ymd("2021-01-06"), ymd_hms("2021-01-06 08:00:00"),
    ymd("2021-01-06"), ymd_hms("2021-01-06 08:00:00"), 120,
    "P01", "ONCE", ymd("2021-01-07"), ymd_hms("2021-01-07 08:00:00"),
    ymd("2021-01-07"), ymd_hms("2021-01-07 08:00:00"), 144,
    "P01", "ONCE", ymd("2021-01-08"), ymd_hms("2021-01-08 08:00:00"),
    ymd("2021-01-08"), ymd_hms("2021-01-08 08:00:00"), 168,
    "P01", "ONCE", ymd("2021-01-09"), ymd_hms("2021-01-09 08:00:00"),
    ymd("2021-01-09"), ymd_hms("2021-01-09 08:00:00"), 192,
    "P01", "ONCE", ymd("2021-01-10"), ymd_hms("2021-01-10 08:00:00"),
    ymd("2021-01-10"), ymd_hms("2021-01-10 08:00:00"), 216,
    "P01", "ONCE", ymd("2021-01-11"), ymd_hms("2021-01-11 08:00:00"),
    ymd("2021-01-11"), ymd_hms("2021-01-11 08:00:00"), 240,
    "P01", "ONCE", ymd("2021-01-12"), ymd_hms("2021-01-12 08:00:00"),
    ymd("2021-01-12"), ymd_hms("2021-01-12 08:00:00"), 264,
    "P01", "ONCE", ymd("2021-01-13"), ymd_hms("2021-01-13 08:00:00"),
    ymd("2021-01-13"), ymd_hms("2021-01-13 08:00:00"), 288,
    "P01", "ONCE", ymd("2021-01-14"), ymd_hms("2021-01-14 08:00:00"),
    ymd("2021-01-14"), ymd_hms("2021-01-14 08:00:00"), 312,
  )

  expect_dfs_equal(
    create_single_dose_dataset(
      dataset = input,
      dose_freq = EXDOSFRQ,
      start_date = ASTDT,
      start_datetime = ASTDTM,
      end_date = AENDT,
      end_datetime = AENDTM,
      lookup_table = dose_freq_lookup,
      lookup_column = CDISC_VALUE,
      nominal_time = NFRLT,
      keep_source_vars = exprs(
        USUBJID, EXDOSFRQ, ASTDT, ASTDTM, AENDT, AENDTM, NFRLT
      )
    ),
    expected_output,
    keys = "ASTDTM"
  )
})

## Test 12: Error if lookup_column contains duplicates ----
test_that("create_single_dose_dataset Test 12: Error if lookup_column contains duplicates", {
  custom_lookup <- tribble(
    ~Value,   ~DOSE_COUNT, ~DOSE_WINDOW, ~CONVERSION_FACTOR,
    "Q30MIN", (1 / 30),    "MINUTE",                      1,
    "Q30MIN", (1 / 30),    "MINUTE",                      1,
    "Q90MIN", (1 / 90),    "MINUTE",                      1
  )

  input <- tribble(
    ~USUBJID, ~EXDOSFRQ, ~ASTDT, ~ASTDTM, ~AENDT, ~AENDTM,
    "P01", "Q30MIN", ymd("2021-01-01"), ymd_hms("2021-01-01T06:00:00"),
    ymd("2021-01-01"), ymd_hms("2021-01-01T07:00:00"),
    "P02", "Q90MIN", ymd("2021-01-01"), ymd_hms("2021-01-01T06:00:00"),
    ymd("2021-01-01"), ymd_hms("2021-01-01T09:00:00")
  )

  expect_error(
    create_single_dose_dataset(input,
      lookup_table = custom_lookup,
      lookup_column = Value,
      start_datetime = ASTDTM,
      end_datetime = AENDTM
    ),
    regexp = paste0(
      "Dataset contains duplicate records with respect to `Value`"
    )
  )
})
