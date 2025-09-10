input <- c(
  "2019-07-18T15:25:40.243",
  "2019-07-18T15:25:40",
  "2019-07-18T15:25",
  "2019-07-18",
  "2019-02",
  "2019",
  "2019---07",
  "2003-12-15T-:15:18",
  "2003-12-15T13:-:19",
  "2020-07--T00:00"
)

# impute_dtc_dtm ----
## Test 1: default: no date imputation, time part set to 00:00:00 ----
test_that("impute_dtc_dtm Test 1: default: no date imputation, time part set to 00:00:00", {
  expected_output <- c(
    "2019-07-18T15:25:40.243",
    "2019-07-18T15:25:40",
    "2019-07-18T15:25:00",
    "2019-07-18T00:00:00",
    NA_character_,
    NA_character_,
    NA_character_,
    "2003-12-15T00:00:00",
    "2003-12-15T13:00:00",
    NA_character_
  )
  expect_equal(impute_dtc_dtm(dtc = input), expected_output)
})

## Test 2: no date imputation, min and sec imputed with 59 ----
test_that("impute_dtc_dtm Test 2: no date imputation, min and sec imputed with 59", {
  expected_output <- c(
    "2019-07-18T15:25:40.243",
    "2019-07-18T15:25:40",
    "2019-07-18T15:25:59",
    NA_character_,
    NA_character_,
    NA_character_,
    NA_character_,
    NA_character_,
    "2003-12-15T13:59:59",
    NA_character_
  )
  expect_equal(
    impute_dtc_dtm(
      dtc = input,
      highest_imputation = "m",
      time_imputation = "23:59:59"
    ),
    expected_output
  )

  expect_equal(
    impute_dtc_dtm(
      dtc = input,
      highest_imputation = "m",
      time_imputation = "LAST"
    ),
    expected_output
  )
})

## Test 3: impute month and day to first, time to 00:00:00 ----
test_that("impute_dtc_dtm Test 3: impute month and day to first, time to 00:00:00", {
  expected_output <- c(
    "2019-07-18T15:25:40.243",
    "2019-07-18T15:25:40",
    "2019-07-18T15:25:00",
    "2019-07-18T00:00:00",
    "2019-02-01T00:00:00",
    "2019-01-01T00:00:00",
    "2019-01-01T00:00:00",
    "2003-12-15T00:00:00",
    "2003-12-15T13:00:00",
    "2020-07-01T00:00:00"
  )

  expect_equal(
    impute_dtc_dtm(
      dtc = input,
      highest_imputation = "M",
      date_imputation = "first"
    ),
    expected_output
  )

  expect_equal(
    impute_dtc_dtm(
      dtc = input,
      highest_imputation = "M",
      date_imputation = "01-01"
    ),
    expected_output
  )
})

## Test 4: impute day to last, time to 23:59:59 ----
test_that("impute_dtc_dtm Test 4: impute day to last, time to 23:59:59", {
  expected_output <- c(
    "2019-07-18T15:25:40.243",
    "2019-07-18T15:25:40",
    "2019-07-18T15:25:59",
    "2019-07-18T23:59:59",
    "2019-02-28T23:59:59",
    NA_character_,
    NA_character_,
    "2003-12-15T23:59:59",
    "2003-12-15T13:59:59",
    "2020-07-31T23:59:59"
  )
  expect_equal(
    impute_dtc_dtm(
      dtc = input,
      highest_imputation = "D",
      date_imputation = "last",
      time_imputation = "last"
    ),
    expected_output
  )
})

## Test 5: impute month, day to last, time to 23:59:59, preserve = TRUE ----
test_that("impute_dtc_dtm Test 5: impute month, day to last, time to 23:59:59, preserve = TRUE", {
  expected_output <- c(
    "2019-07-18T15:25:40.243",
    "2019-07-18T15:25:40",
    "2019-07-18T15:25:59",
    "2019-07-18T23:59:59",
    "2019-02-28T23:59:59",
    "2019-12-31T23:59:59",
    "2019-12-07T23:59:59",
    "2003-12-15T23:15:18",
    "2003-12-15T13:59:19",
    "2020-07-31T00:00:59"
  )
  expect_equal(
    imputes <- impute_dtc_dtm(
      dtc = input,
      highest_imputation = "M",
      date_imputation = "last",
      time_imputation = "last",
      preserve = TRUE
    ),
    expected_output
  )
})

## Test 6: no date imputation, impute second to 59 ----
test_that("impute_dtc_dtm Test 6: no date imputation, impute second to 59", {
  expected_output <- c(
    "2019-07-18T15:25:40.243",
    "2019-07-18T15:25:40",
    "2019-07-18T15:25:59",
    NA_character_,
    NA_character_,
    NA_character_,
    NA_character_,
    NA_character_,
    NA_character_,
    NA_character_
  )

  expect_equal(
    imputes <- impute_dtc_dtm(
      dtc = input,
      highest_imputation = "s",
      time_imputation = "LAST",
      preserve = FALSE
    ),
    expected_output
  )
})

## Test 7: impute month and day to mid, time to first ----
test_that("impute_dtc_dtm Test 7: impute month and day to mid, time to first", {
  expected_output <- c(
    "2019-07-18T15:25:40.243",
    "2019-07-18T15:25:40",
    "2019-07-18T15:25:00",
    "2019-07-18T00:00:00",
    "2019-02-15T00:00:00",
    "2019-06-30T00:00:00",
    "2019-06-30T00:00:00",
    "2003-12-15T00:00:00",
    "2003-12-15T13:00:00",
    "2020-07-15T00:00:00"
  )
  expect_equal(
    imputes <- impute_dtc_dtm(
      dtc = input,
      highest_imputation = "M",
      date_imputation = "mid",
      time_imputation = "first",
      preserve = FALSE
    ),
    expected_output
  )
})

## Test 8: min_dates parameter works ----
test_that("impute_dtc_dtm Test 8: min_dates parameter works", {
  expect_equal(
    impute_dtc_dtm(c("2020-12", "2020-11", NA_character_),
      min_dates = list(
        c(
          ymd_hms("2020-12-06T12:12:12"),
          NA,
          NA
        ),
        c(
          ymd_hms("2020-11-11T11:11:11"),
          ymd_hms("2020-11-11T11:11:11"),
          ymd_hms("2020-11-11T11:11:11")
        )
      ),
      highest_imputation = "Y",
      date_imputation = "first"
    ),
    c("2020-12-06T12:12:12", "2020-11-11T11:11:11", "2020-11-11T11:11:11")
  )
})

## Test 9: max_dates parameter works ----
test_that("impute_dtc_dtm Test 9: max_dates parameter works", {
  expect_equal(
    impute_dtc_dtm(c("2020-12", "2020-11", NA_character_, "2020-02-02"),
      max_dates = list(
        c(ymd_hms("2020-12-06T12:12:12"), NA, ymd_hms("2020-09-13T08:30:00"), NA),
        c(ymd(""), ymd("2020-11-11"), ymd(""), ymd("2020-02-02"))
      ),
      highest_imputation = "Y",
      date_imputation = "last",
      time_imputation = "last"
    ),
    c("2020-12-06T12:12:12", "2020-11-11T23:59:59", "2020-09-13T08:30:00", "2020-02-02T23:59:59")
  )
})

## Test 10: min_dates length mismatch provides error ----
test_that("impute_dtc_dtm Test 10: min_dates length mismatch provides error", {
  expect_snapshot(
    impute_dtc_dtm(
      c("2020-12", NA_character_),
      min_dates = list(
        c(ymd_hms("2020-12-06T12:12:12")),
        c(ymd_hms("2020-11-11T11:11:11"))
      ),
      highest_imputation = "Y"
    ),
    error = TRUE
  )
})

## Test 11: max_dates length mismatch provides error ----
test_that("impute_dtc_dtm Test 11: max_dates length mismatch provides error", {
  expect_snapshot(
    impute_dtc_dtm(
      c("2020-12", NA_character_),
      max_dates = list(
        c(ymd_hms("2020-12-06T12:12:12")),
        c(ymd_hms("2020-11-11T11:11:11"))
      ),
      highest_imputation = "Y",
      date_imputation = "last"
    ),
    error = TRUE
  )
})

## Test 12: Error if null min/max_dates when highest_imputation = Y ----
test_that("impute_dtc_dtm Test 12: Error if null min/max_dates when highest_imputation = Y", {
  expect_snapshot(
    impute_dtc_dtm(
      c("2020-12", NA_character_),
      highest_imputation = "Y"
    ),
    error = TRUE
  )
})

## Test 13: wrong input to `date_imputation` ----
test_that("impute_dtc_dtm Test 13: wrong input to `date_imputation`", {
  # date imputation is not a key when highest_imputation is "D"
  expect_snapshot(
    impute_dtc_dtm(
      dtc = c("2020-12", "2020-11", NA_character_, "2020-02-02"),
      highest_imputation = "D",
      date_imputation = "15",
      time_imputation = "last"
    ),
    error = TRUE
  )

  # wrong format
  expect_snapshot(
    impute_dtc_dtm(
      dtc = c("2020-12", "2020-11", NA_character_, "2020-02-02", "2020"),
      highest_imputation = "M",
      date_imputation = "10:12",
      time_imputation = "last"
    ),
    error = TRUE
  )

  # incomplete date_imputation
  expect_snapshot(
    impute_dtc_dtm(
      dtc = c("2020-12", "2020-11", NA_character_, "2020-02-02", "2020"),
      highest_imputation = "M",
      date_imputation = "10",
      time_imputation = "last"
    ),
    error = TRUE
  )
  # not using key for date_imputation when highest_imputation = "D"
  expect_snapshot(
    impute_dtc_dtm(
      dtc = c("2020-12", "2020-11", NA_character_, "2020-02-02", "2020"),
      highest_imputation = "D",
      date_imputation = "15",
      time_imputation = "last"
    ),
    error = TRUE
  )
  # only first or last is allowed when highest_imputation = "Y"
  expect_snapshot(
    impute_dtc_dtm(
      dtc = c("2020-12", "2020-11", NA_character_, "2020-02-02", "2020"),
      highest_imputation = "Y",
      date_imputation = "2020-02-02",
      time_imputation = "last"
    ),
    error = TRUE
  )

  # wrong time imputation
  expect_snapshot(
    impute_dtc_dtm(
      dtc = c("2020-12", "2020-11", NA_character_, "2020-02-02", "2020"),
      highest_imputation = "M",
      date_imputation = "first",
      time_imputation = "WRONG"
    ),
    error = TRUE
  )

  # impossible time
  expect_snapshot(
    impute_dtc_dtm(
      dtc = c("2020-12", "2020-11", NA_character_, "2020-02-02", "2020"),
      highest_imputation = "M",
      date_imputation = "first",
      time_imputation = "12:12:61"
    ),
    error = TRUE
  )
})

# convert_dtc_to_dtm ----
## Test 14: Convert a complete -- DTC into a date time object ----
test_that("convert_dtc_to_dtm Test 14: Convert a complete -- DTC into a date time object", {
  expect_equal(
    convert_dtc_to_dtm(input[1]),
    ymd_hms("2019-07-18T15:25:40.243")
  )
})

# compute_tmf ----
## Test 15: compute TMF ----
test_that("compute_tmf Test 15: compute TMF", {
  input_dtc <- c(
    "2019-07-18T15:25:52",
    "2019-07-18T15:25",
    "2019-07-18T15",
    "2019-07-18",
    "2019-02",
    "2019",
    "2019---07",
    "2003-12-15T-:15:18",
    "2003-12-15T13:-:19",
    "2020-07--T00:00",
    "2020-07--T00:00:00",
    "2022-05--T00:00",
    "2022-05--T23:00",
    "2022-05--T23:59:00"
  )
  input_dtm <- c(
    ymd_hms("2019-07-18T15:25:52"),
    ymd_hms("2019-07-18T15:25:00"),
    ymd_hms("2019-07-18T15:00:00"),
    as.POSIXct("2019-07-18"),
    as.POSIXct("2019-02-01"),
    as.POSIXct(NA_character_),
    as.POSIXct(NA_character_),
    ymd_hms("2003-12-15T23:15:18"),
    ymd_hms("2003-12-15T13:59:19"),
    ymd_hms("2020-07-31T00:00:59"),
    ymd_hms("2020-07-31T00:00:00"),
    ymd_hms("2022-05-15T23:59:59"),
    ymd_hms("2022-05-15T23:59:59"),
    ymd_hms("2022-05-15T23:59:59")
  )
  expected_output <- c(
    NA_character_,
    "S",
    "M",
    "H",
    "H",
    NA_character_,
    NA_character_,
    "H",
    "M",
    "S",
    NA_character_,
    "H",
    "M",
    "S"
  )

  expect_equal(
    compute_tmf(
      dtc = input_dtc,
      dtm = input_dtm,
      ignore_seconds_flag = FALSE
    ),
    expected_output
  )
})

## Test 16: throws ERROR when ignore_seconds_flag  = T (default) and seconds are present ----
test_that("compute_tmf Test 16: throws ERROR when ignore_seconds_flag  = T and seconds are present", { # nolint
  expect_snapshot(
    compute_tmf(
      dtc = c("2020-11-11T11:11:11", "2020-11-11T11:11"),
      dtm = ymd_hms(c(
        "2020-11-11T11:11:11", "2020-11-11T11:11:00"
      ))
    ),
    error = TRUE
  )
})

## Test 17: ignore_seconds_flag  = TRUE ----
test_that("compute_tmf Test 17: ignore_seconds_flag  = TRUE", {
  expect_equal(
    compute_tmf(
      dtc = c("2020-11-11T11:11", "2020-11-11T11"),
      dtm = ymd_hms(c(
        "2020-11-11T11:11:00", "2020-11-11T11:00:00"
      )),
      ignore_seconds_flag = TRUE
    ),
    c(NA_character_, "M")
  )
})

# derive_vars_dtm ----

input <- tibble::tribble(
  ~XXSTDTC,
  "2019-07-18T15:25:40",
  "2019-07-18T15:25",
  "2019-07-18T15",
  "2019-07-18",
  "2019-02",
  "2019",
  "2019---07"
)

## Test 18: default behavior ----
test_that("derive_vars_dtm Test 18: default behavior", {
  expected_output <- tibble::tribble(
    ~XXSTDTC,              ~ASTDTM,                        ~ASTTMF,
    "2019-07-18T15:25:40", ymd_hms("2019-07-18T15:25:40"), NA_character_,
    "2019-07-18T15:25",    ymd_hms("2019-07-18T15:25:00"), "S",
    "2019-07-18T15",       ymd_hms("2019-07-18T15:00:00"), "M",
    "2019-07-18",          ymd_hms("2019-07-18T00:00:00"), "H",
    "2019-02",             ymd_hms(NA),                    NA_character_,
    "2019",                ymd_hms(NA),                    NA_character_,
    "2019---07",           ymd_hms(NA),                    NA_character_
  )

  actual_output <- derive_vars_dtm(
    highest_imputation = "h",
    input,
    new_vars_prefix = "AST",
    dtc = XXSTDTC,
    ignore_seconds_flag = FALSE
  )

  expect_dfs_equal(
    base = expected_output,
    compare = actual_output,
    keys = "XXSTDTC"
  )
})

## Test 19: date imputed to first, auto DTF/TMF ----
test_that("derive_vars_dtm Test 19: date imputed to first, auto DTF/TMF", {
  expected_output <- tibble::tribble(
    ~XXSTDTC,              ~ASTDTM,                        ~ASTDTF,       ~ASTTMF,
    "2019-07-18T15:25:40", ymd_hms("2019-07-18T15:25:40"), NA_character_, NA_character_,
    "2019-07-18T15:25",    ymd_hms("2019-07-18T15:25:00"), NA_character_, "S",
    "2019-07-18T15",       ymd_hms("2019-07-18T15:00:00"), NA_character_, "M",
    "2019-07-18",          ymd_hms("2019-07-18T00:00:00"), NA_character_, "H",
    "2019-02",             ymd_hms("2019-02-01T00:00:00"), "D",           "H",
    "2019",                ymd_hms("2019-01-01T00:00:00"), "M",           "H",
    "2019---07",           ymd_hms("2019-01-01T00:00:00"), "M",           "H"
  )

  actual_output <- derive_vars_dtm(
    input,
    new_vars_prefix = "AST",
    dtc = XXSTDTC,
    highest_imputation = "M",
    date_imputation = "first",
    ignore_seconds_flag = FALSE
  )

  expect_dfs_equal(
    base = expected_output,
    compare = actual_output,
    keys = "XXSTDTC"
  )
})

## Test 20: date and time imputed to last, no DTF/TMF ----
test_that("derive_vars_dtm Test 20: date and time imputed to last, no DTF/TMF", {
  expected_output <- tibble::tribble(
    ~XXSTDTC,              ~AENDTM,
    "2019-07-18T15:25:40", ymd_hms("2019-07-18T15:25:40"),
    "2019-07-18T15:25",    ymd_hms("2019-07-18T15:25:59"),
    "2019-07-18T15",       ymd_hms("2019-07-18T15:59:59"),
    "2019-07-18",          ymd_hms("2019-07-18T23:59:59"),
    "2019-02",             ymd_hms("2019-02-28T23:59:59"),
    "2019",                ymd_hms("2019-12-31T23:59:59"),
    "2019---07",           ymd_hms("2019-12-31T23:59:59")
  )

  actual_output <- derive_vars_dtm(
    input,
    new_vars_prefix = "AEN",
    dtc = XXSTDTC,
    highest_imputation = "M",
    date_imputation = "LAST",
    time_imputation = "LAST",
    flag_imputation = "none"
  )

  expect_dfs_equal(
    base = expected_output,
    compare = actual_output,
    keys = "XXSTDTC"
  )
})

## Test 21: date and time imputed to last, DTF only ----
test_that("derive_vars_dtm Test 21: date and time imputed to last, DTF only", {
  expected_output <- tibble::tribble(
    ~XXSTDTC,              ~AENDTM,                        ~AENDTF,
    "2019-07-18T15:25:40", ymd_hms("2019-07-18T15:25:40"), NA_character_,
    "2019-07-18T15:25",    ymd_hms("2019-07-18T15:25:59"), NA_character_,
    "2019-07-18T15",       ymd_hms("2019-07-18T15:59:59"), NA_character_,
    "2019-07-18",          ymd_hms("2019-07-18T23:59:59"), NA_character_,
    "2019-02",             ymd_hms("2019-02-28T23:59:59"), "D",
    "2019",                ymd_hms("2019-12-31T23:59:59"), "M",
    "2019---07",           ymd_hms("2019-12-31T23:59:59"), "M"
  )

  actual_output <- derive_vars_dtm(
    input,
    new_vars_prefix = "AEN",
    dtc = XXSTDTC,
    highest_imputation = "M",
    date_imputation = "last",
    time_imputation = "last",
    flag_imputation = "date"
  )

  expect_dfs_equal(
    base = expected_output,
    compare = actual_output,
    keys = "XXSTDTC"
  )
})

## Test 22: date imputed to MID, time to first, TMF only ----
test_that("derive_vars_dtm Test 22: date imputed to MID, time to first, TMF only", {
  expected_output <- tibble::tribble(
    ~XXSTDTC,              ~ASTDTM,                        ~ASTTMF,
    "2019-07-18T15:25:40", ymd_hms("2019-07-18T15:25:40"), NA_character_,
    "2019-07-18T15:25",    ymd_hms("2019-07-18T15:25:00"), "S",
    "2019-07-18T15",       ymd_hms("2019-07-18T15:00:00"), "M",
    "2019-07-18",          ymd_hms("2019-07-18T00:00:00"), "H",
    "2019-02",             ymd_hms("2019-02-15T00:00:00"), "H",
    "2019",                ymd_hms("2019-06-30T00:00:00"), "H",
    "2019---07",           ymd_hms("2019-06-07T00:00:00"), "H"
  )

  actual_output <- derive_vars_dtm(
    input,
    new_vars_prefix = "AST",
    dtc = XXSTDTC,
    highest_imputation = "M",
    date_imputation = "mid",
    flag_imputation = "time",
    preserve = TRUE,
    ignore_seconds_flag = FALSE
  )

  expect_dfs_equal(
    base = expected_output,
    comp = actual_output,
    keys = c("XXSTDTC")
  )
})

## Test 23: No re-derivation is done if --DTF variable already exists ----
test_that("derive_vars_dtm Test 23: No re-derivation is done if --DTF variable already exists", {
  expected_output <- tibble::tribble(
    ~XXSTDTC,              ~ASTDTM,                        ~ASTDTF,       ~ASTTMF,
    "2019-07-18T15:25:40", ymd_hms("2019-07-18T15:25:40"), NA_character_, NA_character_,
    "2019-07-18T15:25",    ymd_hms("2019-07-18T15:25:00"), NA_character_, "S",
    "2019-07-18T15",       ymd_hms("2019-07-18T15:00:00"), NA_character_, "M",
    "2019-07-18",          ymd_hms("2019-07-18T00:00:00"), NA_character_, "H",
    "2019-02",             ymd_hms("2019-02-01T00:00:00"), "D",           "H",
    "2019",                ymd_hms("2019-01-01T00:00:00"), "MD",          "H",
    "2019---07",           ymd_hms("2019-01-01T00:00:00"), "M",           "H"
  ) %>%
    select(XXSTDTC, ASTDTF, everything())

  expect_snapshot(
    actual_output <- derive_vars_dtm(
      mutate(input, ASTDTF = c(NA, NA, NA, NA, "D", "MD", "M")),
      new_vars_prefix = "AST",
      dtc = XXSTDTC,
      highest_imputation = "M",
      date_imputation = "first",
      ignore_seconds_flag = FALSE
    )
  )

  expect_dfs_equal(
    base = expected_output,
    compare = actual_output,
    keys = "XXSTDTC"
  )
})

## Test 24: max_dates parameter works as expected ----
test_that("derive_vars_dtm Test 24: max_dates parameter works as expected", {
  expected_output <- tibble::tribble(
    ~XXSTDTC,    ~ASTDTM,                        ~ASTDTF, ~ASTTMF,
    "2019-02",   ymd_hms("2019-02-10T00:00:00"), "D",     "H",
    "2019",      ymd_hms("2019-02-10T00:00:00"), "M",     "H",
    "2019---07", ymd_hms("2019-02-10T00:00:00"), "M",     "H"
  ) %>%
    mutate(DCUTDT = ymd_hms("2019-02-10T00:00:00"))

  actual_output <- derive_vars_dtm(
    select(expected_output, XXSTDTC, DCUTDT),
    new_vars_prefix = "AST",
    dtc = XXSTDTC,
    highest_imputation = "M",
    date_imputation = "last",
    max_dates = exprs(DCUTDT)
  )

  expect_dfs_equal(
    base = expected_output,
    compare = actual_output,
    keys = c("XXSTDTC")
  )
})

input_secs <- tibble::tribble(
  ~XXSTDTC,
  "2019-07-18T15:25:40",
  "2019-07-18T15:25",
  "2019-07-18T15",
  "2019-07-18",
  "2019-02",
  "2019",
  "2019---07"
)

## Test 25: NA imputation for highest_imputation = Y & max_dates ----
test_that("derive_vars_dtm Test 25: NA imputation for highest_imputation = Y & max_dates", {
  actual <- data.frame(
    AESTDTC = c(NA_character_, NA_character_),
    TRTSDTM = c(ymd_hms("2022-01-01 23:59:59"), NA)
  ) %>%
    mutate(AESTDTC = as.character(AESTDTC)) %>%
    derive_vars_dtm(
      dtc = AESTDTC,
      new_vars_prefix = "AST",
      highest_imputation = "Y",
      date_imputation = "last",
      time_imputation = "last",
      flag_imputation = "both",
      max_dates = exprs(TRTSDTM)
    )

  expected <- data.frame(
    AESTDTC = c(NA_character_, NA_character_),
    TRTSDTM = c(ymd_hms("2022-01-01 23:59:59"), NA),
    ASTDTM  = c(ymd_hms("2022-01-01 23:59:59"), NA),
    ASTDTF  = c("Y", NA),
    ASTTMF  = c("H", NA)
  )

  expect_dfs_equal(actual, expected, keys = c("ASTDTM", "ASTDTF", "ASTTMF"))
})

## Test 26: Error for highest_imputation = Y & max_dates but date_imputation = first ----
test_that("derive_vars_dtm Test 26: Error for highest_imputation = Y & max_dates but date_imputation = first", { # nolint
  expect_snapshot(
    (data.frame(
      AESTDTC = c(NA_character_, NA_character_),
      TRTSDTM = c(ymd_hms("2022-01-01 23:59:59"), NA)
    ) %>%
      mutate(AESTDTC = as.character(AESTDTC)) %>%
      derive_vars_dtm(
        dtc = AESTDTC,
        new_vars_prefix = "AST",
        highest_imputation = "Y",
        date_imputation = "first",
        time_imputation = "first",
        flag_imputation = "both",
        max_dates = exprs(TRTSDTM)
      )),
    error = TRUE
  )
})

## Test 27: NA imputation for highest_imputation = Y & min_dates ----
test_that("derive_vars_dtm Test 27: NA imputation for highest_imputation = Y & min_dates", {
  actual <- data.frame(
    AESTDTC = c(NA_character_, NA_character_),
    TRTSDTM = c(ymd_hms("2022-01-01 23:59:59"), NA)
  ) %>%
    mutate(AESTDTC = as.character(AESTDTC)) %>%
    derive_vars_dtm(
      dtc = AESTDTC,
      new_vars_prefix = "AST",
      highest_imputation = "Y",
      date_imputation = "first",
      time_imputation = "first",
      flag_imputation = "both",
      min_dates = exprs(TRTSDTM)
    )

  expected <- data.frame(
    AESTDTC = c(NA_character_, NA_character_),
    TRTSDTM = c(ymd_hms("2022-01-01 23:59:59"), NA),
    ASTDTM  = c(ymd_hms("2022-01-01 23:59:59"), NA),
    ASTDTF  = c("Y", NA),
    ASTTMF  = c("H", NA)
  )

  expect_dfs_equal(actual, expected, keys = c("ASTDTM", "ASTDTF", "ASTTMF"))
})

## Test 28: Error for highest_imputation = Y & min_dates but date_imputation = last ----
test_that("derive_vars_dtm Test 28: Error for highest_imputation = Y & min_dates but date_imputation = last", { # nolint
  expect_snapshot(
    data.frame(
      AESTDTC = c(NA_character_, NA_character_),
      TRTSDTM = c(ymd_hms("2022-01-01 23:59:59"), NA)
    ) %>%
      mutate(AESTDTC = as.character(AESTDTC)) %>%
      derive_vars_dtm(
        dtc = AESTDTC,
        new_vars_prefix = "AST",
        highest_imputation = "Y",
        date_imputation = "last",
        time_imputation = "last",
        flag_imputation = "both",
        min_dates = exprs(TRTSDTM)
      ),
    error = TRUE
  )
})

## Test 29: NA imputation for highest_imputation = Y but null min/max dates fails ----
test_that("derive_vars_dtm Test 29: NA imputation for highest_imputation = Y but null min/max dates fails", { # nolint
  expect_snapshot(
    data.frame(
      AESTDTC = c(NA_character_, NA_character_),
      TRTSDTM = c(ymd_hms("2022-01-01 23:59:59"), NA)
    ) %>%
      mutate(AESTDTC = as.character(AESTDTC)) %>%
      derive_vars_dtm(
        dtc = AESTDTC,
        new_vars_prefix = "AST",
        highest_imputation = "Y",
        date_imputation = "first",
        time_imputation = "first",
        flag_imputation = "both"
      ),
    error = TRUE
  )
})

## Test 30: Supplying both min/max dates for highest_imputation = Y works ----
test_that("derive_vars_dtm Test 30: Supplying both min/max dates for highest_imputation = Y works", { # nolint
  actual <- data.frame(
    AESTDTC = c(NA_character_, NA_character_),
    TRTSDTM = c(ymd_hms("2022-01-01 23:59:59"), NA),
    TRTEDTM = c(ymd_hms("2022-01-31 23:59:59"), NA)
  ) %>%
    mutate(AESTDTC = as.character(AESTDTC)) %>%
    derive_vars_dtm(
      dtc = AESTDTC,
      new_vars_prefix = "AST",
      highest_imputation = "Y",
      date_imputation = "first",
      time_imputation = "first",
      min_dates = exprs(TRTSDTM),
      max_dates = exprs(TRTEDTM)
    )

  expected <- data.frame(
    AESTDTC = c(NA_character_, NA_character_),
    TRTSDTM = c(ymd_hms("2022-01-01 23:59:59"), NA),
    TRTEDTM = c(ymd_hms("2022-01-31 23:59:59"), NA),
    ASTDTM  = c(ymd_hms("2022-01-01 23:59:59"), NA),
    ASTDTF  = c("Y", NA),
    ASTTMF  = c("H", NA)
  )

  expect_dfs_equal(actual, expected, keys = c("ASTDTM", "ASTDTF", "ASTTMF"))
})

## Test 31: Supplying both min/max dates for highest_imputation = Y works ----
test_that("derive_vars_dtm Test 31: Supplying both min/max dates for highest_imputation = Y works", { # nolint
  actual <- data.frame(
    AESTDTC = c(NA_character_, NA_character_),
    TRTSDTM = c(ymd_hms("2022-01-01 23:59:59"), NA),
    TRTEDTM = c(ymd_hms("2022-01-31 23:59:59"), NA)
  ) %>%
    mutate(AESTDTC = as.character(AESTDTC)) %>%
    derive_vars_dtm(
      dtc = AESTDTC,
      new_vars_prefix = "AEN",
      highest_imputation = "Y",
      date_imputation = "last",
      time_imputation = "last",
      min_dates = exprs(TRTSDTM),
      max_dates = exprs(TRTEDTM)
    )

  expected <- data.frame(
    AESTDTC = c(NA_character_, NA_character_),
    TRTSDTM = c(ymd_hms("2022-01-01 23:59:59"), NA),
    TRTEDTM = c(ymd_hms("2022-01-31 23:59:59"), NA),
    AENDTM  = c(ymd_hms("2022-01-31 23:59:59"), NA),
    AENDTF  = c("Y", NA),
    AENTMF  = c("H", NA)
  )

  expect_dfs_equal(actual, expected, keys = c("AENDTM", "AENDTF", "AENTMF"))
})

## Test 32: catch ignore_seconds_flag error ----
test_that("derive_vars_dtm Test 32: catch ignore_seconds_flag error", {
  expect_snapshot(
    derive_vars_dtm(
      input,
      new_vars_prefix = "AST",
      dtc = XXSTDTC,
      highest_imputation = "M",
      ignore_seconds_flag = TRUE
    ),
    error = TRUE
  )
})

# Test 33: impute_dtc_dtm returns an empty character vector where dtc is empty
test_that("derive_vars_dt Test 30: impute_dtc_dt where dtc is empty", {
  empty_impute <- impute_dtc_dtm(
    dtc = character()
  )

  expect_equal(
    empty_impute,
    character(0)
  )
})

rm(list = c("input"))
