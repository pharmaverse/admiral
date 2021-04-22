context("test-derive_vars_dtm")

date <- tibble::tribble(
  ~XXSTDTC,
  "2019-07-18T15:25:40",
  "2019-07-18T15:25",
  "2019-07-18T15",
  "2019-07-18",
  "2019-02",
  "2019",
  "2019---07"
)

test_that("default: no date imputation, time part set o 00:00:00, add DTF, TMF", {
  expected_output <- tibble::tribble(
    ~XXSTDTC, ~ASTDTM, ~ASTDTF, ~ASTTMF,
    "2019-07-18T15:25:40", ymd_hms("2019-07-18T15:25:40"), "", "",
    "2019-07-18T15:25", ymd_hms("2019-07-18T15:25:00"), "", "S",
    "2019-07-18T15", ymd_hms("2019-07-18T15:00:00"), "", "M",
    "2019-07-18", ymd_hms("2019-07-18T00:00:00"), "", "H",
    "2019-02", ymd_hms(""), "", "",
    "2019", ymd_hms(""), "", "",
    "2019---07", ymd_hms(""), "", "",
  )

  actual_output <- derive_vars_dtm(date,
    new_vars_prefix = "AST",
    dtc = XXSTDTC
  )

  expect_equal(
    expected_output,
    actual_output
  )
})

test_that("Partial date imputed to the first day/month", {
  expected_output <- tibble::tribble(
    ~XXSTDTC, ~ASTDTM, ~ASTDTF, ~ASTTMF,
    "2019-07-18T15:25:40", ymd_hms("2019-07-18T15:25:40"), "", "",
    "2019-07-18T15:25", ymd_hms("2019-07-18T15:25:00"), "", "S",
    "2019-07-18T15", ymd_hms("2019-07-18T15:00:00"), "", "M",
    "2019-07-18", ymd_hms("2019-07-18T00:00:00"), "", "H",
    "2019-02", ymd_hms("2019-02-01T00:00:00"), "D", "H",
    "2019", ymd_hms("2019-01-01T00:00:00"), "M", "H",
    "2019---07", ymd_hms("2019-01-01T00:00:00"), "M", "H"
  )

  actual_output <- derive_vars_dtm(date,
    new_vars_prefix = "AST",
    dtc = XXSTDTC,
    date_imputation = "FIRST"
  )
  actual_output1 <- derive_vars_dtm(date,
    new_vars_prefix = "AST",
    dtc = XXSTDTC,
    date_imputation = "01-01"
  )

  expect_equal(
    expected_output,
    actual_output
  )
  expect_equal(
    expected_output,
    actual_output1
  )
})

test_that("Partial date imputed to the last day/month, Missing time part imputed with 23:59:59", {
  expected_output <- tibble::tribble(
    ~XXSTDTC, ~AENDTM, ~AENDTF, ~AENTMF,
    "2019-07-18T15:25:40", ymd_hms("2019-07-18T15:25:40"), "", "",
    "2019-07-18T15:25", ymd_hms("2019-07-18T15:25:59"), "", "S",
    "2019-07-18T15", ymd_hms("2019-07-18T15:59:59"), "", "M",
    "2019-07-18", ymd_hms("2019-07-18T23:59:59"), "", "H",
    "2019-02", ymd_hms("2019-02-28T23:59:59"), "D", "H",
    "2019", ymd_hms("2019-12-31T23:59:59"), "M", "H",
    "2019---07", ymd_hms("2019-12-31T23:59:59"), "M", "H"
  )

  actual_output <- derive_vars_dtm(date,
    new_vars_prefix = "AEN",
    dtc = XXSTDTC,
    date_imputation = "LAST",
    time_imputation = "LAST"
  )

  actual_output1 <- derive_vars_dtm(date,
    new_vars_prefix = "AEN",
    dtc = XXSTDTC,
    date_imputation = "LAST",
    time_imputation = "23:59:59"
  )

  expect_equal(
    expected_output,
    actual_output
  )
  expect_equal(
    expected_output,
    actual_output1
  )
})

test_that("Partial date imputed to the last day/month, Missing time part imputed with 23:59:59, no imputation flag", { # nolint
  expected_output <- tibble::tribble(
    ~XXSTDTC, ~AENDTM,
    "2019-07-18T15:25:40", ymd_hms("2019-07-18T15:25:40"),
    "2019-07-18T15:25", ymd_hms("2019-07-18T15:25:59"),
    "2019-07-18T15", ymd_hms("2019-07-18T15:59:59"),
    "2019-07-18", ymd_hms("2019-07-18T23:59:59"),
    "2019-02", ymd_hms("2019-02-28T23:59:59"),
    "2019", ymd_hms("2019-12-31T23:59:59"),
    "2019---07", ymd_hms("2019-12-31T23:59:59")
  )

  actual_output <- derive_vars_dtm(date,
    new_vars_prefix = "AEN",
    dtc = XXSTDTC,
    date_imputation = "LAST",
    time_imputation = "LAST",
    flag_imputation = FALSE
  )

  expect_equal(
    expected_output,
    actual_output
  )
})

test_that("Partial date imputed to the mid day/month", {
  expected_output <- tibble::tribble(
    ~XXSTDTC, ~ASTDTM, ~ASTDTF, ~ASTTMF,
    "2019-07-18T15:25:40", ymd_hms("2019-07-18T15:25:40"), "", "",
    "2019-07-18T15:25", ymd_hms("2019-07-18T15:25:00"), "", "S",
    "2019-07-18T15", ymd_hms("2019-07-18T15:00:00"), "", "M",
    "2019-07-18", ymd_hms("2019-07-18T00:00:00"), "", "H",
    "2019-02", ymd_hms("2019-02-15T00:00:00"), "D", "H",
    "2019", ymd_hms("2019-06-15T00:00:00"), "M", "H",
    "2019---07", ymd_hms("2019-06-15T00:00:00"), "M", "H"
  )

  actual_output <- derive_vars_dtm(date,
    new_vars_prefix = "AST",
    dtc = XXSTDTC,
    date_imputation = "MID"
  )
  actual_output1 <- derive_vars_dtm(date,
    new_vars_prefix = "AST",
    dtc = XXSTDTC,
    date_imputation = "15-06"
  )

  expect_equal(
    expected_output,
    actual_output
  )
  expect_equal(
    expected_output,
    actual_output1
  )
})
