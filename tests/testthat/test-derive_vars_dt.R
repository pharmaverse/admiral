context("test-derive_vars_dt")

date <- tibble::tribble(
  ~XXSTDTC,
  "2019-07-18T15:25:40",
  "2019-07-18",
  "2019-02",
  "2019",
  "2019---07"
)

test_that("default: no date imputation, time part set o 00:00:00, add DTF", {
  expected_output <- tibble::tribble(
    ~XXSTDTC, ~ASTDT, ~ASTDTF,
    "2019-07-18T15:25:40", as.Date("2019-07-18"), NA_character_,
    "2019-07-18", as.Date("2019-07-18"), NA_character_,
    "2019-02", as.Date(NA), NA_character_,
    "2019", as.Date(NA), NA_character_,
    "2019---07", as.Date(NA), NA_character_
  )

  actual_output <- derive_vars_dt(date,
    new_vars_prefix = "AST",
    dtc = XXSTDTC
  )

  expect_dfs_equal(
    expected_output,
    actual_output,
    "XXSTDTC"
  )
})

test_that("Partial date imputed to the first day/month", {
  expected_output <- tibble::tribble(
    ~XXSTDTC, ~ASTDT, ~ASTDTF,
    "2019-07-18T15:25:40", as.Date("2019-07-18"), NA_character_,
    "2019-07-18", as.Date("2019-07-18"), NA_character_,
    "2019-02", as.Date("2019-02-01"), "D",
    "2019", as.Date("2019-01-01"), "M",
    "2019---07", as.Date("2019-01-01"), "M"
  )

  actual_output <- derive_vars_dt(date,
    new_vars_prefix = "AST",
    dtc = XXSTDTC,
    date_imputation = "FIRST"
  )
  actual_output1 <- derive_vars_dt(date,
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

test_that("Partial date imputed to the last day/month", {
  expected_output <- tibble::tribble(
    ~XXSTDTC, ~AENDT, ~AENDTF,
    "2019-07-18T15:25:40", as.Date("2019-07-18"), NA_character_,
    "2019-07-18", as.Date("2019-07-18"), NA_character_,
    "2019-02", as.Date("2019-02-28"), "D",
    "2019", as.Date("2019-12-31"), "M",
    "2019---07", as.Date("2019-12-31"), "M"
  )

  actual_output <- derive_vars_dt(date,
    new_vars_prefix = "AEN",
    dtc = XXSTDTC,
    date_imputation = "LAST"
  )

  expect_equal(
    expected_output,
    actual_output
  )
})

test_that("Partial date imputed to the mid day/month", {
  expected_output <- tibble::tribble(
    ~XXSTDTC, ~ASTDT, ~ASTDTF,
    "2019-07-18T15:25:40", as.Date("2019-07-18"), NA_character_,
    "2019-07-18", as.Date("2019-07-18"), NA_character_,
    "2019-02", as.Date("2019-02-15"), "D",
    "2019", as.Date("2019-06-15"), "M",
    "2019---07", as.Date("2019-06-15"), "M"
  )

  actual_output <- derive_vars_dt(date,
    new_vars_prefix = "AST",
    dtc = XXSTDTC,
    date_imputation = "MID"
  )
  actual_output1 <- derive_vars_dt(date,
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


test_that("Partial date imputed to the last day/month", {
  expected_output <- tibble::tribble(
    ~XXSTDTC, ~ASTDT, ~ASTDTF,
    "2019-07-18T15:25:40", as.Date("2019-07-18"), NA_character_,
    "2019-07-18", as.Date("2019-07-18"), NA_character_,
    "2019-02", as.Date("2019-02-28"), "D",
    "2019", as.Date("2019-12-31"), "M",
    "2019---07", as.Date("2019-12-31"), "M"
  )

  actual_output <- derive_vars_dt(date,
    new_vars_prefix = "AST",
    dtc = XXSTDTC,
    date_imputation = "LAST"
  )

  expect_equal(
    expected_output,
    actual_output
  )
})

test_that("Partial date imputed to the mid day/month", {
  expected_output <- tibble::tribble(
    ~XXSTDTC, ~ASTDT, ~ASTDTF,
    "2019-07-18T15:25:40", as.Date("2019-07-18"), NA_character_,
    "2019-07-18", as.Date("2019-07-18"), NA_character_,
    "2019-02", as.Date("2019-02-15"), "D",
    "2019", as.Date("2019-06-15"), "M",
    "2019---07", as.Date("2019-06-15"), "M"
  )

  actual_output <- derive_vars_dt(date,
    new_vars_prefix = "AST",
    dtc = XXSTDTC,
    date_imputation = "MID"
  )
  actual_output1 <- derive_vars_dt(date,
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


test_that("Partial date imputed to the last day/month, no DTF", {
  expected_output <- tibble::tribble(
    ~XXSTDTC, ~ASTDT,
    "2019-07-18T15:25:40", as.Date("2019-07-18"),
    "2019-07-18", as.Date("2019-07-18"),
    "2019-02", as.Date("2019-02-28"),
    "2019", as.Date("2019-12-31"),
    "2019---07", as.Date("2019-12-31")
  )

  actual_output <- derive_vars_dt(date,
    new_vars_prefix = "AST",
    dtc = XXSTDTC,
    date_imputation = "LAST",
    flag_imputation = FALSE
  )

  expect_equal(
    expected_output,
    actual_output
  )
})
