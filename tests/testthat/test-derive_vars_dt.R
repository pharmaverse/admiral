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
    flag_imputation = TRUE,
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
    flag_imputation = TRUE,
    date_imputation = "FIRST"
  )
  actual_output1 <- derive_vars_dt(date,
    new_vars_prefix = "AST",
    dtc = XXSTDTC,
    flag_imputation = TRUE,
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
    flag_imputation = TRUE,
    date_imputation = "LAST"
  )

  expect_equal(
    expected_output,
    actual_output
  )
})


test_that("Partial date imputed to the LAST day/month", {
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
    flag_imputation = TRUE,
    date_imputation = "LAST"
  )

  expect_equal(
    expected_output,
    actual_output
  )
})



test_that("Partial date imputation as MID to the mid day/month", {
  expected_output <- tibble::tribble(
    ~XXSTDTC, ~ASTDT, ~ASTDTF,
    "2019-07-18T15:25:40", as.Date("2019-07-18"), NA_character_,
    "2019-07-18", as.Date("2019-07-18"), NA_character_,
    "2019-02", as.Date("2019-06-30"), "D",
    "2019", as.Date("2019-06-30"), "M",
    "2019---07", as.Date("2019-06-30"), "M"
  )

  actual_output <- derive_vars_dt(
    date,
    new_vars_prefix = "AST",
    dtc = XXSTDTC,
    flag_imputation = TRUE,
    date_imputation = "MID"
  )

  expect_dfs_equal(
    expected_output,
    actual_output,
    keys = c("XXSTDTC")
  )
})

test_that("Partial date imputation as 6-15 to the mid day/month", {
    expected_output <- tibble::tribble(
      ~XXSTDTC, ~ASTDT, ~ASTDTF,
      "2019-07-18T15:25:40", as.Date("2019-07-18"), NA_character_,
      "2019-07-18", as.Date("2019-07-18"), NA_character_,
      "2019-02", as.Date("2019-02-15"), "D",
      "2019", as.Date("2019-06-15"), "M",
      "2019---07", as.Date("2019-06-15"), "M"
    )

  actual_output <- derive_vars_dt(
    date,
    new_vars_prefix = "AST",
    dtc = XXSTDTC,
    flag_imputation = TRUE,
    date_imputation = "06-15"
  )

  expect_dfs_equal(
    expected_output,
    actual_output,
    keys = c("XXSTDTC")
  )
})



date <- tibble::tribble(
  ~XXSTDTC,
  "2019-07-18T15:25:40",
  "2019-07-18",
  "2019-02",
  "2019",
  "2019---07"
)

test_that("Partial date imputation as MID and preserve = TRUE to the mid day/month", { # nolint
  expected_output <- tibble::tribble(
    ~XXSTDTC, ~ASTDT, ~ASTDTF,
    "2019-07-18T15:25:40", as.Date("2019-07-18"), NA_character_,
    "2019-07-18", as.Date("2019-07-18"), NA_character_,
    "2019-02", as.Date("2019-02-15"), "D",
    "2019", as.Date("2019-06-30"), "M",
    "2019---07", as.Date("2019-06-07"), "M"
  )

  actual_output <- derive_vars_dt(
    date,
    new_vars_prefix = "AST",
    dtc = XXSTDTC,
    flag_imputation = TRUE,
    date_imputation = "MID",
    preserve = TRUE
  )

  expect_dfs_equal(
    expected_output,
    actual_output,
    keys = c("XXSTDTC")
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
    date_imputation = "LAST"
  )

  expect_dfs_equal(
    base = expected_output,
    comp = actual_output,
    keys = "XXSTDTC"
  )
})

test_that("Partial date imputed to the last day/month, no DTF and preserve=TRUE", {
  expected_output <- tibble::tribble(
    ~XXSTDTC, ~ASTDT,
    "2019-07-18T15:25:40", as.Date("2019-07-18"),
    "2019-07-18", as.Date("2019-07-18"),
    "2019-02", as.Date("2019-02-28"),
    "2019", as.Date("2019-12-31"),
    "2019---07", as.Date("2019-12-07")
  )

  actual_output <- derive_vars_dt(date,
                                  new_vars_prefix = "AST",
                                  dtc = XXSTDTC,
                                  date_imputation = "LAST",
                                  preserve = TRUE
  )

  expect_dfs_equal(
    base = expected_output,
    comp = actual_output,
    keys = "XXSTDTC"
  )
})

test_that("Partial date imputed to the first day/month, no DTF and preserve=TRUE", {
  expected_output <- tibble::tribble(
    ~XXSTDTC, ~ASTDT,
    "2019-07-18T15:25:40", as.Date("2019-07-18"),
    "2019-07-18", as.Date("2019-07-18"),
    "2019-02", as.Date("2019-02-01"),
    "2019", as.Date("2019-01-01"),
    "2019---07", as.Date("2019-01-07")
  )

  actual_output <- derive_vars_dt(date,
                                  new_vars_prefix = "AST",
                                  dtc = XXSTDTC,
                                  date_imputation = "FIRST",
                                  preserve = TRUE
  )

  expect_dfs_equal(
    base = expected_output,
    comp = actual_output,
    keys = "XXSTDTC"
  )
})
