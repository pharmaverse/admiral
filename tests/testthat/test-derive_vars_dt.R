# impute_dtc_dt ----
input <- c(
  "2019-07-18",
  "2019-02",
  "2019",
  "2019---07"
)

## Test 1: default: no date imputation ----
test_that("impute_dtc_dt Test 1: default: no date imputation", {
  expected_output <- c(
    "2019-07-18",
    NA_character_,
    NA_character_,
    NA_character_
  )
  expect_equal(impute_dtc_dt(dtc = input), expected_output)
})

## Test 2: impute month and day to first ----
test_that("impute_dtc_dt Test 2: impute month and day to first", {
  expected_output <- c(
    "2019-07-18",
    "2019-02-01",
    "2019-01-01",
    "2019-01-01"
  )

  expect_equal(
    impute_dtc_dt(
      dtc = input,
      highest_imputation = "M",
      date_imputation = "first"
    ),
    expected_output
  )

  expect_equal(
    impute_dtc_dt(
      dtc = input,
      highest_imputation = "M",
      date_imputation = "01-01"
    ),
    expected_output
  )
})

## Test 3: impute day to last ----
test_that("impute_dtc_dt Test 3: impute day to last", {
  expected_output <- c(
    "2019-07-18",
    "2019-02-28",
    NA_character_,
    NA_character_
  )
  expect_equal(
    impute_dtc_dt(
      dtc = input,
      highest_imputation = "D",
      date_imputation = "LAST",
      preserve = FALSE
    ),
    expected_output
  )
})

## Test 4: impute month and day to last and preserve = TRUE ----
test_that("impute_dtc_dt Test 4: impute month and day to last and preserve = TRUE", {
  expected_output <- c(
    "2019-07-18",
    "2019-02-28",
    "2019-12-31",
    "2019-12-07"
  )
  expect_equal(
    imputes <- impute_dtc_dt(
      dtc = input,
      highest_imputation = "M",
      date_imputation = "LAST",
      preserve = TRUE
    ),
    expected_output
  )
})


## Test 5: impute month and day to mid ----
test_that("impute_dtc_dt Test 5: impute month and day to mid", {
  expected_output <- c(
    "2019-07-18",
    "2019-02-15",
    "2019-06-30",
    "2019-06-30"
  )
  expect_equal(
    imputes <- impute_dtc_dt(
      dtc = input,
      highest_imputation = "M",
      date_imputation = "mid"
    ),
    expected_output
  )
})

## Test 6: min_dates parameter works ----
test_that("impute_dtc_dt Test 6: min_dates parameter works", {
  expect_equal(
    impute_dtc_dt(
      c("2020-12", "2020-11", NA_character_),
      min_dates = list(
        c(
          ymd("2020-12-06"),
          NA,
          NA
        ),
        c(
          ymd("2020-11-11"),
          ymd("2020-11-11"),
          ymd("2020-11-11")
        )
      ),
      highest_imputation = "Y",
      date_imputation = "first"
    ),
    c("2020-12-06", "2020-11-11", "2020-11-11")
  )
})

## Test 7: max_dates parameter works ----
test_that("impute_dtc_dt Test 7: max_dates parameter works", {
  expect_equal(
    impute_dtc_dt(c("2020-12", "2020-11", NA_character_),
      max_dates = list(
        c(ymd("2020-12-06"), NA, ymd("2020-09-13")),
        c(ymd(""), ymd("2020-11-11"), ymd(""))
      ),
      highest_imputation = "Y",
      date_imputation = "last"
    ),
    c("2020-12-06", "2020-11-11", "2020-09-13")
  )
})


## Test 8: min_dates length mismatch provides error ----
test_that("impute_dtc_dt Test 8: min_dates length mismatch provides error", {
  expect_snapshot(
    impute_dtc_dt(
      c("2020-12", NA_character_),
      min_dates = list(
        c(ymd("2020-12-06")),
        c(ymd("2020-11-11"))
      ),
      highest_imputation = "Y"
    ),
    error = TRUE
  )
})

## Test 9: max_dates length mismatch provides error ----
test_that("impute_dtc_dt Test 9: max_dates length mismatch provides error", {
  expect_snapshot(
    impute_dtc_dt(
      c("2020-12", NA_character_),
      max_dates = list(
        c(ymd("2020-12-06")),
        c(ymd("2020-11-11"))
      ),
      highest_imputation = "Y"
    ),
    error = TRUE
  )
})

## Test 10: Warning if null min/max_dates when highest_imputation = Y ----
test_that("impute_dtc_dt Test 10: Warning if null min/max_dates when highest_imputation = Y", {
  expect_warning(
    impute_dtc_dt(
      c("2020-12", NA_character_),
      highest_imputation = "Y"
    ),
    "If `highest_impuation` = \"Y\" is specified, `min_dates` or `max_dates` should be specified respectively." # nolint
  )
})


# convert_dtc_to_dt ----
inputdtc <- c(
  "2019-07-18T15:25:52",
  "2019-07-18"
)

## Test 11: Convert a complete -- DTC into a date object ----
test_that("convert_dtc_to_dt Test 11: Convert a complete -- DTC into a date object", {
  expected_output <- c(
    as.Date("2019-07-18"),
    as.Date("2019-07-18")
  )
  expect_equal(
    convert_dtc_to_dt(dtc = inputdtc),
    expected_output
  )
})

# compute_dtf ----

inputdtc <- c(
  "2019-07-18",
  "2019-02",
  "2019",
  "2019---07",
  "2019---06T00:00",
  "2019----T00:00",
  "2019-06--T00:00",
  "--06-06T00:00",
  "-----T00:00"
)
inputdt <- c(
  as.Date("2019-07-18"),
  as.Date("2019-02-01"),
  as.Date("2019-01-01"),
  as.Date("2019-01-01"),
  as.Date("2019-06-06"),
  as.Date("2019-06-06"),
  as.Date("2019-06-06"),
  as.Date("2019-06-06"),
  as.Date("2019-06-06")
)

## Test 12: compute DTF ----
test_that("compute_dtf Test 12: compute DTF", {
  expected_output <- c(
    NA_character_,
    "D",
    "M",
    "M",
    "M",
    "M",
    "D",
    "Y",
    "Y"
  )
  expect_equal(
    compute_dtf(
      dtc = inputdtc,
      dt = inputdt
    ),
    expected_output
  )
})

# derive_vars_dt ----

date <- tibble::tribble(
  ~XXSTDTC,
  "2019-07-18T15:25:40",
  "2019-07-18",
  "2019-02",
  "2019",
  "2019---07"
)

## Test 13: default behavior ----
test_that("derive_vars_dt Test 13: default behavior", {
  expected_output <- tibble::tribble(
    ~XXSTDTC,              ~ASTDT,
    "2019-07-18T15:25:40", as.Date("2019-07-18"),
    "2019-07-18",          as.Date("2019-07-18"),
    "2019-02",             as.Date(NA),
    "2019",                as.Date(NA),
    "2019---07",           as.Date(NA)
  )

  actual_output <- derive_vars_dt(
    date,
    new_vars_prefix = "AST",
    dtc = XXSTDTC
  )

  expect_dfs_equal(
    expected_output,
    actual_output,
    "XXSTDTC"
  )
})

## Test 14: no date imputation, add DTF ----
test_that("derive_vars_dt Test 14: no date imputation, add DTF", {
  expected_output <- tibble::tribble(
    ~XXSTDTC,              ~ASTDT,                ~ASTDTF,
    "2019-07-18T15:25:40", as.Date("2019-07-18"), NA_character_,
    "2019-07-18",          as.Date("2019-07-18"), NA_character_,
    "2019-02",             as.Date(NA),           NA_character_,
    "2019",                as.Date(NA),           NA_character_,
    "2019---07",           as.Date(NA),           NA_character_
  )

  actual_output <- derive_vars_dt(
    date,
    new_vars_prefix = "AST",
    flag_imputation = "date",
    dtc = XXSTDTC
  )

  expect_dfs_equal(
    expected_output,
    actual_output,
    "XXSTDTC"
  )
})

## Test 15: date imputed to first, auto DTF ----
test_that("derive_vars_dt Test 15: date imputed to first, auto DTF", {
  expected_output <- tibble::tribble(
    ~XXSTDTC,              ~ASTDT,                ~ASTDTF,
    "2019-07-18T15:25:40", as.Date("2019-07-18"), NA_character_,
    "2019-07-18",          as.Date("2019-07-18"), NA_character_,
    "2019-02",             as.Date("2019-02-01"), "D",
    "2019",                as.Date("2019-01-01"), "M",
    "2019---07",           as.Date("2019-01-01"), "M"
  )

  actual_output <- derive_vars_dt(
    date,
    new_vars_prefix = "AST",
    dtc = XXSTDTC,
    highest_imputation = "M",
    date_imputation = "first"
  )

  expect_dfs_equal(
    base = expected_output,
    compare = actual_output,
    keys = "XXSTDTC"
  )
})

## Test 16: date imputed to last, no DTF ----
test_that("derive_vars_dt Test 16: date imputed to last, no DTF", {
  expected_output <- tibble::tribble(
    ~XXSTDTC,              ~AENDT,
    "2019-07-18T15:25:40", as.Date("2019-07-18"),
    "2019-07-18",          as.Date("2019-07-18"),
    "2019-02",             as.Date("2019-02-28"),
    "2019",                as.Date("2019-12-31"),
    "2019---07",           as.Date("2019-12-31")
  )

  actual_output <- derive_vars_dt(
    date,
    new_vars_prefix = "AEN",
    dtc = XXSTDTC,
    highest_imputation = "M",
    date_imputation = "last",
    flag_imputation = "none"
  )

  expect_dfs_equal(
    base = expected_output,
    compare = actual_output,
    keys = "XXSTDTC"
  )
})

## Test 17: NA imputation for highest_imputation = Y & max_dates ----
test_that("derive_vars_dt Test 17: NA imputation for highest_imputation = Y & max_dates", {
  actual <- data.frame(
    AESTDTC = c(NA_character_, NA_character_),
    TRTSDT = c(ymd("2022-01-01"), NA)
  ) %>%
    mutate(AESTDTC = as.character(AESTDTC)) %>%
    derive_vars_dt(
      dtc = AESTDTC,
      new_vars_prefix = "AST",
      highest_imputation = "Y",
      date_imputation = "last",
      flag_imputation = "auto",
      max_dates = exprs(TRTSDT)
    )

  expected <- data.frame(
    AESTDTC = c(NA_character_, NA_character_),
    TRTSDT = c(ymd("2022-01-01"), NA),
    ASTDT = c(ymd("2022-01-01"), NA),
    ASTDTF = c("Y", NA)
  )

  expect_dfs_equal(actual, expected, keys = c("ASTDT", "ASTDTF"))
})

## Test 18: NA imputation for highest_imputation = Y & max_dates but date_imputation = first ----
test_that("derive_vars_dt Test 18: NA imputation for highest_imputation = Y & max_dates but date_imputation = first", { # nolint
  expect_snapshot(
    data.frame(
      AESTDTC = c(NA_character_, NA_character_),
      TRTSDT = c(ymd("2022-01-01"), NA)
    ) %>%
      mutate(AESTDTC = as.character(AESTDTC)) %>%
      derive_vars_dt(
        dtc = AESTDTC,
        new_vars_prefix = "AST",
        highest_imputation = "Y",
        date_imputation = "first",
        flag_imputation = "auto",
        max_dates = exprs(TRTSDT)
      )
  )
})

## Test 19: NA imputation for highest_imputation = Y & min_dates ----
test_that("derive_vars_dt Test 19: NA imputation for highest_imputation = Y & min_dates", {
  actual <- data.frame(
    AESTDTC = c(NA_character_, NA_character_),
    TRTSDT = c(ymd("2022-01-01"), NA)
  ) %>%
    mutate(AESTDTC = as.character(AESTDTC)) %>%
    derive_vars_dt(
      dtc = AESTDTC,
      new_vars_prefix = "AST",
      highest_imputation = "Y",
      date_imputation = "first",
      flag_imputation = "auto",
      min_dates = exprs(TRTSDT)
    )

  expected <- data.frame(
    AESTDTC = c(NA_character_, NA_character_),
    TRTSDT = c(ymd("2022-01-01"), NA),
    ASTDT = c(ymd("2022-01-01"), NA),
    ASTDTF = c("Y", NA)
  )

  expect_dfs_equal(actual, expected, keys = c("ASTDT", "ASTDTF"))
})

## Test 20: NA imputation for highest_imputation = Y & min_dates but date_imputation = last ----
test_that("derive_vars_dt Test 20: NA imputation for highest_imputation = Y & min_dates but date_imputation = last", { # nolint
  expect_snapshot(
    data.frame(
      AESTDTC = c(NA_character_, NA_character_),
      TRTSDT = c(ymd("2022-01-01"), NA)
    ) %>%
      mutate(AESTDTC = as.character(AESTDTC)) %>%
      derive_vars_dt(
        dtc = AESTDTC,
        new_vars_prefix = "AST",
        highest_imputation = "Y",
        date_imputation = "last",
        flag_imputation = "auto",
        min_dates = exprs(TRTSDT)
      )
  )
})

## Test 21: NA imputation for highest_imputation = Y but null min/max dates fails ----
test_that("derive_vars_dt Test 21: NA imputation for highest_imputation = Y but null min/max dates fails", { # nolint
  expect_snapshot(
    data.frame(
      AESTDTC = c(NA_character_, NA_character_),
      TRTSDT = c(ymd("2022-01-01"), NA)
    ) %>%
      mutate(AESTDTC = as.character(AESTDTC)) %>%
      derive_vars_dt(
        dtc = AESTDTC,
        new_vars_prefix = "AST",
        highest_imputation = "Y",
        date_imputation = "first",
        flag_imputation = "auto"
      ),
    error = TRUE
  )
})

## Test 22: Supplying both min/max dates for highest_imputation = Y works ----
test_that("derive_vars_dt Test 22: Supplying both min/max dates for highest_imputation = Y works", { # nolint
  actual <- data.frame(
    AESTDTC = c(NA_character_, NA_character_),
    TRTSDT = c(ymd("2022-01-01"), NA),
    TRTEDT = c(ymd("2022-01-31"), NA)
  ) %>%
    mutate(AESTDTC = as.character(AESTDTC)) %>%
    derive_vars_dt(
      dtc = AESTDTC,
      new_vars_prefix = "AST",
      highest_imputation = "Y",
      min_dates = exprs(TRTSDT),
      max_dates = exprs(TRTEDT)
    )

  expected <- data.frame(
    AESTDTC = c(NA_character_, NA_character_),
    TRTSDT = c(ymd("2022-01-01"), NA),
    TRTEDT = c(ymd("2022-01-31"), NA),
    ASTDT = c(ymd("2022-01-01"), NA),
    ASTDTF = c("Y", NA)
  )

  expect_dfs_equal(actual, expected, keys = c("ASTDT", "ASTDTF"))
})

## Test 23: Supplying both min/max dates for highest_imputation = Y works ----
test_that("derive_vars_dt Test 23: Supplying both min/max dates for highest_imputation = Y works", { # nolint
  actual <- data.frame(
    AESTDTC = c(NA_character_, NA_character_),
    TRTSDT = c(ymd("2022-01-01"), NA),
    TRTEDT = c(ymd("2022-01-31"), NA)
  ) %>%
    mutate(AESTDTC = as.character(AESTDTC)) %>%
    derive_vars_dt(
      dtc = AESTDTC,
      new_vars_prefix = "AST",
      highest_imputation = "Y",
      date_imputation = "last",
      min_dates = exprs(TRTSDT),
      max_dates = exprs(TRTEDT)
    )

  expected <- data.frame(
    AESTDTC = c(NA_character_, NA_character_),
    TRTSDT = c(ymd("2022-01-01"), NA),
    TRTEDT = c(ymd("2022-01-31"), NA),
    ASTDT = c(ymd("2022-01-31"), NA),
    ASTDTF = c("Y", NA)
  )

  expect_dfs_equal(actual, expected, keys = c("ASTDT", "ASTDTF"))
})
