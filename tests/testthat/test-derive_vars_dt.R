input <- c(
  "2019-07-18", # full date
  "--07-18", # missing year
  "2019", # missing month and day
  "2019-07--", # missing day
  "2019---07" # missing just month
)

input_warnings <- c(
  "", # empty string
  NA_character_, # NA
  "2019/07/18" # inappropriate date format/string
)

## Test 1: default: no date imputation ----
test_that("derive_vars_dt Test 1: default: no date imputation", {
  expected_output <- c(
    "2019-07-18",
    NA_character_,
    NA_character_,
    NA_character_,
    NA_character_
  )
  actual_output <- impute_dtc_dt(dtc = input)

  expect_equal(actual_output, expected_output)
})

## Test 2: impute month and day to first ----
test_that("derive_vars_dt Test 2: impute month and day to first", {
  expected_output <- c(
    "2019-07-18",
    NA_character_,
    "2019-01-01",
    "2019-07-01",
    "2019-01-01"
  )

  actual_output <- impute_dtc_dt(
    dtc = input,
    highest_imputation = "M",
    date_imputation = "first"
  )

  expect_equal(
    actual_output,
    expected_output
  )

  actual_output <- impute_dtc_dt(
    dtc = input,
    highest_imputation = "M",
    date_imputation = "01-01"
  )

  expect_equal(
    actual_output,
    expected_output
  )
})

## Test 3: impute day to last ----
test_that("derive_vars_dt Test 3: impute day to last", {
  expected_output <- c(
    "2019-07-18",
    NA_character_,
    NA_character_,
    "2019-07-31",
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
test_that("derive_vars_dt Test 4: impute month and day to last and preserve = TRUE", {
  expected_output <- c(
    "2019-07-18",
    NA_character_,
    "2019-12-31",
    "2019-07-31",
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
test_that("derive_vars_dt Test 5: impute month and day to mid", {
  expected_output <- c(
    "2019-07-18",
    NA_character_,
    "2019-06-30",
    "2019-07-15",
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
test_that("derive_vars_dt Test 6: min_dates parameter works", {
  expect_equal(
    impute_dtc_dt(
      input,
      min_dates = list(
        c(
          ymd("2019-07-06"),
          ymd("2019-07-06"),
          ymd("2019-07-06"),
          ymd("2019-07-06"),
          ymd("2019-07-06")
        ),
        c(
          ymd("2019-06-06"),
          ymd("2019-06-06"),
          ymd("2019-06-06"),
          ymd("2019-06-06"),
          ymd("2019-06-06")
        )
      ),
      highest_imputation = "Y",
      date_imputation = "first"
    ),
    c("2019-07-18", "2019-07-06", "2019-07-06", "2019-07-06", "2019-07-06")
  )
})

## Test 7: max_dates parameter works ----
test_that("derive_vars_dt Test 7: max_dates parameter works", {
  expect_equal(
    impute_dtc_dt(
      input,
      max_dates = list(
        c(
          ymd("2019-07-06"),
          ymd("2019-07-06"),
          ymd("2019-07-06"),
          ymd("2019-07-06"),
          ymd("2019-07-06")
        ),
        c(
          ymd("2019-06-06"),
          ymd("2019-06-06"),
          ymd("2019-06-06"),
          ymd("2019-06-06"),
          ymd("2019-06-06")
        )
      ),
      highest_imputation = "Y",
      date_imputation = "last"
    ),
    c("2019-07-18", "2019-06-06", "2019-06-06", "2019-07-06", "2019-06-06")
  )
})


## Test 8: min_dates length mismatch provides error ----
test_that("derive_vars_dt Test 8: min_dates length mismatch provides error", {
  expect_snapshot(
    impute_dtc_dt(
      input,
      min_dates = list(
        c(ymd("2019-07-06")),
        c(ymd("2019-06-06"))
      ),
      highest_imputation = "Y",
      date_imputation = "first"
    ),
    error = TRUE
  )
})

## Test 9: max_dates length mismatch provides error ----
test_that("derive_vars_dt Test 9: max_dates length mismatch provides error", {
  expect_snapshot(
    impute_dtc_dt(
      input,
      max_dates = list(
        c(ymd("2019-07-06")),
        c(ymd("2019-06-06"))
      ),
      highest_imputation = "Y",
      date_imputation = "last"
    ),
    error = TRUE
  )
})

## Test 10: Warning if null min/max_dates when highest_imputation = Y ----
test_that("derive_vars_dt Test 10: Warning if null min/max_dates when highest_imputation = Y", {
  expect_warning(
    impute_dtc_dt(
      input,
      highest_imputation = "Y"
    ),
    "If `highest_impuation` = \"Y\" is specified, `min_dates` or `max_dates` should be specified respectively." # nolint
  )
})

## Test 11: appropriate warnings/return object for impute_dtc_dt ----
test_that("derive_vars_dt Test 11: appropriate warnings/return object for impute_dtc_dt", {
  expect_warning(
    impute_dtc_dt(dtc = input_warnings),
    regexp = "incorrect datetime format"
  )

  expect_equal(
    suppressWarnings(impute_dtc_dt(dtc = input_warnings)),
    rep(NA_character_, 3)
  )
})


# convert_dtc_to_dt ----
## Test 12: Convert a complete -- DTC into a date object ----
test_that("convert_dtc_to_dt Test 12: Convert a complete -- DTC into a date object", {
  expected_output <- c(
    as.Date("2019-07-18")
  )
  expect_equal(
    convert_dtc_to_dt(dtc = input[[1]]),
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

## Test 13: compute DTF ----
test_that("compute_dtf Test 13: compute DTF", {
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

# restrict_imputed_dtc_dt ----
## Test 14: restrict_imputed_dtc_dt works as expected ----
test_that("restrict_imputed_dtc_dt Test 14: restrict_imputed_dtc_dt works as expected", {
  imputed_dtc <- impute_dtc_dt(
    input,
    min_dates = list(
      c(
        ymd("2019-07-06"),
        ymd("2019-07-06"),
        ymd("2019-07-06"),
        ymd("2019-07-06"),
        ymd("2019-07-06")
      ),
      c(
        ymd("2019-06-06"),
        ymd("2019-06-06"),
        ymd("2019-06-06"),
        ymd("2019-06-06"),
        ymd("2019-06-06")
      )
    ),
    highest_imputation = "Y",
    date_imputation = "first"
  )
  restricted <- restrict_imputed_dtc_dt(
    input,
    imputed_dtc = imputed_dtc,
    min_dates = list(
      c(
        ymd("2019-07-06"),
        ymd("2019-07-06"),
        ymd("2019-07-06"),
        ymd("2019-07-06"),
        ymd("2019-07-06")
      ),
      c(
        ymd("2019-06-06"),
        ymd("2019-06-06"),
        ymd("2019-06-06"),
        ymd("2019-06-06"),
        ymd("2019-06-06")
      )
    ),
    max_dates = list(
      c(
        ymd("2019-07-06"),
        ymd("2019-07-06"),
        ymd("2019-07-06"),
        ymd("2019-07-06"),
        ymd("2019-07-06")
      ),
      c(
        ymd("2019-06-06"),
        ymd("2019-06-06"),
        ymd("2019-06-06"),
        ymd("2019-06-06"),
        ymd("2019-06-06")
      )
    )
  )
  expect_equal(
    restricted,
    c("2019-07-18", "2019-06-06", "2019-06-06", "2019-07-06", "2019-06-06")
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

## Test 15: default behavior ----
test_that("derive_vars_dt Test 15: default behavior", {
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

## Test 16: no date imputation, add DTF ----
test_that("derive_vars_dt Test 16: no date imputation, add DTF", {
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

## Test 17: date imputed to first, auto DTF ----
test_that("derive_vars_dt Test 17: date imputed to first, auto DTF", {
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

## Test 18: date imputed to last, no DTF ----
test_that("derive_vars_dt Test 18: date imputed to last, no DTF", {
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

## Test 19: NA imputation for highest_imputation = Y & max_dates ----
test_that("derive_vars_dt Test 19: NA imputation for highest_imputation = Y & max_dates", {
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

## Test 20: NA imputation for highest_imputation = Y & max_dates but date_imputation = first ----
test_that("derive_vars_dt Test 20: NA imputation for highest_imputation = Y & max_dates but date_imputation = first", { # nolint
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

## Test 21: NA imputation for highest_imputation = Y & min_dates ----
test_that("derive_vars_dt Test 21: NA imputation for highest_imputation = Y & min_dates", {
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

## Test 22: NA imputation for highest_imputation = Y & min_dates but date_imputation = last ----
test_that("derive_vars_dt Test 22: NA imputation for highest_imputation = Y & min_dates but date_imputation = last", { # nolint
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

## Test 23: NA imputation for highest_imputation = Y but null min/max dates fails ----
test_that("derive_vars_dt Test 23: NA imputation for highest_imputation = Y but null min/max dates fails", { # nolint
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

## Test 24: Supplying both min/max dates for highest_imputation = Y works ----
test_that("derive_vars_dt Test 24: Supplying both min/max dates for highest_imputation = Y works", { # nolint
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

## Test 25: Supplying both min/max dates for highest_imputation = Y works ----
test_that("derive_vars_dt Test 25: Supplying both min/max dates for highest_imputation = Y works", { # nolint
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


## Test 26: no date imputation, DTF present ----
test_that("derive_vars_dt Test 26: no date imputation, DTF present", {
  expected_output <- tibble::tribble(
    ~XXSTDTC,              ~ASTDT,                ~ASTDTF,
    "2019-07-18T15:25:40", as.Date("2019-07-18"), NA_character_,
    "2019-07-18",          as.Date("2019-07-18"), NA_character_,
    "2019-02",             as.Date(NA),           NA_character_,
    "2019",                as.Date(NA),           NA_character_,
    "2019---07",           as.Date(NA),           NA_character_
  )
  date <- select(expected_output, -ASTDT)
  expect_message(
    actual_output <- derive_vars_dt(
      date,
      new_vars_prefix = "AST",
      flag_imputation = "date",
      dtc = XXSTDTC
    ),
    regex =
      "The ASTDTF variable is already present in the input dataset and will not be re-derived."
  )

  expect_dfs_equal(
    expected_output,
    actual_output,
    "XXSTDTC"
  )
})

rm(input, input_warnings, inputdt, inputdtc, date)
