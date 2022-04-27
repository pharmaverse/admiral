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


test_that("default: no date imputation, time part set to 00:00:00, add DTF, TMF", {
  expected_output <- tibble::tribble(
    ~XXSTDTC, ~ASTDTM, ~ASTTMF,
    "2019-07-18T15:25:40", ymd_hms("2019-07-18T15:25:40"),  NA_character_,
    "2019-07-18T15:25", ymd_hms("2019-07-18T15:25:00"), "S",
    "2019-07-18T15", ymd_hms("2019-07-18T15:00:00"), "M",
    "2019-07-18", ymd_hms("2019-07-18T00:00:00"), "H",
    "2019-02", ymd_hms(NA), NA_character_,
    "2019", ymd_hms(NA), NA_character_,
    "2019---07", ymd_hms(NA), NA_character_
  )

  actual_output <- derive_vars_dtm(
    input,
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
    "2019-07-18T15:25:40", ymd_hms("2019-07-18T15:25:40"), NA_character_, NA_character_,
    "2019-07-18T15:25", ymd_hms("2019-07-18T15:25:00"), NA_character_, "S",
    "2019-07-18T15", ymd_hms("2019-07-18T15:00:00"), NA_character_, "M",
    "2019-07-18", ymd_hms("2019-07-18T00:00:00"), NA_character_, "H",
    "2019-02", ymd_hms("2019-02-01T00:00:00"), "D", "H",
    "2019", ymd_hms("2019-01-01T00:00:00"), "M", "H",
    "2019---07", ymd_hms("2019-01-01T00:00:00"), "M", "H"
  )

  actual_output <- derive_vars_dtm(
    input,
    new_vars_prefix = "AST",
    dtc = XXSTDTC,
    date_imputation = "FIRST"
  )
  actual_output1 <- derive_vars_dtm(
    input,
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
    "2019-07-18T15:25:40", ymd_hms("2019-07-18T15:25:40"), NA_character_, NA_character_,
    "2019-07-18T15:25", ymd_hms("2019-07-18T15:25:59"), NA_character_, "S",
    "2019-07-18T15", ymd_hms("2019-07-18T15:59:59"), NA_character_, "M",
    "2019-07-18", ymd_hms("2019-07-18T23:59:59"), NA_character_, "H",
    "2019-02", ymd_hms("2019-02-28T23:59:59"), "D", "H",
    "2019", ymd_hms("2019-12-31T23:59:59"), "M", "H",
    "2019---07", ymd_hms("2019-12-31T23:59:59"), "M", "H"
  )

  actual_output <- derive_vars_dtm(
    input,
    new_vars_prefix = "AEN",
    dtc = XXSTDTC,
    date_imputation = "LAST",
    time_imputation = "LAST"
  )

  actual_output1 <- derive_vars_dtm(
    input,
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

  actual_output <- derive_vars_dtm(
    input,
    new_vars_prefix = "AEN",
    dtc = XXSTDTC,
    date_imputation = "LAST",
    time_imputation = "LAST",
    flag_imputation = "None"
  )

  expect_equal(
    expected_output,
    actual_output
  )
})

test_that("Partial date imputed to the MID day/month", {
  expected_output <- tibble::tribble(
    ~XXSTDTC, ~ASTDTM, ~ASTDTF, ~ASTTMF,
    "2019-07-18T15:25:40", ymd_hms("2019-07-18T15:25:40"), NA_character_, NA_character_,
    "2019-07-18T15:25", ymd_hms("2019-07-18T15:25:00"), NA_character_, "S",
    "2019-07-18T15", ymd_hms("2019-07-18T15:00:00"), NA_character_, "M",
    "2019-07-18", ymd_hms("2019-07-18T00:00:00"), NA_character_, "H",
    "2019-02", ymd_hms("2019-02-15T00:00:00"), "D", "H",
    "2019", ymd_hms("2019-06-30T00:00:00"), "M", "H",
    "2019---07", ymd_hms("2019-06-15T00:00:00"), "M", "H"
  )

  actual_output <- derive_vars_dtm(
    input,
    new_vars_prefix = "AST",
    dtc = XXSTDTC,
    date_imputation = "MID"
  )

  expect_equal(
    expected_output,
    actual_output
  )

})

test_that("Partial date imputed to the 6-15 day/month", {
  expected_output <- tibble::tribble(
    ~XXSTDTC, ~ASTDTM, ~ASTDTF, ~ASTTMF,
    "2019-07-18T15:25:40", ymd_hms("2019-07-18T15:25:40"), NA_character_, NA_character_,
    "2019-07-18T15:25", ymd_hms("2019-07-18T15:25:00"), NA_character_, "S",
    "2019-07-18T15", ymd_hms("2019-07-18T15:00:00"), NA_character_, "M",
    "2019-07-18", ymd_hms("2019-07-18T00:00:00"), NA_character_, "H",
    "2019-02", ymd_hms("2019-02-15T00:00:00"), "D", "H",
    "2019", ymd_hms("2019-06-15T00:00:00"), "M", "H",
    "2019---07", ymd_hms("2019-06-15T00:00:00"), "M", "H"
  )

  actual_output <- derive_vars_dtm(
    input,
    new_vars_prefix = "AST",
    dtc = XXSTDTC,
    date_imputation = "06-15"
  )

  expect_equal(
    expected_output,
    actual_output
  )
})

test_that("No re-derivation is done if --DTF variable already exists", {

  expected_output <- tibble::tribble(
    ~XXSTDTC, ~ASTDTM, ~ASTDTF, ~ASTTMF,
    "2019-07-18T15:25:40", ymd_hms("2019-07-18T15:25:40"), NA_character_, NA_character_,
    "2019-07-18T15:25", ymd_hms("2019-07-18T15:25:00"), NA_character_, "S",
    "2019-07-18T15", ymd_hms("2019-07-18T15:00:00"), NA_character_, "M",
    "2019-07-18", ymd_hms("2019-07-18T00:00:00"), NA_character_, "H",
    "2019-02", ymd_hms("2019-02-01T00:00:00"), "D", "H",
    "2019", ymd_hms("2019-01-01T00:00:00"), "M", "H",
    "2019---07", ymd_hms("2019-01-01T00:00:00"), "M", "H"
  ) %>%
    select(XXSTDTC, ASTDTF, everything())

  actual_output <- expect_message(
    derive_vars_dtm(
      mutate(input, ASTDTF = c(NA, NA, NA, NA, "D", "M", "M")),
      new_vars_prefix = "AST",
      dtc = XXSTDTC,
      date_imputation = "FIRST"
    ),
    regexp = "^The .* variable is already present in the input dataset and will not be re-derived."
  )

  expect_equal(expected_output, actual_output)

})

input_maxed <- input %>%
  filter(!str_detect(XXSTDTC, "18")) %>%
  mutate(DCUTDT = ymd_hms("2019-02-10T00:00:00"))

test_that("max_dates parameter works as expected", {
  expected_output <- tibble::tribble(
    ~XXSTDTC, ~ASTDTM, ~ASTDTF, ~ASTTMF,
    "2019-02", ymd_hms("2019-02-10T00:00:00"), "D", "H",
    "2019", ymd_hms("2019-02-10T00:00:00"), "M", "H",
    "2019---07", ymd_hms("2019-02-10T00:00:00"), "M", "H"
  ) %>%
    mutate(DCUTDT = ymd_hms("2019-02-10T00:00:00"))

  actual_output <- derive_vars_dtm(
    input_maxed,
    new_vars_prefix = "AST",
    dtc = XXSTDTC,
    date_imputation = "LAST",
    max_dates = vars(DCUTDT)
  )

  expect_dfs_equal(expected_output, actual_output, keys = c("XXSTDTC"))

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

test_that("Ignore Seconds Flag is not used when not present in the function call", {

expected_output <- tibble::tribble(
  ~XXSTDTC, ~ASTDTM, ~ASTDTF, ~ASTTMF,
  "2019-07-18T15:25:40", ymd_hms("2019-07-18T15:25:40"), NA_character_, NA_character_,
  "2019-07-18T15:25", ymd_hms("2019-07-18T15:25:00"), NA_character_, "S",
  "2019-07-18T15", ymd_hms("2019-07-18T15:00:00"), NA_character_, "M",
  "2019-07-18", ymd_hms("2019-07-18T00:00:00"), NA_character_, "H",
  "2019-02", ymd_hms("2019-02-01T00:00:00"), "D", "H",
  "2019", ymd_hms("2019-01-01T00:00:00"), "M", "H",
  "2019---07", ymd_hms("2019-01-01T00:00:00"), "M", "H"
)

actual_output <- derive_vars_dtm(
  input_secs,
  new_vars_prefix = "AST",
  dtc = XXSTDTC,
  date_imputation = "FIRST",
  time_imputation = "FIRST"
)

expect_equal(expected_output, actual_output)
})

test_that("Ignore Seconds Flag is not used when set to FALSE in function call", {


  expected_output <- tibble::tribble(
    ~XXSTDTC, ~ASTDTM, ~ASTDTF, ~ASTTMF,
    "2019-07-18T15:25:40", ymd_hms("2019-07-18T15:25:40"), NA_character_, NA_character_,
    "2019-07-18T15:25", ymd_hms("2019-07-18T15:25:00"), NA_character_, "S",
    "2019-07-18T15", ymd_hms("2019-07-18T15:00:00"), NA_character_, "M",
    "2019-07-18", ymd_hms("2019-07-18T00:00:00"), NA_character_, "H",
    "2019-02", ymd_hms("2019-02-01T00:00:00"), "D", "H",
    "2019", ymd_hms("2019-01-01T00:00:00"), "M", "H",
    "2019---07", ymd_hms("2019-01-01T00:00:00"), "M", "H"
  )

  actual_output <- derive_vars_dtm(
    input_secs,
    new_vars_prefix = "AST",
    dtc = XXSTDTC,
    date_imputation = "FIRST",
    time_imputation = "FIRST",
    ignore_seconds_flag = FALSE
  )

  expect_equal(expected_output, actual_output)
})


input_no_s <- tibble::tribble(
  ~XXSTDTC,
  "2019-07-18T15:25",
  "2019-07-18T15:25",
  "2019-07-18T15",
  "2019-07-18",
  "2019-02",
  "2019",
  "2019---07"
)


test_that("Ignore Seconds Flag remove the Seconds Flag, S, from XXDTF variable when set to TRUE", { # nolint

  expected_output <- tibble::tribble(
    ~XXSTDTC, ~ASTDTM, ~ASTDTF, ~ASTTMF,
    "2019-07-18T15:25", ymd_hms("2019-07-18T15:25:00"), NA_character_, NA_character_,
    "2019-07-18T15:25", ymd_hms("2019-07-18T15:25:00"), NA_character_, NA_character_,
    "2019-07-18T15", ymd_hms("2019-07-18T15:00:00"), NA_character_, "M",
    "2019-07-18", ymd_hms("2019-07-18T00:00:00"), NA_character_, "H",
    "2019-02", ymd_hms("2019-02-01T00:00:00"), "D", "H",
    "2019", ymd_hms("2019-01-01T00:00:00"), "M", "H",
    "2019---07", ymd_hms("2019-01-01T00:00:00"), "M", "H"
  )

  actual_output <- derive_vars_dtm(
    input_no_s,
    new_vars_prefix = "AST",
    dtc = XXSTDTC,
    date_imputation = "FIRST",
    time_imputation = "FIRST",
    ignore_seconds_flag = TRUE
  )

  expect_equal(expected_output, actual_output)
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

test_that("Function throws ERROR when Ignore Seconds Flag is invoked and seconds is present in the data ", { # nolint


  expect_error(
    derive_vars_dtm(
      input_secs,
      new_vars_prefix = "AST",
      dtc = XXSTDTC,
      date_imputation = "FIRST",
      time_imputation = "FIRST",
      ignore_seconds_flag = TRUE
    ),
    regexp =  "Seconds detected in data while ignore_seconds_flag is invoked")

})

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

test_that("Partial date imputation as MID and preserve = TRUE to the mid day/month", { # nolint

  expected_output <- tibble::tribble(
      ~XXSTDTC, ~ASTDTM, ~ASTDTF, ~ASTTMF,
      "2019-07-18T15:25:40", ymd_hms("2019-07-18T15:25:40"), NA_character_, NA_character_,
      "2019-07-18T15:25", ymd_hms("2019-07-18T15:25:00"), NA_character_, "S",
      "2019-07-18T15", ymd_hms("2019-07-18T15:00:00"), NA_character_, "M",
      "2019-07-18", ymd_hms("2019-07-18T00:00:00"), NA_character_, "H",
      "2019-02", ymd_hms("2019-02-15T00:00:00"), "D", "H",
      "2019", ymd_hms("2019-06-30T00:00:00"), "M", "H",
      "2019---07", ymd_hms("2019-06-07T00:00:00"), "M", "H"
    )

    actual_output <- derive_vars_dtm(
      input,
      new_vars_prefix = "AST",
      dtc = XXSTDTC,
      date_imputation = "MID",
      preserve = TRUE
    )

    expect_equal(
      expected_output,
      actual_output
    )

  })
