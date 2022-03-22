test_that("Test adding absolute records for each by group", {
  input <- tibble::tribble(
  ~USUBJID, ~LBTESTCD, ~LBSTRESN, ~LBSTRESU, ~VISIT,
  "P01",    "WBC",     33,   "10^9/L",       "CYCLE 1 DAY 1",
  "P01",    "WBC",     38,   "10^9/L",       "CYCLE 2 DAY 1",
  "P01",    "LYMLE",   90,   "fraction of 1", "CYCLE 1 DAY 1",
  "P01",    "LYMLE",   70,   "fraction of 1", "CYCLE 2 DAY 1"
  )
  expected_output <- tibble::tribble(
    ~USUBJID, ~LBTESTCD, ~LBSTRESN, ~LBSTRESU, ~VISIT, ~PARAMCD, ~DTYPE, ~AVAL,
    "P01",    "WBC",     33,   "10^9/L",       "CYCLE 1 DAY 1", NA_character_, NA_character_, NA,
    "P01",    "WBC",     38,   "10^9/L",       "CYCLE 2 DAY 1", NA_character_, NA_character_, NA,
    "P01",    "LYMLE",   90,   "fraction of 1", "CYCLE 1 DAY 1", NA_character_, NA_character_, NA,
    "P01",    "LYMLE",   70,   "fraction of 1", "CYCLE 2 DAY 1", NA_character_, NA_character_, NA,
    "P01",    NA_character_, NA, NA_character_, "CYCLE 1 DAY 1", "LYMPHSI", "CALCULATION", 29.7,
    "P01",    NA_character_, NA, NA_character_, "CYCLE 2 DAY 1", "LYMPHSI", "CALCULATION", 26.6
  )

  expect_equal(
    derive_param_wbc_abs(dataset = input,
                          by_vars = vars(USUBJID, VISIT),
                          set_values_to = vars(PARAMCD = "LYMPHSI",
                                               DTYPE = "CALCULATION"),
                                  filter_diff = LBSTRESU == "fraction of 1",
                                  wbc_code = "WBC",
                                  diff_code = "LYMLE"),
    expected_output
  )
})

test_that("Test when only one of WBC/differential is present", {
  input <- tibble::tribble(
    ~USUBJID, ~LBTESTCD, ~LBSTRESN, ~LBSTRESU, ~VISIT,
    "P01",    "WBC",     33,   "10^9/L",       "CYCLE 1 DAY 1",
    "P01",    "WBC",     38,   "10^9/L",       "CYCLE 2 DAY 1",
    "P01",    "LYMLE",   70,   "fraction of 1", "CYCLE 2 DAY 1",
    "P01",    "LYMLE",   80,   "fraction of 1", "CYCLE 3 DAY 1"
  )
  expected_output <- tibble::tribble(
    ~USUBJID, ~LBTESTCD, ~LBSTRESN, ~LBSTRESU, ~VISIT, ~PARAMCD, ~DTYPE, ~AVAL,
    "P01",    "WBC",     33,   "10^9/L",       "CYCLE 1 DAY 1", NA_character_, NA_character_, NA,
    "P01",    "WBC",     38,   "10^9/L",       "CYCLE 2 DAY 1", NA_character_, NA_character_, NA,
    "P01",    "LYMLE",   70,   "fraction of 1", "CYCLE 2 DAY 1", NA_character_, NA_character_, NA,
    "P01",    "LYMLE",   80,   "fraction of 1", "CYCLE 3 DAY 1", NA_character_, NA_character_, NA,
    "P01",    NA_character_, NA, NA_character_, "CYCLE 2 DAY 1", "LYMPHSI", "CALCULATION", 26.6
  )

  expect_equal(
    derive_param_wbc_abs(dataset = input,
                          by_vars = vars(USUBJID, VISIT),
                          set_values_to = vars(PARAMCD = "LYMPHSI",
                                               DTYPE = "CALCULATION"),
                          filter_diff = LBSTRESU == "fraction of 1",
                          wbc_code = "WBC",
                          diff_code = "LYMLE"),
    expected_output
  )
})

test_that("Test when absolute record already present in source dataset 1", {
  input <- tibble::tribble(
    ~USUBJID, ~LBTESTCD, ~LBSTRESN, ~LBSTRESU, ~VISIT,
    "P01",    "WBC",     33,   "10^9/L",       "CYCLE 1 DAY 1",
    "P01",    "WBC",     38,   "10^9/L",       "CYCLE 2 DAY 1",
    "P01",    "LYMPHSI", 27,   "10^9/L",    "CYCLE 2 DAY 1",
    "P01",    "LYMLE",   70,   "fraction of 1", "CYCLE 2 DAY 1",
    "P01",    "LYMLE",   80,   "fraction of 1", "CYCLE 3 DAY 1"
  )
  expected_output <- tibble::tribble(
    ~USUBJID, ~LBTESTCD, ~LBSTRESN, ~LBSTRESU, ~VISIT,
    "P01",    "WBC",     33,   "10^9/L",       "CYCLE 1 DAY 1",
    "P01",    "WBC",     38,   "10^9/L",       "CYCLE 2 DAY 1",
    "P01",    "LYMPHSI", 27,   "10^9/L",    "CYCLE 2 DAY 1",
    "P01",    "LYMLE",   70,   "fraction of 1", "CYCLE 2 DAY 1",
    "P01",    "LYMLE",   80,   "fraction of 1", "CYCLE 3 DAY 1"
  )

  expect_equal(
    derive_param_wbc_abs(dataset = input,
                          by_vars = vars(USUBJID, VISIT),
                          set_values_to = vars(PARAMCD = "LYMPHSI",
                                               DTYPE = "CALCULATION"),
                          filter_diff = LBSTRESU == "fraction of 1",
                          wbc_code = "WBC",
                          diff_code = "LYMLE"),
    expected_output
  )
})


test_that("Test when absolute record already present in source dataset 2", {
  input <- tibble::tribble(
    ~USUBJID, ~LBTESTCD, ~LBSTRESN, ~LBSTRESU, ~VISIT,
    "P01",    "WBC",     33,   "10^9/L",       "CYCLE 1 DAY 1",
    "P01",    "WBC",     38,   "10^9/L",       "CYCLE 2 DAY 1",
    "P01",    "LYMPHSI", 27,   "10^9/L",    "CYCLE 2 DAY 1",
    "P01",    "LYMLE",   70,   "fraction of 1", "CYCLE 2 DAY 1",
    "P01",    "LYMLE",   80,   "fraction of 1", "CYCLE 1 DAY 1"
  )

  expected_output <- tibble::tribble(
    ~USUBJID, ~LBTESTCD, ~LBSTRESN, ~LBSTRESU, ~VISIT, ~PARAMCD, ~DTYPE, ~AVAL,
    "P01",    "WBC",     33,   "10^9/L",       "CYCLE 1 DAY 1", NA_character_, NA_character_, NA,
    "P01",    "WBC",     38,   "10^9/L",       "CYCLE 2 DAY 1", NA_character_, NA_character_, NA,
    "P01",    "LYMPHSI", 27,   "10^9/L",    "CYCLE 2 DAY 1",  NA_character_, NA_character_, NA,
    "P01",    "LYMLE",   70,   "fraction of 1", "CYCLE 2 DAY 1", NA_character_, NA_character_, NA,
    "P01",    "LYMLE",   80,   "fraction of 1", "CYCLE 1 DAY 1", NA_character_, NA_character_, NA,
    "P01",    NA_character_, NA, NA_character_, "CYCLE 1 DAY 1", "LYMPHSI", "CALCULATION", 26.4
  )

  expect_equal(
    derive_param_wbc_abs(dataset = input,
                          by_vars = vars(USUBJID, VISIT),
                          set_values_to = vars(PARAMCD = "LYMPHSI",
                                               DTYPE = "CALCULATION"),
                          filter_diff = LBSTRESU == "fraction of 1",
                          wbc_code = "WBC",
                          diff_code = "LYMLE"),
    expected_output
  )
})
