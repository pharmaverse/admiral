test_that("Test adding absolute records for each by group", {
  input <- tibble::tribble(
    ~USUBJID, ~PARAMCD, ~AVAL, ~PARAM, ~VISIT,
  "P01",    "WBC",     33,   "Leukocyte Count (10^9/L)",       "CYCLE 1 DAY 1",
  "P01",    "WBC",     38,   "Leukocyte Count (10^9/L)",       "CYCLE 2 DAY 1",
  "P01",    "LYMLE",   90,   "Lymphocytes (fraction of 1)", "CYCLE 1 DAY 1",
  "P01",    "LYMLE",   70,   "Lymphocytes (fraction of 1)", "CYCLE 2 DAY 1"
  )
  expected_output <- tibble::tribble(
    ~USUBJID, ~PARAMCD, ~AVAL, ~PARAM, ~VISIT, ~DTYPE,
    "P01",    "WBC",     33,   "Leukocyte Count (10^9/L)",       "CYCLE 1 DAY 1", NA_character_,
    "P01",    "WBC",     38,   "Leukocyte Count (10^9/L)",       "CYCLE 2 DAY 1", NA_character_,
    "P01",    "LYMLE",   90,   "Lymphocytes (fraction of 1)", "CYCLE 1 DAY 1", NA_character_,
    "P01",    "LYMLE",   70,   "Lymphocytes (fraction of 1)", "CYCLE 2 DAY 1", NA_character_,
    "P01",    "LYMPH", 29.7, "Lymphocytes Abs (10^9/L)", "CYCLE 1 DAY 1", "CALCULATION",
    "P01",    "LYMPH", 26.6, "Lymphocytes Abs (10^9/L)", "CYCLE 2 DAY 1", "CALCULATION"
  )

  expect_equal(
    derive_param_wbc_abs(dataset = input,
                         by_vars = vars(USUBJID, VISIT),
                         set_values_to = vars(PARAMCD = "LYMPH",
                                              PARAM = "Lymphocytes Abs (10^9/L)",
                                              DTYPE = "CALCULATION"),
                         get_unit_expr = extract_unit(PARAM),
                         wbc_code = "WBC",
                         diff_code = "LYMLE"),
    expected_output
  )
})

test_that("Test when only one of WBC/differential is present", {
  input <- tibble::tribble(
    ~USUBJID, ~PARAMCD, ~AVAL, ~PARAM, ~VISIT,
    "P01",    "WBC",     33,   "Leukocyte Count (10^9/L)",       "CYCLE 1 DAY 1",
    "P01",    "WBC",     38,   "Leukocyte Count (10^9/L)",       "CYCLE 2 DAY 1",
    "P01",    "LYMLE",   70,   "Lymphocytes (fraction of 1)", "CYCLE 2 DAY 1",
    "P01",    "LYMLE",   80,   "Lymphocytes (fraction of 1)", "CYCLE 3 DAY 1"
  )
  expected_output <- tibble::tribble(
    ~USUBJID, ~PARAMCD, ~AVAL, ~PARAM, ~VISIT, ~DTYPE,
    "P01",    "WBC",     33,   "Leukocyte Count (10^9/L)",       "CYCLE 1 DAY 1", NA_character_,
    "P01",    "WBC",     38,   "Leukocyte Count (10^9/L)",       "CYCLE 2 DAY 1", NA_character_,
    "P01",    "LYMLE",   70,   "Lymphocytes (fraction of 1)", "CYCLE 2 DAY 1", NA_character_,
    "P01",    "LYMLE",   80,   "Lymphocytes (fraction of 1)", "CYCLE 3 DAY 1", NA_character_,
    "P01",    "LYMPH" , 26.6,  "Lymphocytes Abs (10^9/L)", "CYCLE 2 DAY 1", "CALCULATION"
  )

  expect_equal(
    derive_param_wbc_abs(dataset = input,
                          by_vars = vars(USUBJID, VISIT),
                          set_values_to = vars(PARAMCD = "LYMPH",
                                               PARAM = "Lymphocytes Abs (10^9/L)",
                                               DTYPE = "CALCULATION"),
                          get_unit_expr = extract_unit(PARAM),
                          wbc_code = "WBC",
                          diff_code = "LYMLE"),
    expected_output
  )
})

test_that("Test when absolute record already present in source dataset 1", {
  input <- tibble::tribble(
    ~USUBJID, ~PARAMCD, ~AVAL, ~PARAM, ~VISIT,
    "P01",    "WBC",     33,   "Leukocyte Count (10^9/L)",       "CYCLE 1 DAY 1",
    "P01",    "WBC",     38,   "Leukocyte Count (10^9/L)",       "CYCLE 2 DAY 1",
    "P01",    "LYMPH",   27,   "Lymphocytes Abs (10^9/L)",    "CYCLE 2 DAY 1",
    "P01",    "LYMLE",   70,   "Lymphocytes (fraction of 1)", "CYCLE 2 DAY 1",
    "P01",    "LYMLE",   80,   "Lymphocytes (fraction of 1)", "CYCLE 3 DAY 1"
  )
  expected_output <- tibble::tribble(
    ~USUBJID, ~PARAMCD, ~AVAL, ~PARAM, ~VISIT,
    "P01",    "WBC",     33,   "Leukocyte Count (10^9/L)",       "CYCLE 1 DAY 1",
    "P01",    "WBC",     38,   "Leukocyte Count (10^9/L)",       "CYCLE 2 DAY 1",
    "P01",    "LYMPH",   27,   "Lymphocytes Abs (10^9/L)",    "CYCLE 2 DAY 1",
    "P01",    "LYMLE",   70,   "Lymphocytes (fraction of 1)", "CYCLE 2 DAY 1",
    "P01",    "LYMLE",   80,   "Lymphocytes (fraction of 1)", "CYCLE 3 DAY 1"
  )

  expect_equal(
    derive_param_wbc_abs(dataset = input,
                          by_vars = vars(USUBJID, VISIT),
                          set_values_to = vars(PARAMCD = "LYMPH",
                                               PARAM = "Lymphocytes Abs (10^9/L)",
                                               DTYPE = "CALCULATION"),
                          get_unit_expr = extract_unit(PARAM),
                          wbc_code = "WBC",
                          diff_code = "LYMLE"),
    expected_output
  )
})


test_that("Test when absolute record already present in source dataset 2", {
  input <- tibble::tribble(
    ~USUBJID, ~PARAMCD, ~AVAL, ~PARAM, ~VISIT,
    "P01",    "WBC",     33,   "Leukocyte Count (10^9/L)",       "CYCLE 1 DAY 1",
    "P01",    "WBC",     38,   "Leukocyte Count (10^9/L)",       "CYCLE 2 DAY 1",
    "P01",    "LYMPH",   27,   "Lymphocytes Abs (10^9/L)",    "CYCLE 2 DAY 1",
    "P01",    "LYMLE",   70,   "Lymphocytes (fraction of 1)", "CYCLE 2 DAY 1",
    "P01",    "LYMLE",   80,   "Lymphocytes (fraction of 1)", "CYCLE 1 DAY 1"
  )

  expected_output <- tibble::tribble(
    ~USUBJID, ~PARAMCD, ~AVAL, ~PARAM, ~VISIT, ~DTYPE,
    "P01",    "WBC",     33,   "Leukocyte Count (10^9/L)",       "CYCLE 1 DAY 1", NA_character_,
    "P01",    "WBC",     38,   "Leukocyte Count (10^9/L)",       "CYCLE 2 DAY 1", NA_character_,
    "P01",    "LYMPH", 27,     "Lymphocytes Abs (10^9/L)",    "CYCLE 2 DAY 1",  NA_character_,
    "P01",    "LYMLE",   70,   "Lymphocytes (fraction of 1)", "CYCLE 2 DAY 1", NA_character_,
    "P01",    "LYMLE",   80,   "Lymphocytes (fraction of 1)", "CYCLE 1 DAY 1", NA_character_,
    "P01",    "LYMPH", 26.4,   "Lymphocytes Abs (10^9/L)", "CYCLE 1 DAY 1", "CALCULATION"
  )

  expect_equal(
    derive_param_wbc_abs(dataset = input,
                          by_vars = vars(USUBJID, VISIT),
                          set_values_to = vars(PARAMCD = "LYMPH",
                                               PARAM = "Lymphocytes Abs (10^9/L)",
                                               DTYPE = "CALCULATION"),
                          get_unit_expr = extract_unit(PARAM),
                          wbc_code = "WBC",
                          diff_code = "LYMLE"),
    expected_output
  )
})
