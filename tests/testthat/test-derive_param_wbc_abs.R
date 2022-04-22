test_that("Test 1: Test adding absolute records for each by group", {
  input <- tibble::tribble(
    ~USUBJID, ~PARAMCD, ~AVAL, ~PARAM, ~VISIT,
  "P01",    "WBC",     33,   "Leukocyte Count (10^9/L)",       "CYCLE 1 DAY 1",
  "P01",    "WBC",     38,   "Leukocyte Count (10^9/L)",       "CYCLE 2 DAY 1",
  "P01",    "LYMLE",   0.9,   "Lymphocytes (fraction of 1)", "CYCLE 1 DAY 1",
  "P01",    "LYMLE",   0.6,   "Lymphocytes (fraction of 1)", "CYCLE 2 DAY 1"
  )
  expected_output <- tibble::tribble(
    ~USUBJID, ~PARAMCD, ~AVAL, ~PARAM, ~VISIT, ~DTYPE,
    "P01",    "WBC",     33,   "Leukocyte Count (10^9/L)",       "CYCLE 1 DAY 1", NA_character_,
    "P01",    "WBC",     38,   "Leukocyte Count (10^9/L)",       "CYCLE 2 DAY 1", NA_character_,
    "P01",    "LYMLE",   0.9,   "Lymphocytes (fraction of 1)", "CYCLE 1 DAY 1", NA_character_,
    "P01",    "LYMLE",   0.6,   "Lymphocytes (fraction of 1)", "CYCLE 2 DAY 1", NA_character_,
    "P01",    "LYMPH", 29.7, "Lymphocytes Abs (10^9/L)", "CYCLE 1 DAY 1", "CALCULATION",
    "P01",    "LYMPH", 22.8, "Lymphocytes Abs (10^9/L)", "CYCLE 2 DAY 1", "CALCULATION"
  )

  expect_equal(
    derive_param_wbc_abs(dataset = input,
                         by_vars = vars(USUBJID, VISIT),
                         set_values_to = vars(PARAMCD = "LYMPH",
                                              PARAM = "Lymphocytes Abs (10^9/L)",
                                              DTYPE = "CALCULATION"),
                         get_unit_expr = extract_unit(PARAM),
                         wbc_code = "WBC",
                         diff_code = "LYMLE",
                         diff_type = "fraction"),
    expected_output
  )
})

test_that("Test 2: Test when only one of WBC/differential is present", {
  input <- tibble::tribble(
    ~USUBJID, ~PARAMCD, ~AVAL, ~PARAM, ~VISIT,
    "P01",    "WBC",     33,   "Leukocyte Count (10^9/L)",       "CYCLE 1 DAY 1",
    "P01",    "WBC",     38,   "Leukocyte Count (10^9/L)",       "CYCLE 2 DAY 1",
    "P01",    "LYMLE",   0.9,   "Lymphocytes (fraction of 1)", "CYCLE 2 DAY 1",
    "P01",    "LYMLE",   0.8,   "Lymphocytes (fraction of 1)", "CYCLE 3 DAY 1"
  )
  expected_output <- tibble::tribble(
    ~USUBJID, ~PARAMCD, ~AVAL, ~PARAM, ~VISIT, ~DTYPE,
    "P01",    "WBC",     33,   "Leukocyte Count (10^9/L)",       "CYCLE 1 DAY 1", NA_character_,
    "P01",    "WBC",     38,   "Leukocyte Count (10^9/L)",       "CYCLE 2 DAY 1", NA_character_,
    "P01",    "LYMLE",   0.9,   "Lymphocytes (fraction of 1)", "CYCLE 2 DAY 1", NA_character_,
    "P01",    "LYMLE",   0.8,   "Lymphocytes (fraction of 1)", "CYCLE 3 DAY 1", NA_character_,
    "P01",    "LYMPH", 34.2,  "Lymphocytes Abs (10^9/L)", "CYCLE 2 DAY 1", "CALCULATION"
  )

  expect_equal(
    derive_param_wbc_abs(dataset = input,
                          by_vars = vars(USUBJID, VISIT),
                          set_values_to = vars(PARAMCD = "LYMPH",
                                               PARAM = "Lymphocytes Abs (10^9/L)",
                                               DTYPE = "CALCULATION"),
                          get_unit_expr = extract_unit(PARAM),
                          wbc_code = "WBC",
                          diff_code = "LYMLE",
                          diff_type = "fraction"),
    expected_output
  )
})

test_that("Test 3: Test when absolute record already present in source dataset 1", {
  input <- tibble::tribble(
    ~USUBJID, ~PARAMCD, ~AVAL, ~PARAM, ~VISIT,
    "P01",    "WBC",     33,   "Leukocyte Count (10^9/L)",       "CYCLE 1 DAY 1",
    "P01",    "WBC",     38,   "Leukocyte Count (10^9/L)",       "CYCLE 2 DAY 1",
    "P01",    "LYMPH",   27,   "Lymphocytes Abs (10^9/L)",    "CYCLE 2 DAY 1",
    "P01",    "LYMLE",   0.7,   "Lymphocytes (fraction of 1)", "CYCLE 2 DAY 1",
    "P01",    "LYMLE",   0.8,   "Lymphocytes (fraction of 1)", "CYCLE 3 DAY 1"
  )
  expected_output <- tibble::tribble(
    ~USUBJID, ~PARAMCD, ~AVAL, ~PARAM, ~VISIT,
    "P01",    "WBC",     33,   "Leukocyte Count (10^9/L)",       "CYCLE 1 DAY 1",
    "P01",    "WBC",     38,   "Leukocyte Count (10^9/L)",       "CYCLE 2 DAY 1",
    "P01",    "LYMPH",   27,   "Lymphocytes Abs (10^9/L)",    "CYCLE 2 DAY 1",
    "P01",    "LYMLE",   0.7,   "Lymphocytes (fraction of 1)", "CYCLE 2 DAY 1",
    "P01",    "LYMLE",   0.8,   "Lymphocytes (fraction of 1)", "CYCLE 3 DAY 1"
  )

  expect_equal(
    derive_param_wbc_abs(dataset = input,
                          by_vars = vars(USUBJID, VISIT),
                          set_values_to = vars(PARAMCD = "LYMPH",
                                               PARAM = "Lymphocytes Abs (10^9/L)",
                                               DTYPE = "CALCULATION"),
                          get_unit_expr = extract_unit(PARAM),
                          wbc_code = "WBC",
                          diff_code = "LYMLE",
                          diff_type = "fraction"),
    expected_output
  )
})


test_that("Test 4: Test when absolute record already present in source dataset 2", {
  input <- tibble::tribble(
    ~USUBJID, ~PARAMCD, ~AVAL, ~PARAM, ~VISIT,
    "P01",    "WBC",     33,   "Leukocyte Count (10^9/L)",       "CYCLE 1 DAY 1",
    "P01",    "WBC",     38,   "Leukocyte Count (10^9/L)",       "CYCLE 2 DAY 1",
    "P01",    "LYMPH",   27,   "Lymphocytes Abs (10^9/L)",    "CYCLE 2 DAY 1",
    "P01",    "LYMLE",   0.7,   "Lymphocytes (fraction of 1)", "CYCLE 2 DAY 1",
    "P01",    "LYMLE",   0.9,   "Lymphocytes (fraction of 1)", "CYCLE 1 DAY 1"
  )

  expected_output <- tibble::tribble(
    ~USUBJID, ~PARAMCD, ~AVAL, ~PARAM, ~VISIT, ~DTYPE,
    "P01",    "WBC",     33,   "Leukocyte Count (10^9/L)",       "CYCLE 1 DAY 1", NA_character_,
    "P01",    "WBC",     38,   "Leukocyte Count (10^9/L)",       "CYCLE 2 DAY 1", NA_character_,
    "P01",    "LYMPH", 27,     "Lymphocytes Abs (10^9/L)",    "CYCLE 2 DAY 1",  NA_character_,
    "P01",    "LYMLE",   0.7,   "Lymphocytes (fraction of 1)", "CYCLE 2 DAY 1", NA_character_,
    "P01",    "LYMLE",   0.9,   "Lymphocytes (fraction of 1)", "CYCLE 1 DAY 1", NA_character_,
    "P01",    "LYMPH", 29.7,   "Lymphocytes Abs (10^9/L)", "CYCLE 1 DAY 1", "CALCULATION"
  )

  expect_equal(
    derive_param_wbc_abs(dataset = input,
                          by_vars = vars(USUBJID, VISIT),
                          set_values_to = vars(PARAMCD = "LYMPH",
                                               PARAM = "Lymphocytes Abs (10^9/L)",
                                               DTYPE = "CALCULATION"),
                          get_unit_expr = extract_unit(PARAM),
                          wbc_code = "WBC",
                          diff_code = "LYMLE",
                          diff_type = "fraction"),
    expected_output
  )
})


test_that("Test 5: Test percent differential type", {
  input <- tibble::tribble(
    ~USUBJID, ~PARAMCD, ~AVAL, ~PARAM, ~VISIT,
    "P01",    "WBC",     33,   "Leukocyte Count (10^9/L)",       "CYCLE 1 DAY 1",
    "P01",    "WBC",     38,   "Leukocyte Count (10^9/L)",       "CYCLE 2 DAY 1",
    "P01",    "LYMPH",   27,   "Lymphocytes Abs (10^9/L)",    "CYCLE 2 DAY 1",
    "P01",    "LYMLE",   70,   "Lymphocytes (%)", "CYCLE 2 DAY 1",
    "P01",    "LYMLE",   90,   "Lymphocytes (%)", "CYCLE 1 DAY 1"
  )

  expected_output <- tibble::tribble(
    ~USUBJID, ~PARAMCD, ~AVAL, ~PARAM, ~VISIT, ~DTYPE,
    "P01",    "WBC",     33,   "Leukocyte Count (10^9/L)",       "CYCLE 1 DAY 1", NA_character_,
    "P01",    "WBC",     38,   "Leukocyte Count (10^9/L)",       "CYCLE 2 DAY 1", NA_character_,
    "P01",    "LYMPH", 27,     "Lymphocytes Abs (10^9/L)",    "CYCLE 2 DAY 1",  NA_character_,
    "P01",    "LYMLE",   70,   "Lymphocytes (%)", "CYCLE 2 DAY 1", NA_character_,
    "P01",    "LYMLE",   90,   "Lymphocytes (%)", "CYCLE 1 DAY 1", NA_character_,
    "P01",    "LYMPH", 29.7,   "Lymphocytes Abs (10^9/L)", "CYCLE 1 DAY 1", "CALCULATION"
  )

  expect_equal(
    derive_param_wbc_abs(dataset = input,
                         by_vars = vars(USUBJID, VISIT),
                         set_values_to = vars(PARAMCD = "LYMPH",
                                              PARAM = "Lymphocytes Abs (10^9/L)",
                                              DTYPE = "CALCULATION"),
                         get_unit_expr = extract_unit(PARAM),
                         wbc_code = "WBC",
                         diff_code = "LYMLE",
                         diff_type = "percent"),
    expected_output
  )
})
