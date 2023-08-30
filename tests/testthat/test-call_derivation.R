## Test 1: Test that call_derivation generates expected summary  ----
# ---- call_derivation Test 1:  Test that call_derivation generates expected summary ----
test_that("call_derivation Test 1:  Test that call_derivation generates expected summary", {
  input <- pharmaversesdtm::vs[sample(seq_len(nrow(pharmaversesdtm::vs)), 1000), ]

  expected_output <- input %>%
    derive_summary_records(
      by_vars = exprs(USUBJID, VSTESTCD),
      analysis_var = VSSTRESN,
      summary_fun = function(x) mean(x, na.rm = TRUE),
      set_values_to = exprs(DTYPE = "AVERAGE"),
      filter = dplyr::n() >= 2L
    ) %>%
    derive_summary_records(
      by_vars = exprs(USUBJID, VSTESTCD),
      analysis_var = VSSTRESN,
      summary_fun = function(x) max(x, na.rm = TRUE),
      set_values_to = exprs(DTYPE = "MAXIMUM"),
      filter = dplyr::n() >= 2L
    ) %>%
    derive_summary_records(
      by_vars = exprs(USUBJID, VSTESTCD),
      analysis_var = VSSTRESN,
      summary_fun = function(x) min(x, na.rm = TRUE),
      set_values_to = exprs(DTYPE = "MINIMUM"),
      filter = dplyr::n() >= 2L
    )

  actual_output <- call_derivation(
    dataset = input,
    derivation = derive_summary_records,
    variable_params = list(
      params(
        summary_fun = function(x) mean(x, na.rm = TRUE),
        set_values_to = exprs(DTYPE = "AVERAGE")
      ),
      params(
        summary_fun = function(x) max(x, na.rm = TRUE),
        set_values_to = exprs(DTYPE = "MAXIMUM")
      ),
      params(
        summary_fun = function(x) min(x, na.rm = TRUE),
        set_values_to = exprs(DTYPE = "MINIMUM")
      )
    ),
    by_vars = exprs(USUBJID, VSTESTCD),
    analysis_var = VSSTRESN,
    filter = dplyr::n() >= 2L
  )

  expect_dfs_equal(
    expected_output,
    actual_output,
    keys = c("USUBJID", "VSTESTCD", "VISIT", "DTYPE", "VSSEQ")
  )
})

## Test 2: Test that call_derivation generates expected imputation  ----
# ---- call_derivation Test 2: Test that call_derivation generates expected imputation ----
test_that("call_derivation Test 2: Test that call_derivation generates expected imputation", {
  input <- pharmaversesdtm::ae[sample(seq_len(nrow(pharmaversesdtm::ae)), 1000), ] %>%
    left_join(admiral_adsl, by = "USUBJID")

  expected_output <- input %>%
    derive_vars_dt(
      new_vars_prefix = "AST",
      dtc = AESTDTC,
      date_imputation = "first",
      min_dates = exprs(TRTSDT),
      max_dates = exprs(TRTEDT)
    ) %>%
    derive_vars_dt(
      new_vars_prefix = "AEN",
      dtc = AEENDTC,
      date_imputation = "last",
      min_dates = exprs(TRTSDT),
      max_dates = exprs(TRTEDT)
    )

  actual_output <- call_derivation(
    dataset = input,
    derivation = derive_vars_dt,
    variable_params = list(
      params(dtc = AESTDTC, date_imputation = "first", new_vars_prefix = "AST"),
      params(dtc = AEENDTC, date_imputation = "last", new_vars_prefix = "AEN")
    ),
    min_dates = exprs(TRTSDT),
    max_dates = exprs(TRTEDT)
  )

  expect_dfs_equal(expected_output, actual_output, keys = c("USUBJID", "AESEQ"))
})

## Test 3: Test that Error is thrown if ... has no arguments  ----
# ---- call_derivation Test 3: Test that Error is thrown if ... has no arguments ----
test_that("call_derivation Test 3: Test that Error is thrown if ... has no arguments", {
  input <- pharmaversesdtm::ae[sample(seq_len(nrow(pharmaversesdtm::ae)), 1000), ] %>%
    left_join(admiral_adsl, by = "USUBJID")

  expect_error(
    call_derivation(
      dataset = input,
      derivation = derive_vars_dt,
      variable_params = list(
        params(dtc = AESTDTC, date_imputation = "first", new_vars_prefix = "AST"),
        params(dtc = AEENDTC, date_imputation = "last", new_vars_prefix = "AEN")
      )
    ), "At least one argument must be set inside `...`"
  )
})

## Test 4: Error is thrown if ... arguments are not properly named ----
# ---- call_derivation Test 4: Error is thrown if ... arguments are not properly named ----
test_that("call_derivation Test 4: Error is thrown if ... arguments are not properly named", {
  input <- pharmaversesdtm::ae[sample(seq_len(nrow(pharmaversesdtm::ae)), 1000), ] %>%
    left_join(admiral_adsl, by = "USUBJID")

  expect_error(
    call_derivation(
      dataset = input,
      derivation = derive_vars_dt,
      variable_params = list(
        params(dtc = AESTDTC, date_imputation = "first", new_vars_prefix = "AST"),
        params(dtc = AEENDTC, date_imputation = "last", new_vars_prefix = "AEN")
      ),
      XYZSDT,
      XYZEDT
    )
  )
})

## Test 5: Error is thrown params is empty ----
# ---- call_derivation Test 5: Error is thrown params is empty ----
test_that("call_derivation Test 5: Error is thrown params is empty", {
  input <- pharmaversesdtm::ae[sample(seq_len(nrow(pharmaversesdtm::ae)), 1000), ] %>%
    left_join(admiral_adsl, by = "USUBJID")

  expect_error(
    call_derivation(
      dataset = input,
      derivation = derive_vars_dt,
      variable_params = list(
        params(),
        params()
      ),
      min_dates = exprs(TRTSDT),
      max_dates = exprs(TRTEDT)
    ), "At least one argument must be provided"
  )
})

## Test 6: Error is thrown if passed params are not properly named ----
# ---- call_derivation Test 6: Error is thrown if passed params are not properly named ----
test_that("call_derivation Test 6: Error is thrown if passed params are not properly named", {
  input <- pharmaversesdtm::ae[sample(seq_len(nrow(pharmaversesdtm::ae)), 1000), ] %>%
    left_join(admiral_adsl, by = "USUBJID")

  expect_error(
    call_derivation(
      dataset = input,
      derivation = derive_vars_dt,
      variable_params = list(
        params(XYZ),
        params(XYZ)
      ),
      min_dates = exprs(TRTSDT),
      max_dates = exprs(TRTEDT)
    ), "All arguments passed to `params()` must be named",
    fixed = TRUE
  )
})

## Test 7: Error is thrown if `...` arguments are not properly named ----
# ---- call_derivation Test 7: Error is thrown if `...` arguments are not properly named ----
test_that("call_derivation Test 7: Error is thrown if `...` arguments are not properly named", {
  input <- pharmaversesdtm::ae[sample(seq_len(nrow(pharmaversesdtm::ae)), 1000), ] %>%
    left_join(admiral_adsl, by = "USUBJID")

  expect_error(
    call_derivation(
      dataset = input,
      derivation = derive_vars_dt,
      variable_params = list(
        params(dtc = AESTDTC, date_imputation = "first", new_vars_prefix = "AST"),
        params(dtc = AEENDTC, date_imputation = "last", new_vars_prefix = "AEN")
      ),
      XYZSDT,
      XYZEDT
    )
  )
})

## Test 8: Error is thrown if duplicate parameters ----
# ---- call_derivation Test 8: Error is thrown if duplicate parameters ----
test_that("call_derivation Test 8: Error is thrown if duplicate parameters", {
  expect_error(
    params(dtc = VSDTC, dtc = VSDTC, new_vars_prefix = "A"),
    "The following parameters have been specified more than once: `dtc`",
    fixed = TRUE
  )
})
