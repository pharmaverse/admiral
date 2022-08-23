library(admiral.test)
data(admiral_ae)
data(admiral_vs)

test_that("call_derivation works", {
  library(dplyr)
  input <- admiral_vs[sample(seq_len(nrow(admiral_vs)), 1000), ]

  expected_output <- input %>%
    derive_summary_records(
      by_vars = vars(USUBJID, VSTESTCD),
      analysis_var = VSSTRESN,
      summary_fun = function(x) mean(x, na.rm = TRUE),
      set_values_to = vars(DTYPE = "AVERAGE"),
      filter = n() >= 2L
    ) %>%
    derive_summary_records(
      by_vars = vars(USUBJID, VSTESTCD),
      analysis_var = VSSTRESN,
      summary_fun = function(x) max(x, na.rm = TRUE),
      set_values_to = vars(DTYPE = "MAXIMUM"),
      filter = n() >= 2L
    ) %>%
    derive_summary_records(
      by_vars = vars(USUBJID, VSTESTCD),
      analysis_var = VSSTRESN,
      summary_fun = function(x) min(x, na.rm = TRUE),
      set_values_to = vars(DTYPE = "MINIMUM"),
      filter = n() >= 2L
    )

  actual_output <- call_derivation(
    dataset = input,
    derivation = derive_summary_records,
    variable_params = list(
      params(
        summary_fun = function(x) mean(x, na.rm = TRUE),
        set_values_to = vars(DTYPE = "AVERAGE")
      ),
      params(
        summary_fun = function(x) max(x, na.rm = TRUE),
        set_values_to = vars(DTYPE = "MAXIMUM")
      ),
      params(
        summary_fun = function(x) min(x, na.rm = TRUE),
        set_values_to = vars(DTYPE = "MINIMUM")
      )
    ),
    by_vars = vars(USUBJID, VSTESTCD),
    analysis_var = VSSTRESN,
    filter = n() >= 2L
  )

  expect_dfs_equal(
    expected_output,
    actual_output,
    keys = c("USUBJID", "VSTESTCD", "VISIT", "DTYPE", "VSSEQ")
  )
})

test_that("call_derivation works", {
  input <- admiral_ae[sample(seq_len(nrow(admiral_ae)), 1000), ] %>%
    left_join(admiral_adsl, by = "USUBJID")

  expected_output <- input %>%
    derive_vars_dt(
      new_vars_prefix = "AST",
      dtc = AESTDTC,
      date_imputation = "first",
      min_dates = vars(TRTSDT),
      max_dates = vars(TRTEDT)
    ) %>%
    derive_vars_dt(
      new_vars_prefix = "AEN",
      dtc = AEENDTC,
      date_imputation = "last",
      min_dates = vars(TRTSDT),
      max_dates = vars(TRTEDT)
    )

  actual_output <- call_derivation(
    dataset = input,
    derivation = derive_vars_dt,
    variable_params = list(
      params(dtc = AESTDTC, date_imputation = "first", new_vars_prefix = "AST"),
      params(dtc = AEENDTC, date_imputation = "last", new_vars_prefix = "AEN")
    ),
    min_dates = vars(TRTSDT),
    max_dates = vars(TRTEDT)
  )

  expect_dfs_equal(expected_output, actual_output, keys = c("USUBJID", "AESEQ"))
})


test_that("call_derivation - Error is thrown if ... has no arguments", {
  input <- admiral_ae[sample(seq_len(nrow(admiral_ae)), 1000), ] %>%
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

test_that("call_derivation - Error is thrown if ... arguments are not properly named", {
  input <- admiral_ae[sample(seq_len(nrow(admiral_ae)), 1000), ] %>%
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

test_that("call_derivation - Error is thrown params is empty", {
  input <- admiral_ae[sample(seq_len(nrow(admiral_ae)), 1000), ] %>%
    left_join(admiral_adsl, by = "USUBJID")

  expect_error(
    call_derivation(
      dataset = input,
      derivation = derive_vars_dt,
      variable_params = list(
        params(),
        params()
      ),
      min_dates = vars(TRTSDT),
      max_dates = vars(TRTEDT)
    ), "At least one argument must be provided"
  )
})

test_that("call_derivation - Error is thrown if passed params are not proprely named", {
  input <- admiral_ae[sample(seq_len(nrow(admiral_ae)), 1000), ] %>%
    left_join(admiral_adsl, by = "USUBJID")

  expect_error(
    call_derivation(
      dataset = input,
      derivation = derive_vars_dt,
      variable_params = list(
        params(XYZ),
        params(XYZ)
      ),
      min_dates = vars(TRTSDT),
      max_dates = vars(TRTEDT)
    ), "All arguments passed to `params()` must be named",
    fixed = TRUE
  )
})

test_that("call_derivation - Error is thrown if `...` arguments are not properly named", {
  input <- admiral_ae[sample(seq_len(nrow(admiral_ae)), 1000), ] %>%
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
