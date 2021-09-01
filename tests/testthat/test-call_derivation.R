data(ae)
data(vs)

test_that("call_derivation works", {
  input <- vs[sample(1:nrow(vs), 1000), ]

  expected_output <- input %>%
    derive_summary_records(
      by_vars = vars(USUBJID, VSTESTCD),
      analysis_var = VSSTRESN,
      summary_function = function(x) mean(x, na.rm = TRUE),
      set_values_to = vars(DTYPE = "AVERAGE"),
      filter = dplyr::n() >= 2L
    ) %>%
    derive_summary_records(
      by_vars = vars(USUBJID, VSTESTCD),
      analysis_var = VSSTRESN,
      summary_function = function(x) max(x, na.rm = TRUE),
      set_values_to = vars(DTYPE = "MAXIMUM"),
      filter = dplyr::n() >= 2L
    ) %>%
    derive_summary_records(
      by_vars = vars(USUBJID, VSTESTCD),
      analysis_var = VSSTRESN,
      summary_function = function(x) min(x, na.rm = TRUE),
      set_values_to = vars(DTYPE = "MINIMUM"),
      filter = dplyr::n() >= 2L
    )

  actual_output <- call_derivation(
    dataset = input,
    derivation = derive_summary_records,
    variable_params = list(
      params(analysis_var = VSSTRESN, summary_function = function(x) mean(x, na.rm = TRUE), set_values_to = vars(DTYPE = "AVERAGE")),
      params(analysis_var = VSSTRESN, summary_function = function(x) max(x, na.rm = TRUE), set_values_to = vars(DTYPE = "MAXIMUM")),
      params(analysis_var = VSSTRESN, summary_function = function(x) min(x, na.rm = TRUE), set_values_to = vars(DTYPE = "MINIMUM"))
    ),
    by_vars = vars(USUBJID, VSTESTCD),
    filter = dplyr::n() >= 2L
  )

  expect_dfs_equal(expected_output, actual_output, keys = c("USUBJID", "VSTESTCD", "VISIT", "DTYPE", "VSSEQ"))
})

test_that("call_derivation works", {
  input <- ae[sample(1:nrow(ae), 1000), ] %>%
    left_join(adsl, by = "USUBJID")

  expected_output <- input %>%
    derive_vars_dt(
      new_vars_prefix = "AST",
      dtc = AESTDTC,
      date_imputation = "first",
      min_dates = list(TRTSDT),
      max_dates = list(TRTEDT)
    ) %>%
    derive_vars_dt(
      new_vars_prefix = "AEN",
      dtc = AEENDTC,
      date_imputation = "last",
      min_dates = list(TRTSDT),
      max_dates = list(TRTEDT)
    )

  actual_output <- call_derivation(
    dataset = input,
    derivation = derive_vars_dt,
    variable_params = list(
      params(dtc = AESTDTC, date_imputation = "first", new_vars_prefix = "AST"),
      params(dtc = AEENDTC, date_imputation = "last", new_vars_prefix = "AEN")
    ),
    min_dates = list(TRTSDT),
    max_dates = list(TRTEDT)
  )

  expect_dfs_equal(expected_output, actual_output, keys = c("USUBJID", "AESEQ"))
})
