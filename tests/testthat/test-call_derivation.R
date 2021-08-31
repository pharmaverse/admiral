data(ae)
data(vs)

test_that("call_derivation works", {
  input <- vs[sample(1:nrow(vs), 1000), ]

  expected_output <- input %>%
    derive_summary_records(
      by_vars = vars(USUBJID, VSTESTCD),
      fns = list(VSSTRESN ~ mean(., na.rm = TRUE)),
      set_values_to = vars(DTYPE = "AVERAGE"),
      filter_rows = dplyr::n() >= 2L
    ) %>%
    derive_summary_records(
      by_vars = vars(USUBJID, VSTESTCD),
      fns = list(VSSTRESN ~ max(., na.rm = TRUE)),
      set_values_to = vars(DTYPE = "MAXIMUM"),
      filter_rows = dplyr::n() >= 2L
    ) %>%
    derive_summary_records(
      by_vars = vars(USUBJID, VSTESTCD),
      fns = list(VSSTRESN ~ min(., na.rm = TRUE)),
      set_values_to = vars(DTYPE = "MINIMUM"),
      filter_rows = dplyr::n() >= 2L
    )

  actual_output <- call_derivation(
    dataset = input,
    derivation = derive_summary_records,
    variable_params = list(
      params(fns = list(VSSTRESN ~ mean(., na.rm = TRUE)), set_values_to = vars(DTYPE = "AVERAGE")),
      params(fns = list(VSSTRESN ~ max(., na.rm = TRUE)), set_values_to = vars(DTYPE = "MAXIMUM")),
      params(fns = list(VSSTRESN ~ min(., na.rm = TRUE)), set_values_to = vars(DTYPE = "MINIMUM"))
    ),
    by_vars = vars(USUBJID, VSTESTCD),
    filter_rows = dplyr::n() >= 2L
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
