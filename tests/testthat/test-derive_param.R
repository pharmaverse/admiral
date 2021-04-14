context("test-derive_param")

test_that("mapping works if `by_vars` are all variables in `mapping`", {
  mapping <- tibble::tribble(
    ~VSTESTCD, ~VSPOS,    ~PARAM,                   ~PARAMCD,
    "WEIGHT",  "",        "Weight (kg)",            "WEIGHT",
    "DIABP",   "",        "Dia. BP (mmHG)",         "DIABP",
    "DIABP",   "SITTING", "Sitting Dia. BP (mmHG)", "SDIABP"
  )
  input <- tibble::tribble(
    ~USUBJID, ~VSTESTCD, ~VSPOS,
    "P01",    "WEIGHT",  "",
    "P01",    "DIABP",   "",
    "P02",    "WEIGHT",  "",
    "P02",    "DIABP",   "SITTING"
  )
  new_cols <- tibble::tribble(
    ~PARAM,                   ~PARAMCD,
    "Weight (kg)",            "WEIGHT",
    "Dia. BP (mmHG)",         "DIABP",
    "Weight (kg)",            "WEIGHT",
    "Sitting Dia. BP (mmHG)", "SDIABP"
  )
  expected_output <- dplyr::bind_cols(input, new_cols)

  expect_dfs_equal(
    expected_output,
    derive_param(input, mapping, by_vars = c("VSTESTCD", "VSPOS")),
    keys = c("USUBJID", "VSTESTCD")
  )
})

test_that("mapping works if `by_vars` is a subset of variables in `mapping`", {
  mapping <- tibble::tribble(
    ~VSTESTCD, ~VSPOS,    ~PARAM,                   ~PARAMCD,
    "WEIGHT",  "",        "Weight (kg)",            "WEIGHT",
    "DIABP",   "",        "Dia. BP (mmHG)",         "DIABP",
    "DIABP",   "SITTING", "Sitting Dia. BP (mmHG)", "SDIABP"
  )
  input <- tibble::tribble(
    ~USUBJID, ~VSTESTCD, ~VSPOS,
    "P01",    "WEIGHT",  "",
    "P01",    "DIABP",   "",
    "P02",    "WEIGHT",  "",
    "P02",    "DIABP",   "SITTING"
  )
  new_cols <- tibble::tribble(
    ~PARAM,           ~PARAMCD,
    "Weight (kg)",    "WEIGHT",
    "Dia. BP (mmHG)", "DIABP",
    "Weight (kg)",    "WEIGHT",
    "Dia. BP (mmHG)", "DIABP"
  )
  expected_output <- dplyr::bind_cols(input, new_cols)

  expect_dfs_equal(
    expected_output,
    derive_param(input, mapping, by_vars = "VSTESTCD"),
    keys = c("USUBJID", "VSTESTCD")
  )
})

test_that("a warning is issued when a parameter cannot be mapped", {
  mapping <- tibble::tribble(
    ~VSTESTCD, ~VSPOS, ~PARAM,        ~PARAMCD,
    "WEIGHT",  "",     "Weight (kg)", "WEIGHT",
  )
  input <- tibble::tribble(
    ~USUBJID, ~VSTESTCD, ~VSPOS,
    "P01",    "WEIGHT",  "",
    "P01",    "DIABP",   "",
    "P02",    "WEIGHT",  "",
    "P02",    "DIABP",   "SITTING"
  )

  expect_warning(
    derive_param(input, mapping, by_vars = "VSTESTCD"),
    "Mapping a `PARAM` failed for the following records:"
  )
})
