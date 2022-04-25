dm <- tibble::tribble(
  ~STUDYID, ~USUBJID,
  "TEST01", "PAT01",
  "TEST01", "PAT02",
  "TEST01", "PAT03",
  "TEST01", "PAT04"
)

ds <- tibble::tribble(
  ~STUDYID, ~USUBJID, ~DSCAT, ~DSDECOD, ~DSTERM, ~DSSTDTC,
  "TEST01", "PAT01", "PROTOCOL MILESTONE", "INFORMED CONSENT OBTAINED", "INFORMED CONSENT OBTAINED", "2021-04-01", # nolint
  "TEST01", "PAT01", "PROTOCOL MILESTONE", "RANDOMIZATION", "RANDOMIZED", "2021-04-11",
  "TEST01", "PAT01", "DISPOSITION EVENT", "ADVERSE EVENT", "ADVERSE EVENT", "2021-12-01",
  "TEST01", "PAT01", "OTHER EVENT", "DEATH", "DEATH", "2022-02-01",
  "TEST01", "PAT02", "PROTOCOL MILESTONE", "INFORMED CONSENT OBTAINED", "INFORMED CONSENT OBTAINED", "2021-04-02", # nolint
  "TEST01", "PAT02", "PROTOCOL MILESTONE", "RANDOMIZATION", "RANDOMIZED", "2021-04-11",
  "TEST01", "PAT02", "DISPOSITION EVENT", "COMPLETED", NA_character_, "2021-12-01",
  "TEST01", "PAT02", "OTHER EVENT", "DEATH", "DEATH", "2022-04",
  "TEST01", "PAT03", "PROTOCOL MILESTONE", "INFORMED CONSENT OBTAINED", "INFORMED CONSENT OBTAINED", "2021-04-02", # nolint
  "TEST01", "PAT03", "PROTOCOL MILESTONE", "RANDOMIZATION", "RANDOMIZED", "2021-04-11",
  "TEST01", "PAT03", "DISPOSITION EVENT", "OTHER", "MISTAKE IN CALCULATION", "2021-04-29",
  "TEST01", "PAT03", "OTHER EVENT", "DEATH", "DEATH", "2022-04",
  "TEST01", "PAT04", "PROTOCOL MILESTONE", "INFORMED CONSENT OBTAINED", "INFORMED CONSENT OBTAINED", "2021-04-02", # nolint
  "TEST01", "PAT04", "PROTOCOL MILESTONE", "RANDOMIZATION", "RANDOMIZED", "2021-04-11"
)


test_that("Derive DCSREAS using default mapping", {
  expected_output <- tibble::tribble(
    ~STUDYID, ~USUBJID, ~DCSREAS,
    "TEST01", "PAT01", "ADVERSE EVENT",
    "TEST01", "PAT02", NA_character_,
    "TEST01", "PAT03", "OTHER",
    "TEST01", "PAT04", NA_character_
  )

  actual_output <- derive_vars_disposition_reason(
    dataset = dm,
    dataset_ds = ds,
    new_var = DCSREAS,
    reason_var = DSDECOD,
    filter_ds = DSCAT == "DISPOSITION EVENT"
  )

  expect_dfs_equal(
    expected_output,
    actual_output,
    keys = c("STUDYID", "USUBJID")
  )
})

test_that("Derive DCSREAS DCSREASP using default mapping", {
  expected_output <- tibble::tribble(
    ~STUDYID, ~USUBJID, ~DCSREAS, ~DCTREASP,
    "TEST01", "PAT01", "ADVERSE EVENT", NA_character_,
    "TEST01", "PAT02", NA_character_, NA_character_,
    "TEST01", "PAT03", "OTHER", "MISTAKE IN CALCULATION",
    "TEST01", "PAT04", NA_character_, NA_character_
  )

  actual_output <- derive_vars_disposition_reason(
    dataset = dm,
    dataset_ds = ds,
    new_var = DCSREAS,
    reason_var = DSDECOD,
    new_var_spe = DCTREASP,
    reason_var_spe = DSTERM,
    filter_ds = DSCAT == "DISPOSITION EVENT"
  )

  expect_dfs_equal(
    expected_output,
    actual_output,
    keys = c("STUDYID", "USUBJID")
  )
})


test_that("Derive DCTREAS, DCTREASP using a study specific mapping", {
  format_dctreas <- function(x, y = NULL) {
    if (is.null(y)){
      if_else(x %notin% c("COMPLETED", "SCREEN FAILURE") & !is.na(x), x, NA_character_)
    } else {
      if_else (x == "OTHER", y, NA_character_)
    }
  }
  expected_output <- tibble::tribble(
    ~STUDYID, ~USUBJID, ~DCTREAS, ~DCTREASP,
    "TEST01", "PAT01", "ADVERSE EVENT", NA_character_,
    "TEST01", "PAT02", NA_character_, NA_character_,
    "TEST01", "PAT03", "OTHER", "MISTAKE IN CALCULATION",
    "TEST01", "PAT04", NA_character_, NA_character_
  )

  actual_output <- derive_vars_disposition_reason(
    dataset = dm,
    dataset_ds = ds,
    new_var = DCTREAS,
    reason_var = DSDECOD,
    new_var_spe = DCTREASP,
    reason_var_spe = DSTERM,
    format_new_var = format_dctreas,
    filter_ds = DSCAT == "DISPOSITION EVENT"
  )

  expect_dfs_equal(
    expected_output,
    actual_output,
    keys = c("STUDYID", "USUBJID")
  )
})

test_that("derive_vars_disposition_reason checks new_var_spe and reason_var_spe", {

  expect_error(
    derive_vars_disposition_reason(
      dataset = dm,
      dataset_ds = ds,
      new_var = DCSREAS,
      new_var_spe = DCTREASP,
      reason_var = DSDECOD,
      filter_ds = DSCAT == "DISPOSITION EVENT"
    ),
    regexp = paste(
      "^`new_var_spe` is specified as  .* but `reason_var_spe` is NULL.",
      "Please specify `reason_var_spe` together with `new_var_spe`."
    )
  )

})
