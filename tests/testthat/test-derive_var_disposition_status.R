context("test-derive_var_disposition_status")

dm <- tibble::tribble(
  ~STUDYID, ~USUBJID,
  "TEST01", "PAT01",
  "TEST01", "PAT02",
  "TEST01", "PAT03",
  "TEST01", "PAT04"
)

ds <- tibble::tribble(
  ~STUDYID, ~USUBJID, ~DSCAT, ~DSDECOD, ~DSSTDTC,
  "TEST01", "PAT01", "PROTOCOL MILESTONE", "INFORMED CONSENT OBTAINED", "2021-04-01",
  "TEST01", "PAT01", "PROTOCOL MILESTONE", "RANDOMIZATION", "2021-04-11",
  "TEST01", "PAT01", "DISPOSITION EVENT", "ADVERSE EVENT", "2021-12-01",
  "TEST01", "PAT01", "OTHER EVENT", "DEATH", "2022-02-01",
  "TEST01", "PAT02", "PROTOCOL MILESTONE", "INFORMED CONSENT OBTAINED", "2021-04-02",
  "TEST01", "PAT02", "PROTOCOL MILESTONE", "RANDOMIZATION", "2021-04-11",
  "TEST01", "PAT02", "DISPOSITION EVENT", "COMPLETED", "2021-12-01",
  "TEST01", "PAT02", "OTHER EVENT", "DEATH", "2022-04",
  "TEST01", "PAT03", "PROTOCOL MILESTONE", "INFORMED CONSENT OBTAINED", "2021-04-02",
  "TEST01", "PAT03", "PROTOCOL MILESTONE", "RANDOMIZATION", "2021-04-11",
  "TEST01", "PAT03", "DISPOSITION EVENT", "PROGRESSIVE DISEASE", "2021-05-01",
  "TEST01", "PAT03", "OTHER EVENT", "DEATH", "2022-04",
  "TEST01", "PAT04", "PROTOCOL MILESTONE", "INFORMED CONSENT OBTAINED", "2021-04-02",
  "TEST01", "PAT04", "PROTOCOL MILESTONE", "RANDOMIZATION", "2021-04-11")



test_that("Derive EOSSTT using default mapping", {
  expected_output <- tibble::tribble(
    ~STUDYID, ~USUBJID, ~EOSSTT,
    "TEST01", "PAT01", "DISCONTINUED",
    "TEST01", "PAT02", "COMPLETED",
    "TEST01", "PAT03", "DISCONTINUED",
    "TEST01", "PAT04", "ONGOING"
  )

  actual_output <- derive_var_disposition_status(
    dataset = dm,
    dataset_ds = ds,
    new_var = EOSSTT,
    status_var = DSDECOD,
    filter_ds = DSCAT == "DISPOSITION EVENT"
  )

  expect_dfs_equal(
    expected_output,
    actual_output,
    keys = c("STUDYID", "USUBJID")
  )
})

test_that("Derive EOTSTT using a study specific mapping", {
  format_eosstt <- function(x) {
    case_when(
      x == "COMPLETED" ~ "COMPLETED",
      x == "ADVERSE EVENT" ~ "DISCONTINUED DUE TO AE",
      x %notin% c("ADVERSE EVENT", "COMPLETED") & !is.na(x) ~ "DISCONTINUED NOT DUE TO AE",
      TRUE ~ "ONGOING"
    )
  }
  expected_output <- tibble::tribble(
    ~STUDYID, ~USUBJID, ~EOSSTT,
    "TEST01", "PAT01", "DISCONTINUED DUE TO AE",
    "TEST01", "PAT02", "COMPLETED",
    "TEST01", "PAT03", "DISCONTINUED NOT DUE TO AE",
    "TEST01", "PAT04", "ONGOING"
  )

  actual_output <- derive_var_disposition_status(
    dataset = dm,
    dataset_ds = ds,
    new_var = EOSSTT,
    status_var = DSDECOD,
    format_new_var = format_eosstt,
    filter_ds = DSCAT == "DISPOSITION EVENT"
  )

  expect_dfs_equal(
    expected_output,
    actual_output,
    keys = c("STUDYID", "USUBJID")
  )
})
