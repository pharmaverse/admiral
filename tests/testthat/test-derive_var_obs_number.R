## Test 1: create default new variable ASEQ ----
test_that("derive_var_obs_number Test 1: create default new variable ASEQ", {
  input <- tibble::tribble(
    ~STUDYID, ~USUBJID, ~VSTESTCD, ~VISITNUM, ~VSTPTNUM,
    "TEST01", "PAT01", "DIABP", 1, 1,
    "TEST01", "PAT01", "PULSE", 1, 1,
    "TEST01", "PAT01", "PULSE", 1, 2,
    "TEST01", "PAT01", "PULSE", 2, 1,
    "TEST01", "PAT01", "PULSE", 2, 2,
    "TEST01", "PAT01", "SYSBP", 1, 1,
    "TEST01", "PAT01", "SYSBP", 1, 2,
    "TEST01", "PAT01", "SYSBP", 2, 1,
    "TEST01", "PAT02", "DIABP", 1, 1,
    "TEST01", "PAT02", "DIABP", 2, 1,
    "TEST01", "PAT02", "DIABP", 2, 2
  )

  expected_output <- tibble::tribble(
    ~STUDYID, ~USUBJID, ~VSTESTCD, ~VISITNUM, ~VSTPTNUM, ~ASEQ,
    "TEST01", "PAT01", "DIABP", 1, 1, 1L,
    "TEST01", "PAT01", "PULSE", 1, 1, 1L,
    "TEST01", "PAT01", "PULSE", 1, 2, 2L,
    "TEST01", "PAT01", "PULSE", 2, 1, 3L,
    "TEST01", "PAT01", "PULSE", 2, 2, 4L,
    "TEST01", "PAT01", "SYSBP", 1, 1, 1L,
    "TEST01", "PAT01", "SYSBP", 1, 2, 2L,
    "TEST01", "PAT01", "SYSBP", 2, 1, 3L,
    "TEST01", "PAT02", "DIABP", 1, 1, 1L,
    "TEST01", "PAT02", "DIABP", 2, 1, 2L,
    "TEST01", "PAT02", "DIABP", 2, 2, 3L
  )

  actual_output <- derive_var_obs_number(
    input,
    by_vars = exprs(USUBJID, VSTESTCD),
    order = exprs(VISITNUM, VSTPTNUM)
  )

  expect_dfs_equal(
    base = expected_output,
    compare = actual_output,
    keys = c("USUBJID", "VSTESTCD", "VISITNUM", "VSTPTNUM")
  )
})

## Test 2: sorting missing value  is smallest ----
test_that("derive_var_obs_number Test 2: sorting missing value is smallest", {
  input <- tibble::tribble(
    ~STUDYID, ~USUBJID, ~VSTESTCD, ~VISITNUM, ~VSTPTNUM,
    "TEST01", "PAT01", "DIABP", 1, 1,
    "TEST01", "PAT01", "DIABP", 1, NA,
    "TEST01", "PAT01", "PULSE", 1, 1,
    "TEST01", "PAT01", "PULSE", 1, 2,
    "TEST01", "PAT01", "PULSE", 2, 1,
    "TEST01", "PAT01", "PULSE", 2, 2,
    "TEST01", "PAT01", "SYSBP", 1, 1,
    "TEST01", "PAT01", "SYSBP", 1, 2,
    "TEST01", "PAT01", "SYSBP", 2, 1,
    "TEST01", "PAT02", "DIABP", 1, 1,
    "TEST01", "PAT02", "DIABP", 2, 1,
    "TEST01", "PAT02", "DIABP", 2, 2
  )

  expected_output <- tibble::tribble(
    ~STUDYID, ~USUBJID, ~VSTESTCD, ~VISITNUM, ~VSTPTNUM, ~ASEQ,
    "TEST01", "PAT01", "DIABP", 1, NA, 1L,
    "TEST01", "PAT01", "DIABP", 1, 1, 2L,
    "TEST01", "PAT01", "PULSE", 1, 1, 1L,
    "TEST01", "PAT01", "PULSE", 1, 2, 2L,
    "TEST01", "PAT01", "PULSE", 2, 1, 3L,
    "TEST01", "PAT01", "PULSE", 2, 2, 4L,
    "TEST01", "PAT01", "SYSBP", 1, 1, 1L,
    "TEST01", "PAT01", "SYSBP", 1, 2, 2L,
    "TEST01", "PAT01", "SYSBP", 2, 1, 3L,
    "TEST01", "PAT02", "DIABP", 1, 1, 1L,
    "TEST01", "PAT02", "DIABP", 2, 1, 2L,
    "TEST01", "PAT02", "DIABP", 2, 2, 3L
  )

  actual_output <- derive_var_obs_number(
    input,
    by_vars = exprs(USUBJID, VSTESTCD),
    order = exprs(VISITNUM, !is.na(VSTPTNUM), VSTPTNUM)
  )

  expect_dfs_equal(
    base = expected_output,
    compare = actual_output,
    keys = c("USUBJID", "VSTESTCD", "VISITNUM", "VSTPTNUM")
  )
})


## Test 3: create customized new variable ASEQ1 ----
test_that("derive_var_obs_number Test 3: create customized new variable ASEQ1", {
  input <- tibble::tribble(
    ~STUDYID, ~USUBJID, ~VSTESTCD, ~VISITNUM, ~VSTPTNUM,
    "TEST01", "PAT01", "DIABP", 1, 1,
    "TEST01", "PAT01", "PULSE", 1, 1,
    "TEST01", "PAT01", "PULSE", 1, 2,
    "TEST01", "PAT01", "PULSE", 2, 1,
    "TEST01", "PAT01", "PULSE", 2, 2,
    "TEST01", "PAT01", "SYSBP", 1, 1,
    "TEST01", "PAT01", "SYSBP", 1, 2,
    "TEST01", "PAT01", "SYSBP", 2, 1,
    "TEST01", "PAT02", "DIABP", 1, 1,
    "TEST01", "PAT02", "DIABP", 2, 1,
    "TEST01", "PAT02", "DIABP", 2, 2
  )

  expected_output <- tibble::tribble(
    ~STUDYID, ~USUBJID, ~VSTESTCD, ~VISITNUM, ~VSTPTNUM, ~ASEQ1,
    "TEST01", "PAT01", "DIABP", 1, 1, 1L,
    "TEST01", "PAT01", "PULSE", 1, 1, 1L,
    "TEST01", "PAT01", "PULSE", 1, 2, 2L,
    "TEST01", "PAT01", "PULSE", 2, 1, 3L,
    "TEST01", "PAT01", "PULSE", 2, 2, 4L,
    "TEST01", "PAT01", "SYSBP", 1, 1, 1L,
    "TEST01", "PAT01", "SYSBP", 1, 2, 2L,
    "TEST01", "PAT01", "SYSBP", 2, 1, 3L,
    "TEST01", "PAT02", "DIABP", 1, 1, 1L,
    "TEST01", "PAT02", "DIABP", 2, 1, 2L,
    "TEST01", "PAT02", "DIABP", 2, 2, 3L
  )

  actual_output <- derive_var_obs_number(input,
    by_vars = exprs(USUBJID, VSTESTCD),
    order = exprs(VISITNUM, VSTPTNUM),
    new_var = ASEQ1
  )

  expect_dfs_equal(
    base = expected_output,
    compare = actual_output,
    keys = c("USUBJID", "VSTESTCD", "VISITNUM", "VSTPTNUM")
  )
})

## Test 4: expected condition for `check_type`  ----
test_that("derive_var_obs_number Test 4: expected condition for `check_type` ", {
  input <- tibble::tribble(
    ~STUDYID, ~USUBJID, ~VSTESTCD, ~VISITNUM, ~VSTPTNUM,
    "TEST01", "PAT01", "DIABP", 1, 1,
    "TEST01", "PAT01", "DIABP", 1, 1,
    "TEST01", "PAT01", "DIABP", 2, 1,
    "TEST01", "PAT01", "DIABP", 2, 2
  )

  expect_warning(
    derive_var_obs_number(
      input,
      by_vars = exprs(USUBJID, VSTESTCD),
      order = exprs(VISITNUM, VSTPTNUM),
      check_type = "warning"
    ),
    "duplicate records"
  )

  expect_error(
    derive_var_obs_number(
      input,
      by_vars = exprs(USUBJID, VSTESTCD),
      order = exprs(VISITNUM, VSTPTNUM),
      check_type = "error"
    ),
    "duplicate records"
  )
})
