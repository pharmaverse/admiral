context("test-derive_var_ontrtfl")

test_that("`target` is set to NA when `date` < `ref_start_date`", {
  input <- tibble::tribble(
    ~STUDYID, ~USUBJID, ~ADT,                  ~TRTSDT,
    "TEST01", "PAT01",  as.Date("2021-01-01"), as.Date("2021-01-02")
  )
  expected_output <- tibble::tribble(
    ~STUDYID, ~USUBJID, ~ADT,                  ~TRTSDT,              ~ONTRTFL,
    "TEST01", "PAT01",  as.Date("2021-01-01"), as.Date("2021-01-02"), as.character(NA)
  )

  actual_output <- derive_var_ontrtfl(
    input,
    ADT,
    TRTSDT)

  expect_dfs_equal(
    expected_output,
    actual_output,
    keys = c("STUDYID", "USUBJID", "ADT")
  )
})

test_that("`target` is set to NA when `ref_start_date` is NA", {
  input <- tibble::tribble(
    ~STUDYID, ~USUBJID, ~ADT,                  ~TRTSDT,
    "TEST01", "PAT01",  as.Date("2021-01-01"), as.Date(NA)
  )
  expected_output <- tibble::tribble(
    ~STUDYID, ~USUBJID, ~ADT,                  ~TRTSDT,              ~ONTRTFL,
    "TEST01", "PAT01",  as.Date("2021-01-01"), as.Date(NA),     as.character(NA)
  )

  actual_output <- derive_var_ontrtfl(
    input,
    ADT,
    TRTSDT)

  expect_dfs_equal(
    expected_output,
    actual_output,
    keys = c("STUDYID", "USUBJID", "ADT")
  )
})

test_that("`target` is set to `Y` when `date` is NA", {
  input <- tibble::tribble(
    ~STUDYID, ~USUBJID, ~ADT,                  ~TRTSDT,
    "TEST01", "PAT01",  as.Date(NA), as.Date("2020-01-01"),
    "TEST01", "PAT02",  as.Date(NA), as.Date(NA)
  )
  expected_output <- tibble::tribble(
    ~STUDYID, ~USUBJID, ~ADT,                  ~TRTSDT,               ~ONTRTFL,
    "TEST01", "PAT01",  as.Date(NA),           as.Date("2020-01-01"), "Y",
    "TEST01", "PAT02",  as.Date(NA),           as.Date(NA),           as.character(NA)
  )

  actual_output <- derive_var_ontrtfl(
    input,
    ADT,
    TRTSDT)

  expect_dfs_equal(
    expected_output,
    actual_output,
    keys = c("STUDYID", "USUBJID", "ADT")
  )
})

test_that("`target` is set to `Y` when `date` >= `ref_start_date` and
          `ref_end_date` and `filter_pre_timepoint` are not specified", {
  input <- tibble::tribble(
    ~STUDYID, ~USUBJID, ~ADT,                  ~TRTSDT,
    "TEST01", "PAT01",  as.Date("2020-01-01"), as.Date("2020-01-01"),
    "TEST01", "PAT02",  as.Date("2020-01-02"), as.Date("2020-01-01")
  )
  expected_output <- tibble::tribble(
    ~STUDYID, ~USUBJID, ~ADT,                  ~TRTSDT,               ~ONTRTFL,
    "TEST01", "PAT01",  as.Date("2020-01-01"), as.Date("2020-01-01"), "Y",
    "TEST01", "PAT02",  as.Date("2020-01-02"), as.Date("2020-01-01"), "Y"
  )

  actual_output <- derive_var_ontrtfl(
    input,
    ADT,
    TRTSDT)

  expect_dfs_equal(
    expected_output,
    actual_output,
    keys = c("STUDYID", "USUBJID", "ADT")
  )
})

test_that("`target` is set to `Y` when `filter_pre_timepoint` is not PRE and`date` = `ref_start_date`
          and `ref_end_date` is not specified", {
  input <- tibble::tribble(
    ~STUDYID, ~USUBJID, ~ADT,                  ~TRTSDT,               ~TPT,
    "TEST01", "PAT01",  as.Date("2020-01-01"), as.Date("2020-01-01"), "PRE",
    "TEST01", "PAT02",  as.Date("2020-01-01"), as.Date("2020-01-01"), "POST"
  )
  expected_output <- tibble::tribble(
    ~STUDYID, ~USUBJID, ~ADT,                  ~TRTSDT,               ~TPT,   ~ONTRTFL,
    "TEST01", "PAT01",  as.Date("2020-01-01"), as.Date("2020-01-01"), "PRE",  NA,
    "TEST01", "PAT02",  as.Date("2020-01-01"), as.Date("2020-01-01"), "POST", "Y"
  )

  actual_output <- derive_var_ontrtfl(
    input,
    ADT,
    TRTSDT,
    filter_pre_timepoint = exprs(TPT == "PRE"))

  expect_dfs_equal(
    expected_output,
    actual_output,
    keys = c("STUDYID", "USUBJID", "ADT")
  )
})

test_that("`target` is set to `Y` when `date` >= `ref_start_date` and `date` <=
          `ref_end_date` and no `ref_end_window` is specified, otherwise NA", {
  input <- tibble::tribble(
    ~STUDYID, ~USUBJID, ~ADT,                  ~TRTSDT,               ~TRTEDT,
    "TEST01", "PAT01",  as.Date("2019-12-13"), as.Date("2020-01-01"), as.Date("2020-02-01"),
    "TEST01", "PAT02",  as.Date("2020-01-01"), as.Date("2020-01-01"), as.Date("2020-02-01"),
    "TEST01", "PAT03",  as.Date("2020-01-02"), as.Date("2020-01-01"), as.Date("2020-02-01"),
    "TEST01", "PAT04",  as.Date("2020-02-01"), as.Date("2020-01-01"), as.Date("2020-02-01"),
    "TEST01", "PAT05",  as.Date("2020-02-02"), as.Date("2020-01-01"), as.Date("2020-02-01")
  )
  expected_output <- tibble::tribble(
    ~STUDYID, ~USUBJID, ~ADT,                  ~TRTSDT,               ~TRTEDT,              ~ONTRTFL,
    "TEST01", "PAT01",  as.Date("2019-12-13"), as.Date("2020-01-01"), as.Date("2020-02-01"), NA,
    "TEST01", "PAT02",  as.Date("2020-01-01"), as.Date("2020-01-01"), as.Date("2020-02-01"), "Y",
    "TEST01", "PAT03",  as.Date("2020-01-02"), as.Date("2020-01-01"), as.Date("2020-02-01"), "Y",
    "TEST01", "PAT04",  as.Date("2020-02-01"), as.Date("2020-01-01"), as.Date("2020-02-01"), "Y",
    "TEST01", "PAT05",  as.Date("2020-02-02"), as.Date("2020-01-01"), as.Date("2020-02-01"), NA
  )

  actual_output <- derive_var_ontrtfl(
    input,
    ADT,
    TRTSDT,
    ref_end_date = TRTEDT)

  expect_dfs_equal(
    expected_output,
    actual_output,
    keys = c("STUDYID", "USUBJID", "ADT")
  )
})

test_that("`target` is set to `Y` when `date` >= `ref_start_date` and `date` <=
          `ref_end_date` + 'ref_end_window`", {
  input <- tibble::tribble(
    ~STUDYID, ~USUBJID, ~ADT,                  ~TRTSDT,               ~TRTEDT,
    "TEST01", "PAT01",  as.Date("2020-02-01"), as.Date("2020-01-01"), as.Date("2020-02-01"),
    "TEST01", "PAT02",  as.Date("2020-02-05"), as.Date("2020-01-01"), as.Date("2020-02-01"),
    "TEST01", "PAT03",  as.Date("2020-02-10"), as.Date("2020-01-01"), as.Date("2020-02-01")
  )
  expected_output <- tibble::tribble(
    ~STUDYID, ~USUBJID, ~ADT,                  ~TRTSDT,               ~TRTEDT,              ~ONTRTFL,
    "TEST01", "PAT01",  as.Date("2020-02-01"), as.Date("2020-01-01"), as.Date("2020-02-01"), "Y",
    "TEST01", "PAT02",  as.Date("2020-02-05"), as.Date("2020-01-01"), as.Date("2020-02-01"), "Y",
    "TEST01", "PAT03",  as.Date("2020-02-10"), as.Date("2020-01-01"), as.Date("2020-02-01"), NA
  )

  actual_output <- derive_var_ontrtfl(
    input,
    ADT,
    TRTSDT,
    ref_end_date = TRTEDT,
    ref_end_window = 5)

  expect_dfs_equal(
    expected_output,
    actual_output,
    keys = c("STUDYID", "USUBJID", "ADT")
  )
})
