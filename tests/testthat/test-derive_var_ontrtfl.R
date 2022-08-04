test_that("`target` is set to NA when ` start_date` < `ref_start_date`", {
  library(tibble)

  input <- tribble(
    ~STUDYID, ~USUBJID, ~ADT, ~TRTSDT,
    "TEST01", "PAT01", as.Date("2021-01-01"), as.Date("2021-01-02")
  )
  expected_output <- tribble(
    ~STUDYID, ~USUBJID, ~ADT, ~TRTSDT, ~ONTRTFL,
    "TEST01", "PAT01", as.Date("2021-01-01"), as.Date("2021-01-02"), as.character(NA)
  )

  actual_output <- derive_var_ontrtfl(
    input,
    new_var = ONTRTFL,
    start_date = ADT,
    ref_start_date = TRTSDT
  )

  expect_dfs_equal(
    expected_output,
    actual_output,
    keys = c("STUDYID", "USUBJID", "ADT")
  )
})

test_that("`target` is set to NA when `ref_start_date` is NA", {
  library(tibble)

  input <- tribble(
    ~STUDYID, ~USUBJID, ~ADT, ~TRTSDT,
    "TEST01", "PAT01", as.Date("2021-01-01"), as.Date(NA)
  )
  expected_output <- tribble(
    ~STUDYID, ~USUBJID, ~ADT, ~TRTSDT, ~ONTRTFL,
    "TEST01", "PAT01", as.Date("2021-01-01"), as.Date(NA), as.character(NA)
  )

  actual_output <- derive_var_ontrtfl(
    input,
    new_var = ONTRTFL,
    start_date = ADT,
    ref_start_date = TRTSDT
  )

  expect_dfs_equal(
    expected_output,
    actual_output,
    keys = c("STUDYID", "USUBJID", "ADT")
  )
})

test_that("`target` is set to `Y` when ` start_date` is NA", {
  library(tibble)

  input <- tribble(
    ~STUDYID, ~USUBJID, ~ADT, ~TRTSDT,
    "TEST01", "PAT01", as.Date(NA), as.Date("2020-01-01"),
    "TEST01", "PAT02", as.Date(NA), as.Date(NA)
  )
  expected_output <- tribble(
    ~STUDYID, ~USUBJID, ~ADT, ~TRTSDT, ~ONTRTFL,
    "TEST01", "PAT01", as.Date(NA), as.Date("2020-01-01"), "Y",
    "TEST01", "PAT02", as.Date(NA), as.Date(NA), as.character(NA)
  )

  actual_output <- derive_var_ontrtfl(
    input,
    new_var = ONTRTFL,
    start_date = ADT,
    ref_start_date = TRTSDT
  )

  expect_dfs_equal(
    expected_output,
    actual_output,
    keys = c("STUDYID", "USUBJID", "ADT")
  )
})

test_that("`target` is set to `Y` when ` start_date` >= `ref_start_date` and
            `ref_end_date` and `filter_pre_timepoint` are not specified", {
  library(tibble)

  input <- tribble(
    ~STUDYID, ~USUBJID, ~ADT, ~TRTSDT,
    "TEST01", "PAT01", as.Date("2020-01-01"), as.Date("2020-01-01"),
    "TEST01", "PAT02", as.Date("2020-01-02"), as.Date("2020-01-01")
  )
  expected_output <- tribble(
    ~STUDYID, ~USUBJID, ~ADT, ~TRTSDT, ~ONTRTFL,
    "TEST01", "PAT01", as.Date("2020-01-01"), as.Date("2020-01-01"), "Y",
    "TEST01", "PAT02", as.Date("2020-01-02"), as.Date("2020-01-01"), "Y"
  )

  actual_output <- derive_var_ontrtfl(
    input,
    new_var = ONTRTFL,
    start_date = ADT,
    ref_start_date = TRTSDT
  )

  expect_dfs_equal(
    expected_output,
    actual_output,
    keys = c("STUDYID", "USUBJID", "ADT")
  )
})

test_that("`target` is set to 'Y' when `filter_pre_timepoint` is not 'PRE' and
            ` start_date` = `ref_start_date` and `ref_end_date` is not specified", {
  library(tibble)

  input <- tribble(
    ~STUDYID, ~USUBJID, ~ADT, ~TRTSDT, ~TPT,
    "TEST01", "PAT01", as.Date("2020-01-01"), as.Date("2020-01-01"), "PRE",
    "TEST01", "PAT02", as.Date("2020-01-01"), as.Date("2020-01-01"), "POST"
  )
  expected_output <- tribble(
    ~STUDYID, ~USUBJID, ~ADT, ~TRTSDT, ~TPT, ~ONTRTFL,
    "TEST01", "PAT01", as.Date("2020-01-01"), as.Date("2020-01-01"), "PRE", NA,
    "TEST01", "PAT02", as.Date("2020-01-01"), as.Date("2020-01-01"), "POST", "Y"
  )

  actual_output <- derive_var_ontrtfl(
    input,
    new_var = ONTRTFL,
    start_date = ADT,
    ref_start_date = TRTSDT,
    filter_pre_timepoint = TPT == "PRE"
  )

  expect_dfs_equal(
    expected_output,
    actual_output,
    keys = c("STUDYID", "USUBJID", "ADT")
  )
})

test_that("`target` is set to `Y` when ` start_date` >= `ref_start_date` and ` start_date` <=
            `ref_end_date` and no `ref_end_window` is specified, otherwise NA", {
  library(tibble)

  input <- tribble(
    ~STUDYID, ~USUBJID, ~ADT, ~TRTSDT, ~TRTEDT,
    "TEST01", "PAT01", as.Date("2019-12-13"), as.Date("2020-01-01"), as.Date("2020-02-01"),
    "TEST01", "PAT02", as.Date("2020-01-01"), as.Date("2020-01-01"), as.Date("2020-02-01"),
    "TEST01", "PAT03", as.Date("2020-01-02"), as.Date("2020-01-01"), as.Date("2020-02-01"),
    "TEST01", "PAT04", as.Date("2020-02-01"), as.Date("2020-01-01"), as.Date("2020-02-01"),
    "TEST01", "PAT05", as.Date("2020-02-02"), as.Date("2020-01-01"), as.Date("2020-02-01")
  )
  expected_output <- tribble(
    ~STUDYID, ~USUBJID, ~ADT, ~TRTSDT, ~TRTEDT, ~ONTRTFL, # nolint
    "TEST01", "PAT01", as.Date("2019-12-13"), as.Date("2020-01-01"), as.Date("2020-02-01"), NA,
    "TEST01", "PAT02", as.Date("2020-01-01"), as.Date("2020-01-01"), as.Date("2020-02-01"), "Y",
    "TEST01", "PAT03", as.Date("2020-01-02"), as.Date("2020-01-01"), as.Date("2020-02-01"), "Y",
    "TEST01", "PAT04", as.Date("2020-02-01"), as.Date("2020-01-01"), as.Date("2020-02-01"), "Y",
    "TEST01", "PAT05", as.Date("2020-02-02"), as.Date("2020-01-01"), as.Date("2020-02-01"), NA
  )

  actual_output <- derive_var_ontrtfl(
    input,
    new_var = ONTRTFL,
    start_date = ADT,
    ref_start_date = TRTSDT,
    ref_end_date = TRTEDT
  )

  expect_dfs_equal(
    expected_output,
    actual_output,
    keys = c("STUDYID", "USUBJID", "ADT")
  )
})

test_that("`target` is set to `Y` when ` start_date` >= `ref_start_date` and ` start_date` <=
            `ref_end_date` + `ref_end_window`", {
  library(tibble)

  input <- tribble(
    ~STUDYID, ~USUBJID, ~ADT, ~TRTSDT, ~TRTEDT,
    "TEST01", "PAT01", as.Date("2020-02-01"), as.Date("2020-01-01"), as.Date("2020-02-01"),
    "TEST01", "PAT02", as.Date("2020-02-05"), as.Date("2020-01-01"), as.Date("2020-02-01"),
    "TEST01", "PAT03", as.Date("2020-02-10"), as.Date("2020-01-01"), as.Date("2020-02-01")
  )
  expected_output <- tribble(
    ~STUDYID, ~USUBJID, ~ADT, ~TRTSDT, ~TRTEDT, ~ONTRTFL, # nolint
    "TEST01", "PAT01", as.Date("2020-02-01"), as.Date("2020-01-01"), as.Date("2020-02-01"), "Y",
    "TEST01", "PAT02", as.Date("2020-02-05"), as.Date("2020-01-01"), as.Date("2020-02-01"), "Y",
    "TEST01", "PAT03", as.Date("2020-02-10"), as.Date("2020-01-01"), as.Date("2020-02-01"), NA
  )

  actual_output <- derive_var_ontrtfl(
    input,
    new_var = ONTRTFL,
    start_date = ADT,
    ref_start_date = TRTSDT,
    ref_end_date = TRTEDT,
    ref_end_window = 5
  )

  expect_dfs_equal(
    expected_output,
    actual_output,
    keys = c("STUDYID", "USUBJID", "ADT")
  )
})



test_that("`target` is set to  NA when `end_date`<`ref_start_date`
          regradless of start_date being NA", {
  library(tibble)

  input <- tribble(
    ~USUBJID, ~ASTDT, ~TRTSDT, ~TRTEDT, ~AENDT,
    "PAT01", ymd("2019-04-30"), ymd("2020-01-01"), ymd("2020-03-01"), ymd("2019-03-15"),
    "PAT01", NA, ymd("2020-01-01"), ymd("2020-03-01"), ymd("2019-03-15"),
  )
  expected_output <- tribble(
    ~USUBJID, ~ASTDT, ~TRTSDT, ~TRTEDT, ~AENDT, ~ONTRTFL,
    "PAT01", ymd("2019-04-30"), ymd("2020-01-01"), ymd("2020-03-01"), ymd("2019-03-15"), NA_character_, # nolint
    "PAT01", NA, ymd("2020-01-01"), ymd("2020-03-01"), ymd("2019-03-15"), NA_character_,
  )

  actual_output <- derive_var_ontrtfl(
    input,
    new_var = ONTRTFL,
    start_date = ASTDT,
    end_date = AENDT,
    ref_start_date = TRTSDT,
    ref_end_date = TRTEDT,
    ref_end_window = 60
  )

  expect_dfs_equal(
    expected_output,
    actual_output,
    keys = c("USUBJID", "ASTDT")
  )
})

test_that("`target` is set to  `Y` when `end_date`>`ref_start_date` when `start_date` is missing", {
  library(tibble)

  input <- tribble(
    ~STUDYID, ~USUBJID, ~ASTDT, ~TRTSDT, ~TRTEDT, ~AENDT,
    "TEST01", "PAT01", NA, ymd("2020-01-01"), ymd("2020-03-01"), ymd("2021-03-15"),
  )
  expected_output <- tribble(
    ~STUDYID, ~USUBJID, ~ASTDT, ~TRTSDT, ~TRTEDT, ~AENDT, ~ONTRTFL,
    "TEST01", "PAT01", NA, ymd("2020-01-01"), ymd("2020-03-01"), ymd("2021-03-15"), "Y"
  )

  actual_output <- derive_var_ontrtfl(
    input,
    new_var = ONTRTFL,
    start_date = ASTDT,
    end_date = AENDT,
    ref_start_date = TRTSDT,
    ref_end_date = TRTEDT,
    ref_end_window = 60
  )

  expect_dfs_equal(
    expected_output,
    actual_output,
    keys = c("STUDYID", "USUBJID", "ASTDT")
  )
})

test_that("`target` is set to  NA when `end_date` is missing and
            `start_date` is before `ref_start_date`  a la Roche", {
  library(tibble)

  input <- tribble(
    ~STUDYID, ~USUBJID, ~ASTDT, ~TRTSDT, ~TRTEDT, ~AENDT,
    "TEST01", "PAT01", ymd("2019-04-30"), ymd("2020-01-01"), ymd("2020-03-01"), NA,
  )
  expected_output <- tribble(
    ~STUDYID, ~USUBJID, ~ASTDT, ~TRTSDT, ~TRTEDT, ~AENDT, ~ONTRTFL,
    "TEST01", "PAT01", ymd("2019-04-30"), ymd("2020-01-01"), ymd("2020-03-01"), NA, NA_character_,
  )

  actual_output <- derive_var_ontrtfl(
    input,
    new_var = ONTRTFL,
    start_date = ASTDT,
    end_date = AENDT,
    ref_start_date = TRTSDT,
    ref_end_date = TRTEDT,
    ref_end_window = 60
  )

  expect_dfs_equal(
    expected_output,
    actual_output,
    keys = c("STUDYID", "USUBJID", "ASTDT")
  )
})

test_that("`target` is set to  Y when `end_date` is missing and
            `start_date` is before `ref_start_date`  a la GSK", {
  library(tibble)

  input <- tribble(
    ~STUDYID, ~USUBJID, ~ASTDT, ~TRTSDT, ~TRTEDT, ~AENDT,
    "TEST01", "PAT01", ymd("2019-04-30"), ymd("2020-01-01"), ymd("2020-03-01"), NA,
  )
  expected_output <- tribble(
    ~STUDYID, ~USUBJID, ~ASTDT, ~TRTSDT, ~TRTEDT, ~AENDT, ~ONTRTFL,
    "TEST01", "PAT01", ymd("2019-04-30"), ymd("2020-01-01"), ymd("2020-03-01"), NA, "Y",
  )

  actual_output <- derive_var_ontrtfl(
    input,
    new_var = ONTRTFL,
    start_date = ASTDT,
    end_date = AENDT,
    ref_start_date = TRTSDT,
    ref_end_date = TRTEDT,
    ref_end_window = 60,
    span_period = "Y"
  )

  expect_dfs_equal(
    expected_output,
    actual_output,
    keys = c("STUDYID", "USUBJID", "ASTDT")
  )
})

test_that("`target` is set to  Y when `end_date` is missing and
            `start_date` is before `ref_start_date`  a la GSK", {
  library(tibble)

  input <- tribble(
    ~STUDYID, ~USUBJID, ~ASTDT, ~TRTSDT, ~TRTEDT, ~AENDT,
    "TEST01", "PAT01", ymd("2019-04-30"), ymd("2020-01-01"), ymd("2020-03-01"), NA,
  )
  expected_output <- tribble(
    ~STUDYID, ~USUBJID, ~ASTDT, ~TRTSDT, ~TRTEDT, ~AENDT, ~ONTRTFL,
    "TEST01", "PAT01", ymd("2019-04-30"), ymd("2020-01-01"), ymd("2020-03-01"), NA, "Y",
  )

  actual_output <- derive_var_ontrtfl(
    input,
    new_var = ONTRTFL,
    start_date = ASTDT,
    end_date = AENDT,
    ref_start_date = TRTSDT,
    ref_end_date = TRTEDT,
    ref_end_window = 60,
    span_period = "Y"
  )

  expect_dfs_equal(
    expected_output,
    actual_output,
    keys = c("STUDYID", "USUBJID", "ASTDT")
  )
})

test_that("`target` is set to Y when `start_date` is before `ref_start_date` and
            `end_date` is before `ref_end_date` for Period 01", {
  library(tibble)

  input <- tribble(
    ~STUDYID, ~USUBJID, ~ASTDT, ~AP01SDT, ~AP01EDT, ~AENDT,
    "TEST01", "PAT01", ymd("2019-04-30"), ymd("2020-01-01"), ymd("2020-03-01"), ymd("2020-03-15")
  )

  expected_output <- tribble(
    ~STUDYID, ~USUBJID, ~ASTDT, ~AP01SDT, ~AP01EDT, ~AENDT, ~ONTR01FL,
    "TEST01", "PAT01", ymd("2019-04-30"), ymd("2020-01-01"),
    ymd("2020-03-01"), ymd("2020-03-15"), "Y",
  )

  actual_output <- derive_var_ontrtfl(
    input,
    new_var = ONTR01FL,
    start_date = ASTDT,
    end_date = AENDT,
    ref_start_date = AP01SDT,
    ref_end_date = AP01EDT,
    span_period = "Y"
  )

  expect_dfs_equal(
    expected_output,
    actual_output,
    keys = c("STUDYID", "USUBJID", "ASTDT")
  )
})
