context("test-derive_merged_vars")


test_that("variable is added from the first observation in each group", {
  adsl <- tibble::tibble(STUDYID = "STUDY", USUBJID = 1:3)
  adlb <- tibble::tribble(
    ~USUBJID, ~AVISITN, ~AVAL,
    1L, 1, 12,
    1L, 3, 9,
    2L, 2, 42,
    3L, 3, 14,
    3L, 3, 10)

  expected_output <- mutate(adsl, FIRSTVISN = c(1, 2, 3), FIRSTAVAL = c(12, 42, 10))

  actual_output <- derive_merged_vars(adsl,
                                      dataset_add = adlb,
                                      new_vars = exprs(FIRSTVISN = AVISITN, FIRSTAVAL = AVAL),
                                      filter_order = exprs(AVISITN, AVAL),
                                      filter_mode = "first",
                                      by_vars = exprs(USUBJID))

  expect_dfs_equal(base = expected_output,
                   compare = actual_output,
                   keys = c("USUBJID"))
})

test_that("filter_add parameter works, all observations from dataset are kept", {
  adsl <- tibble::tibble(STUDYID = "STUDY", USUBJID = 1:3)
  adlb <- tibble::tribble(
    ~USUBJID, ~AVISITN, ~AVAL,
    1L, 1, 12,
    1L, 3, 9,
    2L, 2, 42,
    3L, 3, 14,
    3L, 3, 10,
    4L, 1, 1)

  expected_output <- mutate(adsl, FIRSTAVAL = c(9, NA, 14))

  actual_output <- derive_merged_vars(
    adsl,
    dataset_add = adlb,
    filter_add = exprs(AVISITN == 3),
    new_vars = exprs(FIRSTAVAL = AVAL),
    filter_order = exprs(desc(AVAL)),
    filter_mode = "first",
    by_vars = exprs(USUBJID)
  )

  expect_dfs_equal(base = expected_output,
                   compare = actual_output,
                   keys = c("USUBJID"))
})
