context("test-derive_var_trtedtm")

test_that("TRTEDTM variable is added", {
  adsl <- tibble::tibble(STUDYID = "STUDY", USUBJID = 1:3)
  ex <- tibble::tribble(
    ~USUBJID, ~EXENDTC, ~EXSEQ, ~EXDOSE, ~EXTRT,
    1L, "2020-01-01", 1, 12, "ACTIVE",
    1L, "2020-02-03", 2, 9, "ACTIVE",
    2L, "2020-01-02", 1, 0, "PLACEBO",
    3L, "2020-03-13", 1, 14, "ACTIVE",
    3L, "2020-03-21", 2, 0, "ACTIVE")

  expected_output <- mutate(adsl, TRTEDTM = c(ymd_hms("2020-02-03T23:59:59"),
                                              ymd_hms("2020-01-02T23:59:59"),
                                              ymd_hms("2020-03-13T23:59:59")))

  actual_output <- derive_var_trtedtm(adsl,
                                      dataset_ex = ex)
  expect_dfs_equal(base = expected_output,
                   compare = actual_output,
                   keys = c("USUBJID"))
})
