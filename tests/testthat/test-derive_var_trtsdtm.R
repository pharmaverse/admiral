test_that("TRTSDTM variable is added", {
  adsl <- tibble::tibble(STUDYID = "STUDY", USUBJID = 1:3)
  ex <- tibble::tribble(
    ~USUBJID, ~EXSTDTC, ~EXSEQ, ~EXDOSE, ~EXTRT,
    1L, "2020-01-01", 1, 12, "ACTIVE",
    1L, "2020-02-03", 2, 9, "ACTIVE",
    2L, "2020-01-02", 1, 0, "PLACEBO",
    3L, "2020-03-13", 1, 14, "ACTIVE",
    3L, "2020-03-21", 2, 0, "ACTIVE"
  )

  expected_output <- mutate(
    adsl,
    TRTSDTM = ymd_hms(c(
      "2020-01-01T00:00:00",
      "2020-01-02T00:00:00",
      "2020-03-13T00:00:00"
    ))
  )

  actual_output <- derive_var_trtsdtm(adsl, dataset_ex = ex, subject_keys = vars(USUBJID))

  expect_dfs_equal(
    base = expected_output,
    compare = actual_output,
    keys = "USUBJID"
  )
})
