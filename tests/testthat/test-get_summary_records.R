## Test 1: Summarize the average of the triplicate ECG interval values (AVAL) ----
test_that("get_summary_records Test 1: Summarize the average of the triplicate ECG interval values (AVAL)", {
  input <- tibble::tribble(
    ~USUBJID,   ~EGSEQ, ~PARAM,             ~AVISIT,    ~EGDTC,             ~AVAL, ~TRTA,
    "XYZ-1001", 1,      "QTcF Int. (msec)", "Baseline", "2016-02-24T07:50", 385,   NA_character_,
    "XYZ-1001", 2,      "QTcF Int. (msec)", "Baseline", "2016-02-24T07:52", 399,   NA_character_,
    "XYZ-1001", 3,      "QTcF Int. (msec)", "Baseline", "2016-02-24T07:56", 396,   NA_character_,
    "XYZ-1001", 4,      "QTcF Int. (msec)", "Visit 2",  "2016-03-08T09:45", 384,   "Placebo",
    "XYZ-1001", 5,      "QTcF Int. (msec)", "Visit 2",  "2016-03-08T09:48", 393,   "Placebo",
    "XYZ-1001", 6,      "QTcF Int. (msec)", "Visit 2",  "2016-03-08T09:51", 388,   "Placebo",
    "XYZ-1001", 7,      "QTcF Int. (msec)", "Visit 3",  "2016-03-22T10:45", 385,   "Placebo",
    "XYZ-1001", 8,      "QTcF Int. (msec)", "Visit 3",  "2016-03-22T10:48", 394,   "Placebo",
    "XYZ-1001", 9,      "QTcF Int. (msec)", "Visit 3",  "2016-03-22T10:51", 402,   "Placebo",
    "XYZ-1002", 1,      "QTcF Int. (msec)", "Baseline", "2016-02-22T07:58", 399,   NA_character_,
    "XYZ-1002", 2,      "QTcF Int. (msec)", "Baseline", "2016-02-22T07:58", 410,   NA_character_,
    "XYZ-1002", 3,      "QTcF Int. (msec)", "Baseline", "2016-02-22T08:01", 392,   NA_character_,
    "XYZ-1002", 4,      "QTcF Int. (msec)", "Visit 2",  "2016-03-06T09:50", 401,   "Active 20mg",
    "XYZ-1002", 5,      "QTcF Int. (msec)", "Visit 2",  "2016-03-06T09:53", 407,   "Active 20mg",
    "XYZ-1002", 6,      "QTcF Int. (msec)", "Visit 2",  "2016-03-06T09:56", 400,   "Active 20mg",
    "XYZ-1002", 7,      "QTcF Int. (msec)", "Visit 3",  "2016-03-24T10:50", 412,   "Active 20mg",
    "XYZ-1002", 8,      "QTcF Int. (msec)", "Visit 3",  "2016-03-24T10:53", 414,   "Active 20mg",
    "XYZ-1002", 9,      "QTcF Int. (msec)", "Visit 3",  "2016-03-24T10:56", 402,   "Active 20mg"
  )

  actual_output <- input %>%
    get_summary_records(
      by_vars = exprs(USUBJID, PARAM, AVISIT),
      set_values_to = exprs(
        AVAL = mean(AVAL, na.rm = TRUE),
        DTYPE = "AVERAGE"
      )
    ) %>%
    dplyr::mutate(AVAL = round(AVAL))

  expected_output <- tibble::tribble(
    ~USUBJID,   ~PARAM,             ~AVISIT,    ~AVAL, ~DTYPE,
    "XYZ-1001", "QTcF Int. (msec)", "Baseline",  393, "AVERAGE",
    "XYZ-1001", "QTcF Int. (msec)", "Visit 2",   388, "AVERAGE",
    "XYZ-1001", "QTcF Int. (msec)", "Visit 3",   394, "AVERAGE",
    "XYZ-1002", "QTcF Int. (msec)", "Baseline",  400, "AVERAGE",
    "XYZ-1002", "QTcF Int. (msec)", "Visit 2",   403, "AVERAGE",
    "XYZ-1002", "QTcF Int. (msec)", "Visit 3",   409, "AVERAGE"
  )

  expect_dfs_equal(
    base = expected_output,
    compare = actual_output,
    keys = c("USUBJID", "PARAM", "AVISIT")
  )
})
