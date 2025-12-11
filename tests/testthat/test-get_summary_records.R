## Test 1: deprecation message if function is called ----
test_that("get_summary_records Test 1: deprecation message if function is called", {
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

  expect_snapshot(
    df <- input %>%
      get_summary_records(
        by_vars = exprs(USUBJID, PARAM, AVISIT),
        set_values_to = exprs(
          AVAL = mean(AVAL, na.rm = TRUE),
          DTYPE = "AVERAGE"
        )
      ) %>%
      dplyr::mutate(AVAL = round(AVAL))
  )
})

## Test 2: Summarize average of triplicate ECG interval values ----
test_that("get_summary_records Test 2: Summarize average of triplicate ECG interval values", {
  # Suppress lifecycle messages within the test environment
  withr::local_options(list(lifecycle_verbosity = "quiet"))

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
    ~USUBJID, ~PARAM, ~AVISIT, ~AVAL, ~DTYPE,
    "XYZ-1001", "QTcF Int. (msec)", "Baseline", 393, "AVERAGE",
    "XYZ-1001", "QTcF Int. (msec)", "Visit 2", 388, "AVERAGE",
    "XYZ-1001", "QTcF Int. (msec)", "Visit 3", 394, "AVERAGE",
    "XYZ-1002", "QTcF Int. (msec)", "Baseline", 400, "AVERAGE",
    "XYZ-1002", "QTcF Int. (msec)", "Visit 2", 403, "AVERAGE",
    "XYZ-1002", "QTcF Int. (msec)", "Visit 3", 409, "AVERAGE"
  )

  expect_dfs_equal(
    base = expected_output,
    compare = actual_output,
    keys = c("USUBJID", "PARAM", "AVISIT")
  )
})

## Test 3: Derive more than one summary variable ----
test_that("get_summary_records Test 3: Derive more than one summary variable", {
  # Suppress lifecycle messages within the test environment
  withr::local_options(list(lifecycle_verbosity = "quiet"))

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
        AVAL = mean(AVAL),
        ASTDTM = min(convert_dtc_to_dtm(EGDTC)),
        AENDTM = max(convert_dtc_to_dtm(EGDTC)),
        DTYPE = "AVERAGE"
      )
    ) %>%
    dplyr::mutate(AVAL = round(AVAL))

  expected_output <- tibble::tribble(
    ~USUBJID, ~PARAM, ~AVISIT, ~AVAL, ~ASTDTM, ~AENDTM, ~DTYPE,
    # nolint start: line_len
    "XYZ-1001", "QTcF Int. (msec)", "Baseline", 393., "2016-02-24 07:50:00", "2016-02-24 07:56 :00", "AVERAGE",
    "XYZ-1001", "QTcF Int. (msec)", "Visit 2", 388., "2016-03-08 09:45:00", "2016-03-08 09:51:00", "AVERAGE",
    "XYZ-1001", "QTcF Int. (msec)", "Visit 3", 394., "2016-03-22 10:45:00", "2016-03-22 10:51:00", "AVERAGE",
    "XYZ-1002", "QTcF Int. (msec)", "Baseline", 400., "2016-02-22 07:58:00", "2016-02-22 08:01:00", "AVERAGE",
    "XYZ-1002", "QTcF Int. (msec)", "Visit 2", 403., "2016-03-06 09:50:00", "2016-03-06 09:56:00", "AVERAGE",
    "XYZ-1002", "QTcF Int. (msec)", "Visit 3", 409., "2016-03-24 10:50:00", "2016-03-24 10:56:00", "AVERAGE"
    # nolint end
  ) %>%
    dplyr::mutate(
      ASTDTM = as.POSIXct(ASTDTM, tz = "UTC"),
      AENDTM = as.POSIXct(AENDTM, tz = "UTC")
    )

  expect_dfs_equal(
    base = expected_output,
    compare = actual_output,
    keys = c("USUBJID", "PARAM", "AVISIT")
  )
})

## Test 4: Compute avg AVAL only if >2 records within by group ----
test_that("get_summary_records Test 4: Compute avg AVAL only if >2 records within by group", {
  # Suppress lifecycle messages within the test environment
  withr::local_options(list(lifecycle_verbosity = "quiet"))

  input <- tibble::tribble(
    ~USUBJID,   ~EGSEQ, ~PARAM,             ~AVISIT,    ~EGDTC,             ~AVAL, ~TRTA,
    "XYZ-1001", 1,      "QTcF Int. (msec)", "Baseline", "2016-02-24T07:50", 385,   NA_character_,
    "XYZ-1001", 2,      "QTcF Int. (msec)", "Baseline", "2016-02-24T07:52", 399,   NA_character_,
    "XYZ-1001", 3,      "QTcF Int. (msec)", "Baseline", "2016-02-24T07:56", 396,   NA_character_,
    "XYZ-1001", 4,      "QTcF Int. (msec)", "Visit 2",  "2016-03-08T09:48", 393,   "Placebo",
    "XYZ-1001", 5,      "QTcF Int. (msec)", "Visit 2",  "2016-03-08T09:51", 388,   "Placebo",
    "XYZ-1001", 6,      "QTcF Int. (msec)", "Visit 3",  "2016-03-22T10:48", 394,   "Placebo",
    "XYZ-1001", 7,      "QTcF Int. (msec)", "Visit 3",  "2016-03-22T10:51", 402,   "Placebo",
    "XYZ-1002", 1,      "QTcF Int. (msec)", "Baseline", "2016-02-22T07:58", 399,   NA_character_,
    "XYZ-1002", 2,      "QTcF Int. (msec)", "Baseline", "2016-02-22T07:58", 410,   NA_character_,
    "XYZ-1002", 3,      "QTcF Int. (msec)", "Baseline", "2016-02-22T08:01", 392,   NA_character_,
    "XYZ-1002", 4,      "QTcF Int. (msec)", "Visit 2",  "2016-03-06T09:53", 407,   "Active 20mg",
    "XYZ-1002", 5,      "QTcF Int. (msec)", "Visit 2",  "2016-03-06T09:56", 400,   "Active 20mg",
    "XYZ-1002", 6,      "QTcF Int. (msec)", "Visit 3",  "2016-03-24T10:53", 414,   "Active 20mg",
    "XYZ-1002", 7,      "QTcF Int. (msec)", "Visit 3",  "2016-03-24T10:56", 402,   "Active 20mg"
  )

  actual_output <- input %>%
    get_summary_records(
      by_vars = exprs(USUBJID, PARAM, AVISIT),
      filter = n() > 2,
      set_values_to = exprs(
        AVAL = mean(AVAL, na.rm = TRUE),
        DTYPE = "AVERAGE"
      )
    ) %>%
    dplyr::mutate(AVAL = round(AVAL))

  expected_output <- tibble::tribble(
    ~USUBJID,   ~PARAM,             ~AVISIT,    ~AVAL, ~DTYPE,
    "XYZ-1001", "QTcF Int. (msec)", "Baseline", 393,   "AVERAGE",
    "XYZ-1002", "QTcF Int. (msec)", "Baseline", 400,   "AVERAGE"
  )

  expect_dfs_equal(
    base = expected_output,
    compare = actual_output,
    keys = c("USUBJID", "PARAM", "AVISIT")
  )
})
