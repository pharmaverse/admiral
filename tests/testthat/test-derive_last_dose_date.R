context("test-derive_last_dose_date.R")

input_ae <- tibble::tribble(
  ~STUDYID, ~USUBJID, ~AESEQ, ~AESTDTC,
  "my_study", "subject1", 1, "2020-01-02",
  "my_study", "subject1", 2, "2020-08-31",
  "my_study", "subject1", 3, "2020-10-10",
  "my_study", "subject2", 1, "2019-05-15",
  "my_study", "subject2", 2, "2020-02-20",
  "my_study", "subject3", 1, "2020-03-02",
  "my_study", "subject4", 1, "2020-11-02"
)

input_ex <- tibble::tribble(
  ~STUDYID,   ~USUBJID,   ~EXSTDTC,     ~EXENDTC,    ~EXSEQ, ~EXDOSE, ~EXTRT,
  "my_study", "subject1", "2020-01-01", "2020-01-01", 1,     10,      "treatment",
  "my_study", "subject1", "2020-08-29", "2020-08-29", 2,     10,      "treatment",
  "my_study", "subject1", "2020-09-02", "2020-09-02", 3,     10,      "treatment",
  "my_study", "subject1", "2020-10-20", "2020-10-20", 4,     10,      "treatment",
  "my_study", "subject2", "2019-05-25", "2019-05-25", 1,      0,      "placebo",
  "my_study", "subject2", "2020-01-20", "2020-01-20", 2,      0,      "placebo",
  "my_study", "subject3", "2020-03-15", "2020-03-15", 1,     10,      "treatment"
) %>%
  mutate(EXSTDTC = as.Date(EXSTDTC), EXENDTC = as.Date(EXENDTC))


test_that("derive_last_dose works as expected", {

  expected_output <- mutate(
    input_ae,
    LDOSEDTM = as.POSIXct(
      c("2020-01-01", "2020-08-29", "2020-09-02", NA, "2020-01-20", NA, NA), tz = "UTC"
    )
  )

  res <- derive_last_dose_date(
    input_ae,
    input_ex,
    filter_ex = (EXDOSE > 0) | (EXDOSE == 0 & EXTRT == "placebo"),
    by_vars = vars(STUDYID, USUBJID),
    dose_start = EXSTDTC,
    dose_end = EXENDTC,
    analysis_date = AESTDTC,
    dataset_seq_var = AESEQ,
    new_var = LDOSEDTM,
    output_datetime = TRUE,
    check_dates_only = FALSE,
    traceability_vars = NULL
  )

  expect_dfs_equal(expected_output, res, keys = c("STUDYID", "USUBJID", "AESEQ", "AESTDTC"))

})
