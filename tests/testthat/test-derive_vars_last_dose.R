input_ae <- tibble::tribble(
  ~STUDYID,   ~USUBJID,   ~AESEQ, ~AESTDTC,
  "my_study", "subject1",      1, "2020-01-02",
  "my_study", "subject1",      2, "2020-08-31",
  "my_study", "subject1",      3, "2020-10-10",
  "my_study", "subject2",      1, "2019-05-15",
  "my_study", "subject2",      2, "2020-02-20",
  "my_study", "subject3",      1, "2020-03-02",
  "my_study", "subject4",      1, "2020-11-02"
) %>%
  mutate(
    AESTDT = ymd(AESTDTC)
  )

input_ex <- tibble::tribble(
  ~STUDYID,   ~USUBJID,   ~EXSTDTC,     ~EXENDTC,    ~EXSEQ, ~EXDOSE, ~EXTRT,
  "my_study", "subject1", "2020-01-01", "2020-01-01",     1,      10, "treatment",
  "my_study", "subject1", "2020-08-29", "2020-08-29",     2,      10, "treatment",
  "my_study", "subject1", "2020-09-02", "2020-09-02",     3,      10, "treatment",
  "my_study", "subject1", "2020-10-20", "2020-10-20",     4,      10, "treatment",
  "my_study", "subject2", "2019-05-25", "2019-05-25",     1,       0, "placebo",
  "my_study", "subject2", "2020-01-20", "2020-01-20",     2,       0, "placebo",
  "my_study", "subject3", "2020-03-15", "2020-03-15",     1,      10, "treatment"
) %>%
  mutate(EXSTDT = as.Date(EXSTDTC), EXENDT = as.Date(EXENDTC))

# derive_vars_last_dose ----
## Test 1: function works as expected and returns an error message ----
test_that("derive_vars_last_dose Test 1: function works as expected and returns an error message", {
  expected_output <- mutate(
    input_ae,
    EXSTDT = as.Date(c("2020-01-01", "2020-08-29", "2020-09-02", NA, "2020-01-20", NA, NA)),
    EXENDT = as.Date(c("2020-01-01", "2020-08-29", "2020-09-02", NA, "2020-01-20", NA, NA)),
    EXSEQ = c(1, 2, 3, NA, 2, NA, NA),
    EXDOSE = c(10, 10, 10, NA, 0, NA, NA),
    EXTRT = c("treatment", "treatment", "treatment", NA, "placebo", NA, NA)
  )
  expect_error(
    derive_vars_last_dose(
      input_ae,
      input_ex,
      filter_ex = (EXDOSE > 0) | (EXDOSE == 0 & EXTRT == "placebo"),
      by_vars = exprs(STUDYID, USUBJID),
      dose_date = EXENDT,
      new_vars = exprs(EXDOSE, EXTRT, EXSEQ, EXENDT, EXSTDT),
      analysis_date = AESTDT,
      single_dose_condition = (EXSTDTC == EXENDTC),
      traceability_vars = NULL
    ),
    class = "lifecycle_error_deprecated"
  )
})
