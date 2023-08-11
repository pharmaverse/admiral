input_ae <- tibble::tribble(
  ~STUDYID, ~USUBJID, ~AESEQ, ~AESTDTC,
  "my_study", "subject1", 1, "2020-01-02",
  "my_study", "subject1", 2, "2020-08-31",
  "my_study", "subject1", 3, "2020-10-10",
  "my_study", "subject2", 1, "2019-05-15",
  "my_study", "subject2", 2, "2020-02-20",
  "my_study", "subject3", 1, "2020-03-02",
  "my_study", "subject4", 1, "2020-11-02"
) %>% mutate(
  AESTDT = ymd(AESTDTC)
)

input_ex <- tibble::tribble(
  ~STUDYID, ~USUBJID, ~EXSTDTC, ~EXENDTC, ~EXSEQ, ~EXDOSE, ~EXTRT,
  "my_study", "subject1", "2020-01-01", "2020-01-01", 1, 10, "treatment",
  "my_study", "subject1", "2020-08-29", "2020-08-29", 2, 10, "treatment",
  "my_study", "subject1", "2020-09-02", "2020-09-02", 3, 10, "treatment",
  "my_study", "subject1", "2020-10-20", "2020-10-20", 4, 10, "treatment",
  "my_study", "subject2", "2019-05-25", "2019-05-25", 1, 0, "placebo",
  "my_study", "subject2", "2020-01-20", "2020-01-20", 2, 0, "placebo",
  "my_study", "subject3", "2020-03-15", "2020-03-15", 1, 10, "treatment"
) %>%
  mutate(EXSTDT = as.Date(EXSTDTC), EXENDT = as.Date(EXENDTC))

# derive_var_last_dose_date ----
## Test 1: works as expected output_datetime = FALSE ----
test_that("derive_var_last_dose_date Test 1: works as expected output_datetime = FALSE", {
  expected_output <- tibble::tribble(
    ~STUDYID, ~USUBJID, ~AESEQ, ~AESTDTC, ~LDOSEDTM,
    "my_study", "subject1", 1, "2020-01-02", "2020-01-01",
    "my_study", "subject1", 2, "2020-08-31", "2020-08-29",
    "my_study", "subject1", 3, "2020-10-10", "2020-09-02",
    "my_study", "subject2", 1, "2019-05-15", NA_character_,
    "my_study", "subject2", 2, "2020-02-20", "2020-01-20",
    "my_study", "subject3", 1, "2020-03-02", NA_character_,
    "my_study", "subject4", 1, "2020-11-02", NA_character_
  ) %>%
    mutate(
      LDOSEDTM = as.Date(LDOSEDTM),
      AESTDT = ymd(AESTDTC)
    )
  expect_error(
    derive_var_last_dose_date(
      input_ae,
      input_ex,
      filter_ex = (EXDOSE > 0) | (EXDOSE == 0 & EXTRT == "placebo"),
      by_vars = exprs(STUDYID, USUBJID),
      dose_date = EXENDT,
      analysis_date = AESTDT,
      new_var = LDOSEDTM,
      single_dose_condition = (EXSTDTC == EXENDTC),
      output_datetime = FALSE,
      traceability_vars = NULL
    ),
    class = "lifecycle_error_deprecated"
  )
})
