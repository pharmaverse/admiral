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


test_that("derive_var_last_dose_amt works as expected", {

  expected_output <- mutate(
    input_ae,
    LDOSE = c(10, 10, 10, NA, 0, NA, NA)
  )

  res <- derive_var_last_dose_amt(
    input_ae,
    input_ex,
    filter_ex = (EXDOSE > 0) | (EXDOSE == 0 & EXTRT == "placebo"),
    by_vars = vars(STUDYID, USUBJID),
    dose_date = EXENDTC,
    analysis_date = AESTDTC,
    new_var = LDOSE,
    dose_var = EXDOSE,
    single_dose_condition = (EXSTDTC == EXENDTC),
    traceability_vars = NULL
  )

  expect_dfs_equal(expected_output, res, keys = c("STUDYID", "USUBJID", "AESEQ", "AESTDTC"))

})

test_that("derive_var_last_dose_amt returns traceability vars", {

  expected_output <- mutate(
    input_ae,
    LDOSEDOM = c("EX", "EX", "EX", NA, "EX", NA, NA),
    LDOSESEQ = c(1, 2, 3, NA, 2, NA, NA),
    LDOSEVAR = c("EXSTDTC", "EXSTDTC", "EXSTDTC", NA, "EXSTDTC", NA, NA),
    LDOSE = c(10, 10, 10, NA, 0, NA, NA)
  )

  res <- derive_var_last_dose_amt(
    input_ae,
    input_ex,
    filter_ex = (EXDOSE > 0) | (EXDOSE == 0 & EXTRT == "placebo"),
    by_vars = vars(STUDYID, USUBJID),
    dose_date = EXENDTC,
    analysis_date = AESTDTC,
    new_var = LDOSE,
    dose_var = EXDOSE,
    single_dose_condition = (EXSTDTC == EXENDTC),
    traceability_vars = dplyr::vars(LDOSEDOM = "EX", LDOSESEQ = EXSEQ, LDOSEVAR = "EXSTDTC")
  )

  expect_dfs_equal(expected_output, res, keys = c("STUDYID", "USUBJID", "AESEQ", "AESTDTC"))

})
