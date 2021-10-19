context("test-get_last_dose.R")

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


test_that("get_last_dose works as expected", {

  expected_output <- mutate(
    input_ae,
    EXSTDTC = as.Date(c("2020-01-01", "2020-08-29", "2020-09-02", NA, "2020-01-20", NA, NA)),
    EXENDTC = as.Date(c("2020-01-01", "2020-08-29", "2020-09-02", NA, "2020-01-20", NA, NA)),
    EXSEQ = c(1, 2, 3, NA, 2, NA, NA),
    EXDOSE = c(10, 10, 10, NA, 0, NA, NA),
    EXTRT = c("treatment", "treatment", "treatment", NA, "placebo", NA, NA)
  )

  res <- get_last_dose(
    input_ae,
    input_ex,
    filter_ex = (EXDOSE > 0) | (EXDOSE == 0 & EXTRT == "placebo"),
    by_vars = vars(STUDYID, USUBJID),
    dose_start = EXSTDTC,
    dose_end = EXENDTC,
    analysis_date = AESTDTC,
    dataset_seq_var = AESEQ,
    check_dates_only = FALSE,
    traceability_vars = NULL
  )

  expect_dfs_equal(expected_output, res, keys = c("STUDYID", "USUBJID", "AESEQ", "AESTDTC"))

})


test_that("get_last_dose checks validity of start and end dose inputs", {

  input_ex_wrong <- dplyr::bind_rows(
    input_ex,
    tibble::tribble(
      ~STUDYID,   ~USUBJID,   ~EXSTDTC,              ~EXENDTC,        ~EXSEQ, ~EXDOSE, ~EXTRT,
      "my_study", "subject4", as.Date("2020-11-05"), as.Date("2020-11-06"), 1, 10, "treatment")
  )

  expect_error(
    get_last_dose(
      input_ae,
      input_ex_wrong,
      filter_ex = (EXDOSE > 0) | (EXDOSE == 0 & EXTRT == "placebo"),
      by_vars = vars(STUDYID, USUBJID),
      dose_start = EXSTDTC,
      dose_end = EXENDTC,
      analysis_date = AESTDTC,
      dataset_seq_var = AESEQ,
      check_dates_only = FALSE,
      traceability_vars = NULL
    ),
    regexp = "Not all values of EXSTDTC are equal to EXENDTC"
  )

})


test_that(paste("get_last_dose checks validity of start and end dose inputs",
                "- time component (check_dates_only = FALSE)"), {

                  input_ex_wrong <- dplyr::bind_rows(
                    mutate_at(input_ex, c("EXSTDTC", "EXENDTC"), as.POSIXct),
                    tibble::tribble(
                      ~STUDYID,   ~USUBJID,   ~EXSTDTC, ~EXENDTC, ~EXSEQ, ~EXDOSE, ~EXTRT,
                      "my_study", "subject4", as.POSIXct("2020-11-06 00:00:00"),
                      as.POSIXct("2020-11-06 00:00:01"), 1, 10, "treatment")
                  )

                  expect_error(
                    get_last_dose(
                      input_ae,
                      input_ex_wrong,
                      filter_ex = (EXDOSE > 0) | (EXDOSE == 0 & EXTRT == "placebo"),
                      by_vars = vars(STUDYID, USUBJID),
                      dose_start = EXSTDTC,
                      dose_end = EXENDTC,
                      analysis_date = AESTDTC,
                      dataset_seq_var = AESEQ,
                      check_dates_only = FALSE,
                      traceability_vars = NULL
                    ),
                    regexp = "Not all values of EXSTDTC are equal to EXENDTC"
                  )

})



test_that(
  paste("get_last_dose checks validity of start and end dose inputs",
        "- time component (check_dates_only = TRUE)"), {

          input_ex_wrong <- dplyr::bind_rows(
            mutate_at(input_ex, c("EXSTDTC", "EXENDTC"), ~as.POSIXct(.x, tz = "UTC")),
            tibble::tribble(
              ~STUDYID,   ~USUBJID,   ~EXSTDTC, ~EXENDTC, ~EXSEQ, ~EXDOSE, ~EXTRT,
              "my_study", "subject4", as.POSIXct("2020-11-06T00:00:00", tz = "UTC"),
              as.POSIXct("2020-11-06T00:00:01", tz = "UTC"), 1, 10, "treatment")
          )

          expected_output <- mutate(
            input_ae,
            EXSTDTC = as.POSIXct(
              c("2020-01-01", "2020-08-29", "2020-09-02", NA, "2020-01-20", NA, NA), tz = "UTC"),
            EXENDTC = as.POSIXct(
              c("2020-01-01", "2020-08-29", "2020-09-02", NA, "2020-01-20", NA, NA), tz = "UTC"),
            EXSEQ = c(1, 2, 3, NA, 2, NA, NA),
            EXDOSE = c(10, 10, 10, NA, 0, NA, NA),
            EXTRT = c("treatment", "treatment", "treatment", NA, "placebo", NA, NA)
          )

          res <- get_last_dose(
            input_ae,
            input_ex_wrong,
            filter_ex = (EXDOSE > 0) | (EXDOSE == 0 & EXTRT == "placebo"),
            by_vars = vars(STUDYID, USUBJID),
            dose_start = EXSTDTC,
            dose_end = EXENDTC,
            analysis_date = AESTDTC,
            dataset_seq_var = AESEQ,
            check_dates_only = TRUE,
            traceability_vars = NULL
          )

          expect_dfs_equal(expected_output, res, keys = c("STUDYID", "USUBJID", "AESEQ", "AESTDTC"))

})


test_that("derive_last_dose returns traceability vars", {

  expected_output <- mutate(
    input_ae,
    EXSTDTC = as.Date(c("2020-01-01", "2020-08-29", "2020-09-02", NA, "2020-01-20", NA, NA)),
    EXENDTC = as.Date(c("2020-01-01", "2020-08-29", "2020-09-02", NA, "2020-01-20", NA, NA)),
    EXSEQ = c(1, 2, 3, NA, 2, NA, NA),
    EXDOSE = c(10, 10, 10, NA, 0, NA, NA),
    EXTRT = c("treatment", "treatment", "treatment", NA, "placebo", NA, NA),
    LDOSEDOM = c("EX", "EX", "EX", NA, "EX", NA, NA),
    LDOSESEQ = c(1, 2, 3, NA, 2, NA, NA),
    LDOSEVAR = c("EXSTDTC", "EXSTDTC", "EXSTDTC", NA, "EXSTDTC", NA, NA)
    )

  res <- get_last_dose(
    input_ae,
    input_ex,
    filter_ex = (EXDOSE > 0) | (EXDOSE == 0 & EXTRT == "placebo"),
    by_vars = vars(STUDYID, USUBJID),
    dose_start = EXSTDTC,
    dose_end = EXENDTC,
    analysis_date = AESTDTC,
    dataset_seq_var = AESEQ,
    check_dates_only = FALSE,
    traceability_vars = dplyr::vars(LDOSEDOM = "EX", LDOSESEQ = EXSEQ, LDOSEVAR = "EXSTDTC")
  )

  expect_dfs_equal(expected_output, res, keys = c("STUDYID", "USUBJID", "AESEQ", "AESTDTC"))

})
