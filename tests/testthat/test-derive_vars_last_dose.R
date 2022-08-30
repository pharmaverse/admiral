library(tibble)
library(dplyr)
library(lubridate)
input_ae <- tribble(
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

input_ex <- tribble(
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
## Test 1: function works as expected ----
test_that("derive_vars_last_dose Test 1: function works as expected", {
  expected_output <- mutate(
    input_ae,
    EXSTDT = as.Date(c("2020-01-01", "2020-08-29", "2020-09-02", NA, "2020-01-20", NA, NA)),
    EXENDT = as.Date(c("2020-01-01", "2020-08-29", "2020-09-02", NA, "2020-01-20", NA, NA)),
    EXSEQ = c(1, 2, 3, NA, 2, NA, NA),
    EXDOSE = c(10, 10, 10, NA, 0, NA, NA),
    EXTRT = c("treatment", "treatment", "treatment", NA, "placebo", NA, NA)
  )

  res <- derive_vars_last_dose(
    input_ae,
    input_ex,
    filter_ex = (EXDOSE > 0) | (EXDOSE == 0 & EXTRT == "placebo"),
    by_vars = vars(STUDYID, USUBJID),
    dose_date = EXENDT,
    new_vars = vars(EXDOSE, EXTRT, EXSEQ, EXENDT, EXSTDT),
    analysis_date = AESTDT,
    single_dose_condition = (EXSTDTC == EXENDTC),
    traceability_vars = NULL
  )

  expect_dfs_equal(expected_output, res, keys = c("STUDYID", "USUBJID", "AESEQ", "AESTDTC"))
})

## Test 2: function checks validity of start and end dose inputs ----
test_that("derive_vars_last_dose Test 2: function checks validity of start and end dose inputs", {
  input_ex_wrong <- bind_rows(
    input_ex,
    tribble(
      ~STUDYID, ~USUBJID, ~EXSTDTC, ~EXENDTC, ~EXSEQ, ~EXDOSE, ~EXTRT,
      "my_study", "subject4", "2020-11-05", "2020-11-06", 1, 10, "treatment"
    ) %>%
      mutate(
        EXENDT = ymd(EXENDTC),
        EXSTDT = ymd(EXSTDTC)
      )
  )

  expect_error(
    derive_vars_last_dose(
      input_ae,
      input_ex_wrong,
      filter_ex = (EXDOSE > 0) | (EXDOSE == 0 & EXTRT == "placebo"),
      by_vars = vars(STUDYID, USUBJID),
      dose_date = EXENDT,
      analysis_date = AESTDT,
      single_dose_condition = (EXSTDTC == EXENDTC),
      traceability_vars = NULL
    ),
    regexp = "Specified `single_dose_condition` is not satisfied."
  )
})

## Test 3: function returns traceability vars ----
test_that("derive_vars_last_dose Test 3: function returns traceability vars", {
  expected_output <- mutate(
    input_ae,
    EXSTDTC = c("2020-01-01", "2020-08-29", "2020-09-02", NA, "2020-01-20", NA, NA),
    EXENDTC = c("2020-01-01", "2020-08-29", "2020-09-02", NA, "2020-01-20", NA, NA),
    EXENDT = ymd(EXENDTC),
    EXSTDT = ymd(EXSTDTC),
    EXSEQ = c(1, 2, 3, NA, 2, NA, NA),
    EXDOSE = c(10, 10, 10, NA, 0, NA, NA),
    EXTRT = c("treatment", "treatment", "treatment", NA, "placebo", NA, NA),
    LDOSEDOM = c("EX", "EX", "EX", NA, "EX", NA, NA),
    LDOSESEQ = c(1, 2, 3, NA, 2, NA, NA),
    LDOSEVAR = c("EXSTDTC", "EXSTDTC", "EXSTDTC", NA, "EXSTDTC", NA, NA)
  )

  res <- derive_vars_last_dose(
    input_ae,
    input_ex,
    filter_ex = (EXDOSE > 0) | (EXDOSE == 0 & EXTRT == "placebo"),
    by_vars = vars(STUDYID, USUBJID),
    dose_date = EXENDT,
    analysis_date = AESTDT,
    single_dose_condition = (EXSTDTC == EXENDTC),
    traceability_vars = vars(LDOSEDOM = "EX", LDOSESEQ = EXSEQ, LDOSEVAR = "EXSTDTC")
  )

  expect_dfs_equal(expected_output, res, keys = c("STUDYID", "USUBJID", "AESEQ", "AESTDTC"))
})

## Test 4: function errors when multiple doses are on same date ----
test_that("derive_vars_last_dose Test 4: function errors when multiple doses are on same date", {
  input_ex_dup <- bind_rows(
    input_ex,
    tribble(
      ~STUDYID, ~USUBJID, ~EXSTDTC, ~EXENDTC, ~EXSEQ, ~EXDOSE, ~EXTRT,
      "my_study", "subject2", "2020-01-20", "2020-01-20", 3, 0, "placebo"
    ) %>%
      mutate(
        EXSTDT = ymd(EXSTDTC),
        EXENDT = ymd(EXENDTC)
      )
  )

  expected_output <- mutate(
    input_ae,
    EXSTDTC = c("2020-01-01", "2020-08-29", "2020-09-02", NA, "2020-01-20", NA, NA),
    EXENDTC = c("2020-01-01", "2020-08-29", "2020-09-02", NA, "2020-01-20", NA, NA),
    EXSTDT = ymd(EXSTDTC),
    EXENDT = ymd(EXENDTC),
    EXSEQ = c(1, 2, 3, NA, 3, NA, NA),
    EXDOSE = c(10, 10, 10, NA, 0, NA, NA),
    EXTRT = c("treatment", "treatment", "treatment", NA, "placebo", NA, NA)
  )

  expect_error(
    derive_vars_last_dose(
      input_ae,
      input_ex_dup,
      filter_ex = (EXDOSE > 0) | (EXDOSE == 0 & EXTRT == "placebo"),
      by_vars = vars(STUDYID, USUBJID),
      dose_date = EXENDT,
      analysis_date = AESTDT,
      single_dose_condition = (EXSTDTC == EXENDTC),
      traceability_vars = NULL
    ),
    regexp = "Multiple doses exist for the same `dose_date`. Update `dose_id` to identify unique doses." # nolint
  )
})

## Test 5: multiple doses on same date - dose_id supplied ----
test_that("derive_vars_last_dose Test 5: multiple doses on same date - dose_id supplied", {
  input_ex_dup <- bind_rows(
    input_ex,
    tribble(
      ~STUDYID, ~USUBJID, ~EXSTDTC, ~EXENDTC, ~EXSEQ, ~EXDOSE, ~EXTRT,
      "my_study", "subject2", "2020-01-20", "2020-01-20", 3, 0, "placebo"
    ) %>% mutate(
      EXSTDT = ymd(EXSTDTC),
      EXENDT = ymd(EXENDTC)
    )
  )

  expected_output <- mutate(
    input_ae,
    EXSTDT = ymd(c("2020-01-01", "2020-08-29", "2020-09-02", NA, "2020-01-20", NA, NA)),
    EXENDT = ymd(c("2020-01-01", "2020-08-29", "2020-09-02", NA, "2020-01-20", NA, NA)),
    EXSEQ = c(1, 2, 3, NA, 3, NA, NA),
    EXDOSE = c(10, 10, 10, NA, 0, NA, NA),
    EXTRT = c("treatment", "treatment", "treatment", NA, "placebo", NA, NA)
  )

  res <- derive_vars_last_dose(
    input_ae,
    input_ex_dup,
    filter_ex = (EXDOSE > 0) | (EXDOSE == 0 & EXTRT == "placebo"),
    by_vars = vars(STUDYID, USUBJID),
    dose_date = EXENDT,
    dose_id = vars(EXSEQ),
    new_vars = vars(EXDOSE, EXTRT, EXSEQ, EXSTDT, EXENDT),
    analysis_date = AESTDT,
    single_dose_condition = (EXSTDTC == EXENDTC),
    traceability_vars = NULL
  )

  expect_dfs_equal(expected_output, res, keys = c("STUDYID", "USUBJID", "AESEQ", "AESTDTC"))
})

## Test 6: error is issued if same variable is found in both input datasets ----
test_that("derive_vars_last_dose Test 6: error is issued if same variable is found in both input datasets", { # nolint
  input_ae <- tribble(
    ~STUDYID,   ~USUBJID,   ~AESEQ, ~EXSTDTC,
    "my_study", "subject1",      1, "2020-01-02",
    "my_study", "subject1",      2, "2020-08-31",
    "my_study", "subject1",      3, "2020-10-10",
    "my_study", "subject2",      1, "2019-05-15",
    "my_study", "subject2",      2, "2020-02-20",
    "my_study", "subject3",      1, "2020-03-02",
    "my_study", "subject4",      1, "2020-11-02"
  ) %>%
    mutate(
      EXSTDT = ymd(EXSTDTC)
    )

  input_ex <- tribble(
    ~STUDYID, ~USUBJID, ~EXSTDTC, ~EXENDTC, ~EXSEQ, ~EXDOSE, ~EXTRT,
    "my_study", "subject1", "2020-01-01", "2020-01-01", 1, 10, "treatment",
    "my_study", "subject1", "2020-08-29", "2020-08-29", 2, 10, "treatment",
    "my_study", "subject1", "2020-09-02", "2020-09-02", 3, 10, "treatment",
    "my_study", "subject1", "2020-10-20", "2020-10-20", 4, 10, "treatment",
    "my_study", "subject2", "2019-05-25", "2019-05-25", 1, 0, "placebo",
    "my_study", "subject2", "2020-01-20", "2020-01-20", 2, 0, "placebo",
    "my_study", "subject3", "2020-03-15", "2020-03-15", 1, 10, "treatment"
  ) %>%
    mutate(
      EXSTDT = as.Date(EXSTDTC),
      EXENDT = as.Date(EXENDTC)
    )

  expect_error(
    derive_vars_last_dose(
      input_ae,
      input_ex,
      filter_ex = (EXDOSE > 0) | (EXDOSE == 0 & EXTRT == "placebo"),
      by_vars = vars(STUDYID, USUBJID),
      dose_date = EXENDT,
      new_vars = vars(EXDOSE, EXTRT, EXSEQ, EXENDT, EXSTDT),
      analysis_date = EXSTDT,
      single_dose_condition = (EXSTDTC == EXENDTC),
      traceability_vars = NULL
    ),
    "Variable(s) `EXSTDT` found in both datasets, cannot perform join",
    fixed = TRUE
  )
})

## Test 7: no error is raised when setting `dose_date` to a renamed variable ----
test_that("derive_vars_last_dose Test 7: no error is raised when setting `dose_date` to a renamed variable", { # nolint
  adae <- tribble(
    ~USUBJID, ~AESTDTC, ~AENDTC, ~ASTDT, ~AENDT, ~AEDECOD,
    "P01", "2022-01-10", "2022-01-12", ymd("2022-01-10"), ymd("2022-01-12"), "Nausea",
    "P02", "2022-01-31", "2022-01-31", ymd("2022-01-31"), ymd("2022-01-31"), "Vomitting",
    "P02", "2022-02-02", "2022-02-04", ymd("2022-02-02"), ymd("2022-02-04"), "Vomitting"
  )

  adex <- tribble(
    ~USUBJID, ~EXTRT, ~EXDOSFRQ, ~EXSTDTC, ~EXENDTC, ~ASTDT, ~AENDT, ~ASTDTM, ~AENDTM
    "P01", "Drug A", "QD", "2022-01-09", "2022-01-12", ymd("2022-01-09"), ymd("2022-01-12"),
    ymd_hms("2022-01-09 09:30:00"), ymd_hms("2022-01-12 09:30:00"),
    "P02", "Drug A", "QD", "2022-02-01", "2022-02-04", ymd("2022-02-01"), ymd("2022-02-04"),
    ymd_hms("2022-02-01 10:00:00"), ymd_hms("2022-02-04 10:00:00")
  )

  (adex_single <- create_single_dose_dataset(adex))

  expect_error(
    derive_vars_last_dose(
      adae,
      adex_single,
      by_vars = vars(USUBJID),
      dose_date = EXSTDT,
      analysis_date = ASTDT,
      new_vars = vars(EXSTDT = ASTDT)
    ),
    NA
  )
})
