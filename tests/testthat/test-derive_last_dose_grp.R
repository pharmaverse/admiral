input_ae <- tibble::tribble(
  ~STUDYID, ~USUBJID, ~AESEQ, ~AESTDTC,
  "my_study", "subject1", 1, "2020-01-02",
  "my_study", "subject1", 2, "2020-08-31",
  "my_study", "subject1", 3, "2020-10-10",
  "my_study", "subject2", 2, "2020-02-20",
  "my_study", "subject3", 1, "2020-03-02",
  "my_study", "subject4", 1, "2020-11-02"
)

input_ex <- tibble::tribble(
  ~STUDYID,   ~USUBJID,   ~EXSTDTC,     ~EXENDTC,    ~EXSEQ, ~EXDOSE, ~EXTRT,
  "my_study", "subject1", "2020-01-01", "2020-01-01", 1,     1,      "treatment",
  "my_study", "subject1", "2020-08-29", "2020-08-29", 2,     3,      "treatment",
  "my_study", "subject1", "2020-09-02", "2020-09-02", 3,     4,      "treatment",
  "my_study", "subject1", "2020-10-20", "2020-10-20", 4,     4,      "treatment",
  "my_study", "subject2", "2019-05-25", "2019-05-25", 1,     6,      "placebo",
  "my_study", "subject3", "2020-01-20", "2020-01-20", 2,     7,      "placebo",
  "my_study", "subject4", "2020-03-15", "2020-03-15", 1,     13,     "treatment"
) %>%
  mutate(EXSTDTC = as.Date(EXSTDTC), EXENDTC = as.Date(EXENDTC))


derive_last_dose_grp(input_ae,
                     input_ex,
                     filter_ex = (EXDOSE > 0) | (EXDOSE == 0 & EXTRT == "placebo"),
                     by_vars = vars(STUDYID, USUBJID),
                     dose_start = EXSTDTC,
                     dose_end = EXENDTC,
                     grp_var = LDGRP,
                     grp_brks = c(1, 5, 10, 15),
                     grp_lbls = c("G1", "G2", "G3"),
                     dose_var = EXDOSE,
                     analysis_date = AESTDTC,
                     dataset_seq_var = AESEQ,
                     check_dates_only = FALSE,
                     traceability_vars = NULL)
