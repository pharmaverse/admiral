context("test-derive_param_tte")


test_that("new observations with analysis date are derived correctly", {
  adsl <- tibble::tribble(
    ~USUBJID, ~DTHFL, ~DTHDT,            ~LSTALVDT,         ~TRTSDT,           ~TRTSDTF,
    "03",     "Y",    ymd("2021-08-21"), ymd("2021-08-21"), ymd("2021-08-10"), NA,
    "04",     "N",    NA,                ymd("2021-05-24"), ymd("2021-02-03"), NA) %>%
    mutate(STUDYID = "AB42")

  death <- tte_source(
    dataset_name = "adsl",
    filter = DTHFL == "Y",
    date = DTHDT,
    set_values_to = vars(
      EVENTDESC = "DEATH",
      SRCDOM = "ADSL",
      SRCVAR = "DTHDT"))

  lstalv <- tte_source(
    dataset_name = "adsl",
    date = LSTALVDT,
    censor = 1,
    set_values_to = vars(
      EVENTDESC = "LAST KNOWN ALIVE DATE",
      SRCDOM = "ADSL",
      SRCVAR = "LSTALVDT"))

  expected_output <- tibble::tribble(
    ~USUBJID, ~ADT,              ~CNSR, ~EVENTDESC,              ~SRCDOM, ~SRCVAR,
    "03",     ymd("2021-08-21"), 0L,    "DEATH",                 "ADSL", "DTHDT",
    "04",     ymd("2021-05-24"), 1L,    "LAST KNOWN ALIVE DATE", "ADSL", "LSTALVDT") %>%
    mutate(STUDYID = "AB42",
           PARAMCD = "OS",
           PARAM = "Overall Survival") %>%
    left_join(adsl %>% select(USUBJID, STARTDT = TRTSDT, STARTDTF = TRTSDTF),
              by = "USUBJID")

  expect_dfs_equal(
    derive_param_tte(
      dataset_adsl = adsl,
      start_date = TRTSDT,
      start_date_imputation_flag = TRTSDTF,
      event_conditions = list(death),
      censor_conditions = list(lstalv),
      source_datasets = list(adsl = adsl),
      set_values_to = vars(
        PARAMCD = "OS",
        PARAM = "Overall Survival"
      )
    )
    ,
    expected_output,
    keys = c("USUBJID", "PARAMCD")
  )
})
test_that("new observations with analysis datetime are derived correctly", {
  adsl <- tibble::tribble(
    ~USUBJID, ~DTHFL, ~DTHDT,            ~TRTSDTM,                       ~TRTSDTF, ~TRTSTMF,
    "01",     "Y",    ymd("2021-06-12"), ymd_hms("2021-01-01 00:00:00"), "M",      "H",
    "02",     "N",    NA,                ymd_hms("2021-02-03 10:24:00"), NA,       NA,
    "03",     "Y",    ymd("2021-08-21"), ymd_hms("2021-08-10 00:00:00"), NA,       "H",
    "04",     "N",    NA,                ymd_hms("2021-02-03 10:24:00"), NA,       NA,
    "05",     "N",    NA,                ymd_hms("2021-04-05 11:22:33"), NA,       NA) %>%
    mutate(STUDYID = "AB42")

  adrs <- tibble::tribble(
    ~USUBJID, ~AVALC, ~ADTM,                          ~ASEQ,
    "01",     "SD",   ymd_hms("2021-01-03 10:56:00"), 1,
    "01",     "PR",   ymd_hms("2021-03-04 11:13:00"), 2,
    "01",     "PD",   ymd_hms("2021-05-05 12:02:00"), 3,
    "02",     "PD",   ymd_hms("2021-02-03 10:56:00"), 1,
    "04",     "SD",   ymd_hms("2021-02-13 10:56:00"), 1,
    "04",     "PR",   ymd_hms("2021-04-14 11:13:00"), 2,
    "04",     "CR",   ymd_hms("2021-05-15 12:02:00"), 3) %>%
    mutate(STUDYID = "AB42",
           PARAMCD = "OVR")

  pd <- tte_source(
    dataset_name = "adrs",
    filter = AVALC == "PD",
    date = ADTM,
    set_values_to = vars(
      EVENTDESC = "PD",
      SRCDOM = "ADRS",
      SRCVAR = "ADTM",
      SRCSEQ = ASEQ))

  death <- tte_source(
    dataset_name = "adsl",
    filter = DTHFL == "Y",
    date = DTHDT,
    set_values_to = vars(
      EVENTDESC = "DEATH",
      SRCDOM = "ADSL",
      SRCVAR = "DTHDT"))

  lastvisit <- tte_source(
    dataset_name = "adrs",
    date = ADTM,
    censor = 1,
    set_values_to = vars(
      EVENTDESC = "LAST TUMOR ASSESSMENT",
      SRCDOM = "ADRS",
      SRCVAR = "ADTM"))

  start <- tte_source(
    dataset_name = "adsl",
    date = TRTSDTM,
    censor = 1,
    set_values_to = vars(
      EVENTDESC = "TREATMENT START",
      SRCDOM = "ADSL",
      SRCVAR = "TRTSDTM"))

  expected_output <- tibble::tribble(
    ~USUBJID, ~ADTM,                          ~CNSR, ~EVENTDESC,              ~SRCDOM, ~SRCVAR,   ~SRCSEQ,
    "01",     ymd_hms("2021-05-05 12:02:00"), 0L,    "PD",                    "ADRS",  "ADTM",    3,
    "02",     ymd_hms("2021-02-03 10:56:00"), 0L,    "PD",                    "ADRS",  "ADTM",    1,
    "03",     as_datetime(ymd("2021-08-21")), 0L,    "DEATH",                 "ADSL",  "DTHDT",   NA,
    "04",     ymd_hms("2021-05-15 12:02:00"), 1L,    "LAST TUMOR ASSESSMENT", "ADRS",  "ADTM",    NA,
    "05",     ymd_hms("2021-04-05 11:22:33"), 1L,    "TREATMENT START",       "ADSL",  "TRTSDTM", NA) %>%
    mutate(STUDYID = "AB42",
           PARAMCD = "PFS",
           PARAM = "Progression Free Survival") %>%
    left_join(adsl %>% select(USUBJID, STARTDTM = TRTSDTM, STARTDTF = TRTSDTF, STARTTMF = TRTSTMF),
              by = "USUBJID")

  expect_dfs_equal(
    derive_param_tte(
      dataset_adsl = adsl,
      start_date = TRTSDTM,
      start_date_imputation_flag = TRTSDTF,
      start_time_imputation_flag = TRTSTMF,
      event_conditions = list(pd, death),
      censor_conditions = list(lastvisit, start),
      source_datasets = list(adsl = adsl, adrs = adrs),
      create_datetime = TRUE,
      set_values_to = vars(
        PARAMCD = "PFS",
        PARAM = "Progression Free Survival"
      )
    )
    ,
    expected_output,
    keys = c("USUBJID", "PARAMCD")
  )
})

test_that("new observations based on DTC variables are derived correctly", {
  adsl <- tibble::tribble(
    ~USUBJID, ~TRTSDT,           ~EOSDT,
    "01",     ymd("2020-12-06"), ymd("2021-03-06"),
    "02",     ymd("2021-01-16"), ymd("2021-02-03")) %>%
    mutate(STUDYID = "AB42")

  ae <- tibble::tribble(
    ~USUBJID, ~AESTDTC,           ~AESEQ,
    "01",     "2021-01-03T10:56", 1,
    "01",     "2021-03-04",       2,
    "01",     "2021",             3) %>%
    mutate(STUDYID = "AB42")

  ttae <- tte_source(
    dataset_name = "ae",
    date = AESTDTC,
    set_values_to = vars(
      EVENTDESC = "AE",
      SRCDOM = "AE",
      SRCVAR = "AESTDTC",
      SRCSEQ = AESEQ))

  eos <- tte_source(
    dataset_name = "adsl",
    date = EOSDT,
    censor = 1,
    set_values_to = vars(
      EVENTDESC = "END OF STUDY",
      SRCDOM = "ADSL",
      SRCVAR = "EOSDT"))

  expected_output <- tibble::tribble(
    ~USUBJID, ~ADT,              ~CNSR, ~EVENTDESC,     ~SRCDOM, ~SRCVAR,   ~SRCSEQ,
    "01",     ymd("2021-01-01"), 0L,    "AE",           "AE",    "AESTDTC", 3,
    "02",     ymd("2021-02-03"), 1L,    "END OF STUDY", "ADSL",  "EOSDT",   NA) %>%
    mutate(STUDYID = "AB42",
           PARAMCD = "TTAE",
           PARAM = "Time to First Adverse Event") %>%
    left_join(adsl %>% select(USUBJID, STARTDT = TRTSDT),
              by = "USUBJID")

  expect_dfs_equal(
    derive_param_tte(
      dataset_adsl = adsl,
      start_date = TRTSDT,
      event_conditions = list(ttae),
      censor_conditions = list(eos),
      source_datasets = list(adsl = adsl, ae = ae),
      set_values_to = vars(
        PARAMCD = "TTAE",
        PARAM = "Time to First Adverse Event"
      )
    )
    ,
    expected_output,
    keys = c("USUBJID", "PARAMCD")
  )
})

test_that("by_vars parameter works correctly", {
  adsl <- tibble::tribble(
    ~USUBJID, ~TRTSDT,           ~EOSDT,
    "01",     ymd("2020-12-06"), ymd("2021-03-06"),
    "02",     ymd("2021-01-16"), ymd("2021-02-03")) %>%
    mutate(STUDYID = "AB42")

  ae <- tibble::tribble(
    ~USUBJID, ~AESTDTC,           ~AESEQ, ~AEDECOD,
    "01",     "2021-01-03T10:56", 1,      "Flu",
    "01",     "2021-03-04",       2,      "Cough",
    "01",     "2021",             3,      "Flu") %>%
    mutate(STUDYID = "AB42")

  ttae <- tte_source(
    dataset_name = "ae",
    date = AESTDTC,
    set_values_to = vars(
      EVENTDESC = "AE",
      SRCDOM = "AE",
      SRCVAR = "AESTDTC",
      SRCSEQ = AESEQ))

  eos <- tte_source(
    dataset_name = "adsl",
    date = EOSDT,
    censor = 1,
    set_values_to = vars(
      EVENTDESC = "END OF STUDY",
      SRCDOM = "ADSL",
      SRCVAR = "EOSDT"))

  expected_output <- tibble::tribble(
    ~USUBJID, ~ADT,              ~CNSR, ~EVENTDESC,     ~SRCDOM, ~SRCVAR,   ~SRCSEQ, ~PARCAT2, ~PARAMCD,
    "01",     ymd("2021-01-01"), 0L,    "AE",           "AE",    "AESTDTC", 3,       "Flu",    "TTAE2",
    "02",     ymd("2021-02-03"), 1L,    "END OF STUDY", "ADSL",  "EOSDT",   NA,      "Flu",    "TTAE2",
    "01",     ymd("2021-03-04"), 0L,    "AE",           "AE",    "AESTDTC", 2,       "Cough",  "TTAE1",
    "02",     ymd("2021-02-03"), 1L,    "END OF STUDY", "ADSL",  "EOSDT",   NA,      "Cough",  "TTAE1") %>%
    mutate(STUDYID = "AB42",
           PARCAT1 = "TTAE",
           PARAM = paste("Time to First", PARCAT2, "Adverse Event")) %>%
    left_join(adsl %>% select(USUBJID, STARTDT = TRTSDT),
              by = "USUBJID")

  expect_dfs_equal(
    derive_param_tte(
      dataset_adsl = adsl,
      by_vars = vars(AEDECOD),
      start_date = TRTSDT,
      event_conditions = list(ttae),
      censor_conditions = list(eos),
      source_datasets = list(adsl = adsl, ae = ae),
      set_values_to = vars(
        PARAMCD = paste0("TTAE", as.numeric(as.factor(AEDECOD))),
        PARAM = paste("Time to First", AEDECOD, "Adverse Event"),
        PARCAT1 = "TTAE",
        PARCAT2 = AEDECOD
      )
    )
    ,
    expected_output,
    keys = c("USUBJID", "PARAMCD")
  )
})

test_that("an error is issued if some of the by variables are missing", {
  adsl <- tibble::tribble(
    ~USUBJID, ~TRTSDT,           ~EOSDT,
    "01",     ymd("2020-12-06"), ymd("2021-03-06"),
    "02",     ymd("2021-01-16"), ymd("2021-02-03")) %>%
    mutate(STUDYID = "AB42")

  ae <- tibble::tribble(
    ~USUBJID, ~AESTDTC,           ~AESEQ, ~AEDECOD,
    "01",     "2021-01-03T10:56", 1,      "Flu",
    "01",     "2021-03-04",       2,      "Cough",
    "01",     "2021",             3,      "Flu") %>%
    mutate(STUDYID = "AB42")

  ttae <- tte_source(
    dataset_name = "ae",
    date = AESTDTC,
    set_values_to = vars(
      EVENTDESC = "AE",
      SRCDOM = "AE",
      SRCVAR = "AESTDTC",
      SRCSEQ = AESEQ))

  eos <- tte_source(
    dataset_name = "adsl",
    date = EOSDT,
    censor = 1,
    set_values_to = vars(
      EVENTDESC = "END OF STUDY",
      SRCDOM = "ADSL",
      SRCVAR = "EOSDT"))

  expect_error(
    derive_param_tte(
      dataset_adsl = adsl,
      by_vars = vars(AEBODSYS, AEDECOD),
      start_date = TRTSDT,
      event_conditions = list(ttae),
      censor_conditions = list(eos),
      source_datasets = list(adsl = adsl, ae = ae),
      set_values_to = vars(
        PARAMCD = paste0("TTAE", as.numeric(as.factor(AEDECOD))),
        PARAM = paste("Time to First", AEDECOD, "Adverse Event"),
        PARCAT1 = "TTAE",
        PARCAT2 = AEDECOD
      )
    )
    ,
    regexp = "^Only AEDECOD are included in source dataset.*"
  )
})
test_that("an error is issued all by variables are missing in all source datasets", {
  adsl <- tibble::tribble(
    ~USUBJID, ~TRTSDT,           ~EOSDT,
    "01",     ymd("2020-12-06"), ymd("2021-03-06"),
    "02",     ymd("2021-01-16"), ymd("2021-02-03")) %>%
    mutate(STUDYID = "AB42")

  ae <- tibble::tribble(
    ~USUBJID, ~AESTDTC,           ~AESEQ, ~AEDECOD,
    "01",     "2021-01-03T10:56", 1,      "Flu",
    "01",     "2021-03-04",       2,      "Cough",
    "01",     "2021",             3,      "Flu") %>%
    mutate(STUDYID = "AB42")

  ttae <- tte_source(
    dataset_name = "ae",
    date = AESTDTC,
    set_values_to = vars(
      EVENTDESC = "AE",
      SRCDOM = "AE",
      SRCVAR = "AESTDTC",
      SRCSEQ = AESEQ))

  eos <- tte_source(
    dataset_name = "adsl",
    date = EOSDT,
    censor = 1,
    set_values_to = vars(
      EVENTDESC = "END OF STUDY",
      SRCDOM = "ADSL",
      SRCVAR = "EOSDT"))

  expect_error(
    derive_param_tte(
      dataset_adsl = adsl,
      by_vars = vars(AEBODSYS),
      start_date = TRTSDT,
      event_conditions = list(ttae),
      censor_conditions = list(eos),
      source_datasets = list(adsl = adsl, ae = ae),
      set_values_to = vars(
        PARAMCD = paste0("TTAE", as.numeric(as.factor(AEDECOD))),
        PARAM = paste("Time to First", AEDECOD, "Adverse Event"),
        PARCAT1 = "TTAE",
        PARCAT2 = AEDECOD
      )
    )
    ,
    regexp = "The by variables (AEBODSYS) are not contained in any of the source datasets.",
    fixed = TRUE
  )
})

test_that("an error is issued if there is no one to one mapping between PARAMCD and by_vars", {
  adsl <- tibble::tribble(
    ~USUBJID, ~TRTSDT,           ~EOSDT,
    "01",     ymd("2020-12-06"), ymd("2021-03-06"),
    "02",     ymd("2021-01-16"), ymd("2021-02-03")) %>%
    mutate(STUDYID = "AB42")

  ae <- tibble::tribble(
    ~USUBJID, ~AESTDTC,           ~AESEQ, ~AEDECOD,
    "01",     "2021-01-03T10:56", 1,      "Flu",
    "01",     "2021-03-04",       2,      "Cough",
    "01",     "2021",             3,      "Flu") %>%
    mutate(STUDYID = "AB42")

  ttae <- tte_source(
    dataset_name = "ae",
    date = AESTDTC,
    set_values_to = vars(
      EVENTDESC = "AE",
      SRCDOM = "AE",
      SRCVAR = "AESTDTC",
      SRCSEQ = AESEQ))

  eos <- tte_source(
    dataset_name = "adsl",
    date = EOSDT,
    censor = 1,
    set_values_to = vars(
      EVENTDESC = "END OF STUDY",
      SRCDOM = "ADSL",
      SRCVAR = "EOSDT"))

  expect_error(
    derive_param_tte(
      dataset_adsl = adsl,
      by_vars = vars(AEDECOD),
      start_date = TRTSDT,
      event_conditions = list(ttae),
      censor_conditions = list(eos),
      source_datasets = list(adsl = adsl, ae = ae),
      set_values_to = vars(
        PARAMCD = "TTAE",
        PARCAT2 = AEDECOD
      )
    )
    ,
    regexp = paste0("For some values of PARAMCD there is more than one value of AEDECOD.\n",
                    "Call `get_one_to_many_dataset()` to get all one to many values."),
    fixed = TRUE
  )
})

test_that("an error if issued set_values_to contains invalid expressions", {
  adsl <- tibble::tribble(
    ~USUBJID, ~TRTSDT,           ~EOSDT,
    "01",     ymd("2020-12-06"), ymd("2021-03-06"),
    "02",     ymd("2021-01-16"), ymd("2021-02-03")) %>%
    mutate(STUDYID = "AB42")

  ae <- tibble::tribble(
    ~USUBJID, ~AESTDTC,           ~AESEQ, ~AEDECOD,
    "01",     "2021-01-03T10:56", 1,      "Flu",
    "01",     "2021-03-04",       2,      "Cough",
    "01",     "2021",             3,      "Flu") %>%
    mutate(STUDYID = "AB42")

  ttae <- tte_source(
    dataset_name = "ae",
    date = AESTDTC,
    set_values_to = vars(
      EVENTDESC = "AE",
      SRCDOM = "AE",
      SRCVAR = "AESTDTC",
      SRCSEQ = AESEQ))

  eos <- tte_source(
    dataset_name = "adsl",
    date = EOSDT,
    censor = 1,
    set_values_to = vars(
      EVENTDESC = "END OF STUDY",
      SRCDOM = "ADSL",
      SRCVAR = "EOSDT"))

  expect_error(
    derive_param_tte(
      dataset_adsl = adsl,
      by_vars = vars(AEDECOD),
      start_date = TRTSDT,
      event_conditions = list(ttae),
      censor_conditions = list(eos),
      source_datasets = list(adsl = adsl, ae = ae),
      set_values_to = vars(
        PARAMCD = paste0("TTAE", as.numeric(as.factor(AEDECOD))),
        PARAM = past("Time to First", AEDECOD, "Adverse Event"),
        PARCAT1 = "TTAE",
        PARCAT2 = AEDECOD
      )
    )
    ,
    regexp = paste0("Assigning new variables failed!\n",
                    "set_values_to = \\(\n",
                    "  PARAMCD = paste0\\(\"TTAE\", as.numeric\\(as.factor\\(AEDECOD\\)\\)\\)\n",
                    "  PARAM = past\\(\"Time to First\", AEDECOD, \"Adverse Event\"\\)\n",
                    "  PARCAT1 = TTAE\n",
                    "  PARCAT2 = AEDECOD\n",
                    "\\)\n",
                    "Error message:\n",
                    "  .*")
  )
})

test_that("new observations analysis datetime based on DTC variables are derived correctly", {
  adsl <- tibble::tribble(
    ~USUBJID, ~TRTSDTM,                       ~EOSDT,
    "01",     ymd_hms("2020-12-06 14:23:00"), ymd("2021-03-06"),
    "02",     ymd_hms("2021-01-16 13:09:00"), ymd("2021-02-03")) %>%
    mutate(STUDYID = "AB42")

  ae <- tibble::tribble(
    ~USUBJID, ~AESTDTC,           ~AESEQ,
    "01",     "2021-01-03T10:56", 1,
    "01",     "2021-03-04",       2,
    "01",     "2021-03",          3) %>%
    mutate(STUDYID = "AB42")

  ttae <- tte_source(
    dataset_name = "ae",
    date = AESTDTC,
    set_values_to = vars(
      EVENTDESC = "AE",
      SRCDOM = "AE",
      SRCVAR = "AESTDTC",
      SRCSEQ = AESEQ))

  eos <- tte_source(
    dataset_name = "adsl",
    date = EOSDT,
    censor = 1,
    set_values_to = vars(
      EVENTDESC = "END OF STUDY",
      SRCDOM = "ADSL",
      SRCVAR = "EOSDT"))

  expected_output <- tibble::tribble(
    ~USUBJID, ~ADTM,                          ~CNSR, ~EVENTDESC,     ~SRCDOM, ~SRCVAR,   ~SRCSEQ,
    "01",     ymd_hms("2021-01-03 10:56:00"), 0L,    "AE",           "AE",    "AESTDTC", 1,
    "02",     ymd_hms("2021-02-03 00:00:00"), 1L,    "END OF STUDY", "ADSL",  "EOSDT",   NA) %>%
    mutate(STUDYID = "AB42",
           PARAMCD = "TTAE",
           PARAM = "Time to First Adverse Event") %>%
    left_join(adsl %>% select(USUBJID, STARTDTM = TRTSDTM),
              by = "USUBJID")

  expect_dfs_equal(
    derive_param_tte(
      dataset_adsl = adsl,
      start_date = TRTSDTM,
      event_conditions = list(ttae),
      censor_conditions = list(eos),
      source_datasets = list(adsl = adsl, ae = ae),
      create_datetime = TRUE,
      set_values_to = vars(
        PARAMCD = "TTAE",
        PARAM = "Time to First Adverse Event"
      )
    )
    ,
    expected_output,
    keys = c("USUBJID", "PARAMCD")
  )
})

test_that("error is issued if parameter code already exists", {
  adsl <- tibble::tribble(
    ~USUBJID, ~TRTSDT,           ~EOSDT,
    "01",     ymd("2020-12-06"), ymd("2021-03-06"),
    "02",     ymd("2021-01-16"), ymd("2021-02-03")) %>%
    mutate(STUDYID = "AB42")

  ae <- tibble::tribble(
    ~USUBJID, ~AESTDTC,           ~AESEQ,
    "01",     "2021-01-03T10:56", 1,
    "01",     "2021-03-04",       2,
    "01",     "2021",             3) %>%
    mutate(STUDYID = "AB42")

  ttae <- tte_source(
    dataset_name = "ae",
    date = AESTDTC,
    set_values_to = vars(
      EVENTDESC = "AE",
      SRCDOM = "AE",
      SRCVAR = "AESTDTC",
      SRCSEQ = AESEQ))

  eos <- tte_source(
    dataset_name = "adsl",
    date = EOSDT,
    censor = 1,
    set_values_to = vars(
      EVENTDESC = "END OF STUDY",
      SRCDOM = "ADSL",
      SRCVAR = "EOSDT"))

  expected_output <- tibble::tribble(
    ~USUBJID, ~ADT,              ~CNSR, ~EVENTDESC,     ~SRCDOM, ~SRCVAR,   ~SRCSEQ,
    "01",     ymd("2021-01-01"), 0L,    "AE",           "AE",    "AESTDTC", 3,
    "02",     ymd("2021-02-03"), 1L,    "END OF STUDY", "ADSL",  "EOSDT",   NA) %>%
    mutate(STUDYID = "AB42",
           PARAMCD = "TTAE",
           PARAM = "Time to First Adverse Event") %>%
    left_join(adsl %>% select(USUBJID, STARTDT = TRTSDT),
              by = "USUBJID")

  expect_error(
    derive_param_tte(
      expected_output,
      dataset_adsl = adsl,
      start_date = TRTSDT,
      event_conditions = list(ttae),
      censor_conditions = list(eos),
      source_datasets = list(adsl = adsl, ae = ae),
      set_values_to = vars(
        PARAMCD = "TTAE",
        PARAM = "Time to First Adverse Event"
      )
    )
    ,
    regexp = "^The parameter code TTAE does already exist in `dataset`.$"
  )
})

test_that("error is issued if events with censor != 0 exists", {
  adsl <- tibble::tribble(
    ~USUBJID, ~TRTSDT,           ~EOSDT,
    "01",     ymd("2020-12-06"), ymd("2021-03-06"),
    "02",     ymd("2021-01-16"), ymd("2021-02-03")) %>%
    mutate(STUDYID = "AB42")

  ae <- tibble::tribble(
    ~USUBJID, ~AESTDTC,           ~AESEQ,
    "01",     "2021-01-03T10:56", 1,
    "01",     "2021-03-04",       2,
    "01",     "2021",             3) %>%
    mutate(STUDYID = "AB42")

  ttae <- tte_source(
    dataset_name = "ae",
    date = AESTDTC,
    censor = 1,
    set_values_to = vars(
      EVENTDESC = "AE",
      SRCDOM = "AE",
      SRCVAR = "AESTDTC",
      SRCSEQ = AESEQ))

  eos <- tte_source(
    dataset_name = "adsl",
    date = EOSDT,
    censor = 1,
    set_values_to = vars(
      EVENTDESC = "END OF STUDY",
      SRCDOM = "ADSL",
      SRCVAR = "EOSDT"))

  expected_output <- tibble::tribble(
    ~USUBJID, ~ADT,              ~CNSR, ~EVENTDESC,     ~SRCDOM, ~SRCVAR,   ~SRCSEQ,
    "01",     ymd("2021-01-01"), 0L,    "AE",           "AE",    "AESTDTC", 3,
    "02",     ymd("2021-02-03"), 1L,    "END OF STUDY", "ADSL",  "EOSDT",   NA) %>%
    mutate(STUDYID = "AB42",
           PARAMCD = "TTAE",
           PARAM = "Time to First Adverse Event") %>%
    left_join(adsl %>% select(USUBJID, STARTDT = TRTSDT),
              by = "USUBJID")

  expect_error(
    derive_param_tte(
      expected_output,
      dataset_adsl = adsl,
      start_date = TRTSDT,
      event_conditions = list(ttae),
      censor_conditions = list(eos),
      source_datasets = list(adsl = adsl, ae = ae),
      set_values_to = vars(
        PARAMCD = "TTAE",
        PARAM = "Time to First Adverse Event"
      )
    )
    ,
    regexp = "^The censor value of events must be zero..*"
  )
})
test_that("error is issued if events with censor != 0 exists", {
  adsl <- tibble::tribble(
    ~USUBJID, ~TRTSDT,           ~EOSDT,
    "01",     ymd("2020-12-06"), ymd("2021-03-06"),
    "02",     ymd("2021-01-16"), ymd("2021-02-03")) %>%
    mutate(STUDYID = "AB42")

  ae <- tibble::tribble(
    ~USUBJID, ~AESTDTC,           ~AESEQ,
    "01",     "2021-01-03T10:56", 1,
    "01",     "2021-03-04",       2,
    "01",     "2021",             3) %>%
    mutate(STUDYID = "AB42")

  ttae <- tte_source(
    dataset_name = "ae",
    date = AESTDTC,
    set_values_to = vars(
      EVENTDESC = "AE",
      SRCDOM = "AE",
      SRCVAR = "AESTDTC",
      SRCSEQ = AESEQ))

  eos <- tte_source(
    dataset_name = "adsl",
    date = EOSDT,
    censor = 0,
    set_values_to = vars(
      EVENTDESC = "END OF STUDY",
      SRCDOM = "ADSL",
      SRCVAR = "EOSDT"))

  expected_output <- tibble::tribble(
    ~USUBJID, ~ADT,              ~CNSR, ~EVENTDESC,     ~SRCDOM, ~SRCVAR,   ~SRCSEQ,
    "01",     ymd("2021-01-01"), 0L,    "AE",           "AE",    "AESTDTC", 3,
    "02",     ymd("2021-02-03"), 1L,    "END OF STUDY", "ADSL",  "EOSDT",   NA) %>%
    mutate(STUDYID = "AB42",
           PARAMCD = "TTAE",
           PARAM = "Time to First Adverse Event") %>%
    left_join(adsl %>% select(USUBJID, STARTDT = TRTSDT),
              by = "USUBJID")

  expect_error(
    derive_param_tte(
      expected_output,
      dataset_adsl = adsl,
      start_date = TRTSDT,
      event_conditions = list(ttae),
      censor_conditions = list(eos),
      source_datasets = list(adsl = adsl, ae = ae),
      set_values_to = vars(
        PARAMCD = "TTAE",
        PARAM = "Time to First Adverse Event"
      )
    )
    ,
    regexp = "^The censor value of censorings must be positive..*"
  )
})
