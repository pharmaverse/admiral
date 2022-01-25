library(testthat)
context("Running test: derive_var_adae_trtemfl")

# Create test data 

adae_toxgr <- tibble::tribble(
  ~STUDYID, ~USUBJID, ~AEGRPID, ~AEDECOD, ~TRTSDTM, ~ASTDTM, ~AENDTM, ~AETOXGR, ~AEITOXGR, ~EXP_TRTEMFL,
  # AE after Treatment Start - expected to be flagged with 'Y'
  "AB99876", "AB99876-001", "AE1", "HEADACHE", ymd_hms("2020-06-30T12:12:12"), ymd_hms("2020-07-08T12:12:12"), ymd_hms("2020-07-08T22:10:30"), "1", "1", "Y",
  
  # AE starts before treatment, continues after, has severity changes in FAAE after treatment start - expected to be flagged with 'Y'
  "AB99876", "AB99876-002", "AE1", "COMMON COLD", ymd_hms("2020-12-06T12:12:12"), ymd_hms("2020-12-05T23:00:00"), ymd_hms("2020-12-08T22:10:30"), NA_character_, "1", "Y",
  
  # AE starts before treatment, continues after, has NO severity changes in FAAE - expected not to be flagged
  "AB99876", "AB99876-002", "AE2", "RASH", ymd_hms("2020-12-06T12:12:12"), ymd_hms("2020-12-05T23:00:00"), ymd_hms("2020-12-08T22:10:30"), "1", "1", NA_character_,
  
  # AE starts before treatment, continues after, no severity changes in FAAE, but missing AETOXGR
  # expected to be flagged with 'Y'
  "AB99876", "AB99876-003", "AE1", "COMMON COLD", ymd_hms("2020-12-06T12:12:12"), ymd_hms("2020-12-05T23:00:00"), ymd_hms("2020-12-08T22:10:30"), NA_character_, "1", "Y",
  
  # AE starts and ends before treatment - expected not to be flagged
  "AB99876", "AB99876-004", "AE1", "ALLERGIC REACTION", ymd_hms("2020-12-06T12:12:12"), ymd_hms("2020-12-04T23:00:00"), ymd_hms("2020-12-05T22:10:30"), "3", "1", NA_character_,
  
  # AE starts before treatment, continues after, has severity changes in FAAE before treatment start 
  # expected to be flagged with 'Y'
  "AB99876", "AB99876-004", "AE2", "VOMITING", ymd_hms("2020-12-06T12:12:12"), ymd_hms("2020-12-05T23:00:00"), ymd_hms("2020-12-08T22:10:30"), "3", NA_character_, "Y",
  
  # Missing TRTSDTM, but a valid dose exists - expected to be flagged with 'Y'
  "AB99876", "AB99876-005", "AE1", "DIARRHOEA", ymd_hms(""), ymd_hms("2020-12-05T23:00:00"), ymd_hms("2020-12-08T22:10:30"), "3", "1", "Y",
  
  # Missing TRTSDTM and no valid dose exists - expected not to be flagged
  "AB99876", "AB99876-006", "AE1", "DIARRHOEA", ymd_hms(""), ymd_hms("2020-12-05T23:00:00"), ymd_hms("2020-12-08T22:10:30"), "3", "1", NA_character_
  )

adae_toxgr2 <- tibble::tribble(
  ~STUDYID, ~USUBJID, ~AEGRPID, ~AEDECOD, ~TRTSDTM, ~ASTDTM, ~AENDTM, ~AETOXGR, ~AEITOXGR, ~EXP_TRTEMFL,
  # AE after Treatment Start - expected to be flagged with 'Y'
  "AB99876", "AB99876-001", "AE1", "HEADACHE", ymd_hms("2020-06-30T12:12:12"), ymd_hms("2020-07-08T12:12:12"), ymd_hms("2020-07-08T22:10:30"), "1", "1", "Y",
  
  # AE starts before treatment, continues after, no FAAE data, but tox change - expected to be flagged with 'Y'
  "AB99876", "AB99876-002", "AE1", "COMMON COLD", ymd_hms("2020-12-06T12:12:12"), ymd_hms("2020-12-05T23:00:00"), ymd_hms("2020-12-08T22:10:30"), "2", "1", "Y",
  
  # AE starts before treatment, continues after, no FAAE data - expected not to be flagged
  "AB99876", "AB99876-002", "AE2", "RASH", ymd_hms("2020-12-06T12:12:12"), ymd_hms("2020-12-05T23:00:00"), ymd_hms("2020-12-08T22:10:30"), "1", "1", NA_character_,
  
  # AE starts before treatment, continues after, no FAAE data, but missing AETOXGR
  # expected to be flagged with 'Y'
  "AB99876", "AB99876-003", "AE1", "COMMON COLD", ymd_hms("2020-12-06T12:12:12"), ymd_hms("2020-12-05T23:00:00"), ymd_hms("2020-12-08T22:10:30"), NA_character_, "1", "Y",
  
  # AE starts and ends before treatment - expected not to be flagged
  "AB99876", "AB99876-004", "AE1", "ALLERGIC REACTION", ymd_hms("2020-12-06T12:12:12"), ymd_hms("2020-12-04T23:00:00"), ymd_hms("2020-12-05T22:10:30"), "3", "1", NA_character_,
  
  # AE starts before treatment, continues after, no FAAE data, but missing initial tox 
  # expected to be flagged with 'Y'
  "AB99876", "AB99876-004", "AE2", "VOMITING", ymd_hms("2020-12-06T12:12:12"), ymd_hms("2020-12-05T23:00:00"), ymd_hms("2020-12-08T22:10:30"), "3", NA_character_, "Y",
  
  # Missing TRTSDTM, but a valid dose exists - expected to be flagged with 'Y'
  "AB99876", "AB99876-005", "AE1", "DIARRHOEA", ymd_hms(""), ymd_hms("2020-12-05T23:00:00"), ymd_hms("2020-12-08T22:10:30"), "3", "1", "Y",
  
  # Missing TRTSDTM and no valid dose exists - expected not to be flagged
  "AB99876", "AB99876-006", "AE1", "DIARRHOEA", ymd_hms(""), ymd_hms("2020-12-05T23:00:00"), ymd_hms("2020-12-08T22:10:30"), "3", "1", NA_character_
)

adae_sev <- tibble::tribble(
  ~STUDYID, ~USUBJID, ~AEGRPID, ~AEDECOD, ~TRTSDTM, ~ASTDTM, ~AENDTM, ~AESEV, ~AESEVIN, ~EXP_TRTEMFL,
  # AE after Treatment Start - expected to be flagged with 'Y'
  "AB99876", "AB99876-001", "AE1", "HEADACHE", ymd_hms("2020-06-30T12:12:12"), ymd_hms("2020-07-08T12:12:12"), ymd_hms("2020-07-08T22:10:30"), "MILD", "MILD", "Y",
  
  # AE starts before treatment, continues after, has severity changes in FAAE after treatment start - expected to be flagged with 'Y'
  "AB99876", "AB99876-002", "AE1", "COMMON COLD", ymd_hms("2020-12-06T12:12:12"), ymd_hms("2020-12-05T23:00:00"), ymd_hms("2020-12-08T22:10:30"), NA_character_, "MILD", "Y",
  
  # AE starts before treatment, continues after, has NO severity changes in FAAE - expected not to be flagged
  "AB99876", "AB99876-002", "AE2", "RASH", ymd_hms("2020-12-06T12:12:12"), ymd_hms("2020-12-05T23:00:00"), ymd_hms("2020-12-08T22:10:30"), "MILD", "MILD", NA_character_,
  
  # AE starts before treatment, continues after, no severity changes in FAAE, but missing AETOXGR
  # expected to be flagged with 'Y'
  "AB99876", "AB99876-003", "AE1", "COMMON COLD", ymd_hms("2020-12-06T12:12:12"), ymd_hms("2020-12-05T23:00:00"), ymd_hms("2020-12-08T22:10:30"), NA_character_, "MILD", "Y",
  
  # AE starts and ends before treatment - expected not to be flagged
  "AB99876", "AB99876-004", "AE1", "ALLERGIC REACTION", ymd_hms("2020-12-06T12:12:12"), ymd_hms("2020-12-04T23:00:00"), ymd_hms("2020-12-05T22:10:30"), "SEVERE", "MILD", NA_character_,
  
  # AE starts before treatment, continues after, has severity changes in FAAE before treatment start 
  # expected to be flagged with 'Y'
  "AB99876", "AB99876-004", "AE2", "VOMITING", ymd_hms("2020-12-06T12:12:12"), ymd_hms("2020-12-05T23:00:00"), ymd_hms("2020-12-08T22:10:30"), "MODERATE", NA_character_, "Y",
  
  # Missing TRTSDTM, but a valid dose exists - expected to be flagged with 'Y'
  "AB99876", "AB99876-005", "AE1", "DIARRHOEA", ymd_hms(""), ymd_hms("2020-12-05T23:00:00"), ymd_hms("2020-12-08T22:10:30"), "MODERATE", "MILD", "Y",
  
  # Missing TRTSDTM and no valid dose exists - expected not to be flagged
  "AB99876", "AB99876-006", "AE1", "DIARRHOEA", ymd_hms(""), ymd_hms("2020-12-05T23:00:00"), ymd_hms("2020-12-08T22:10:30"), "MODERATE", "MILD", NA_character_
)

adae_sev2 <- tibble::tribble(
  ~STUDYID, ~USUBJID, ~AEGRPID, ~AEDECOD, ~TRTSDTM, ~ASTDTM, ~AENDTM, ~AESEV, ~AESEVIN, ~EXP_TRTEMFL,
  # AE after Treatment Start - expected to be flagged with 'Y'
  "AB99876", "AB99876-001", "AE1", "HEADACHE", ymd_hms("2020-06-30T12:12:12"), ymd_hms("2020-07-08T12:12:12"), ymd_hms("2020-07-08T22:10:30"), "MILD", "MILD", "Y",
  
  # AE starts before treatment, continues after, no FAAE data, but sev change - expected to be flagged with 'Y'
  "AB99876", "AB99876-002", "AE1", "COMMON COLD", ymd_hms("2020-12-06T12:12:12"), ymd_hms("2020-12-05T23:00:00"), ymd_hms("2020-12-08T22:10:30"), "MODERATE", "MILD", "Y",
  
  # AE starts before treatment, continues after, no FAAE data - expected not to be flagged
  "AB99876", "AB99876-002", "AE2", "RASH", ymd_hms("2020-12-06T12:12:12"), ymd_hms("2020-12-05T23:00:00"), ymd_hms("2020-12-08T22:10:30"), "MILD", "MILD", NA_character_,
  
  # AE starts before treatment, continues after, no FAAE data, but missing AESEV
  # expected to be flagged with 'Y'
  "AB99876", "AB99876-003", "AE1", "COMMON COLD", ymd_hms("2020-12-06T12:12:12"), ymd_hms("2020-12-05T23:00:00"), ymd_hms("2020-12-08T22:10:30"), NA_character_, "MILD", "Y",
  
  # AE starts and ends before treatment - expected not to be flagged
  "AB99876", "AB99876-004", "AE1", "ALLERGIC REACTION", ymd_hms("2020-12-06T12:12:12"), ymd_hms("2020-12-04T23:00:00"), ymd_hms("2020-12-05T22:10:30"), "SEVERE", "MILD", NA_character_,
  
  # AE starts before treatment, continues after, no FAAE data, but missing initial sev 
  # expected to be flagged with 'Y'
  "AB99876", "AB99876-004", "AE2", "VOMITING", ymd_hms("2020-12-06T12:12:12"), ymd_hms("2020-12-05T23:00:00"), ymd_hms("2020-12-08T22:10:30"), "SEVERE", NA_character_, "Y",
  
  # Missing TRTSDTM, but a valid dose exists - expected to be flagged with 'Y'
  "AB99876", "AB99876-005", "AE1", "DIARRHOEA", ymd_hms(""), ymd_hms("2020-12-05T23:00:00"), ymd_hms("2020-12-08T22:10:30"), "SEVERE", "MILD", "Y",
  
  # Missing TRTSDTM and no valid dose exists - expected not to be flagged
  "AB99876", "AB99876-006", "AE1", "DIARRHOEA", ymd_hms(""), ymd_hms("2020-12-05T23:00:00"), ymd_hms("2020-12-08T22:10:30"), "SEVERE", "MILD", NA_character_
)

ex <- tibble::tribble(
  ~STUDYID, ~USUBJID, ~EXDOSE, ~EXTRT, ~EXOCCUR, ~EXPILTK1, ~EXPILTK2,
  "AB99876", "AB99876-001", 100, "TREATMENT1", NA_character_, "", "",
  "AB99876", "AB99876-002", 0, "PLACEBO", "Y", "", "",
  "AB99876", "AB99876-003", 20, "TREATMENT1", "", "1", "0",
  "AB99876", "AB99876-004", 50, "TREATMENT1", "", "", "",
  "AB99876", "AB99876-005", 10, "TREATMENT1", "Y", "", ""
)

faae_toxgr <- tibble::tribble(
  ~STUDYID, ~USUBJID, ~FAGRPID, ~FATESTCD, ~FAORRES, ~FADTC,
  "AB99876", "AB99876-002", "AE1", "GRADE", "2", "2020-12-07",
  "AB99876", "AB99876-002", "AE1", "GRADE", "1", "2020-12-07",
  "AB99876", "AB99876-004", "AE2", "TOXICITY", "2", "2020-12-06T08:00"
)

faae_sev <- tibble::tribble(
  ~STUDYID, ~USUBJID, ~FAGRPID, ~FATESTCD, ~FAORRES, ~FADTC,
  "AB99876", "AB99876-002", "AE1", "SEVERITY", "MODERATE", "2020-12-07",
  "AB99876", "AB99876-002", "AE1", "SEVERITY", "MILD", "2020-12-07",
  "AB99876", "AB99876-004", "AE2", "SEV", "MILD", "2020-12-06T08:00"
)


test_that("derive_var_adae_trtemfl runs successfully", {
  # Toxicity with FAAE data
  result_toxgr <- derive_var_adae_trtemfl(dataset = adae_toxgr, dataset_faae = faae_toxgr, 
                                    dataset_ex = ex, sevtox = "AETOXGR")
  expect_equal(result_toxgr$EXP_TRTEMFL, result_toxgr$TRTEMFL)
  # Toxicity without FAAE data
  result_toxgr_nofaae <- derive_var_adae_trtemfl(dataset = adae_toxgr2, dataset_ex = ex, sevtox = "AETOXGR")
  expect_equal(result_toxgr_nofaae$EXP_TRTEMFL, result_toxgr_nofaae$TRTEMFL)
  
  # Severity with FAAE data
  result_sev <- derive_var_adae_trtemfl(dataset = adae_sev, dataset_faae = faae_sev, 
                                          dataset_ex = ex, sevtox = "AESEV")
  expect_equal(result_sev$EXP_TRTEMFL, result_sev$TRTEMFL)
  # Severity without FAAE data
  result_sev_nofaae <- derive_var_adae_trtemfl(dataset = adae_sev2, dataset_ex = ex, sevtox = "AESEV")
  expect_equal(result_sev_nofaae$EXP_TRTEMFL, result_sev_nofaae$TRTEMFL)
})

test_that("derive_var_adae_trtemfl with TRTEMFL already present", {
  adae0 <- tibble::tribble(
    ~STUDYID, ~USUBJID, ~AEGRPID, ~AEDECOD, ~TRTSDTM, ~ASTDTM, ~AENDTM, ~AETOXGR, ~AEITOXGR, ~TRTEMFL,
    "AB99876", "AB99876-001", "AE1", "HEADACHE", ymd_hms("2020-06-30T12:12:12"), ymd_hms("2020-07-08T12:12:12"), ymd_hms("2020-07-08T22:10:30"), "1", "1", "Y",
  )  
  expect_warning(derive_var_adae_trtemfl(dataset=adae0, dataset_ex = ex, sevtox = "AETOXGR"), "Variable `TRTEMFL` already exists in the dataset.*")
})

test_that("derive_var_adae_trtemfl with missing input variables", {
  adae1 <- tibble::tribble(
    ~STUDYID, ~USUBJID, ~AEGRPID, ~AEDECOD, ~TRTSDTM, ~ASTDTM, ~AENDTM, ~AETOXGR,
    "AB99876", "AB99876-001", "AE1", "HEADACHE", ymd_hms("2020-06-30T12:12:12"), ymd_hms("2020-07-08T12:12:12"), ymd_hms("2020-07-08T22:10:30"), "1"
  )  
  expect_error(derive_var_adae_trtemfl(dataset=adae1, dataset_ex = ex, sevtox = "AETOXGR"), "Required variable `AEITOXGR` is missing")
  expect_error(derive_var_adae_trtemfl(dataset=adae1,sevtox = "AETOXGR"), 'argument "dataset_ex" is missing, with no default')
  })

test_that("derive_var_adae_trtemfl receives bad input ", {
  expect_error(derive_var_adae_trtemfl("df", dataset_ex = ex), ".*`dataset` must be a data frame.*")
  expect_error(derive_var_adae_trtemfl(dataset=adae_toxgr2, dataset_ex = ex, sevtox = "AETOXGR", subject_keys = c("USUBJID")), "`subject_keys` must be a list of unquoted variable names.*")
  expect_error(derive_var_adae_trtemfl(dataset=adae_toxgr, dataset_ex = ex, sevtox = "AETOX"), "`sevtox` must be one of 'AESEV' or 'AETOXGR'.*")
})

