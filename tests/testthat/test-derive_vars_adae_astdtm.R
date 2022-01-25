library(testthat)
context("Running test: derive_vars_adae_astdtm")

# Create a test data set that covers all potential cases and also add the expected result
# for comparison
adae <- tibble::tribble(
  ~AESTDTC, ~TRTSDTM,  ~AENDTM, ~EXP_ASTDTM, ~EXP_ASTDTF, ~EXP_ASTTMF,
  # Set to Start Date/Time of Adverse Event [ADAE.AESTDTC] where Start Date/Time of
  # Adverse Event [ADAE.AESTDTC] is complete.
  "2020-07-02T22:10:30", ymd_hms("2020-12-06T12:12:12"), ymd_hms("2020-07-08T12:12:12"), ymd_hms("2020-07-02T22:10:30"), NA_character_, NA_character_,

  # Missing Time parts are imputed with 23:59:59
  "2020-07-02T22:10", ymd_hms("2020-12-06T12:12:12"), ymd_hms("2020-07-08T12:12:12"), ymd_hms("2020-07-02T22:10:59"), NA_character_, "S",
  "2020-07-02T22", ymd_hms("2020-12-06T12:12:12"), ymd_hms("2020-07-08T12:12:12"), ymd_hms("2020-07-02T22:59:59"), NA_character_, "M",
  "2020-07-02", ymd_hms("2020-12-06T12:12:12"), ymd_hms("2020-07-08T12:12:12"), ymd_hms("2020-07-02T23:59:59"), NA_character_, "H",

  # If imputation of missing Time parts leads to a date/time after AENDTM, AENDTM is used
  # instead
  "2020-07-08", ymd_hms("2020-12-06T12:12:12"), ymd_hms("2020-07-08T12:12:12"), ymd_hms("2020-07-08T12:12:12"), NA_character_, "H",

  ##### Missing Day ######
  # Imputed as Datetime of First Exposure to Treatment [ADSL.TRTSDTM]
  # (with time replaced as 23:59:59) if Datetime of First Exposure to Treatment
  # [ADSL.TRTSDTM] is in same month and year as start date of adverse event [ADAE.AESTDTC].
  "2020-07", ymd_hms("2020-07-06T12:30:00"), ymd_hms("2020-07-08T12:12:12"), ymd_hms("2020-07-06T23:59:59"), "D", "H",

  # Else set missing day to '01' and time to '23:59:59'.
  "2020-07", ymd_hms("2019-01-06T12:12:12"), ymd_hms("2020-07-08T12:12:12"), ymd_hms("2020-07-01T23:59:59"), "D", "H",
  "2020-07", ymd_hms(""), ymd_hms("2020-07-08T12:12:12"), ymd_hms("2020-07-01T23:59:59"), "D", "H",
  "2020-07", ymd_hms("2019-01-06T12:12:12"), ymd_hms(""), ymd_hms("2020-07-01T23:59:59"), "D", "H",

  # Set missing day to '01' and time to '23:59:59' if Analysis End Date/Time [ADAE.AENDTM]
  # is before Datetime of First Exposure to Treatment [ADSL.TRTSDTM].
  "2020-07", ymd_hms("2020-07-26T12:12:12"), ymd_hms("2020-07-08T12:12:12"), ymd_hms("2020-07-01T23:59:59"), "D", "H",

  ##### Missing Month ######
  # Set to Datetime of First Exposure to Treatment [ADSL.TRTSDTM] (with time replaced
  # as 23:59:59) if Datetime of First Exposure to Treatment [ADSL.TRTSDTM] is not missing
  # and has matching year part to Start Date/Time of Adverse Event [ADAE.AESTDTC].
  "2020", ymd_hms("2020-07-06T12:12:12"), ymd_hms("2020-07-08T12:12:12"), ymd_hms("2020-07-06T23:59:59"), "M", "H",
  "2020---03", ymd_hms("2020-07-06T12:12:12"), ymd_hms("2020-07-08T12:12:12"), ymd_hms("2020-07-06T23:59:59"), "M", "H",

  # Else set to Analysis End Date/Time [ADAE.AENDTM] if Datetime of First Exposure to
  # Treatment [ADSL.TRTSDTM] is missing or has a year part greater than the year part of
  # Start Date/Time of Adverse Event [ADAE.AESTDTC], and Analysis End Date/Time
  # [ADAE.AENDTM] has matching year part to Start Date/Time of Adverse Event [ADAE.AESTDTC].
  "2020", ymd_hms(""), ymd_hms("2020-07-08T12:12:12"), ymd_hms("2020-07-08T12:12:12"),  "M", "H",
  "2020", ymd_hms("2021-01-06T12:12:12"), ymd_hms("2020-07-08T12:12:12"), ymd_hms("2020-07-08T12:12:12"), "M", "H",

  # If start date imputation takes date past Analysis End Date/Time [ADAE.AENDTM],
  # then set to Analysis End Date/Time [ADAE.AENDTM].
  "2020", ymd_hms("2020-12-06T12:12:12"), ymd_hms("2020-07-08T12:12:12"), ymd_hms("2020-07-08T12:12:12"), "M", "H",

  ##### Missing Year ######
  # Set to earliest of non-missing Datetime of First Exposure to Treatment [ADSL.TRTSDTM]
  # or the (imputed) Analysis End Date/Time [ADAE.AENDTM].
  "", ymd_hms("2020-12-06T12:12:12"), ymd_hms("2020-05-08T12:12:12"),  ymd_hms("2020-05-08T12:12:12"), "Y", "H",
  "", ymd_hms("2020-01-06T12:12:12"), ymd_hms("2020-05-08T12:12:12"), ymd_hms("2020-01-06T12:12:12"), "Y", "H",
  "", ymd_hms("2020-12-06T12:12:12"), ymd_hms(""), ymd_hms("2020-12-06T12:12:12"), "Y", "H",
  "", ymd_hms(""), ymd_hms("2020-07-08T12:12:12"), ymd_hms("2020-07-08T12:12:12"), "Y", "H",

  # Set to null, if Datetime of First Exposure to Treatment [ADSL.TRTSDTM] and Analysis
  # End Date/Time [ADAE.AENDTM] are is missing.
  "", ymd_hms(""), ymd_hms(""), ymd_hms(""), NA_character_, NA_character_
)


test_that("derive_vars_adae_astdtm runs successfully", {
  result <- derive_vars_adae_astdtm(adae)
  expect_identical(result$EXP_ASTDTM, result$ASTDTM)
  expect_identical(result$EXP_ASTDTF, result$ASTDTF)
  expect_identical(result$EXP_ASTTMF, result$ASTTMF)
})

test_that("derive_vars_adae_astdtm drops empty flags", {
  adae0 <- tibble::tribble(
    ~AESTDTC, ~TRTSDTM, ~AENDTM,
    "2020-07-02T22:10:30", ymd_hms("2020-12-06T12:12:12"), ymd_hms("2020-07-08T12:12:12"),
    "2020-07-02T22:10", ymd_hms("2020-12-06T12:12:12"), ymd_hms("2020-07-08T12:12:12")
  )
  result0 <- derive_vars_adae_astdtm(adae0, drop_empty_flag = T)
  expect_equal("ASTDTF" %in% colnames(result0), FALSE)
  expect_equal("ASTTMF" %in% colnames(result0), TRUE)
})

test_that("derive_vars_adae_astdtm with missing input variables", {
  adae1 <- tibble::tribble(
    ~AESTDTC, ~TRTSDTM,
    "2020-07-02T22:10:30", ymd_hms("2020-12-06T12:12:12")
  )
  expect_error(derive_vars_adae_astdtm(adae1), "Required variable `AENDTM` is missing")
  adae2 <- tibble::tribble(
    ~AESTDTC, ~AENDTM,
    "2020-07-02T22:10:30", ymd_hms("2020-12-06T12:12:12")
  )
  expect_error(derive_vars_adae_astdtm(adae2), "Required variable `TRTSDTM` is missing")
})

test_that("derive_vars_adae_astdtm receives bad input ", {
  expect_error(derive_vars_adae_astdtm("df"), ".*`dataset` must be a data frame.*")
  expect_error(derive_vars_adae_astdtm(adae, drop_empty_flag = "Y"), "`drop_empty_flag` must be either `TRUE` or `FALSE`.*")
})

