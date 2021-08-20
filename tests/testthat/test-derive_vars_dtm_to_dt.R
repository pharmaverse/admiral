context("test-derive_vars_dt")

adcm <- tibble::tribble(
  ~STUDYID, ~USUBJID, ~TRTSDTM, ~ASTDTM, ~AENDTM,
  "TEST01", "PAT01", "2012-02-25 23:00:00", "2012-02-28 19:00:00", "2012-02-25 23:00:00",
  "TEST01", "PAT01", "", "2012-02-28 19:00:00", "",
  "TEST01", "PAT01", "2017-02-25 23:00:00", "2013-02-25 19:00:00", "2014-02-25 19:00:00",
  "TEST01", "PAT01", "2017-02-25 16:00:00", "2017-02-25 14:00:00", "2017-03-25 23:00:00",
  "TEST01", "PAT01", "2017-02-25 16:00:00", "2017-02-25 14:00:00", "2017-04-29 14:00:00",
) %>% dplyr::mutate(
  TRTSDTM = lubridate::as_datetime(TRTSDTM),
  ASTDTM = lubridate::as_datetime(ASTDTM),
  AENDTM = lubridate::as_datetime(AENDTM)
)

test_that("Convert a complete -- DTM into a date object", {
  expected_output <- tibble::tribble(
    ~STUDYID, ~USUBJID, ~TRTSDTM,               ~ASTDTM,               ~AENDTM,             ~TRTSDT, ~ASTDT,~AENDT,
    "TEST01", "PAT01",  "2012-02-25 23:00:00", "2012-03-28 19:00:00", "2012-05-25 23:00:00", "2012-02-25", "2012-02-28","2012-02-25",
    "TEST01", "PAT01",  "",                    "2012-02-28 19:00:00", ""                   , ""          , "2012-02-28", "",
    "TEST01", "PAT01",  "2017-02-25 23:00:00", "2013-02-25 19:00:00", "2014-02-25 19:00:00", "2017-02-25", "2013-02-25", "2014-02-25",
    "TEST01", "PAT01",  "2017-02-25 23:00:00", "2017-02-25 14:00:00", "2017-03-25 23:00:00", "2017-02-25", "2017-02-25", "2017-03-25",
    "TEST01", "PAT01",  "2017-02-25 23:00:00", "2017-01-25 14:00:00", "2018-04-29 14:00:00", "2017-02-25", "2017-02-25", "2017-04-29",
  )%>% dplyr::mutate(TRTSDTM = lubridate::as_datetime(TRTSDTM),
                     ASTDTM = lubridate::as_datetime(ASTDTM),
                     AENDTM = lubridate::as_datetime(AENDTM),
                     TRTSDT = lubridate::as_date(TRTSDT),
                     ASTDT = lubridate::as_date(ASTDT),
                     AENDT = lubridate::as_date(AENDT))

  actual_output <- derive_vars_dtm_to_dt(dataset = adcm,
                   source_vars = vars(TRTSDTM,ASTDTM,AENDTM))
  expect_dfs_equal(
    expected_output,
    actual_output,
    keys = c("STUDYID","USUBJID","TRTSDTM")
  )

})
