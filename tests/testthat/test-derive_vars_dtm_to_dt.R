#test 1 - multiple vars
test_that("Convert a complete -- DTM into a date object", {
  input <- tibble::tribble(
    ~STUDYID, ~USUBJID, ~TRTSDTM,             ~ASTDTM,               ~AENDTM,
    "TEST01", "PAT01", "2012-02-25 23:00:00", "2012-02-28 19:00:00", "2013-02-25 23:00:00",
    "TEST01", "PAT01", "",                    "2012-02-28 19:00:00", "",
    "TEST01", "PAT01", "2017-02-25 23:00:00", "2013-02-25 19:00:00", "2014-02-25 19:00:00",
    "TEST01", "PAT01", "2017-02-25 16:00:00", "2017-02-25 14:00:00", "2017-03-25 23:00:00",
    "TEST01", "PAT01", "2017-02-25 16:00:00", "2017-02-25 14:00:00", "2018-04-29 14:00:00",
  ) %>% mutate(
    TRTSDTM = lubridate::as_datetime(TRTSDTM),
    ASTDTM = lubridate::as_datetime(ASTDTM),
    AENDTM = lubridate::as_datetime(AENDTM)
  )

  # nolint start
  expected_output <- tibble::tribble(
    ~STUDYID, ~USUBJID, ~TRTSDTM,               ~ASTDTM,               ~AENDTM,             ~TRTSDT,       ~ASTDT,      ~AENDT,
    "TEST01", "PAT01",  "2012-02-25 23:00:00", "2012-02-28 19:00:00", "2013-02-25 23:00:00", "2012-02-25", "2012-02-28","2013-02-25",
    "TEST01", "PAT01",  "",                    "2012-02-28 19:00:00", ""                   , ""          , "2012-02-28", "",
    "TEST01", "PAT01",  "2017-02-25 23:00:00", "2013-02-25 19:00:00", "2014-02-25 19:00:00", "2017-02-25", "2013-02-25", "2014-02-25",
    "TEST01", "PAT01",  "2017-02-25 16:00:00", "2017-02-25 14:00:00", "2017-03-25 23:00:00", "2017-02-25", "2017-02-25", "2017-03-25",
    "TEST01", "PAT01",  "2017-02-25 16:00:00", "2017-02-25 14:00:00", "2018-04-29 14:00:00", "2017-02-25", "2017-02-25", "2018-04-29",
  ) %>% mutate(
    TRTSDTM = lubridate::as_datetime(TRTSDTM),
    ASTDTM = lubridate::as_datetime(ASTDTM),
    AENDTM = lubridate::as_datetime(AENDTM),
    TRTSDT = lubridate::as_date(TRTSDT),
    ASTDT = lubridate::as_date(ASTDT),
    AENDT = lubridate::as_date(AENDT)
  )
  # nolint end

  actual_output <- derive_vars_dtm_to_dt(input, vars(TRTSDTM, ASTDTM, AENDTM))

  expect_dfs_equal(
    expected_output,
    actual_output,
    keys = c("STUDYID", "USUBJID", "TRTSDTM", "ASTDTM", "AENDTM")
  )
})

#test 2 - single var
test_that("Convert a complete -- DTM into a date object", {
  input <- tibble::tribble(
    ~STUDYID, ~USUBJID, ~TRTSDTM,             ~ASTDTM,               ~AENDTM,
    "TEST01", "PAT01", "2012-02-25 23:00:00", "2012-02-28 19:00:00", "2013-02-25 23:00:00",
    "TEST01", "PAT01", "",                    "2012-02-28 19:00:00", "",
    "TEST01", "PAT01", "2017-02-25 23:00:00", "2013-02-25 19:00:00", "2014-02-25 19:00:00",
    "TEST01", "PAT01", "2017-02-25 16:00:00", "2017-02-25 14:00:00", "2017-03-25 23:00:00",
    "TEST01", "PAT01", "2017-02-25 16:00:00", "2017-02-25 14:00:00", "2018-04-29 14:00:00",
  ) %>% mutate(
    TRTSDTM = lubridate::as_datetime(TRTSDTM),
    ASTDTM = lubridate::as_datetime(ASTDTM),
    AENDTM = lubridate::as_datetime(AENDTM)
  )

  # nolint start
  expected_output <- tibble::tribble(
    ~STUDYID, ~USUBJID, ~TRTSDTM,               ~ASTDTM,               ~AENDTM,             ~TRTSDT,
    "TEST01", "PAT01",  "2012-02-25 23:00:00", "2012-02-28 19:00:00", "2013-02-25 23:00:00", "2012-02-25",
    "TEST01", "PAT01",  "",                    "2012-02-28 19:00:00", ""                   , ""          ,
    "TEST01", "PAT01",  "2017-02-25 23:00:00", "2013-02-25 19:00:00", "2014-02-25 19:00:00", "2017-02-25",
    "TEST01", "PAT01",  "2017-02-25 16:00:00", "2017-02-25 14:00:00", "2017-03-25 23:00:00", "2017-02-25",
    "TEST01", "PAT01",  "2017-02-25 16:00:00", "2017-02-25 14:00:00", "2018-04-29 14:00:00", "2017-02-25",
  ) %>% mutate(
    TRTSDTM = lubridate::as_datetime(TRTSDTM),
    ASTDTM = lubridate::as_datetime(ASTDTM),
    AENDTM = lubridate::as_datetime(AENDTM),
    TRTSDT = lubridate::as_date(TRTSDT),
  )
  # nolint end

  actual_output <- derive_vars_dtm_to_dt(input, vars(TRTSDTM))

  expect_dfs_equal(
    expected_output,
    actual_output,
    keys = c("STUDYID", "USUBJID", "TRTSDTM", "ASTDTM", "AENDTM")
  )
})
