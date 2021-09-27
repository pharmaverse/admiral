context("test-derive_vars_dtm_to_tm")

input <- tibble::tribble(
  ~STUDYID, ~USUBJID, ~TRTSDTM,               ~ASTDTM,               ~AENDTM,
  "TEST01", "PAT01",  "2012-02-25 23:41:10", "2012-02-28 19:03:00", "2013-02-25 23:32:16",
  "TEST01", "PAT01",  "",                    "2012-02-28 19:00:00", "",
  "TEST01", "PAT01",  "2017-02-25 23:00:02", "2013-02-25 19:00:15", "2014-02-25 19:00:56",
  "TEST01", "PAT01",  "2017-02-25 16:00:00", "2017-02-25 14:25:00", "2017-03-25 23:00:00",
  "TEST01", "PAT01",  "2017-02-25 16:05:17", "2017-02-25 14:20:00", "2018-04-29 14:06:45",
) %>% mutate(
  TRTSDTM = lubridate::as_datetime(TRTSDTM),
  ASTDTM = lubridate::as_datetime(ASTDTM),
  AENDTM = lubridate::as_datetime(AENDTM)
)

# nolint start
#expected output, displaying HH:MM:SS
expected_output2 <- tibble::tribble(
  ~STUDYID, ~USUBJID, ~TRTSDTM,               ~ASTDTM,               ~AENDTM,             ~TRTSTM,   ~ASTTM, ~AENTM,
  "TEST01", "PAT01",  "2012-02-25 23:41:10", "2012-02-28 19:03:00", "2013-02-25 23:32:16", "23:41:10", "19:03:00", "23:32:16",
  "TEST01", "PAT01",  "",                    "2012-02-28 19:00:00", "",                 NA_character_, "19:00:00", NA_character_ ,
  "TEST01", "PAT01",  "2017-02-25 23:00:02", "2013-02-25 19:00:15", "2014-02-25 19:00:56", "23:00:02", "19:00:15", "19:00:56",
  "TEST01", "PAT01",  "2017-02-25 16:00:00", "2017-02-25 14:25:00", "2017-03-25 23:00:00", "16:00:00", "14:25:00", "23:00:00",
  "TEST01", "PAT01",  "2017-02-25 16:05:17", "2017-02-25 14:20:00", "2018-04-29 14:06:45", "16:05:17", "14:20:00", "14:06:45",
) %>% mutate(
  TRTSDTM = lubridate::as_datetime(TRTSDTM),
  ASTDTM = lubridate::as_datetime(ASTDTM),
  AENDTM = lubridate::as_datetime(AENDTM),
  TRTSTM = hms::as_hms(TRTSTM),
  ASTTM = hms::as_hms(ASTTM),
  AENTM = hms::as_hms(AENTM)

)
# nolint end

##############Test
test_that("Convert a complete -- DTM into --TM, TM out is HH:MM:SS", {
  actual_output <- derive_vars_dtm_to_tm(input,
                                         vars(TRTSDTM, ASTDTM, AENDTM))
  expect_dfs_equal(
    expected_output2,
    actual_output,
    keys = c("STUDYID", "USUBJID", "TRTSDTM", "ASTDTM", "AENDTM")
  )
})
