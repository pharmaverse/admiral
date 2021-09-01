context("test-derive_var_atirel")

test_that("Derive ATIREL", {
  input <- tibble::tribble(
    ~STUDYID, ~USUBJID, ~TRTSDTM,               ~ASTDTM,               ~AENDTM,             ~ASTTMF,
    "TEST01", "PAT01",  "2012-02-25 23:00:00", "2012-03-28 19:00:00", "2012-05-25 23:00:00", "",
    "TEST01", "PAT01",  "",                    "2012-02-28 19:00:00", "",                    "",
    "TEST01", "PAT01",  "2017-02-25 23:00:00", "2013-02-25 19:00:00", "2014-02-25 19:00:00", "",
    "TEST01", "PAT01",  "2017-02-25 23:00:00", "2017-02-25 14:00:00", "2017-03-25 23:00:00", "m",
    "TEST01", "PAT01",  "2017-02-25 23:00:00", "2017-01-25 14:00:00", "2018-04-29 14:00:00", "",
  ) %>% mutate(
    TRTSDTM = lubridate::as_datetime(TRTSDTM),
    ASTDTM = lubridate::as_datetime(ASTDTM),
    AENDTM = lubridate::as_datetime(AENDTM)
  )

  # nolint start
  expected_output <- tibble::tribble(
    ~STUDYID, ~USUBJID, ~TRTSDTM,               ~ASTDTM,              ~AENDTM,               ~ASTTMF, ~ATIREL,
    "TEST01", "PAT01",  "2012-02-25 23:00:00", "2012-03-28 19:00:00", "2012-05-25 23:00:00", "",      "CONCOMITANT",
    "TEST01", "PAT01",  "",                    "2012-02-28 19:00:00", "",                    "",      NA_character_,
    "TEST01", "PAT01",  "2017-02-25 23:00:00", "2013-02-25 19:00:00", "2014-02-25 19:00:00", "",      "PRIOR",
    "TEST01", "PAT01",  "2017-02-25 23:00:00", "2017-02-25 14:00:00", "2017-03-25 23:00:00", "m",     "CONCOMITANT",
    "TEST01", "PAT01",  "2017-02-25 23:00:00", "2017-01-25 14:00:00", "2018-04-29 14:00:00", "",      "PRIOR_CONCOMITANT"
  ) %>% mutate(
    TRTSDTM = lubridate::as_datetime(TRTSDTM),
    ASTDTM = lubridate::as_datetime(ASTDTM),
    AENDTM = lubridate::as_datetime(AENDTM)
  )
  # nolint end

  actual_output <- derive_var_atirel(
    dataset = input,
    flag_var = ASTTMF,
    new_var = ATIREL
  )

  expect_dfs_equal(
    expected_output,
    actual_output,
    keys = c("STUDYID", "USUBJID", "ASTDTM")
  )
})
