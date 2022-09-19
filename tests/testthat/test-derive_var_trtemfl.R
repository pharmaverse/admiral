library(tibble)
library(lubridate)

expected <- tribble(
  ~USUBJID, ~ASTDTM,            ~AENDTM,            ~AEITOXGR, ~AETOXGR, ~TRTEMFL, ~TRTEM2FL,
  # before treatment
  "1",      "2021-12-13T20:15", "2021-12-15T12:45", "1",       "1",      NA,       NA,
  "1",      "2021-12-14T20:15", "2021-12-14T22:00", "1",       "3",      NA,       NA,
  # starting before treatment and ending during treatment
  "1",      "2021-12-30T20:00", "2022-01-14T11:00", "1",       "3",      NA,       "Y",
  "1",      "2021-12-31T20:15", "2022-01-01T01:23", "1",       "1",      NA,       NA,
  # starting during treatment
  "1",      "2022-01-01T12:00", "2022-01-02T23:25", "3",       "4",      "Y",      "Y",
  # after treatment
  "1",      "2022-05-10T11:00", "2022-05-10T13:05", "2",       "2",      "Y",      "Y",
  "1",      "2022-05-11T11:00", "2022-05-11T13:05", "2",       "2",      "Y",      NA,
  # missing dates
  "1",      "",                 "",                 "3",       "4",      "Y",      "Y",
  "1",      "2021-12-30T09:00", "",                 "3",       "4",      NA,       "Y",
  "1",      "2021-12-30T11:00", "",                 "3",       "3",      NA,       NA,
  "1",      "",                 "2022-01-04T09:00", "3",       "4",      "Y",      "Y",
  "1",      "",                 "2021-12-24T19:00", "3",       "4",      NA,       NA,
  "1",      "",                 "2022-06-04T09:00", "3",       "4",      "Y",      "Y",

  # without treatment
  "2",      "",                 "2021-12-03T12:00", "1",       "2",      NA,       NA,
  "2",      "2021-12-01T12:00", "2021-12-03T12:00", "1",       "2",      NA,       NA,
  "2",      "2021-12-06T18:00", "",                 "1",       "2",      NA,       NA
) %>%
  mutate(
    ASTDTM = ymd_hm(ASTDTM),
    AENDTM = ymd_hm(AENDTM),
    TRTSDTM = if_else(USUBJID == "1", ymd_hm("2022-01-01T01:01"), ymd_hms("")),
    TRTEDTM = if_else(USUBJID == "1", ymd_hm("2022-04-30T23:59"), ymd_hms(""))
  )

adae <- select(expected, -starts_with("TRTEM"))

## Test 1: without end_window and worsening ----
test_that("derive_var_trtemfl Test 1: without end_window and worsening", {
  expect_dfs_equal(
    base = select(expected, -TRTEM2FL),
    comp = derive_var_trtemfl(adae),
    keys = c("USUBJID", "ASTDTM", "AENDTM"))
})

## Test 2: with end_window and worsening ----
test_that("derive_var_trtemfl Test 2: with end_window and worsening", {
  expect_dfs_equal(
    base = select(expected, -TRTEMFL),
    comp = derive_var_trtemfl(
      adae,
      new_var = TRTEM2FL,
      trt_end_date = TRTEDTM,
      end_window = 10,
      initial_intensity = AEITOXGR,
      intensity = AETOXGR
      ),
    keys = c("USUBJID", "ASTDTM", "AENDTM"))
})

## Test 3: error if `end_window` without `trt_end_date` ----
test_that("derive_var_trtemfl Test 3: error if `end_window` without `trt_end_date`", {
  expect_error(
    derive_var_trtemfl(
      adae,
      end_window = 10
    ),
    paste(
      "`end_window` argument was specified but not `trt_end_date`",
      "Either both or none of them must be specified.",
      sep = "\n"),
    fixed = TRUE
  )
})

## Test 4: error if `initial_intensity` without `intensity` ----
test_that("derive_var_trtemfl Test 4: error if `initial_intensity` without `intensity`", {
  expect_error(
    derive_var_trtemfl(
      adae,
      initial_intensity = AEITOXGR
    ),
    paste(
      "`initial_intensity` argument was specified but not `intensity`",
      "Either both or none of them must be specified.",
      sep = "\n"),
    fixed = TRUE
  )
})

## Test 5: error if `intensity` without `initial_intensity` ----
test_that("derive_var_trtemfl Test 5: error if `intensity` without `initial_intensity`", {
  expect_error(
    derive_var_trtemfl(
      adae,
      intensity = AETOXGR
    ),
    paste(
      "`intensity` argument was specified but not `initial_intensity`",
      "Either both or none of them must be specified.",
      sep = "\n"),
    fixed = TRUE
  )
})
