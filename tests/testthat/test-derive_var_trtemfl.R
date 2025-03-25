expected <- tibble::tribble(
  ~USUBJID, ~ASTDTM,            ~AENDTM,            ~AEITOXGR, ~AETOXGR, ~TRTEMFL, ~TRTEM2FL, ~TRTEM3FL, # nolint
  # before treatment
  "1",      "2021-12-13T20:15", "2021-12-15T12:45", "1",       "1",      NA,       NA,        NA,
  "1",      "2021-12-14T20:15", "2021-12-14T22:00", "1",       "3",      NA,       NA,        NA,
  # starting before treatment and ending during treatment
  "1",      "2021-12-30T20:00", "2022-01-14T11:00", "1",       "3",      NA,       "Y",       "Y",
  "1",      "2021-12-31T20:15", "2022-01-01T01:23", "1",       "1",      NA,       NA,        NA,
  # starting during treatment
  "1",      "2022-01-01T12:00", "2022-01-02T23:25", "3",       "4",      "Y",      "Y",       "Y",
  # after treatment
  "1",      "2022-05-10T11:00", "2022-05-10T13:05", "2",       "2",      "Y",      "Y",       "Y",
  "1",      "2022-05-10T12:00", "2022-05-10T13:05", "2",       "2",      "Y",      "Y",       NA,
  "1",      "2022-05-11T11:00", "2022-05-11T13:05", "2",       "2",      "Y",      NA,        NA,
  # missing dates
  "1",      "",                 "",                 "3",       "4",      "Y",      "Y",       "Y",
  "1",      "2021-12-30T09:00", "",                 "3",       "4",      NA,       "Y",       "Y",
  "1",      "2021-12-30T11:00", "",                 "3",       "3",      NA,       NA,        NA,
  "1",      "",                 "2022-01-04T09:00", "3",       "4",      "Y",      "Y",       "Y",
  "1",      "",                 "2021-12-24T19:00", "3",       "4",      NA,       NA,        NA,
  "1",      "",                 "2022-06-04T09:00", "3",       "4",      "Y",      "Y",       "Y",

  # without treatment
  "2",      "",                 "2021-12-03T12:00", "1",       "2",      NA,       NA,        NA,
  "2",      "2021-12-01T12:00", "2021-12-03T12:00", "1",       "2",      NA,       NA,        NA,
  "2",      "2021-12-06T18:00", "",                 "1",       "2",      NA,       NA,        NA
) %>%
  mutate(
    ASTDTM = lubridate::ymd_hm(ASTDTM),
    AENDTM = lubridate::ymd_hm(AENDTM),
    TRTSDTM = if_else(USUBJID == "1", lubridate::ymd_hm("2022-01-01T01:01"), ymd_hms("")),
    TRTEDTM = if_else(USUBJID == "1", lubridate::ymd_hm("2022-04-30T11:30"), ymd_hms(""))
  )

adae <- select(expected, -starts_with("TRTEM"))

expected2 <- tribble(
  ~USUBJID, ~ASTDTM, ~AENDTM, ~AEITOXGR, ~AETOXGR, ~AEGRPID, ~TRTEMFL, ~TRTEM2FL,
  # before treatment
  "1", "2021-12-13T20:15", "2021-12-15T12:45", "1", "1", "1", NA, NA,
  "1", "2021-12-14T20:15", "2021-12-14T22:00", "1", "3", "1", NA, NA,
  # starting before treatment and ending during treatment
  "1", "2021-12-30T20:15", "2022-01-14T01:23", "3", "3", "2", NA, NA,
  "1", "2022-01-05T20:00", "2022-06-01T11:00", "3", "1", "2", NA, NA,
  "1", "2022-01-10T20:15", "2022-01-11T01:23", "3", "2", "2", "Y", "Y",
  "1", "2022-01-13T20:15", "2022-03-01T01:23", "3", "1", "2", "Y", "Y",
  # starting during treatment
  "1", "2022-01-01T12:00", "2022-01-02T23:25", "4", "4", "3", "Y", "Y",

  # after treatment
  "1", "2022-05-10T11:00", "2022-05-10T13:05", "2", "2", "4", "Y", "Y",
  "1", "2022-05-10T12:00", "2022-05-10T13:05", "2", "2", "4", "Y", "Y",
  "1", "2022-05-11T11:00", "2022-05-11T13:05", "2", "2", "4", NA, NA,
  # missing dates
  "1", "", "", "3", "4", "5", "Y", "Y",
  "1", "2021-12-30T09:00", "", "3", "4", "5", NA, NA,
  "1", "2021-12-30T11:00", "", "3", "3", "5", NA, NA,
  "1", "", "2022-01-04T09:00", "3", "4", "5", "Y", "Y",
  "1", "", "2021-12-24T19:00", "3", "4", "5", NA, NA,
  "1", "", "2022-06-04T09:00", "3", "4", "5", "Y", "Y",
  # without treatment
  "2", "", "2021-12-03T12:00", "1", "2", "1", NA, NA,
  "2", "2021-12-01T12:00", "2021-12-03T12:00", "1", "2", "2", NA, NA,
  "2", "2021-12-06T18:00", "", "1", "2", "3", NA, NA
) %>%
  mutate(
    STUDYID = "ABC12345",
    ASTDTM = lubridate::ymd_hm(ASTDTM),
    AENDTM = lubridate::ymd_hm(AENDTM),
    TRTSDTM = if_else(USUBJID != "2", lubridate::ymd_hm("2022-01-01T01:01"), ymd_hms("")),
    TRTEDTM = if_else(USUBJID != "2", lubridate::ymd_hm("2022-04-30T23:59"), ymd_hms(""))
  )

adae2 <- select(expected2, -starts_with("TRTEM"))

## Test 1: end_window and worsening parameters not specfied ----
test_that("derive_var_trtemfl Test 1: end_window and worsening parameters not specfied", {
  expect_dfs_equal(
    base = select(expected, -TRTEM2FL, -TRTEM3FL),
    comp = derive_var_trtemfl(adae),
    keys = c("USUBJID", "ASTDTM", "AENDTM")
  )
})

## Test 2: with end_window and worsening ----
test_that("derive_var_trtemfl Test 2: with end_window and worsening", {
  expect_dfs_equal(
    base = select(expected, -TRTEMFL, -TRTEM3FL),
    comp = derive_var_trtemfl(
      adae,
      new_var = TRTEM2FL,
      trt_end_date = TRTEDTM,
      end_window = 10,
      initial_intensity = AEITOXGR,
      intensity = AETOXGR
    ),
    keys = c("USUBJID", "ASTDTM", "AENDTM")
  )
})

## Test 3: with end_window and worsening within grouping variable----
test_that("derive_var_trtemfl Test 3: with end_window and worsening within grouping variable", {
  expect_dfs_equal(
    base = select(expected2, -TRTEM2FL),
    comp = derive_var_trtemfl(
      adae2,
      new_var = TRTEMFL,
      trt_end_date = TRTEDTM,
      end_window = 10,
      intensity = AETOXGR,
      group_var = AEGRPID
    ),
    keys = c("USUBJID", "ASTDTM", "AENDTM")
  )
})


## Test 4: considering trt end time without grouping variable----
test_that("derive_var_trtemfl Test 4: considering trt end time", {
  expect_dfs_equal(
    base = select(expected, -TRTEMFL, -TRTEM2FL),
    comp = derive_var_trtemfl(
      adae,
      new_var = TRTEM3FL,
      trt_end_date = TRTEDTM,
      end_window = 10,
      ignore_time_for_trt_end = FALSE,
      initial_intensity = AEITOXGR,
      intensity = AETOXGR
    ),
    keys = c("USUBJID", "ASTDTM", "AENDTM")
  )
})

## Test 5: considering trt end time with grouping variable----
test_that("derive_var_trtemfl Test 5: considering trt end time", {
  expect_dfs_equal(
    base = select(expected2, -TRTEMFL),
    comp = derive_var_trtemfl(
      adae2,
      new_var = TRTEM2FL,
      trt_end_date = TRTEDTM,
      end_window = 10,
      ignore_time_for_trt_end = FALSE,
      intensity = AETOXGR,
      group_var = AEGRPID
    ),
    keys = c("USUBJID", "ASTDTM", "AENDTM")
  )
})

## Test 6: error if `end_window` without `trt_end_date` ----
test_that("derive_var_trtemfl Test 6: error if `end_window` without `trt_end_date`", {
  expect_snapshot(
    derive_var_trtemfl(
      adae,
      end_window = 10
    ),
    error = TRUE
  )
})

## Test 7: error if `initial_intensity` without `intensity` ----
test_that("derive_var_trtemfl Test 7: error if `initial_intensity` without `intensity`", {
  expect_snapshot(
    derive_var_trtemfl(
      adae,
      initial_intensity = AEITOXGR
    ),
    error = TRUE
  )
})

## Test 8: error if `intensity` without `initial_intensity` ----
test_that("derive_var_trtemfl Test 8: error if `intensity` without `initial_intensity`", {
  expect_snapshot(
    derive_var_trtemfl(
      adae,
      intensity = AETOXGR
    ),
    error = TRUE
  )
})

## Test 9: error if `intensity` without `initial_intensity` ----
test_that("derive_var_trtemfl Test 9: error if `intensity` without `initial_intensity`", {
  expect_snapshot(
    derive_var_trtemfl(
      adae2,
      intensity = AETOXGR
    ),
    error = TRUE
  )
})

## Test 10: warning if both `initial_intensity` and `group_var` are specified ----
test_that("derive_var_trtemfl Test 10: error if `intensity` without `initial_intensity`", {
  expect_warning(
    derive_var_trtemfl(
      adae2,
      initial_intensity = AETOXGR,
      group_var = AEGRPID
    ),
    "`initial_intensity` argument is ignored when `group_var` is specified"
  )
})

## Test 11: error if `group_var` specified without `subject_keys` ----
test_that("derive_var_trtemfl Test 11: error if `group_var` specified without `subject_keys`", {
  expect_snapshot(
    derive_var_trtemfl(
      adae2,
      group_var = AEGRPID,
      subject_keys = NULL
    ),
    error = TRUE
  )
})

# nolint start

## Test 12: Checking test cases from PHUSE White Paper on Treatment Emergent AEs  ----
test_that("derive_var_trtemfl Test 12: Checking test cases from PHUSE White Paper on Treatment Emergent AEs", {
  adae_phuse <- tribble(
    ~USUBJID, ~TRTSDTM,           ~TRTEDTM,           ~ASTDTM,            ~AENDTM,            ~AEITOXGR, ~AETOXGR, ~TRTEMFL,
    # Patient 1: Pre-treatment AE
    "1",      "2021-01-01T00:01", "2021-12-31T23:59", "2020-12-20T00:00", "2020-12-21T00:00", "2",       "2",      NA_character_,
    # Patient 2: On-treatment AE
    "2",      "2021-01-01T00:01", "2021-12-31T23:59", "2021-12-20T00:00", "2021-12-21T00:00", "2",       "2",      "Y",
    # Patient 3: Pre-treatment AE, then on-treatment AE at same intensity
    "3",      "2021-01-01T00:01", "2021-12-31T23:59", "2020-12-20T00:00", "2020-12-21T00:00", "2",       "2",      NA_character_,
    "3",      "2021-01-01T00:01", "2021-12-31T23:59", "2021-12-20T00:00", "2021-12-21T00:00", "2",       "2",      "Y",
    # Patient 4: Pre-treatment AE, then on-treatment AE at worsened intensity
    "4",      "2021-01-01T00:01", "2021-12-31T23:59", "2020-12-20T00:00", "2020-12-21T00:00", "2",       "2",      NA_character_,
    "4",      "2021-01-01T00:01", "2021-12-31T23:59", "2021-12-20T00:00", "2021-12-21T00:00", "2",       "3",      "Y",
    # Patient 5: Pre-treatment AE, then on-treatment AE at improved intensity
    "5",      "2021-01-01T00:01", "2021-12-31T23:59", "2020-12-20T00:00", "2020-12-21T00:00", "2",       "2",      NA_character_,
    "5",      "2021-01-01T00:01", "2021-12-31T23:59", "2021-12-20T00:00", "2021-12-21T00:00", "2",       "1",      "Y",
    # Patient 6: AE starting pre-treatment and continuing on-treatment, then second AE at same intensity
    "6",      "2021-01-01T00:01", "2021-12-31T23:59", "2020-12-23T00:00", "2021-01-21T00:00", "2",       "2",      NA_character_,
    "6",      "2021-01-01T00:01", "2021-12-31T23:59", "2021-12-20T00:00", "2021-12-21T00:00", "2",       "2",      "Y",
    # Patient 7: AE starting pre-treatment and continuing on-treatment, then second AE at worsened intensity
    "7",      "2021-01-01T00:01", "2021-12-31T23:59", "2020-12-23T00:00", "2021-01-21T00:00", "2",       "2",      NA_character_,
    "7",      "2021-01-01T00:01", "2021-12-31T23:59", "2021-12-20T00:00", "2021-12-21T00:00", "2",       "3",      "Y",
    # Patient 8: AE starting pre-treatment and continuing on-treatment, then second AE at improved intensity
    "8",      "2021-01-01T00:01", "2021-12-31T23:59", "2020-12-23T00:00", "2021-01-21T00:00", "2",       "2",      NA_character_,
    "8",      "2021-01-01T00:01", "2021-12-31T23:59", "2021-12-20T00:00", "2021-12-21T00:00", "2",       "1",      "Y",
    # Patient 9: AE starting pre-treatment and continuing on-treatment, and no change in intensity
    "9",      "2021-01-01T00:01", "2021-12-31T23:59", "2020-12-23T00:00", "2021-01-21T00:00", "2",       "2",      NA_character_,
    # Patient 10: AE starting pre-treatment and continuing on-treatment, and worsening intensity
    "10",     "2021-01-01T00:01", "2021-12-31T23:59", "2020-12-23T00:00", "2021-01-21T00:00", "2",       "4",      "Y",
    # Patient 11: AE starting pre-treatment and continuing on-treatment, and improving intensity
    "11",     "2021-01-01T00:01", "2021-12-31T23:59", "2020-12-23T00:00", "2021-01-21T00:00", "2",       "1",      NA_character_,
    # Patient 12: AE starting pre-treatment, worsening, then improving
    "12",     "2021-01-01T00:01", "2021-12-31T23:59", "2020-12-23T00:00", "2021-01-21T00:00", "3",       "2",      NA_character_,
    # Patient 13: AE starting pre-treatment, improving, then worsening
    "13",     "2021-01-01T00:01", "2021-12-31T23:59", "2020-12-23T00:00", "2021-01-21T00:00", "1",       "2",      "Y",
  ) %>%
    mutate(
      ASTDTM = lubridate::ymd_hm(ASTDTM),
      AENDTM = lubridate::ymd_hm(AENDTM),
      TRTSDTM = lubridate::ymd_hm(TRTSDTM),
      TRTEDTM = lubridate::ymd_hm(TRTEDTM),
    )

  # nolint end

  expect_dfs_equal(
    base = adae_phuse,
    comp = derive_var_trtemfl(
      select(adae_phuse, -TRTEMFL),
      new_var = TRTEMFL,
      trt_end_date = TRTEDTM,
      end_window = 0,
      initial_intensity = AEITOXGR,
      intensity = AETOXGR,
      subject_keys = exprs(USUBJID)
    ),
    keys = c("USUBJID", "ASTDTM", "AENDTM")
  )
})
