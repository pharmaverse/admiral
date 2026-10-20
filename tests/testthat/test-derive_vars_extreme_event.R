## Test 1: derive_vars_extreme_event ----
test_that("derive_vars_extreme_event Test 1: derive_vars_extreme_event", {
  adsl <- tribble(
    ~STUDYID, ~USUBJID, ~TRTEDT, ~DTHDT,
    "PILOT01", "01-1130", ymd("2014-08-16"), ymd("2014-09-13"),
    "PILOT01", "01-1133", ymd("2013-04-28"), ymd(""),
    "PILOT01", "01-1211", ymd("2013-01-12"), ymd(""),
    "PILOT01", "09-1081", ymd("2014-04-27"), ymd(""),
    "PILOT01", "09-1088", ymd("2014-10-09"), ymd("2014-11-01"),
  )

  actual <- derive_vars_extreme_event(
    adsl,
    by_vars = exprs(STUDYID, USUBJID),
    events = list(
      event(
        dataset_name = "adsl",
        condition = !is.na(DTHDT),
        set_values_to = exprs(LSTALVDT = DTHDT, DTHFL = "Y")
      ),
      event(
        dataset_name = "adsl",
        condition = !is.na(TRTEDT),
        set_values_to = exprs(LSTALVDT = TRTEDT, DTHFL = "N")
      )
    ),
    source_datasets = list(adsl = adsl),
    order = exprs(LSTALVDT),
    mode = "last",
    new_vars = exprs(LSTALVDT = LSTALVDT, DTHFL = DTHFL)
  )

  expected <- tribble(
    ~STUDYID, ~USUBJID, ~TRTEDT, ~DTHDT, ~LSTALVDT, ~DTHFL,
    "PILOT01", "01-1130", ymd("2014-08-16"), ymd("2014-09-13"), ymd("2014-09-13"), "Y",
    "PILOT01", "01-1133", ymd("2013-04-28"), ymd(""), ymd("2013-04-28"), "N",
    "PILOT01", "01-1211", ymd("2013-01-12"), ymd(""), ymd("2013-01-12"), "N",
    "PILOT01", "09-1081", ymd("2014-04-27"), ymd(""), ymd("2014-04-27"), "N",
    "PILOT01", "09-1088", ymd("2014-10-09"), ymd("2014-11-01"), ymd("2014-11-01"), "Y",
  )

  expect_dfs_equal(
    base = expected,
    compare = actual,
    keys = c("STUDYID", "USUBJID")
  )
  expect_s3_class(actual, "admiral_df")
})

## Test 2: source precedence, study keys, and unmatched subjects ----
test_that("derive_vars_extreme_event Test 2: source precedence and subject retention", {
  adsl <- tibble::tribble(
    ~STUDYID, ~USUBJID,
    "A",      "1",
    "A",      "2",
    "B",      "1"
  )
  reports <- tibble::tribble(
    ~STUDYID, ~USUBJID, ~SOURCE, ~CAUSE,      ~ADT,
    "A",      "1",      "DS",    "DS cause",  as.Date("2024-01-01"),
    "A",      "1",      "AE",    "AE cause",  as.Date("2024-01-01"),
    "B",      "1",      "DS",    "B cause",   as.Date("2024-02-01")
  )

  actual <- derive_vars_extreme_event(
    adsl,
    by_vars = exprs(STUDYID, USUBJID),
    events = list(
      event(dataset_name = "reports", condition = SOURCE == "AE"),
      event(dataset_name = "reports", condition = SOURCE == "DS")
    ),
    source_datasets = list(reports = reports),
    tmp_event_nr_var = event_nr,
    order = exprs(ADT, event_nr),
    mode = "first",
    check_type = "error",
    new_vars = exprs(DTHCAUS = CAUSE, DTHDT = ADT)
  )
  expected <- dplyr::mutate(
    adsl,
    DTHCAUS = c("AE cause", NA_character_, "B cause"),
    DTHDT = as.Date(c("2024-01-01", NA_character_, "2024-02-01"))
  )
  expect_dfs_equal(actual, expected, keys = c("STUDYID", "USUBJID"))
})

## Test 3: ambiguous reports stop derivation ----
test_that("derive_vars_extreme_event Test 3: ambiguous reports stop derivation", {
  adsl <- tibble::tibble(USUBJID = "1")
  reports <- tibble::tribble(
    ~USUBJID, ~ADT,                  ~CAUSE,
    "1",      as.Date("2024-01-01"), "Cause A",
    "1",      as.Date("2024-01-01"), "Cause B"
  )

  expect_error(
    derive_vars_extreme_event(
      adsl,
      by_vars = exprs(USUBJID),
      events = list(event(dataset_name = "reports")),
      source_datasets = list(reports = reports),
      order = exprs(ADT),
      mode = "first",
      check_type = "error",
      new_vars = exprs(DTHCAUS = CAUSE, DTHDT = ADT)
    ),
    class = "duplicate_records"
  )
})

## Test 4: confirmation requires consecutive days above the threshold ----
test_that("derive_vars_extreme_event Test 4: consecutive-day confirmation", {
  adsl <- tibble::tibble(USUBJID = c("1", "2", "3", "4", "5", "6"))
  adlb <- tibble::tribble(
    ~USUBJID, ~ADT, ~AVAL, ~ANRHI,
    "1", as.Date("2024-01-01"), 160, 40,
    "1", as.Date("2024-01-02"), 200, 40,
    "1", as.Date("2024-01-03"), 240, 40,
    "2", as.Date("2024-01-01"), 160, 40,
    "2", as.Date("2024-01-03"), 200, 40,
    "3", as.Date("2024-01-01"), 120, 40,
    "3", as.Date("2024-01-02"), 200, 40,
    "4", as.Date("2024-01-01"), 160, NA_real_,
    "4", as.Date("2024-01-02"), 200, 40,
    "6", as.Date("2024-01-01"), 160, 40,
    "6", as.Date("2024-01-02"), NA_real_, 40
  )

  actual <- derive_vars_extreme_event(
    adsl,
    by_vars = exprs(USUBJID),
    events = list(event_joined(
      dataset_name = "adlb",
      order = exprs(ADT),
      join_vars = exprs(ADT, AVAL, ANRHI),
      join_type = "before",
      condition = AVAL > 3 * ANRHI & AVAL.join > 3 * ANRHI.join &
        ADT == ADT.join + 1
    )),
    source_datasets = list(adlb = adlb),
    order = exprs(ADT),
    mode = "first",
    new_vars = exprs(ALTCDT = ADT)
  )
  expected <- dplyr::mutate(
    adsl,
    ALTCDT = as.Date(c("2024-01-02", rep(NA_character_, 5)))
  )
  expect_dfs_equal(actual, expected, keys = "USUBJID")
})
