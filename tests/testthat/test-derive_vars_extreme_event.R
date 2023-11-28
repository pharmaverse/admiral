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
})
