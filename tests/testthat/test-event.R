## Test 1: Cover event$order ----
test_that("event Test 1: Cover event$order", {
  adsl_ext <- tribble(
    ~STUDYID, ~USUBJID, ~TRTEDT, ~DTHDT,
    "PILOT01", "01-1130", ymd("2014-08-16"), ymd("2014-09-13"),
    "PILOT01", "01-1133", ymd("2013-04-28"), ymd(""),
    "PILOT01", "01-1211", ymd("2013-01-12"), ymd(""),
    "PILOT01", "09-1081", ymd("2014-04-27"), ymd(""),
    "PILOT01", "09-1088", ymd("2014-10-09"), ymd("2014-11-01"),
  )

  lb_ext <- tribble(
    ~STUDYID,  ~DOMAIN,  ~USUBJID, ~LBSEQ,             ~LBDTC,
    "PILOT01",    "LB", "01-1130",    219, "2014-06-07T13:20",
    "PILOT01",    "LB", "01-1130",    322, "2014-08-16T13:10",
    "PILOT01",    "LB", "01-1133",    268, "2013-04-18T15:30",
    "PILOT01",    "LB", "01-1133",    304, "2013-05-01T10:13",
    "PILOT01",    "LB", "01-1211",      8, "2012-10-30T14:26",
    "PILOT01",    "LB", "01-1211",    162, "2013-01-08T12:13",
    "PILOT01",    "LB", "09-1081",     47, "2014-02-01T10:55",
    "PILOT01",    "LB", "09-1081",    219, "2014-05-10T11:15",
    "PILOT01",    "LB", "09-1088",    283, "2014-09-27T12:13",
    "PILOT01",    "LB", "09-1088",    322, "2014-10-09T13:25"
  ) %>%
    mutate(
      ADT = convert_dtc_to_dt(LBDTC)
    )

  expected_output <- tribble(
    ~STUDYID, ~USUBJID, ~TRTEDT, ~DTHDT, ~LSTALVDT, ~DTHFL,
    "PILOT01", "01-1130", ymd("2014-08-16"), ymd("2014-09-13"), ymd("2014-09-13"), "Y",
    "PILOT01", "01-1133", ymd("2013-04-28"), ymd(""), ymd("2013-05-01"), "N",
    "PILOT01", "01-1211", ymd("2013-01-12"), ymd(""), ymd("2013-01-12"), "N",
    "PILOT01", "09-1081", ymd("2014-04-27"), ymd(""), ymd("2014-05-10"), "N",
    "PILOT01", "09-1088", ymd("2014-10-09"), ymd("2014-11-01"), ymd("2014-11-01"), "Y",
  )

  actual_output <- derive_vars_extreme_event(
    adsl_ext,
    by_vars = exprs(STUDYID, USUBJID),
    events = list(
      event(
        dataset_name = "adsl_ext",
        condition = !is.na(DTHDT),
        set_values_to = exprs(LSTALVDT = DTHDT, DTHFL = "Y")
      ),
      event(
        dataset_name = "lb_ext",
        condition = !is.na(ADT),
        order = exprs(ADT, LBSEQ),
        mode = "last",
        set_values_to = exprs(LSTALVDT = ADT, DTHFL = "N")
      ),
      event(
        dataset_name = "adsl_ext",
        condition = !is.na(TRTEDT),
        order = exprs(TRTEDT),
        mode = "last",
        set_values_to = exprs(LSTALVDT = TRTEDT, DTHFL = "N")
      )
    ),
    source_datasets = list(adsl_ext = adsl_ext, lb_ext = lb_ext),
    tmp_event_nr_var = event_nr,
    order = exprs(LSTALVDT, event_nr),
    mode = "last",
    new_vars = exprs(LSTALVDT = LSTALVDT, DTHFL = DTHFL)
  )


  expect_dfs_equal(
    base = expected_output,
    compare = actual_output,
    keys = c("USUBJID")
  )
})
