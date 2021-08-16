context("test-derive_param_qtcb")


test_that("new observations are derived correctly", {
  input <- tibble::tribble(
    ~USUBJID,      ~PARAMCD, ~PARAM,        ~AVAL, ~AVALU,      ~VISIT,
    "01-701-1015", "HR",     "Heart Rate",  70.14, "beats/min", "BASELINE",
    "01-701-1015", "QT",     "QT Duration", 370,   "msec",      "WEEK 2",
    "01-701-1015", "HR",     "Heart Rate",  62.66, "beats/min", "WEEK 1",
    "01-701-1015", "RR",     "RR Duration", 710,   "msec",      "WEEK 2",
    "01-701-1028", "HR",     "Heart Rate",  85.45, "beats/min", "BASELINE",
    "01-701-1028", "QT",     "QT Duration", 480,   "msec",      "WEEK 2",
    "01-701-1028", "QT",     "QT Duration", 350,   "msec",      "WEEK 3",
    "01-701-1028", "HR",     "Heart Rate",  56.54, "beats/min", "WEEK 3",
    "01-701-1028", "RR",     "RR Duration", 842,   "msec",      "WEEK 2",
  )

  new_obs <- input %>%
    filter(PARAMCD == "HR") %>%
    select(USUBJID, VISIT, AVAL) %>%
    mutate(
      AVAL = 60000 / AVAL,
      PARAMCD = "RRR",
      PARAM = "RR Duration Rederived (msec)",
      AVALU = "msec"
    )
  expected_output <- bind_rows(input, new_obs)

  expect_dfs_equal(
    derive_param_rr(
      input,
      by_vars = vars(USUBJID, VISIT),
      set_values_to = vars(PARAMCD = "RRR", PARAM = "RR Duration Rederived (msec)"),
      unit_var = AVALU
    ),
    expected_output,
    keys = c("USUBJID", "PARAMCD", "VISIT")
  )
})

test_that("new observations are derived correctly without unit", {
  input <- tibble::tribble(
    ~USUBJID,      ~PARAMCD, ~PARAM,        ~AVAL, ~VISIT,
    "01-701-1015", "HR",     "Heart Rate",  70.14, "BASELINE",
    "01-701-1015", "HR",     "Heart Rate",  62.66, "WEEK 1",
    "01-701-1028", "HR",     "Heart Rate",  85.45, "BASELINE",
    "01-701-1028", "HR",     "Heart Rate",  56.54, "WEEK 3",
  )

  new_obs <- input %>%
    filter(PARAMCD == "HR") %>%
    select(USUBJID, VISIT, AVAL) %>%
    mutate(
      AVAL = 60000 / AVAL,
      PARAMCD = "RRR",
      PARAM = "RR Duration Rederived (msec)"
    )
  expected_output <- bind_rows(input, new_obs)

  expect_dfs_equal(
    derive_param_rr(
      input,
      by_vars = vars(USUBJID, VISIT),
      set_values_to = vars(PARAMCD = "RRR", PARAM = "RR Duration Rederived (msec)")
    ),
    expected_output,
    keys = c("USUBJID", "PARAMCD", "VISIT")
  )
})
