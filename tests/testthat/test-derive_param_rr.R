# derive_param_rr ----

## Test 1: new observations are derived correctly ----
test_that("derive_param_rr Test 1: new observations are derived correctly", {
  input <- tibble::tribble(
    ~USUBJID,      ~PARAMCD, ~PARAM,        ~AVAL, ~AVALU,      ~VISIT,
    "01-701-1015", "HR",     "Heart Rate",  70.14, "beats/min", "BASELINE",
    "01-701-1015", "HR",     "Heart Rate",  62.66, "beats/min", "WEEK 1",
    "01-701-1015", "RR",     "RR Duration", 710,   "ms",        "WEEK 2",
    "01-701-1028", "HR",     "Heart Rate",  85.45, "beats/min", "BASELINE",
    "01-701-1028", "HR",     "Heart Rate",  56.54, "beats/min", "WEEK 3",
    "01-701-1028", "RR",     "RR Duration", 842,   "ms",        "WEEK 2",
  )

  new_obs <- input %>%
    filter(PARAMCD == "HR") %>%
    select(USUBJID, VISIT, AVAL) %>%
    mutate(
      AVAL = 60000 / AVAL,
      PARAMCD = "RRR",
      PARAM = "RR Duration Rederived (ms)",
      AVALU = "ms"
    )
  expected_output <- bind_rows(input, new_obs)

  expect_dfs_equal(
    derive_param_rr(
      input,
      by_vars = exprs(USUBJID, VISIT),
      set_values_to = exprs(
        PARAMCD = "RRR",
        PARAM = "RR Duration Rederived (ms)",
        AVALU = "ms"
      ),
      get_unit_expr = AVALU
    ),
    expected_output,
    keys = c("USUBJID", "PARAMCD", "VISIT")
  )
})

## Test 2: Message if no new records ----
test_that("derive_param_rr Test 2: Message if no new records", {
  input <- tibble::tribble(
    ~USUBJID,      ~PARAMCD, ~PARAM,        ~AVAL, ~AVALU,      ~VISIT,
    "01-701-1015", "HR",     "Heart Rate",     NA, "beats/min", "BASELINE",
    "01-701-1015", "HR",     "Heart Rate",     NA, "beats/min", "WEEK 1",
    "01-701-1015", "RR",     "RR Duration",   710, "ms",        "WEEK 2"
  )


  expect_snapshot(
    actual <- derive_param_rr(
      input,
      by_vars = exprs(USUBJID, VISIT),
      set_values_to = exprs(
        PARAMCD = "RRR",
        PARAM = "RR Duration Rederived (ms)",
        AVALU = "ms"
      ),
      get_unit_expr = AVALU
    )
  )

  expect_dfs_equal(
    base = input,
    compare = actual,
    keys = c("USUBJID", "PARAMCD", "VISIT")
  )
})
