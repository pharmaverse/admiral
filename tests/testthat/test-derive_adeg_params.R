# derive_param_qtc ----
## Test 1: Bazett's method ----
test_that("derive_param_qtc Test 1: Bazett's method", {
  input <- tibble::tribble(
    ~USUBJID,      ~PARAMCD, ~PARAM,        ~AVAL, ~AVALU, ~VISIT,
    "01-701-1015", "QT",     "QT Duration",   370, "ms",   "WEEK 2",
    "01-701-1015", "RR",     "RR Duration",   710, "ms",   "WEEK 2",
    "01-701-1028", "QT",     "QT Duration",   480, "ms",   "WEEK 2",
    "01-701-1028", "QT",     "QT Duration",   350, "ms",   "WEEK 3",
    "01-701-1028", "RR",     "RR Duration",   842, "ms",   "WEEK 2",
  )

  expected <- bind_rows(
    input,
    tibble::tribble(
      ~USUBJID,      ~PARAMCD, ~VISIT,   ~AVAL,
      "01-701-1015", "QTCBR",  "WEEK 2", 370 / sqrt(710 / 1000),
      "01-701-1028", "QTCBR",  "WEEK 2", 480 / sqrt(842 / 1000)
    )
  )
  actual <- derive_param_qtc(
    input,
    by_vars = exprs(USUBJID, VISIT),
    method = "Bazett",
    get_unit_expr = AVALU
  )

  expect_dfs_equal(
    base = expected,
    compare = actual,
    keys = c("USUBJID", "PARAMCD", "VISIT")
  )
})

## Test 2: Fridericia's method ----
test_that("derive_param_qtc Test 2: Fridericia's method", {
  input <- tibble::tribble(
    ~USUBJID,      ~PARAMCD, ~PARAM,        ~AVAL, ~AVALU, ~VISIT,
    "01-701-1015", "QT",     "QT Duration",   370, "msec", "WEEK 2",
    "01-701-1015", "RR",     "RR Duration",   710, "msec", "WEEK 2",
    "01-701-1028", "QT",     "QT Duration",   480, "msec", "WEEK 2",
    "01-701-1028", "QT",     "QT Duration",   350, "msec", "WEEK 3",
    "01-701-1028", "RR",     "RR Duration",   842, "msec", "WEEK 2",
  )

  expected <- bind_rows(
    input,
    tibble::tribble(
      ~USUBJID,      ~PARAMCD, ~VISIT,   ~AVAL,
      "01-701-1015", "QTCFR",  "WEEK 2", 370 / 710^(1 / 3) * 10,
      "01-701-1028", "QTCFR",  "WEEK 2", 480 / 842^(1 / 3) * 10
    )
  )
  actual <- derive_param_qtc(
    input,
    by_vars = exprs(USUBJID, VISIT),
    method = "Fridericia",
    get_unit_expr = AVALU
  )

  expect_dfs_equal(
    base = expected,
    compare = actual,
    keys = c("USUBJID", "PARAMCD", "VISIT")
  )
})

## Test 3: Sagie's method ----
test_that("derive_param_qtc Test 3: Sagie's method", {
  input <- tibble::tribble(
    ~USUBJID,      ~PARAMCD, ~PARAM,        ~AVAL, ~AVALU, ~VISIT,
    "01-701-1015", "QT",     "QT Duration",   370, "ms",   "WEEK 2",
    "01-701-1015", "RR",     "RR Duration",   710, "ms",   "WEEK 2",
    "01-701-1028", "QT",     "QT Duration",   480, "ms",   "WEEK 2",
    "01-701-1028", "QT",     "QT Duration",   350, "ms",   "WEEK 3",
    "01-701-1028", "RR",     "RR Duration",   842, "ms",   "WEEK 2",
  )

  expected <- bind_rows(
    input,
    tibble::tribble(
      ~USUBJID,      ~PARAMCD, ~VISIT,   ~AVAL,
      "01-701-1015", "QTLCR",  "WEEK 2", 370 + 154 * (1 - 710 / 1000),
      "01-701-1028", "QTLCR",  "WEEK 2", 480 + 154 * (1 - 842 / 1000)
    )
  )
  actual <- derive_param_qtc(
    input,
    by_vars = exprs(USUBJID, VISIT),
    method = "Sagie",
    get_unit_expr = AVALU
  )

  expect_dfs_equal(
    base = expected,
    compare = actual,
    keys = c("USUBJID", "PARAMCD", "VISIT")
  )
})

## Test 4: Message if no new records ----
test_that("derive_param_qtc Test 4: Message if no new records", {
  input <- tibble::tribble(
    ~USUBJID,      ~PARAMCD, ~PARAM,        ~AVAL, ~AVALU, ~VISIT,
    "01-701-1015", "QT",     "QT Duration",   370, "ms",   "WEEK 2",
    "01-701-1015", "RR",     "RR Duration",    NA, "ms",   "WEEK 2",
    "01-701-1028", "QT",     "QT Duration",   480, "ms",   "WEEK 2",
    "01-701-1028", "QT",     "QT Duration",   350, "ms",   "WEEK 3",
    "01-701-1028", "RR",     "RR Duration",    NA, "ms",   "WEEK 2",
  )

  expect_snapshot(
  actual <- derive_param_qtc(
    input,
    by_vars = exprs(USUBJID, VISIT),
    method = "Bazett",
    get_unit_expr = AVALU
  )
  )

  expect_dfs_equal(
    base = input,
    compare = actual,
    keys = c("USUBJID", "PARAMCD", "VISIT")
  )
})

# derive_param_rr ----

## Test 5: new observations are derived correctly ----
test_that("derive_param_rr Test 5: new observations are derived correctly", {
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

## Test 6: Message if no new records ----
test_that("derive_param_rr Test 6: Message if no new records", {
  input <- tibble::tribble(
    ~USUBJID,      ~PARAMCD, ~PARAM,        ~AVAL, ~AVALU,      ~VISIT,
    "01-701-1015", "HR",     "Heart Rate",  NA   , "beats/min", "BASELINE",
    "01-701-1015", "HR",     "Heart Rate",  NA   , "beats/min", "WEEK 1",
    "01-701-1015", "RR",     "RR Duration", 710,   "ms",        "WEEK 2"
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
