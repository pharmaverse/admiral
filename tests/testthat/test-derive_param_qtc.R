test_that("new observations are derived correctly", {
  input <- tibble::tribble(
    ~USUBJID,      ~PARAMCD, ~PARAM,        ~AVAL,  ~AVALU,      ~VISIT,
    "01-701-1015", "HR",     "Heart Rate",   70.14, "beats/min", "BASELINE",
    "01-701-1015", "QT",     "QT Duration", 370,    "msec",      "WEEK 2",
    "01-701-1015", "HR",     "Heart Rate",   62.66, "beats/min", "WEEK 1",
    "01-701-1015", "RR",     "RR Duration", 710,    "msec",      "WEEK 2",
    "01-701-1028", "HR",     "Heart Rate",   85.45, "beats/min", "BASELINE",
    "01-701-1028", "QT",     "QT Duration", 480,    "msec",      "WEEK 2",
    "01-701-1028", "QT",     "QT Duration", 350,    "msec",      "WEEK 3",
    "01-701-1028", "HR",     "Heart Rate",   56.54, "beats/min", "WEEK 3",
    "01-701-1028", "RR",     "RR Duration", 842,    "msec",      "WEEK 2",
  )
  by <- vars(USUBJID, VISIT)
  output_bazett <- derive_param_qtc(input, by_vars = by, method = "Bazett")
  output_fridericia <- derive_param_qtc(input, by_vars = by, method = "Fridericia")
  output_sagie <- derive_param_qtc(input, by_vars = by, method = "Sagie")

  expect_identical(nrow(output_bazett), nrow(input) + 2L)
  expect_identical(nrow(output_fridericia), nrow(input) + 2L)
  expect_identical(nrow(output_sagie), nrow(input) + 2L)

  expect_identical(
    output_bazett %>% filter(PARAMCD == "QTCBR") %>% pull(AVAL),
    compute_qtcb(c(370, 480), c(710, 842))
  )
  expect_identical(
    output_fridericia %>% filter(PARAMCD == "QTCFR") %>% pull(AVAL),
    compute_qtcf(c(370, 480), c(710, 842))
  )
  expect_identical(
    output_sagie %>% filter(PARAMCD == "QTLCR") %>% pull(AVAL),
    compute_qtlc(c(370, 480), c(710, 842))
  )
})
