test_that("new observations are derived correctly", {
  input <- tibble::tribble(
    ~USUBJID, ~PARAMCD, ~PARAM, ~AVAL, ~AVALU, ~VISIT,
    "01-701-1015", "HR", "Heart Rate", 70.14, "beats/min", "BASELINE",
    "01-701-1015", "QT", "QT Duration", 370, "msec", "WEEK 2",
    "01-701-1015", "HR", "Heart Rate", 62.66, "beats/min", "WEEK 1",
    "01-701-1015", "RR", "RR Duration", 710, "msec", "WEEK 2",
    "01-701-1028", "HR", "Heart Rate", 85.45, "beats/min", "BASELINE",
    "01-701-1028", "QT", "QT Duration", 480, "msec", "WEEK 2",
    "01-701-1028", "QT", "QT Duration", 350, "msec", "WEEK 3",
    "01-701-1028", "HR", "Heart Rate", 56.54, "beats/min", "WEEK 3",
    "01-701-1028", "RR", "RR Duration", 842, "msec", "WEEK 2",
  )
  methods <- c("Bazett", "Fridericia", "Sagie")
  outputs <- lapply(methods, function(method) {
    derive_param_qtc(
      input,
      by_vars = vars(USUBJID, VISIT),
      method = method,
      get_unit_expr = AVALU
    )
  })
  names(outputs) <- methods

  expect_identical(nrow(outputs$Bazett), nrow(input) + 2L)
  expect_identical(nrow(outputs$Fridericia), nrow(input) + 2L)
  expect_identical(nrow(outputs$Sagie), nrow(input) + 2L)

  expect_identical(
    outputs$Bazett %>% filter(PARAMCD == "QTCBR") %>% pull(AVAL),
    compute_qtc(c(370, 480), c(710, 842), method = "Bazett")
  )
  expect_identical(
    outputs$Fridericia %>% filter(PARAMCD == "QTCFR") %>% pull(AVAL),
    compute_qtc(c(370, 480), c(710, 842), method = "Fridericia")
  )
  expect_identical(
    outputs$Sagie %>% filter(PARAMCD == "QTLCR") %>% pull(AVAL),
    compute_qtc(c(370, 480), c(710, 842), method = "Sagie")
  )
})
