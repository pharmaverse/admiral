context("test-derive_param_qtlc")


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

  new_obs <-
    inner_join(input %>% filter(PARAMCD == "QT") %>% select(USUBJID, VISIT, AVAL),
               input %>% filter(PARAMCD == "RR") %>% select(USUBJID, VISIT, AVAL),
               by = c("USUBJID", "VISIT"),
               suffix = c(".QT", ".RR")) %>%
    mutate(AVAL = 1000 * (AVAL.QT / 1000 + 0.154 * (1 - AVAL.RR / 1000)),
           PARAMCD = "QTLCR",
           AVALU = "msec") %>%
    select(-AVAL.QT, -AVAL.RR)
  expected_output <- bind_rows(input, new_obs)

  expect_dfs_equal(derive_param_qtlc(input,
                                     unit_var = AVALU,
                                     by_vars = vars(USUBJID, VISIT)),
                   expected_output,
                   keys = c("USUBJID", "PARAMCD", "VISIT"))
})

test_that("new observations are derived correctly without unit", {
  input <- tibble::tribble(
    ~USUBJID,      ~PARAMCD, ~PARAM,        ~AVAL, ~VISIT,
    "01-701-1015", "HR",     "Heart Rate",  70.14, "BASELINE",
    "01-701-1015", "QT",     "QT Duration", 370,   "WEEK 2",
    "01-701-1015", "HR",     "Heart Rate",  62.66, "WEEK 1",
    "01-701-1015", "RR",     "RR Duration", 710,   "WEEK 2",
    "01-701-1028", "HR",     "Heart Rate",  85.45, "BASELINE",
    "01-701-1028", "QT",     "QT Duration", 480,   "WEEK 2",
    "01-701-1028", "QT",     "QT Duration", 350,   "WEEK 3",
    "01-701-1028", "HR",     "Heart Rate",  56.54, "WEEK 3",
    "01-701-1028", "RR",     "RR Duration", 842,   "WEEK 2",
  )

  new_obs <-
    inner_join(input %>% filter(PARAMCD == "QT") %>% select(USUBJID, VISIT, AVAL),
               input %>% filter(PARAMCD == "RR") %>% select(USUBJID, VISIT, AVAL),
               by = c("USUBJID", "VISIT"),
               suffix = c(".QT", ".RR")) %>%
    mutate(AVAL = 1000 * (AVAL.QT / 1000 + 0.154 * (1 - AVAL.RR / 1000)),
           PARAMCD = "QTLCR") %>%
    select(-AVAL.QT, -AVAL.RR)
  expected_output <- bind_rows(input, new_obs)

  expect_dfs_equal(derive_param_qtlc(input,
                                     by_vars = vars(USUBJID, VISIT)),
                   expected_output,
                   keys = c("USUBJID", "PARAMCD", "VISIT"))
})
