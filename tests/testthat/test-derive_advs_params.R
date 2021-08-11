context("test-derive_advs_params")

input <- tibble::tribble(
  ~USUBJID,      ~PARAMCD,     ~PARAM,        ~AVAL, ~AVALU,  ~VISIT,
  "01-701-1015", "HEIGHT",     "Height (cm)", 170,   "cm",   "BASELINE",
  "01-701-1015", "WEIGHT",     "Weight (kg)",  75,   "kg",   "BASELINE",
  "01-701-1015", "WEIGHT",     "Weight (kg)",  78,   "kg",   "MONTH 1",
  "01-701-1015", "WEIGHT",     "Weight (kg)",  80,   "kg",   "MONTH 2",
  "01-701-1028", "HEIGHT",     "Height (cm)", 185,   "cm",   "BASELINE",
  "01-701-1028", "WEIGHT",     "Weight (kg)",  90,   "kg",   "BASELINE",
  "01-701-1028", "WEIGHT",     "Weight (kg)",  88,   "kg",   "MONTH 1",
  "01-701-1028", "WEIGHT",     "Weight (kg)",  85,   "kg",   "MONTH 2",
)

test_that("new observations are derived correctly with Mosteller method", {

  expected_output <- tibble::tribble(
    ~USUBJID,      ~PARAMCD,  ~PARAM,              ~AVAL, ~AVALU,  ~VISIT,
    "01-701-1015", "HEIGHT",  "Height (cm)",       170,   "cm",   "BASELINE",
    "01-701-1015", "WEIGHT",  "Weight (kg)",        75,   "kg",   "BASELINE",
    "01-701-1015", "WEIGHT",  "Weight (kg)",        78,   "kg",   "MONTH 1",
    "01-701-1015", "WEIGHT",  "Weight (kg)",        80,   "kg",   "MONTH 2",
    "01-701-1028", "HEIGHT",  "Height (cm)",       185,   "cm",   "BASELINE",
    "01-701-1028", "WEIGHT",  "Weight (kg)",        90,   "kg",   "BASELINE",
    "01-701-1028", "WEIGHT",  "Weight (kg)",        88,   "kg",   "MONTH 1",
    "01-701-1028", "WEIGHT",  "Weight (kg)",        85,   "kg",   "MONTH 2",
    "01-701-1015", "BSA",     "Body Surface Area",  1.88, "m^2",  "BASELINE",
    "01-701-1028", "BSA",     "Body Surface Area",  2.15, "m^2",  "BASELINE",
  )

  expect_dfs_equal(
    derive_param_bsa(input, by_vars = vars(USUBJID, VISIT), method = "Mosteller"),
    expected_output,
    keys = c("USUBJID", "PARAMCD", "VISIT")
  )
})

test_that("new observations are derived correctly with DuBois & DuBois method", {

  expected_output <- tibble::tribble(
    ~USUBJID,      ~PARAMCD,  ~PARAM,              ~AVAL, ~AVALU,  ~VISIT,
    "01-701-1015", "HEIGHT",  "Height (cm)",       170,   "cm",   "BASELINE",
    "01-701-1015", "WEIGHT",  "Weight (kg)",        75,   "kg",   "BASELINE",
    "01-701-1015", "WEIGHT",  "Weight (kg)",        78,   "kg",   "MONTH 1",
    "01-701-1015", "WEIGHT",  "Weight (kg)",        80,   "kg",   "MONTH 2",
    "01-701-1028", "HEIGHT",  "Height (cm)",       185,   "cm",   "BASELINE",
    "01-701-1028", "WEIGHT",  "Weight (kg)",        90,   "kg",   "BASELINE",
    "01-701-1028", "WEIGHT",  "Weight (kg)",        88,   "kg",   "MONTH 1",
    "01-701-1028", "WEIGHT",  "Weight (kg)",        85,   "kg",   "MONTH 2",
    "01-701-1015", "BSA",     "Body Surface Area",  1.86, "m^2",  "BASELINE",
    "01-701-1028", "BSA",     "Body Surface Area",  2.14, "m^2",  "BASELINE",
  )

  expect_dfs_equal(
    derive_param_bsa(input, by_vars = vars(USUBJID, VISIT), method = "DuBois-DuBois"),
    expected_output,
    keys = c("USUBJID", "PARAMCD", "VISIT")
  )
})
