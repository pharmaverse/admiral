# filter_exist ----
## Test 1: filter_exist() works as expected ----
test_that("filter_exist Test 1: filter_exist() works as expected", {
  input_dataset <- tibble::tribble(
    ~USUBJID,     ~AGE, ~SEX,
    "01-701-1015",  63,  "F",
    "01-701-1034",  77,  "F",
    "01-701-1115",  84,  "M",
    "01-701-1444",  63,  "M"
  )

  input_dataset_add <- tibble::tribble(
    ~USUBJID,                         ~AEDECOD,     ~AESTDTC,
    "01-701-1015",                 "DIARRHOEA", "2014-01-09",
    "01-701-1034", "APPLICATION SITE PRURITUS", "2014-08-27",
    "01-701-1034",                   "FATIGUE", "2014-11-02",
    "01-701-1115",                   "FATIGUE", "2013-01-14"
  )

  expected_output <- input_dataset %>%
    filter(USUBJID %in% c("01-701-1034", "01-701-1115"))

  expect_equal(
    filter_exist(
      dataset = input_dataset,
      dataset_add = input_dataset_add,
      by_vars = exprs(USUBJID),
      filter_add = AEDECOD == "FATIGUE"
    ),
    expected_output
  )
})

# filter_not_exist ----
## Test 2: filter_not_exist() works as expected ----
test_that("filter_not_exist Test 2: filter_not_exist() works as expected", {
  input_dataset <- tibble::tribble(
    ~USUBJID,     ~AGE, ~SEX,
    "01-701-1015",  63,  "F",
    "01-701-1034",  77,  "F",
    "01-701-1115",  84,  "M",
    "01-701-1444",  63,  "M"
  )

  input_dataset_add <- tibble::tribble(
    ~USUBJID,                         ~AEDECOD,     ~AESTDTC,
    "01-701-1015",                 "DIARRHOEA", "2014-01-09",
    "01-701-1034", "APPLICATION SITE PRURITUS", "2014-08-27",
    "01-701-1034",                   "FATIGUE", "2014-11-02",
    "01-701-1115",                   "FATIGUE", "2013-01-14"
  )

  expected_output <- input_dataset %>%
    filter(USUBJID %in% c("01-701-1015", "01-701-1444"))

  expect_equal(
    filter_not_exist(
      dataset = input_dataset,
      dataset_add = input_dataset_add,
      by_vars = exprs(USUBJID),
      filter_add = AEDECOD == "FATIGUE"
    ),
    expected_output
  )
})
