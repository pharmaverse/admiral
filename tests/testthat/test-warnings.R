library(admiral.test)

test_that("A warning is issued when a variable to be derived already exists in the input dataset", {
  data(dm)

  expect_warning(
    warn_if_vars_exist(dm, "AGE"),
    "Variable `AGE` already exists in the dataset"
  )
  expect_warning(
    warn_if_vars_exist(dm, c("AGE", "AGEU", "ARM")),
    "Variables `AGE`, `AGEU` and `ARM` already exist in the dataset"
  )
  expect_warning(
    warn_if_vars_exist(dm, c("AAGE", "AGEU", "ARM")),
    "Variables `AGEU` and `ARM` already exist in the dataset"
  )
  expect_warning(
    warn_if_vars_exist(dm, "AAGE"),
    NA
  )
})

test_that("A warning is issued when a vector contain unknown datetime format", {
  expect_warning(
    warn_if_invalid_dtc(dtc = "2021-04-06T-:30:30")
  )
})

test_that("A warning is issued when a vector contain an incomplete dtc", {
  expect_warning(
    warn_if_incomplete_dtc("2021-04-06", n = 19)
  )
})

test_that("A warning is issued when two lists are inconsistent", {
  expect_warning(
    warn_if_inconsistent_list(
      base = vars(DTHDOM = "DM", DTHSEQ = DMSEQ, DTHVAR = "text"),
      compare = vars(DTHDOM = "DM", DTHSEQ = DMSEQ),
      list_name = "Test"
    )
  )
})
