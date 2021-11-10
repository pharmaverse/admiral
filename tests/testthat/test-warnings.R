library(cdiscpilot)

test_that("a warning is issued when a variable to be derived already exists in the input dataset", {
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
