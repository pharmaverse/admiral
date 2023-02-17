# warn_if_vars_exist ----
## Test 1: warning if a variable already exists in the input dataset ----
test_that("warn_if_vars_exist Test 1: warning if a variable already exists in the input dataset", {
  dm <- admiral.test::admiral_dm

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

# warn_if_invalud_dtc ----
## Test 2: Warning if vector contains unknown datetime format ----
test_that("warn_if_invalud_dtc Test 2: Warning if vector contains unknown datetime format", {
  expect_warning(
    warn_if_invalid_dtc(dtc = "20210406T12:30:30")
  )
})

# warn_if_inclomplete_dtc ----
## Test 3: Warning if vector contains an incomplete dtc ----
test_that("warn_if_inclomplete_dtc Test 3: Warning if vector contains an incomplete dtc", {
  expect_warning(
    warn_if_incomplete_dtc("2021-04-06", n = 19)
  )
})

# warn_if_inconsistent_list ----
## Test 4: Warning if two lists are inconsistent ----
test_that("warn_if_inconsistent_list Test 4: Warning if two lists are inconsistent", {
  expect_warning(
    warn_if_inconsistent_list(
      base = exprs(DTHDOM = "DM", DTHSEQ = DMSEQ, DTHVAR = "text"),
      compare = exprs(DTHDOM = "DM", DTHSEQ = DMSEQ),
      list_name = "Test"
    )
  )
})
