# warn_if_vars_exist ----
## Test 1: warning if a variable already exists in the input dataset ----
test_that("warn_if_vars_exist Test 1: warning if a variable already exists in the input dataset", {
  dm <- dplyr::tribble(
    ~USUBJID,      ~AGE,   ~AGEU,      ~ARM,
    "01-701-1015",   25, "YEARS", "Placebo",
    "01-701-1016",   27, "YEARS", "Placebo",
  )

  expect_snapshot(
    warn_if_vars_exist(dm, "AGE")
  )
  expect_snapshot(
    warn_if_vars_exist(dm, c("AGE", "AGEU", "ARM"))
  )
  expect_snapshot(
    warn_if_vars_exist(dm, c("AAGE", "AGEU", "ARM"))
  )
  expect_warning(
    warn_if_vars_exist(dm, "AAGE"),
    NA
  )
})

# warn_if_invalid_dtc ----
## Test 2: Warning if vector contains unknown datetime format ----
test_that("warn_if_invalid_dtc Test 2: Warning if vector contains unknown datetime format", {
  expect_snapshot(
    warn_if_invalid_dtc(dtc = "20210406T12:30:30")
  )
})

# warn_if_inclomplete_dtc ----
## Test 3: Warning if vector contains an incomplete dtc ----
test_that("warn_if_inclomplete_dtc Test 3: Warning if vector contains an incomplete date", {
  expect_warning(
    warn_if_incomplete_dtc("2021-04", n = 10)
  )
})

# warn_if_inclomplete_dtc ----
## Test 4: Warning if vector contains an incomplete dtc ----
test_that("warn_if_inclomplete_dtc Test 4: Warning if vector contains an incomplete datetime", {
  expect_warning(
    warn_if_incomplete_dtc("2021-04-06T12:30", n = 19)
  )
})

# warn_if_inconsistent_list ----
## Test 5: Warning if two lists are inconsistent ----
test_that("warn_if_inconsistent_list Test 5: Warning if two lists are inconsistent", {
  expect_warning(
    warn_if_inconsistent_list(
      base = exprs(DTHDOM = "DM", DTHSEQ = DMSEQ, DTHVAR = "text"),
      compare = exprs(DTHDOM = "DM", DTHSEQ = DMSEQ),
      list_name = "Test"
    )
  )
})

# suppress_warning ----
## Test 6: Suppress certain warnings issued by an expression ----
test_that("suppress_warning Test 6: Suppress certain warnings issued by an expression", {
  # Verify if warning is issued
  expect_warning(
    suppress_warning(as.numeric("fun"), "x")
  )

  # Actual result
  actual_result <- as.numeric(NA)

  # Call the suppress_warning() to suppress warning messages containing "NAs introduced by coercion"
  expected_result <- suppress_warning(as.numeric("fun"), "coercion")

  # Expect that the warning message has been suppressed and that the result is NA
  expect_equal(expected_result, actual_result)
})
