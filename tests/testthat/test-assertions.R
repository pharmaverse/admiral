library(admiral.test)

# assert_has_variables ----
## Test 1: error if a required variable is missing ----
test_that("assert_has_variables Test 1: error if a required variable is missing", {
  data(admiral_dm)

  expect_error(
    assert_has_variables(admiral_dm, "TRT01P"),
    "Required variable `TRT01P` is missing."
  )
})

## Test 2: no error if a required variable exists ----
test_that("assert_has_variables Test 2: no error if a required variable exists", {
  data(admiral_dm)

  expect_error(assert_has_variables(admiral_dm, "USUBJID"), NA)
})

# assert_filter_cond ----
## Test 3: `assert_filter_cond` works as expected ----
test_that("assert_filter_cond Test 3: `assert_filter_cond` works as expected", {
  fc <- quo(AGE == 64)
  expect_identical(
    assert_filter_cond(fc),
    fc
  )

  fc <- quo()
  expect_error(
    assert_filter_cond(arg = fc),
    "Argument `fc` is missing, with no default"
  )

  expect_identical(
    assert_filter_cond(arg = fc, optional = TRUE),
    fc
  )

  fc <- quo("string")
  expect_error(
    assert_filter_cond(arg = fc),
    "`fc` must be a filter condition but is `\"string\"`"
  )
})

# is_valid_sec_min ----
## Test 4: is_valid_sec_min works as expected ----
test_that("is_valid_sec_min Test 4: is_valid_sec_min works as expected", {
  expect_true(is_valid_sec_min(59))
})

# is_valid_hour ----
## Test 5: is_valid_hour works as expected ----
test_that("is_valid_hour Test 5: is_valid_hour works as expected", {
  expect_true(is_valid_hour(23))
})

# assert_data_fram ----
## Test 6: error if not a dataframe ----
test_that("assert_data_fram Test 6: error if not a dataframe", {
  example_fun <- function(dataset) {
    assert_data_frame(dataset, required_vars = vars(STUDYID, USUBJID))
  }
  expect_error(
    example_fun(c(1, 2, 3))
  )
})

## Test 7: error if dataframe is grouped ----
test_that("assert_data_fram Test 7: error if dataframe is grouped", {
  example_fun <- function(dataset) {
    assert_data_frame(dataset, required_vars = vars(STUDYID, USUBJID))
  }

  admiral_dm <- admiral_dm %>% group_by(ARMCD)

  expect_error(
    example_fun(admiral_dm)
  )
})

# assert_character_scalar ----
## Test 8: error if not a character scaler string ----
test_that("assert_character_scalar Test 8: error if not a character scaler string", {
  example_fun2 <- function(msg_type) {
    msg_type <- assert_character_scalar(msg_type,
      values = c("warning", "error"), case_sensitive = FALSE
    )

    if (msg_type == "warning") {
      print("A warning was requested.")
    }
  }
  expect_error(example_fun2(2))
})

## Test 9: error if input is a vector ----
test_that("assert_character_scalar Test 9: error if input is a vector", {
  example_fun2 <- function(msg_type) {
    msg_type <- assert_character_scalar(msg_type,
      values = c("warning", "error"), case_sensitive = FALSE
    )

    if (msg_type == "warning") {
      print("A warning was requested.")
    }
  }
  expect_error(example_fun2(c("admiral", "admiralonco")))
})

# assert_vars ----
test_that("no error if expected input", {
  expect_invisible(assert_vars(vars(USUBJID, PARAMCD)))
  expect_invisible(assert_vars(
    vars(APERSDT = APxxSDT, APEREDT = APxxEDT),
    expect_names = TRUE))
})

test_that("error if unexpected input", {
  expect_error(assert_vars(AVAL + 1))
  expect_error(assert_vars(rlang::exprs(USUBJID, PARAMCD)))
  expect_error(assert_vars(c("USUBJID", "PARAMCD", "VISIT")))
  expect_error(assert_vars(vars(USUBJID, AVAL + 2)))
  expect_error(assert_vars(vars(APERSDT = APxxSDT, APxxEDT), expect_names = TRUE))
}
)

# assert_order_vars ----
## Test 10: returns invisible if used correctly ----
test_that("assert_order_vars Test 10: returns invisible if used correctly", {
  expect_invisible(assert_order_vars(vars(USUBJID, PARAMCD, desc(AVISITN))))
})

## Test 11: returns errors if used incorrectly ----
test_that("assert_order_vars Test 11: returns errors if used incorrectly", {
  expect_error(assert_order_vars(rlang::exprs(USUBJID, PARAMCD)))
  expect_error(assert_order_vars(c("USUBJID", "PARAMCD", "VISIT")))
  expect_error(assert_order_vars(vars(USUBJID, toupper(PARAMCD), -AVAL)))
})
