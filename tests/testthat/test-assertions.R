library(admiral.test)

test_that("Test 1 : `assert_has_variables` an error is thrown if a required
           variable is missing", {
  data(admiral_dm)

  expect_error(
    assert_has_variables(admiral_dm, "TRT01P"),
    "Required variable `TRT01P` is missing."
  )
})

test_that("Test 2 : `assert_has_variables` no error is thrown if a required
           variable exists", {
  data(admiral_dm)

  expect_error(assert_has_variables(admiral_dm, "USUBJID"), NA)
})

test_that("Test 3 : `assert_filter_cond` works as expected", {
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

test_that("Test 4 : is_valid_sec_min works as expected", {
  expect_true(is_valid_sec_min(59))
})

test_that("Test 5 : is_valid_hour works as expected", {
  expect_true(is_valid_hour(23))
})

test_that("Test 6 : `assert_data_frame` throws an error if not a dataframe", {
  example_fun <- function(dataset) {
    assert_data_frame(dataset, required_vars = vars(STUDYID, USUBJID))
  }
  expect_error(
    example_fun(c(1, 2, 3))
  )
})

test_that("Test 7 : `assert_data_frame` throws an error if dataframe is grouped", {
  example_fun <- function(dataset) {
    assert_data_frame(dataset, required_vars = vars(STUDYID, USUBJID))
  }

  admiral_dm <- admiral_dm %>% group_by(ARMCD)

  expect_error(
    example_fun(admiral_dm)
  )
})

test_that("Test 8 : `assert_character_scalar` throws an error if not a character scaler string", {
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

test_that("Test 9 : `assert_character_scalar` throws an error if input is a vector", {
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

test_that("Test 10 : `assert_order_vars` returns invisible if used correctly", {
  expect_invisible(assert_order_vars(vars(USUBJID, PARAMCD, desc(AVISITN))))
})

test_that("Test 11 : `assert_order_vars` returns errors if used incorrectly", {
  expect_error(assert_order_vars(rlang::exprs(USUBJID, PARAMCD)))
  expect_error(assert_order_vars(c("USUBJID", "PARAMCD", "VISIT")))
  expect_error(assert_order_vars(vars(USUBJID, toupper(PARAMCD), -AVAL)))
})

test_that("Test 12 : `assert_data_frame` doesn't throw an error if optional is TRUE and arg is NULL", {
  example_fun <- function(dataset) {
    assert_data_frame(dataset, required_vars = vars(STUDYID, USUBJID), optional=TRUE)
  }

  expect_invisible(
    example_fun(NULL)
  )
})

test_that("Test 13 : `assert_data_frame` throws an error if required variables are missing", {
  example_fun <- function(dataset) {
    assert_data_frame(dataset, required_vars = vars(STUDYID, USUBJID))
  }

  admiral_dm <- admiral_dm %>% select(-c(STUDYID,USUBJID))

  expect_error(
    example_fun(admiral_dm)
  )
})

test_that("Test 14 : `assert_character_scalar` doesn't throw an error if optional is TRUE and arg is NULL", {
  example_fun <- function(character) {
    assert_character_scalar(character, optional=TRUE)
  }

  expect_invisible(
    example_fun(NULL)
  )
})

test_that("Test 15 : `assert_character_scalar` doesn't throw an error if case_sensitive is FALSE", {
  example_fun <- function(character) {
    assert_character_scalar(character, values=c("test"), case_sensitive=FALSE)
  }

  expect_invisible(
    example_fun(character="TEST")
  )
})

test_that("Test 16 : `assert_character_scalar` throws an error if values are not NULL and arg not in values", {
  example_fun <- function(character) {
    assert_character_scalar(character, values=c("test"))
  }

  expect_error(
    example_fun(character="oak")
  )
})

test_that("Test 17 : `assert_character_vector` throws an error if arg is not a character vector", {

  expect_error(
    assert_character_vector(c(1,2,3))
  )
})

test_that("Test 18 : `assert_character_vector` throws an error if values are not NULL and arg not in values", {
  example_fun <- function(character) {
    assert_character_vector(character, values=c("test","oak"))
  }

  expect_error(
    example_fun(character=c("oak","mint"))
  )
})

test_that("Test 19 : `assert_logical_scalar` doesn't throw an error if optional is TRUE and arg is NULL", {
  example_fun <- function(arg) {
    assert_logical_scalar(arg, optional=TRUE)
  }

  expect_invisible(
    example_fun(NULL)
  )
})
