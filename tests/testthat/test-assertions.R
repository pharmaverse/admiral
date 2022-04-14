library(admiraltest)

test_that("Test 1 : `assert_has_variables` an error is thrown if a required
           variable is missing", {
  data(dm)

  expect_error(
    assert_has_variables(dm, "TRT01P"),
    "Required variable `TRT01P` is missing."
  )
})

test_that("Test 2 : `assert_has_variables` no error is thrown if a required
           variable exists", {
  data(dm)

  expect_error(assert_has_variables(dm, "USUBJID"), NA)
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

  dm <- dm %>% group_by(ARMCD)

  expect_error(
    example_fun(dm)
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

test_that("Test 10 : `quo_not_missing` returns TRUE if no missing argument", {
  test_fun <- function(x) {
    x <- rlang::enquo(x)
    assertthat::assert_that(quo_not_missing(x))
  }
  expect_true(test_fun(my_variable))
})

test_that("Test 11 : `quo_not_missing` throws and Error if missing argument", {
  test_fun <- function(x) {
    x <- rlang::enquo(x)
    assertthat::assert_that(quo_not_missing(x))
  }
  expect_error(test_fun()) # missing argument -> throws error
})
