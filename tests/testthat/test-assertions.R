# assert_has_variables ----
## Test 1: error if a required variable is missing ----
test_that("assert_has_variables Test 1: error if a required variable is missing", {
  data <- tibble::tribble(
    ~USUBJID,
    "1"
  )

  expect_error(
    assert_has_variables(data, "TRT01P"),
    "Required variable `TRT01P` is missing."
  )

  expect_error(
    assert_has_variables(admiral.test::admiral_dm, c("TRT01P", "AVAL"))
  )
})

## Test 2: no error if a required variable exists ----
test_that("assert_has_variables Test 2: no error if a required variable exists", {
  data <- tibble::tribble(
    ~USUBJID,
    "1"
  )

  expect_error(assert_has_variables(data, "USUBJID"), NA)
})

# assert_filter_cond ----
## Test 3: `assert_filter_cond` works as expected ----
test_that("assert_filter_cond Test 3: `assert_filter_cond` works as expected", {
  fc <- expr(AGE == 64)
  expect_identical(
    assert_filter_cond(fc),
    fc
  )

  expect_identical(
    assert_filter_cond(arg = fc, optional = TRUE),
    fc
  )

  fc <- expr("string")
  expect_error(
    assert_filter_cond(arg = fc),
    "`fc` must be a filter condition but is `\"string\"`"
  )

  vals <- c("A", "B")
  fc <- expr(VAR %in% !!vals)
  expect_identical(
    assert_filter_cond(arg = fc),
    fc
  )
})

# assert_data_frame ----
## Test 4: error if not a dataframe ----
test_that("assert_data_frame Test 4: error if not a dataframe", {
  example_fun <- function(dataset) {
    assert_data_frame(dataset, required_vars = exprs(STUDYID, USUBJID))
  }
  expect_error(
    example_fun(c(1, 2, 3))
  )
})

## Test 5: error if dataframe is grouped ----
test_that("assert_data_frame Test 5: error if dataframe is grouped", {
  example_fun <- function(dataset) {
    assert_data_frame(dataset, required_vars = exprs(STUDYID, USUBJID))
  }

  data <- tibble::tribble(
    ~STUDYID, ~USUBJID, ~ARMCD,
    "xyz",    "1",      "PLACEBO",
    "xyz",    "2",      "ACTIVE"
  ) %>% group_by(ARMCD)

  expect_error(
    example_fun(data)
  )
})

# assert_character_scalar ----
## Test 6: error if not a character scaler string ----
test_that("assert_character_scalar Test 6: error if not a character scaler string", {
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

## Test 7: error if input is a vector ----
test_that("assert_character_scalar Test 7: error if input is a vector", {
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
## Test 8: no error if expected input ----
test_that("assert_vars Test 8: no error if expected input", {
  expect_invisible(assert_vars(exprs(USUBJID, PARAMCD)))
  expect_invisible(assert_vars(
    exprs(APERSDT = APxxSDT, APEREDT = APxxEDT),
    expect_names = TRUE
  ))
})

## Test 9: error if unexpected input ----
test_that("assert_vars Test 9: error if unexpected input", {
  expect_error(assert_vars(AVAL + 1))
  expect_error(assert_vars(rlang::quos(USUBJID, PARAMCD)))
  expect_error(assert_vars(c("USUBJID", "PARAMCD", "VISIT")))
  expect_error(assert_vars(exprs(USUBJID, AVAL + 2)))
  expect_error(assert_vars(exprs(APERSDT = APxxSDT, APxxEDT), expect_names = TRUE))
})

# assert_data_frame ----
## Test 10: no error if optional is TRUE and `arg` is NULL ----
test_that("assert_data_frame Test 10: no error if optional is TRUE and `arg` is NULL", {
  example_fun <- function(dataset) {
    assert_data_frame(dataset, required_vars = exprs(STUDYID, USUBJID), optional = TRUE)
  }

  expect_invisible(
    example_fun(NULL)
  )
})

## Test 11: error if required variables are missing ----
test_that("assert_data_frame Test 11: error if required variables are missing", {
  example_fun <- function(dataset) {
    assert_data_frame(dataset, required_vars = exprs(STUDYID, USUBJID))
  }

  admiral_dm <- admiral.test::admiral_dm %>% select(-c(STUDYID, USUBJID))

  expect_error(
    example_fun(admiral_dm)
  )
})

## Test 12: error if required variable is missing ----
test_that("assert_data_frame Test 12: error if required variable is missing", {
  example_fun <- function(dataset) {
    assert_data_frame(dataset, required_vars = exprs(STUDYID, USUBJID))
  }

  admiral_dm <- admiral.test::admiral_dm %>% select(-c(USUBJID))

  expect_error(
    example_fun(admiral_dm)
  )
})

# assert_character_scalar ----
## Test 13: no error if optional is TRUE and `arg` is NULL ----
test_that("assert_character_scalar Test 13: no error if optional is TRUE and `arg` is NULL", {
  example_fun <- function(character) {
    assert_character_scalar(character, optional = TRUE)
  }

  expect_invisible(
    example_fun(NULL)
  )
})

## Test 14: no error, case_sensitive = FALSE ----
test_that("assert_character_scalar Test 14: no error, case_sensitive = FALSE", {
  example_fun <- function(character) {
    assert_character_scalar(character, values = c("test"), case_sensitive = FALSE)
  }

  out <- expect_invisible(example_fun(character = "TEST"))
  expect_equal(out, "test")

  check_unit <- function(duration_unit) {
    assert_character_scalar(
      duration_unit,
      values = c("years", "months", "weeks", "days", "hours", "minutes", "seconds"),
      case_sensitive <- FALSE
    )
  }

  out <- expect_invisible(check_unit("months"))
  expect_equal(out, "months")

  out <- expect_invisible(check_unit("MONTHS"))
  expect_equal(out, "months")

  check_unit2 <- function(duration_unit) {
    assert_character_scalar(
      duration_unit,
      values = c("YEARS", "MONTHS", "WEEKS", "DAYS", "HOURS", "MINUTES", "SECONDS"),
      case_sensitive <- FALSE
    )
  }

  out <- expect_invisible(check_unit2("months"))
  expect_equal(out, "months")

  out <- expect_invisible(check_unit2("MONTHS"))
  expect_equal(out, "months")
})

## Test 15: error if `arg` not in values ----
test_that("assert_character_scalar Test 15: error if `arg` not in values", {
  example_fun <- function(character) {
    assert_character_scalar(character, values = c("test"))
  }

  expect_error(
    example_fun(character = "oak")
  )

  check_unit <- function(duration_unit) {
    assert_character_scalar(
      duration_unit,
      values = c("years", "months", "weeks", "days", "hours", "minutes", "seconds"),
      case_sensitive <- FALSE
    )
  }

  expect_error(
    check_unit("month"),
    paste0(
      "`duration_unit` must be one of 'years', 'months', 'weeks', 'days', ",
      "'hours', 'minutes' or 'seconds' but is 'month'"
    )
  )

  expect_error(
    check_unit("MONTH"),
    paste0(
      "`duration_unit` must be one of 'years', 'months', 'weeks', 'days', ",
      "'hours', 'minutes' or 'seconds' but is 'MONTH'"
    )
  )

  check_unit2 <- function(duration_unit) {
    assert_character_scalar(
      duration_unit,
      values = c("YEARS", "MONTHS", "WEEKS", "DAYS", "HOURS", "MINUTES", "SECONDS"),
      case_sensitive <- FALSE
    )
  }

  expect_error(
    check_unit2("month"),
    paste0(
      "`duration_unit` must be one of 'YEARS', 'MONTHS', 'WEEKS', 'DAYS', ",
      "'HOURS', 'MINUTES' or 'SECONDS' but is 'month'"
    )
  )

  expect_error(
    check_unit2("MONTH"),
    paste0(
      "`duration_unit` must be one of 'YEARS', 'MONTHS', 'WEEKS', 'DAYS', ",
      "'HOURS', 'MINUTES' or 'SECONDS' but is 'MONTH'"
    )
  )
})

## Test 16: error if `arg` not a character vector ----
test_that("assert_character_scalar Test 16: error if `arg` not a character vector", {
  arg <- c(1, 2, 3)

  expect_error(
    assert_character_vector(arg)
  )
})

## Test 17: error if `arg` is not in values ----
test_that("assert_character_scalar Test 17: error if `arg` is not in values", {
  example_fun <- function(character) {
    assert_character_vector(character, values = c("test", "oak"))
  }

  expect_error(
    example_fun(character = c("oak", "mint"))
  )
})

# assert_logical_scalar ----
## Test 18: no error if optional is TRUE and `arg` is NULL ----
test_that("assert_logical_scalar Test 18: no error if optional is TRUE and `arg` is NULL", {
  example_fun <- function(arg) {
    assert_logical_scalar(arg, optional = TRUE)
  }

  expect_invisible(
    example_fun(NULL)
  )
})

## Test 19: error if `arg` is not TRUE or FALSE ----
test_that("assert_logical_scalar Test 19: error if `arg` is not TRUE or FALSE", {
  example_fun <- function(arg) {
    assert_logical_scalar(arg)
  }
  arg <- c()
  expect_error(example_fun(NA))
  expect_error(example_fun(arg))
  expect_error(example_fun("test"))
})

# assert_symbol ----
## Test 20: no error if optional = TRUE and `arg` = NULL ----
test_that("assert_symbol Test 20: no error if optional = TRUE and `arg` = NULL", {
  f <- function(var) {
    v <- enexpr(var)
  }

  example_fun <- function(arg) {
    assert_symbol(arg, optional = TRUE)
  }

  expect_invisible(
    example_fun(
      f(NULL)
    )
  )
})

## Test 21: `assert_symbol` throws an error if `arg` is missing ----
test_that("assert_symbol Test 21: `assert_symbol` throws an error if `arg` is missing", {
  f <- function(var) {
    v <- enexpr(var)
  }

  example_fun <- function(arg) {
    assert_symbol(arg)
  }

  expect_error(
    example_fun(
      f()
    )
  )
})

## Test 22: `assert_symbol` throws an error if `arg` is not a symbol ----
test_that("assert_symbol Test 22: `assert_symbol` throws an error if `arg` is not a symbol", {
  f <- function(var) {
    v <- enexpr(var)
  }

  example_fun <- function(arg) {
    assert_symbol(arg)
  }

  expect_error(
    example_fun(
      f(NULL)
    )
  )
})

## Test 23: `assert_symbol` does not throw an error if `arg` is a symbol ----
test_that("assert_symbol Test 23: `assert_symbol` does not throw an error if `arg` is a symbol", {
  f <- function(var) {
    v <- enexpr(var)
  }

  admiral_dm <- admiral.test::admiral_dm

  example_fun <- function(arg) {
    assert_symbol(arg)
  }

  expect_invisible(
    example_fun(
      f(admiral_dm)
    )
  )
})

# assert_expr ----
## Test 24: `assert_expr` does not throw an error if `arg` is an expression ----
test_that("assert_expr Test 24: `assert_expr` does not throw an error if `arg` is an expression", {
  f <- function(var) {
    v <- enexpr(var)
  }

  example_fun <- function(arg) {
    assert_expr(arg)
  }

  expect_invisible(
    example_fun(
      f(admiral.test::admiral_dm)
    )
  )
})

## Test 25: no error if optional is TRUE and `arg` is NULL ----
test_that("assert_expr Test 25: no error if optional is TRUE and `arg` is NULL", {
  f <- function(var) {
    v <- enexpr(var)
  }

  example_fun <- function(arg) {
    assert_expr(arg, optional = TRUE)
  }

  expect_invisible(
    example_fun(
      f(NULL)
    )
  )
})

## Test 26: `assert_expr` throws an error if `arg` is missing ----
test_that("assert_expr Test 26: `assert_expr` throws an error if `arg` is missing", {
  f <- function(var) {
    v <- enexpr(var)
  }

  example_fun <- function(arg) {
    assert_expr(arg)
  }

  expect_error(
    example_fun(f())
  )
})

## Test 27: `assert_expr` throws an error if `arg` is not an expression ----
test_that("assert_expr Test 27: `assert_expr` throws an error if `arg` is not an expression", {
  f <- function(var) {
    v <- enexpr(var)
  }

  example_fun <- function(arg) {
    assert_expr(arg)
  }

  expect_error(
    example_fun(
      f(NULL)
    )
  )
})

# assert_vars ----
## Test 28: no error if `arg` is not a list of unquoted variable names ----
test_that("assert_vars Test 28: no error if `arg` is not a list of unquoted variable names", {
  example_fun <- function(arg) {
    assert_vars(arg)
  }

  expect_error(
    example_fun(c("USUBJID", "PARAMCD", "VISIT"))
  )
  expect_error(
    example_fun(TRUE)
  )
})

## Test 29: error if some elements of `arg` are not unquoted variable names ----
test_that("assert_vars Test 29: error if some elements of `arg` are not unquoted variable names", {
  example_fun <- function(arg) {
    assert_vars(arg)
  }

  expect_error(
    example_fun(exprs(USUBJID, PARAMCD, NULL))
  )
})

# assert_order_vars ----
## Test 30: error if `arg` is not a list variable names or `desc()` ----
test_that("assert_order_vars Test 30: error if `arg` is not a list variable names or `desc()`", {
  example_fun <- function(arg) {
    assert_order_vars(arg)
  }

  expect_error(
    example_fun(TRUE)
  )
  expect_error(
    example_fun(1)
  )
})

## Test 31: no error if `arg` is NULL and optional is TRUE ----
test_that("assert_order_vars Test 31: no error if `arg` is NULL and optional is TRUE", {
  example_fun <- function(arg) {
    assert_order_vars(arg, optional = TRUE)
  }

  expect_invisible(
    example_fun(NULL)
  )
})

# assert_integer_scalar ----
## Test 32: no error if `arg` is NULL and optional is TRUE ----
test_that("assert_integer_scalar Test 32: no error if `arg` is NULL and optional is TRUE", {
  example_fun <- function(arg) {
    assert_integer_scalar(arg, optional = TRUE)
  }

  expect_invisible(
    example_fun(NULL)
  )
})

## Test 33: error if chosen subset not in subsets ----
test_that("assert_integer_scalar Test 33: error if chosen subset not in subsets", {
  example_fun <- function(arg) {
    assert_integer_scalar(arg, subset = "infinity")
  }

  expect_error(
    example_fun(1)
  )
})

## Test 34: no error if `arg` is in selected subset ----
test_that("assert_integer_scalar Test 34: no error if `arg` is in selected subset", {
  example_fun <- function(arg) {
    assert_integer_scalar(arg, subset = "positive")
  }

  expect_invisible(
    example_fun(1)
  )
})

## Test 35: error if `arg` is not an integer scalar ----
test_that("assert_integer_scalar Test 35: error if `arg` is not an integer scalar", {
  example_fun <- function(arg) {
    assert_integer_scalar(arg)
  }

  arg <- c()

  expect_error(example_fun(TRUE))
  expect_error(example_fun(arg))
  expect_error(example_fun(Inf))
  expect_error(example_fun(1.5))
})

# assert_numeric_vector ----
## Test 36: no error if `arg` is NULL and optional is TRUE ----
test_that("assert_numeric_vector Test 36: no error if `arg` is NULL and optional is TRUE", {
  example_fun <- function(arg) {
    assert_numeric_vector(arg, optional = TRUE)
  }

  expect_invisible(
    example_fun(NULL)
  )
})

# assert_integer_scalar ----
## Test 37: error if `arg` is not an integer scalar ----
test_that("assert_integer_scalar Test 37: error if `arg` is not an integer scalar", {
  example_fun <- function(arg) {
    assert_numeric_vector(arg)
  }

  arg <- c()

  expect_error(example_fun(TRUE))
  expect_error(example_fun(arg))
  expect_error(example_fun("1.5"))
})

# assert_s3_class ----
## Test 38: error if `arg` is not an object of a specific class S3 ----
test_that("assert_s3_class Test 38: error if `arg` is not an object of a specific class S3", {
  example_fun <- function(arg) {
    assert_s3_class(arg, "factor")
  }

  expect_error(example_fun("test"))
})

## Test 39: no error if `arg` is NULL and optional is TRUE ----
test_that("assert_s3_class Test 39: no error if `arg` is NULL and optional is TRUE", {
  example_fun <- function(arg) {
    assert_s3_class(arg, class = "factor", optional = TRUE)
  }

  expect_invisible(
    example_fun(NULL)
  )
})

## Test 40: error if `arg` is NULL and optional is FALSE ----
test_that("assert_s3_class Test 40: error if `arg` is NULL and optional is FALSE", {
  example_fun <- function(arg) {
    assert_s3_class(arg, class = "factor", optional = FALSE)
  }

  expect_error(
    example_fun(NULL),
    "`arg` must be an object of class 'factor' but is `NULL`",
    fixed = TRUE
  )
})

## Test 41: no error if `arg` is an object of a specific class S3 ----
test_that("assert_s3_class Test 41: no error if `arg` is an object of a specific class S3", {
  example_fun <- function(arg) {
    assert_s3_class(arg, "factor")
  }

  expect_invisible(example_fun(as.factor("test")))
})

# assert_list_of ----
## Test 42: error if `arg` is not a list of specific class S3 objects ----
test_that("assert_list_of Test 42: error if `arg` is not a list of specific class S3 objects", {
  example_fun <- function(arg) {
    assert_list_of(arg, "factor")
  }

  expect_error(example_fun(list("test")))
})

## Test 43: no error if `arg` is NULL and optional is TRUE ----
test_that("assert_list_of Test 43: no error if `arg` is NULL and optional is TRUE", {
  example_fun <- function(arg) {
    assert_list_of(arg, class = "factor", optional = TRUE)
  }

  expect_invisible(
    example_fun(NULL)
  )
})

## Test 44: error if `arg` is NULL and optional is FALSE ----
test_that("assert_list_of Test 44: error if `arg` is NULL and optional is FALSE", {
  example_fun <- function(arg) {
    assert_list_of(arg, class = "factor", optional = FALSE)
  }

  expect_error(
    example_fun(NULL),
    "`arg` must be an object of class 'list' but is `NULL`",
    fixed = TRUE
  )
})

## Test 45: no error if `arg` is a list of specific class S3 objects ----
test_that("assert_list_of Test 45: no error if `arg` is a list of specific class S3 objects", {
  example_fun <- function(arg) {
    assert_list_of(arg, "factor")
  }

  expect_invisible(
    example_fun(
      list(as.factor("test"), as.factor(1))
    )
  )
})

## Test 44: error if `arg` is not a named list (no elements named) ----
test_that("assert_list_of Test 44: error if `arg` is not a named list (no elements named)", {
  expect_error(
    assert_list_of(mylist <- list(1, 2, 3), class = "numeric", named = TRUE),
    paste(
      "All elements of mylist must be named.",
      "No element is named.",
      sep = "\n"
    )
  )
})

## Test 45: error if `arg` is not a named list (some elements named) ----
test_that("assert_list_of Test 45: error if `arg` is not a named list (some elements named)", {
  expect_error(
    assert_list_of(mylist <- list(1, 2, 3, d = 4), class = "numeric", named = TRUE),
    paste(
      "All elements of mylist must be named.",
      "The following elements are not named: 1, 2 and 3",
      sep = "\n"
    )
  )
})

## Test 46: no error if `arg` is a named list ----
test_that("assert_list_of Test 46: no error if `arg` is a named list", {
  expect_invisible(
    assert_list_of(mylist <- list(a = 1, b = 2, c = 3), class = "numeric", named = TRUE)
  )
})

# assert_named_exprs ----
## Test 46: error if `arg` is not a named list of expressions ----
test_that("assert_named_exprs Test 46: error if `arg` is not a named list of expressions", {
  example_fun <- function(arg) {
    assert_named_exprs(arg)
  }

  arg <- list("test")
  names(arg) <- c("")

  expect_error(example_fun(5))
  expect_error(example_fun(admiral.test::admiral_dm))
  expect_error(example_fun(list(1, 2, TRUE)))
  expect_error(example_fun(arg))
})

## Test 47: no error if `arg` is NULL and optional is TRUE ----
test_that("assert_named_exprs Test 47: no error if `arg` is NULL and optional is TRUE", {
  example_fun <- function(arg) {
    assert_named_exprs(arg, optional = TRUE)
  }

  expect_invisible(
    example_fun(NULL)
  )
})

## Test 48: no error if `arg` is a named list of expressions ----
test_that("assert_named_exprs Test 48: no error if `arg` is a named list of expressions", {
  example_fun <- function(arg) {
    assert_named_exprs(arg)
  }

  expect_invisible(
    example_fun(
      rlang::exprs()
    )
  )
})

# assert_function ----
## Test 49: error if `arg` is not a function ----
test_that("assert_function Test 49: error if `arg` is not a function", {
  example_fun <- function(arg) {
    assert_function(arg)
  }

  expect_error(example_fun(5))
  expect_error(example_fun())
})

## Test 50: no error if `arg` is NULL and optional is TRUE ----
test_that("assert_function Test 50: no error if `arg` is NULL and optional is TRUE", {
  example_fun <- function(arg) {
    assert_function(arg, optional = TRUE)
  }

  expect_invisible(
    example_fun(NULL)
  )
})

## Test 51: no error if `arg` is a function with all parameters defined ----
test_that("assert_function Test 51: no error if `arg` is a function with all parameters defined", {
  example_fun <- function(arg) {
    assert_function(arg, params = c("x"))
  }

  expect_invisible(example_fun(mean))
})

## Test 52: error if  `params`  is missing with no default ----
test_that("assert_function Test 52: error if  `params`  is missing with no default", {
  example_fun <- function(arg) {
    assert_function(arg, params = c("x"))
  }

  expect_error(example_fun(sum))

  example_fun <- function(arg) {
    assert_function(arg, params = c("x", "y"))
  }

  expect_error(example_fun(sum))
})


# assert_function_param ----
## Test 53: no error if `arg` is a parameter of a function ----
test_that("assert_function_param Test 53: no error if `arg` is a parameter of a function", {
  hello <- function(name) {
    print(sprintf("Hello %s", name))
  }

  expect_invisible(assert_function_param("hello", "name"))
})

## Test 54: error if expected function parameters are missing ----
test_that("assert_function_param Test 54: error if expected function parameters are missing", {
  hello <- function(name) {
    print(sprintf("Hello %s", name))
  }

  expect_error(assert_function_param("hello", "surname"))
  expect_error(assert_function_param("hello", params = c("surname", "firstname")))
})

# assert_unit ----
## Test 55: no error if the parameter is provided in the expected unit ----
test_that("assert_unit Test 55: no error if the parameter is provided in the expected unit", {
  advs <- tibble::tribble(
    ~USUBJID, ~VSTESTCD, ~VSTRESN, ~VSSTRESU, ~PARAMCD, ~AVAL,
    "P01",    "WEIGHT",      80.1, "kg",      "WEIGHT",  80.1,
    "P02",    "WEIGHT",      85.7, "kg",      "WEIGHT",  85.7
  )

  expect_invisible(
    assert_unit(advs, param = "WEIGHT", required_unit = "kg", get_unit_expr = VSSTRESU)
  )
})

## Test 56: error if there are multiple units in the input dataset ----
test_that("assert_unit Test 56: error if there are multiple units in the input dataset", {
  advs <- tibble::tribble(
    ~USUBJID, ~VSTESTCD, ~VSTRESN, ~VSSTRESU, ~PARAMCD, ~AVAL,
    "P01",    "WEIGHT",      80.1, "kg",      "WEIGHT",  80.1,
    "P02",    "WEIGHT",      85.7, "lb",      "WEIGHT",  85.7
  )

  expect_error(
    assert_unit(advs, param = "WEIGHT", required_unit = "kg", get_unit_expr = VSSTRESU)
  )
})

## Test 57: error if unexpected unit in the input dataset ----
test_that("assert_unit Test 57: error if unexpected unit in the input dataset", {
  advs <- tibble::tribble(
    ~USUBJID, ~VSTESTCD, ~VSTRESN, ~VSSTRESU, ~PARAMCD, ~AVAL,
    "P01",    "WEIGHT",      80.1, "kg",      "WEIGHT",  80.1,
    "P02",    "WEIGHT",      85.7, "kg",      "WEIGHT",  85.7
  )

  expect_error(
    assert_unit(advs, param = "WEIGHT", required_unit = "lb", get_unit_expr = VSSTRESU)
  )
})

# assert_param_does_not_exist ----
## Test 58: error if parameter exists in the input dataset ----
test_that("assert_param_does_not_exist Test 58: error if parameter exists in the input dataset", {
  advs <- tibble::tribble(
    ~USUBJID, ~VSTESTCD, ~VSTRESN, ~VSSTRESU, ~PARAMCD, ~AVAL,
    "P01",    "WEIGHT",      80.1, "kg",      "WEIGHT",  80.1,
    "P02",    "WEIGHT",      85.7, "kg",      "WEIGHT",  85.7
  )

  expect_error(
    assert_param_does_not_exist(advs, param = "WEIGHT")
  )
})

## Test 59: no error if the parameter exists in the dataset ----
test_that("assert_param_does_not_exist Test 59: no error if the parameter exists in the dataset", {
  advs <- tibble::tribble(
    ~USUBJID, ~VSTESTCD, ~VSTRESN, ~VSSTRESU, ~PARAMCD, ~AVAL,
    "P01",    "WEIGHT",      80.1, "kg",      "WEIGHT",  80.1,
    "P02",    "WEIGHT",      85.7, "kg",      "WEIGHT",  85.7
  )

  expect_invisible(
    assert_param_does_not_exist(advs, param = "HR")
  )
})

# assert_varval_list ----
## Test 60: error if `arg` is not a list of var-value expressions ----
test_that("assert_varval_list Test 60: error if `arg` is not a list of var-value expressions", {
  example_fun <- function(arg) {
    assert_varval_list(arg, accept_var = FALSE)
  }

  expect_error(
    example_fun(c("USUBJID", "PARAMCD", "VISIT"))
  )
})

## Test 61: error if `arg` is not a list of var-value expressions ----
test_that("assert_varval_list Test 61: error if `arg` is not a list of var-value expressions", {
  example_fun <- function(arg) {
    assert_varval_list(arg, accept_var = TRUE)
  }

  expect_error(
    example_fun(exprs(USUBJID, PARAMCD, NULL))
  )
})

## Test 62: error if `required_elements` are missing from `arg` ----
test_that("assert_varval_list Test 62: error if `required_elements` are missing from `arg`", {
  example_fun <- function(arg) {
    assert_varval_list(arg, required_elements = "DTHDOM")
  }

  expect_error(
    example_fun(exprs(DTHSEQ = AESEQ))
  )
})

## Test 63: no error if `arg` is NULL and optional is TRUE ----
test_that("assert_varval_list Test 63: no error if `arg` is NULL and optional is TRUE", {
  example_fun <- function(arg) {
    assert_varval_list(arg, optional = TRUE)
  }

  expect_invisible(
    example_fun(NULL)
  )
})

## Test 64: error if `accept_expr` is TRUE and value is invalid ----
test_that("assert_varval_list Test 64: error if `accept_expr` is TRUE and value is invalid", {
  example_fun <- function(arg) {
    assert_varval_list(arg, accept_expr = TRUE)
  }

  expect_error(
    example_fun(exprs(DTHSEQ = TRUE))
  )
})

## Test 65: error if `accept_expr` is FALSE and value is invalid ----
test_that("assert_varval_list Test 65: error if `accept_expr` is FALSE and value is invalid", {
  example_fun <- function(arg) {
    assert_varval_list(arg, accept_expr = FALSE)
  }

  expect_error(
    example_fun(exprs(DTHSEQ = exprs()))
  )
})

## Test 66: no error if an argument is a variable-value list ----
test_that("assert_varval_list Test 66: no error if an argument is a variable-value list", {
  example_fun <- function(arg) {
    assert_varval_list(arg)
  }

  expect_invisible(
    example_fun(exprs(DTHDOM = "AE", DTHSEQ = AESEQ))
  )
})

# assert_list_element ----
## Test 67: no error if the elements fulfill a certain condition ----
test_that("assert_list_element Test 67: no error if the elements fulfill a certain condition", {
  expect_invisible(
    assert_list_element(
      list(
        list(var = expr(DTHDT), val = 1),
        list(var = expr(EOSDT), val = 0)
      ),
      element = "val",
      condition = val >= 0,
      message_text = ""
    )
  )
})

## Test 68: error if the elements do not fulfill the condition ----
test_that("assert_list_element Test 68: error if the elements do not fulfill the condition", {
  expect_error(
    assert_list_element(
      input <- list(
        list(var = expr(DTHDT), val = 1),
        list(var = expr(EOSDT), val = -1)
      ),
      element = "val",
      condition = val >= 0,
      message_text = "Invalid value for `val`:"
    ),
    "Invalid value for `val`:\ninput[[2]]$val = -1",
    fixed = TRUE
  )
})

# assert_one_to_one ----
## Test 69: error if there is a one to many mapping ----
test_that("assert_one_to_one Test 69: error if there is a one to many mapping", {
  expect_error(
    assert_one_to_one(admiral.test::admiral_dm, exprs(DOMAIN), exprs(USUBJID))
  )
})

## Test 70: error if there is a many to one mapping ----
test_that("assert_one_to_one Test 70: error if there is a many to one mapping", {
  expect_error(
    assert_one_to_one(admiral.test::admiral_dm, exprs(USUBJID), exprs(DOMAIN))
  )
})

# assert_date_var ----
## Test 71: error if variable is not a date or datetime variable ----
test_that("assert_date_var Test 71: error if variable is not a date or datetime variable", {
  example_fun <- function(dataset, var) {
    var <- assert_symbol(enexpr(var))
    assert_date_var(dataset = dataset, var = !!var)
  }

  my_data <- tibble::tribble(
    ~USUBJID, ~ADT,
    "1",      ymd("2020-12-06"),
    "2",      ymd("")
  )

  expect_error(
    example_fun(
      dataset = my_data,
      var = USUBJID
    )
  )
})

# assert_date_vector ----
## Test 72: returns error if input vector is not a date formatted ----
test_that("assert_date_vector Test 72: returns error if input vector is not a date formatted", {
  expect_error(assert_date_vector("2018-08-23"))
})

## Test 73: returns invisible if input is date formatted ----
test_that("assert_date_vector Test 73: returns invisible if input is date formatted", {
  expect_invisible(assert_date_vector(as.Date("2022-10-25")))
})

## Test 74: no error if `arg` is NULL and optional is TRUE ----
test_that("assert_date_vector Test 74: no error if `arg` is NULL and optional is TRUE", {
  example_fun <- function(arg) {
    assert_date_vector(arg, optional = TRUE)
  }

  expect_invisible(
    example_fun(NULL)
  )
})

## Test 75: error if `arg` is NULL and optional is FALSE ----
test_that("assert_date_vector Test 75: error if `arg` is NULL and optional is FALSE", {
  example_fun <- function(arg) {
    assert_date_vector(arg, optional = FALSE)
  }

  expect_error(
    example_fun(NULL),
    "`arg` must be a date or datetime variable but it's `NULL`",
    fixed = TRUE
  )
})


# assert_atomic_vector ----
## Test 76: error if input is not atomic vector ----
test_that("assert_atomic_vector Test 76: error if input is not atomic vector", {
  x <- list("a", "a", "b", "c", "d", "d", 1, 1, 4)
  expect_error(assert_atomic_vector(x))
})

# assert_same_type ----
## Test 77: no error if same type ----
test_that("assert_same_type Test 77: no error if same type", {
  true_value <- "Y"
  false_value <- "N"
  expect_invisible(assert_same_type(true_value, false_value))
})

## Test 78: error if different type ----
test_that("assert_same_type Test 78: error if different type", {
  true_value <- "Y"
  false_value <- "N"
  missing_value <- 0
  expect_error(
    assert_same_type(true_value, false_value, missing_value),
    regexp = paste(
      "All arguments must be of the same type.",
      "Argument: Type",
      "--------------",
      "true_value: character",
      "false_value: character",
      "missing_value: double",
      sep = "\n"
    ),
    fixed = TRUE
  )
})
