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
## Test 10: no error if expected input ----
test_that("assert_vars Test 10: no error if expected input", {
  expect_invisible(assert_vars(vars(USUBJID, PARAMCD)))
  expect_invisible(assert_vars(
    vars(APERSDT = APxxSDT, APEREDT = APxxEDT),
    expect_names = TRUE
  ))
})

## Test 11: error if unexpected input ----
test_that("assert_vars Test 11: error if unexpected input", {
  expect_error(assert_vars(AVAL + 1))
  expect_error(assert_vars(rlang::exprs(USUBJID, PARAMCD)))
  expect_error(assert_vars(c("USUBJID", "PARAMCD", "VISIT")))
  expect_error(assert_vars(vars(USUBJID, AVAL + 2)))
  expect_error(assert_vars(vars(APERSDT = APxxSDT, APxxEDT), expect_names = TRUE))
})

# assert_data_frame ----
test_that("Test 12 : `assert_data_frame` does not throw an error
          if optional is TRUE and `arg` is NULL", {
  example_fun <- function(dataset) {
    assert_data_frame(dataset, required_vars = vars(STUDYID, USUBJID), optional = TRUE)
  }

  expect_invisible(
    example_fun(NULL)
  )
})

test_that("Test 13 : `assert_data_frame` throws an error if required variables are missing", {
  example_fun <- function(dataset) {
    assert_data_frame(dataset, required_vars = vars(STUDYID, USUBJID))
  }

  admiral_dm <- admiral.test::admiral_dm %>% select(-c(STUDYID, USUBJID))

  expect_error(
    example_fun(admiral_dm)
  )
})

test_that("Test 14 : `assert_data_frame` throws an error if required variable is missing", {
  example_fun <- function(dataset) {
    assert_data_frame(dataset, required_vars = vars(STUDYID, USUBJID))
  }

  admiral_dm <- admiral.test::admiral_dm %>% select(-c(USUBJID))

  expect_error(
    example_fun(admiral_dm)
  )
})

# assert_character_scalar ----
test_that("Test 15 : `assert_character_scalar` doesn't throw an error if
          optional is TRUE and `arg` is NULL", {
  example_fun <- function(character) {
    assert_character_scalar(character, optional = TRUE)
  }

  expect_invisible(
    example_fun(NULL)
  )
})

test_that("Test 16 : `assert_character_scalar` does not throw an error if case_sensitive is FALSE", {
  example_fun <- function(character) {
    assert_character_scalar(character, values = c("test"), case_sensitive = FALSE)
  }

  expect_invisible(
    example_fun(character = "TEST")
  )
})

test_that("Test 17 : `assert_character_scalar` throws an error if
          values are not NULL and `arg` not in values", {
  example_fun <- function(character) {
    assert_character_scalar(character, values = c("test"))
  }

  expect_error(
    example_fun(character = "oak")
  )
})

test_that("Test 18 : `assert_character_vector` throws an error if `arg` not a character vector", {
  arg <- c(1, 2, 3)

  expect_error(
    assert_character_vector(arg)
  )
})

test_that("Test 19 : `assert_character_vector` throws an error
          if values are not NULL and `arg` not in values", {
  example_fun <- function(character) {
    assert_character_vector(character, values = c("test", "oak"))
  }

  expect_error(
    example_fun(character = c("oak", "mint"))
  )
})

# assert_logical_scalar ----
test_that("Test 20 : `assert_logical_scalar` doesn't throw an error
          if optional is TRUE and `arg` is NULL", {
  example_fun <- function(arg) {
    assert_logical_scalar(arg, optional = TRUE)
  }

  expect_invisible(
    example_fun(NULL)
  )
})

test_that("Test 21 : `assert_logical_scalar` throws an error if `arg` is not TRUE or FALSE", {
  example_fun <- function(arg) {
    assert_logical_scalar(arg)
  }
  arg <- c()
  expect_error(example_fun(NA))
  expect_error(example_fun(arg))
  expect_error(example_fun("test"))
})

# assert_symbol ----
test_that("Test 22 : `assert_symbol` doesn't throw an error if optional = TRUE and `arg` = NULL", {
  f <- function(var) {
    v <- enquo(var)
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

test_that("Test 23 : `assert_symbol` throws an error if `arg` is missing", {
  f <- function(var) {
    v <- enquo(var)
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

test_that("Test 24 : `assert_symbol` throws an error if `arg` is not a symbol", {
  f <- function(var) {
    v <- enquo(var)
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

test_that("Test 25 : `assert_symbol` doesn't throw an error if `arg` is a symbol", {
  f <- function(var) {
    v <- enquo(var)
  }

  example_fun <- function(arg) {
    assert_symbol(arg)
  }

  library(admiral.test)

  expect_invisible(
    example_fun(
      f(admiral_dm)
    )
  )
})

# assert_expr ----
test_that("Test 26 : `assert_expr` doesn't throw an error if `arg` is an expression", {
  f <- function(var) {
    v <- enquo(var)
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

test_that("Test 27 : `assert_expr` doesn't throw an error if
          optional is TRUE and `arg` is NULL", {
  f <- function(var) {
    v <- enquo(var)
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

test_that("Test 28 : `assert_expr` throws an error if `arg` is missing", {
  f <- function(var) {
    v <- enquo(var)
  }

  example_fun <- function(arg) {
    assert_expr(arg)
  }

  expect_error(
    example_fun(f())
  )
})

test_that("Test 29 : `assert_expr` throws an error if `arg` is not an expression", {
  f <- function(var) {
    v <- enquo(var)
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
test_that("Test 30 : `assert_vars` throws an error
          if `arg` is not a list of unquoted variable names", {
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

test_that("Test 31 : `assert_vars` throws an error
          if some elements of `arg` are not unquoted variable names", {
  example_fun <- function(arg) {
    assert_vars(arg)
  }

  expect_error(
    example_fun(vars(USUBJID, PARAMCD, NULL))
  )
})

# assert_order_vars ----
test_that("Test 32 : `assert_order_vars` throws an error
          if `arg` is not a list of unquoted variable names or `desc()` calls", {
  example_fun <- function(arg) {
    assert_order_vars(arg)
  }

  expect_error(
    example_fun(TRUE)
  )
})

test_that("Test 33 : `assert_order_vars` doesn't throw an error
          if `arg` is NULL and optional is TRUE", {
  example_fun <- function(arg) {
    assert_order_vars(arg, optional = TRUE)
  }

  expect_invisible(
    example_fun(NULL)
  )
})

# assert_integer_scalar ----
test_that("Test 34 : `assert_integer_scalar` doesn't throw an error
          if `arg` is NULL and optional is TRUE", {
  example_fun <- function(arg) {
    assert_integer_scalar(arg, optional = TRUE)
  }

  expect_invisible(
    example_fun(NULL)
  )
})

test_that("Test 35 : `assert_integer_scalar` throws an error if chosen subset not in subsets", {
  example_fun <- function(arg) {
    assert_integer_scalar(arg, subset = "infinity")
  }

  expect_error(
    example_fun(1)
  )
})

test_that("Test 36 : `assert_integer_scalar` doesn't throw an error
          if `arg` is an integer scalar in selected subset", {
  example_fun <- function(arg) {
    assert_integer_scalar(arg, subset = "positive")
  }

  expect_invisible(
    example_fun(1)
  )
})

test_that("Test 37 : `assert_integer_scalar` throws an error
          if `arg` is not an integer scalar", {
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
test_that("Test 38 : `assert_numeric_vector` doesn't throw an error
          if `arg` is NULL and optional is TRUE", {
  example_fun <- function(arg) {
    assert_numeric_vector(arg, optional = TRUE)
  }

  expect_invisible(
    example_fun(NULL)
  )
})

# assert_integer_scalar ----
test_that("Test 39 : `assert_integer_scalar` throws an error if `arg` is not an integer scalar", {
  example_fun <- function(arg) {
    assert_numeric_vector(arg)
  }

  arg <- c()

  expect_error(example_fun(TRUE))
  expect_error(example_fun(arg))
  expect_error(example_fun("1.5"))
})

# assert_s3_class ----
test_that("Test 40 : `assert_s3_class` throws an error
          if `arg` is not an object of a specific class S3", {
  example_fun <- function(arg) {
    assert_s3_class(arg, "factor")
  }

  expect_error(example_fun("test"))
})

test_that("Test 41 : `assert_s3_class` doesn't throw an error
          if `arg` is NULL and optional is TRUE", {
  example_fun <- function(arg) {
    assert_s3_class(arg, class = "factor", optional = TRUE)
  }

  expect_invisible(
    example_fun(NULL)
  )
})

test_that("Test 42 : `assert_s3_class` doesn't throw an error
          if `arg` is an object of a specific class S3", {
  example_fun <- function(arg) {
    assert_s3_class(arg, "factor")
  }

  expect_invisible(example_fun(as.factor("test")))
})

# assert_list_of ----
test_that("Test 43 : `assert_list_of` throws an error
          if `arg` is not a list of objects of a specific class S3", {
  example_fun <- function(arg) {
    assert_list_of(arg, "factor")
  }

  expect_error(example_fun(list("test")))
})

test_that("Test 44 : `assert_list_of` doesn't throw an error
          if `arg` is NULL and optional is TRUE", {
  example_fun <- function(arg) {
    assert_list_of(arg, class = "factor", optional = TRUE)
  }

  expect_invisible(
    example_fun(NULL)
  )
})

test_that("Test 45 : `assert_list_of` doesn't throw an error
          if `arg` is a list of objects of a specific class S3", {
  example_fun <- function(arg) {
    assert_list_of(arg, "factor")
  }

  expect_invisible(
    example_fun(
      list(as.factor("test"), as.factor(1))
    )
  )
})

# assert_named_exprs ----
test_that("Test 46 : `assert_named_exprs` throws an error
          if `arg` is not a named list of expressions", {
  example_fun <- function(arg) {
    assert_named_exprs(arg)
  }

  arg <- list("test")
  names(arg) <- c("")

  expect_error(example_fun(5))
  expect_error(example_fun(admiral_df))
  expect_error(example_fun(list(1, 2, TRUE)))
  expect_error(example_fun(arg))
})

test_that("Test 47 : `assert_named_exprs` doesn't throw an error
          if `arg` is NULL and optional is TRUE", {
  example_fun <- function(arg) {
    assert_named_exprs(arg, optional = TRUE)
  }

  expect_invisible(
    example_fun(NULL)
  )
})

test_that("Test 48 : `assert_named_exprs` doesn't throw an error
          if `arg` is a named list of expressions", {
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
test_that("Test 49 : `assert_function` throws an error
          if `arg` is not a function", {
  example_fun <- function(arg) {
    assert_function(arg)
  }

  expect_error(example_fun(5))
  expect_error(example_fun())
})

test_that("Test 50 : `assert_function` doesn't throw an error
          if `arg` is NULL and optional is TRUE", {
  example_fun <- function(arg) {
    assert_function(arg, optional = TRUE)
  }

  expect_invisible(
    example_fun(NULL)
  )
})

test_that("Test 51 : `assert_function` doesn't throw an error
          if `arg` is a function with all parameters defined", {
  example_fun <- function(arg) {
    assert_function(arg, params = c("x"))
  }

  expect_invisible(example_fun(mean))
})

test_that("Test 52 : `assert_function` throws an error
          if  `params`  is missing with no default", {
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
test_that("Test 53 : `assert_function_param` doesn't throw an error
          if `arg` is a parameter of a function", {
  hello <- function(name) {
    print(sprintf("Hello %s", name))
  }

  expect_invisible(assert_function_param("hello", "name"))
})

test_that("Test 54 : `assert_function_param` throws an error
          if any elements of `params` is not a parameter of the function given by `arg`", {
  hello <- function(name) {
    print(sprintf("Hello %s", name))
  }

  expect_error(assert_function_param("hello", "surname"))
  expect_error(assert_function_param("hello", params = c("surname", "firstname")))
})

# assert_unit ----
test_that("Test 55 : `assert_unit` doesn't throw an error
          if the parameter is provided in the expected unit", {
  advs <- tibble::tribble(
    ~USUBJID, ~VSTESTCD, ~VSTRESN, ~VSSTRESU, ~PARAMCD, ~AVAL,
    "P01",    "WEIGHT",      80.1, "kg",      "WEIGHT",  80.1,
    "P02",    "WEIGHT",      85.7, "kg",      "WEIGHT",  85.7
  )

  expect_invisible(
    assert_unit(advs, param = "WEIGHT", required_unit = "kg", get_unit_expr = VSSTRESU)
  )
})

test_that("Test 56 : `assert_unit` throws an error
          if there are multiple units in the input dataset", {
  advs <- tibble::tribble(
    ~USUBJID, ~VSTESTCD, ~VSTRESN, ~VSSTRESU, ~PARAMCD, ~AVAL,
    "P01",    "WEIGHT",      80.1, "kg",      "WEIGHT",  80.1,
    "P02",    "WEIGHT",      85.7, "lb",      "WEIGHT",  85.7
  )

  expect_error(
    assert_unit(advs, param = "WEIGHT", required_unit = "kg", get_unit_expr = VSSTRESU)
  )
})

test_that("Test 57 : `assert_unit` throws an error if the unit variable differs from
          the unit for any observation of the parameter in the input dataset", {
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
test_that("Test 58 : `assert_param_does_not_exist` throws an error
          if the parameter exists in the input dataset", {
  advs <- tibble::tribble(
    ~USUBJID, ~VSTESTCD, ~VSTRESN, ~VSSTRESU, ~PARAMCD, ~AVAL,
    "P01",    "WEIGHT",      80.1, "kg",      "WEIGHT",  80.1,
    "P02",    "WEIGHT",      85.7, "kg",      "WEIGHT",  85.7
  )

  expect_error(
    assert_param_does_not_exist(advs, param = "WEIGHT")
  )
})

test_that("Test 59 : `assert_param_does_not_exist` doesn't throw an error
          if the parameter exists in the input dataset", {
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
test_that("Test 60 : `assert_varval_list` throws an error
          if `arg` is not a list of variable-value expressions", {
  example_fun <- function(arg) {
    assert_varval_list(arg, accept_var = FALSE)
  }

  expect_error(
    example_fun(c("USUBJID", "PARAMCD", "VISIT"))
  )
})

test_that("Test 61 : `assert_varval_list` throws an error
          if `arg` is not a list of variable-value expressions", {
  example_fun <- function(arg) {
    assert_varval_list(arg, accept_var = TRUE)
  }

  expect_error(
    example_fun(vars(USUBJID, PARAMCD, NULL))
  )
})

test_that("Test 62 : `assert_varval_list` throws an error
          if `required_elements` are missing from `arg`", {
  example_fun <- function(arg) {
    assert_varval_list(arg, required_elements = "DTHDOM")
  }

  expect_error(
    example_fun(vars(DTHSEQ = AESEQ))
  )
})

test_that("Test 63 : `assert_varval_list` doesn't throw an error
          if `arg` is NULL and optional is TRUE", {
  example_fun <- function(arg) {
    assert_varval_list(arg, optional = TRUE)
  }

  expect_invisible(
    example_fun(NULL)
  )
})

test_that("Test 64 : `assert_varval_list` throws an error if `accept_expr` is TRUE
          and the expressions in `arg` are not: a symbol, a string, a numeric , an
          `NA` or an expression", {
  example_fun <- function(arg) {
    assert_varval_list(arg, accept_expr = TRUE)
  }

  expect_error(
    example_fun(vars(DTHSEQ = TRUE))
  )
})

test_that("Test 65 : `assert_varval_list` throws an error if `accept_expr` is FALSE and
          the expressions in `arg` are not: a symbol, a string, a numeric or `NA`", {
  example_fun <- function(arg) {
    assert_varval_list(arg, accept_expr = FALSE)
  }

  expect_error(
    example_fun(vars(DTHSEQ = vars()))
  )
})

test_that("Test 66 : `assert_varval_list` doesn't throw an error
          if an argument is a variable-value list", {
  example_fun <- function(arg) {
    assert_varval_list(arg)
  }

  expect_invisible(
    example_fun(vars(DTHDOM = "AE", DTHSEQ = AESEQ))
  )
})

# assert_list_element ----
test_that("Test 67 : `assert_list_element` doesn't throw an error
          if the elements of a list of named lists/classes fulfill a certain condition", {
  library(admiral.test)
  expect_invisible(
    assert_list_element(vars(DTHDOM = "AE", DTHSEQ = AESEQ), "DTHSEQ", TRUE, message_text = "")
  )
  expect_invisible(
    assert_list_element(vars(DTHDOM = "AE", DTHSEQ = admiral.test::admiral_dm), "DTHSEQ",
      admiral_dm$DOMAIN == "DM",
      message_text = ""
    )
  )
})

test_that("Test 68 : `assert_list_element` throws an error
          if the elements of a list of named lists/classes don't fulfill a certain condition", {
  expect_error(
    assert_list_element(vars(DTHDOM = "AE", DTHSEQ = admiral.test::admiral_dm), "DTHSEQ",
      admiral_dm$DOMAIN == "GM",
      message_text = ""
    )
  )
})

# assert_one_to_one ----
test_that("Test 69 : `assert_one_to_one` throws an error
          if there is a one to many mapping between two lists of variables", {
  expect_error(
    assert_one_to_one(admiral.test::admiral_dm, vars(DOMAIN), vars(USUBJID))
  )
})

test_that("Test 70 : `assert_one_to_one` throws an error
          if there is a many to one mapping between two lists of variables", {
  expect_error(
    assert_one_to_one(admiral.test::admiral_dm, vars(USUBJID), vars(DOMAIN))
  )
})

# assert_date_var ----
test_that("Test 71 : `assert_date_var` throws an error
          if a variable in a dataset is not a date or datetime variable", {
  example_fun <- function(dataset, var) {
    var <- assert_symbol(enquo(var))
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
