#' Is an Argument a Data Frame?
#'
#' Checks if an argument is a data frame and (optionally) whether is contains
#' a set of required variables
#'
#' @param arg A function argument to be checked
#' @param required_vars A list of variables created using `vars()`
#' @param optional Is the checked parameter optional? If set to `FALSE` and `arg`
#' is `NULL` then an error is thrown
#'
#' @author Thomas Neitmann
#'
#' @return
#' The function throws an error if `arg` is not a data frame or if `arg`
#' is a data frame but misses any variable specified in `required_vars`. Otherwise,
#' the input is returned invisibly.
#'
#' @export
#'
#' @keywords assertion
#'
#' @examples
#' data(dm)
#'
#' example_fun <- function(dataset) {
#'   assert_data_frame(dataset, required_vars = vars(STUDYID, USUBJID))
#' }
#'
#' example_fun(dm)
#'
#' try(example_fun(dplyr::select(dm, -STUDYID)))
#'
#' try(example_fun("Not a dataset"))
assert_data_frame <- function(arg, required_vars = NULL, optional = FALSE) {
  assert_vars(required_vars, optional = TRUE)
  assert_logical_scalar(optional)

  if (optional && is.null(arg)) {
    return(invisible(arg))
  }

  if (!is.data.frame(arg)) {
    err_msg <- sprintf(
      "`%s` must be a data frame but is %s",
      arg_name(substitute(arg)),
      what_is_it(arg)
    )
    abort(err_msg)
  }

  if (dplyr::is_grouped_df(arg)) {
    err_msg <- sprintf(
      "`%s` is a grouped data frame, please `ungroup()` it first",
      arg_name(substitute(arg))
    )
    abort(err_msg)
  }

  if (!is.null(required_vars)) {
    required_vars <- vars2chr(required_vars)
    is_missing <- !required_vars %in% colnames(arg)
    if (any(is_missing)) {
      missing_vars <- required_vars[is_missing]
      if (length(missing_vars) == 1L) {
        err_msg <- sprintf("Required variable `%s` is missing", missing_vars)
      } else {
        err_msg <- sprintf("Required variables %s are missing", enumerate(missing_vars))
      }
      abort(err_msg)
    }
  }

  invisible(arg)
}

#' Is an Argument a Character Scalar (String)?
#'
#' Checks if an argument is a character scalar and (optionally) whether it matches
#' one of the provided `values`.
#'
#' @param arg A function argument to be checked
#' @param values A `character` vector of valid values for `arg`
#' @param optional Is the checked parameter optional? If set to `FALSE` and `arg`
#' is `NULL` then an error is thrown
#'
#' @author Thomas Neitmann
#'
#' @return
#' The function throws an error if `arg` is not a character vector or if `arg`
#' is a character vector but of length > 1 or if its value is not one of the `values`
#' specified. Otherwise, the input is returned invisibly.
#'
#' @export
#'
#' @keywords assertion
#'
#' @examples
#' example_fun <- function(msg_type) {
#'   assert_character_scalar(msg_type, values = c("warning", "error"))
#' }
#'
#' example_fun("warning")
#'
#' try(example_fun("message"))
#'
#' try(example_fun(TRUE))
assert_character_scalar <- function(arg, values = NULL, optional = FALSE) {
  assert_character_vector(values, optional = TRUE)
  assert_logical_scalar(optional)

  if (optional && is.null(arg)) {
    return(invisible(arg))
  }

  if (!is.character(arg)) {
    err_msg <- sprintf(
      "`%s` must be a character scalar but is %s",
      arg_name(substitute(arg)),
      what_is_it(arg)
    )
    abort(err_msg)
  }

  if (length(arg) != 1L) {
    err_msg <- sprintf(
      "`%s` must be a character scalar but is a character vector of length %d",
      arg_name(substitute(arg)),
      length(arg)
    )
    abort(err_msg)
  }

  if (!is.null(values) && arg %notin% values) {
    err_msg <- sprintf(
      "`%s` must be one of %s but is '%s'",
      arg_name(substitute(arg)),
      enumerate(values, quote_fun = squote, conjunction = "or"),
      arg
    )
    abort(err_msg)
  }

  invisible(arg)
}

#' Is an Argument a Character Vector?
#'
#' Checks if an argument is a character vector
#'
#' @param arg A function argument to be checked
#' @param values A `character` vector of valid values for `arg`
#' @param optional Is the checked parameter optional? If set to `FALSE` and `arg`
#' is `NULL` then an error is thrown
#'
#' @author Thomas Neitmann
#'
#' @return The function throws an error if `arg` is not a character vector or if
#' any element is not included in the list of valid values. Otherwise, the input
#' is returned invisibly.
#'
#' @export
#'
#' @keywords assertion
#'
#' @examples
#' example_fun <- function(chr) {
#'   assert_character_vector(chr)
#' }
#'
#' example_fun(letters)
#'
#' try(example_fun(1:10))
assert_character_vector <- function(arg, values = NULL, optional = FALSE) {
  assert_logical_scalar(optional)

  if (optional && is.null(arg)) {
    return(invisible(arg))
  }

  if (!is.character(arg)) {
    err_msg <- sprintf(
      "`%s` must be a character vector but is %s",
      arg_name(substitute(arg)),
      what_is_it(arg)
    )
    abort(err_msg)
  }

  assert_character_vector(values, optional = TRUE)
  if (!is.null(values)) {
    mismatches <- unique(arg[!map_lgl(arg, `%in%`, values)])
    if (length(mismatches) > 0) {
      abort(paste0("`", arg_name(substitute(arg)),
                   "` contains invalid values:\n",
                   enumerate(mismatches), "\n",
                   "Valid values:\n",
                   enumerate(values)))
    }
  }
}

#' Is an Argument a Logical Scalar (Boolean)?
#'
#' Checks if an argument is a logical scalar
#'
#' @param arg A function argument to be checked
#'
#' @author Thomas Neitmann
#'
#' @return
#' The function throws an error if `arg` is neither `TRUE` or `FALSE`. Otherwise,
#' the input is returned invisibly.
#'
#' @export
#'
#' @keywords assertion
#'
#' @examples
#' example_fun <- function(flag) {
#'   assert_logical_scalar(flag)
#' }
#'
#' example_fun(FALSE)
#'
#' try(example_fun(NA))
#'
#' try(example_fun(c(TRUE, FALSE, FALSE)))
#'
#' try(example_fun(1:10))
assert_logical_scalar <- function(arg) {
  if (!is.logical(arg) || length(arg) != 1L || is.na(arg)) {
    err_msg <- sprintf(
      "`%s` must be either `TRUE` or `FALSE` but is %s",
      arg_name(substitute(arg)),
      what_is_it(arg)
    )
    abort(err_msg)
  }

  invisible(arg)
}

#' Is an Argument a Symbol?
#'
#' Checks if an argument is a symbol
#'
#' @param arg A function argument to be checked. Must be a `quosure`. See examples.
#' @param optional Is the checked parameter optional? If set to `FALSE` and `arg`
#' is `NULL` then an error is thrown
#'
#' @author Thomas Neitmann
#'
#' @return
#' The function throws an error if `arg` is not a symbol and returns the input
#' invisibly otherwise.
#'
#' @export
#'
#' @keywords assertion
#'
#' @examples
#' data(dm)
#'
#' example_fun <- function(dat, var) {
#'   var <- assert_symbol(rlang::enquo(var))
#'   dplyr::select(dat, !!var)
#' }
#'
#' example_fun(dm, USUBJID)
#'
#' try(example_fun(dm))
#'
#' try(example_fun(dm, "USUBJID"))
#'
#' try(example_fun(dm, toupper(PARAMCD)))
assert_symbol <- function(arg, optional = FALSE) {
  assert_logical_scalar(optional)

  if (optional && quo_is_null(arg)) {
    return(invisible(arg))
  }

  if (quo_is_missing(arg)) {
    err_msg <- sprintf("Argument `%s` missing, with no default", arg_name(substitute(arg)))
    abort(err_msg)
  }

  if (!quo_is_symbol(arg)) {
    err_msg <- sprintf(
      "`%s` must be a symbol but is %s",
      arg_name(substitute(arg)),
      what_is_it(quo_get_expr(arg))
    )
    abort(err_msg)
  }

  invisible(arg)
}

#' Is an Argument a Filter Condition?
#'
#' @param arg Quosure - filtering condition.
#' @param optional Logical - is the argument optional? Defaults to `FALSE`.
#'
#' @details Check if `arg` is a suitable filtering condition to be used in
#' functions like `subset` or `dplyr::filter`.
#'
#' @return Performs necessary checks and returns `arg` if all pass.
#' Otherwise throws an informative error.
#'
#' @export
#' @keywords assertion
#' @author Ondrej Slama
#'
#' @examples
#' data(dm)
#'
#' # typical usage in a function as a parameter check
#' example_fun <- function(dat, x) {
#'   x <- assert_filter_cond(rlang::enquo(x))
#'   dplyr::filter(dat, !!x)
#' }
#'
#' example_fun(dm, AGE == 64)
#'
#' try(example_fun(dm, USUBJID))
assert_filter_cond <- function(arg, optional = FALSE) {
  stopifnot(is_quosure(arg))
  assert_logical_scalar(optional)

  if (optional && quo_is_null(arg)) {
    return(invisible(arg))
  }

  provided <- quo_not_missing(arg)
  if (!provided & !optional) {
    err_msg <- sprintf("Argument `%s` is missing, with no default", arg_name(substitute(arg)))
    abort(err_msg)
  }

  if (provided & !quo_is_call(arg)) {
    err_msg <- sprintf(
      "`%s` must be a filter condition but is %s",
      arg_name(substitute(arg)),
      what_is_it(quo_get_expr(arg))
    )
    abort(err_msg)
  }

  invisible(arg)
}

#' Is an Argument a List of Variables?
#'
#' Checks if an argument is a valid list of variables created using `vars()`
#'
#' @param arg A function argument to be checked
#' @param optional Is the checked parameter optional? If set to `FALSE` and `arg`
#' is `NULL` then an error is thrown
#'
#' @author Samia Kabi
#'
#' @return
#' The function throws an error if `arg` is not a list of variables created using `vars()`
#' and returns the input invisibly otherwise.
#'
#' @export
#'
#' @keywords assertion
#'
#' @examples
#' example_fun <- function(by_vars) {
#'   assert_vars(by_vars)
#' }
#'
#' example_fun(vars(USUBJID, PARAMCD))
#'
#' try(example_fun(exprs(USUBJID, PARAMCD)))
#'
#' try(example_fun(c("USUBJID", "PARAMCD", "VISIT")))
#'
#' try(example_fun(vars(USUBJID, toupper(PARAMCD), desc(AVAL))))
assert_vars <- function(arg, optional = FALSE) {
  assert_logical_scalar(optional)

  default_err_msg <- sprintf(
    "`%s` must be a a list of unquoted variable names, e.g. `vars(USUBJID, VISIT)`",
    arg_name(substitute(arg))
  )

  if (isTRUE(tryCatch(force(arg), error = function(e) TRUE))) {
    abort(default_err_msg)
  }

  if (optional && is.null(arg)) {
    return(invisible(arg))
  }

  if (!inherits(arg, "quosures")) {
    abort(default_err_msg)
  }

  is_symbol <- map_lgl(arg, quo_is_symbol)
  if (!all(is_symbol)) {
    expr_list <- map_chr(arg, quo_text)
    err_msg <- paste0(
      default_err_msg,
      ", but the following elements are not: ",
      enumerate(expr_list[!is_symbol])
    )
    abort(err_msg)
  }

  invisible(arg)
}

#' Is an Argument a List of Order Variables?
#'
#' Checks if an argument is a valid list of order variables created using `vars()`
#'
#' @param arg A function argument to be checked
#' @param optional Is the checked parameter optional? If set to `FALSE` and `arg`
#' is `NULL` then an error is thrown
#'
#' @author Stefan Bundfuss
#'
#' @return
#' The function throws an error if `arg` is not a list of variables or `desc()`
#' calls created using `vars()` and returns the input invisibly otherwise.
#'
#' @export
#'
#' @keywords assertion
#'
#' @examples
#' example_fun <- function(by_vars) {
#'   assert_order_vars(by_vars)
#' }
#'
#' example_fun(vars(USUBJID, PARAMCD, desc(AVISITN)))
#'
#' try(example_fun(exprs(USUBJID, PARAMCD)))
#'
#' try(example_fun(c("USUBJID", "PARAMCD", "VISIT")))
#'
#' try(example_fun(vars(USUBJID, toupper(PARAMCD), -AVAL)))
assert_order_vars <- function(arg, optional = FALSE) {
  assert_logical_scalar(optional)

  default_err_msg <- paste(
    backquote(arg_name(substitute(arg))),
    "must be a a list of unquoted variable names or `desc()` calls,",
    "e.g. `vars(USUBJID, desc(VISITNUM))`"
  )

  if (isTRUE(tryCatch(force(arg), error = function(e) TRUE))) {
    abort(default_err_msg)
  }

  if (optional && is.null(arg)) {
    return(invisible(arg))
  }

  if (!inherits(arg, "quosures")) {
    abort(default_err_msg)
  }

  assert_that(is_order_vars(arg))

  invisible(arg)
}

#' Is an Argument an Integer Scalar?
#'
#' Checks if an argument is an integer scalar
#'
#' @param arg A function argument to be checked
#' @param subset A subset of integers that `arg` should be part of. Should be one
#'   of `"none"` (the default), `"positive"`, `"non-negative"` or `"negative"`.
#' @param optional Is the checked parameter optional? If set to `FALSE` and `arg`
#'   is `NULL` then an error is thrown
#'
#' @author Thomas Neitmann
#'
#' @return
#' The function throws an error if `arg` is not an integer belonging to the
#' specified `subset`.
#'
#' @export
#'
#' @keywords assertion
#'
#' @examples
#' example_fun <- function(num1, num2) {
#'   assert_integer_scalar(num1, subset = "positive")
#'   assert_integer_scalar(num2, subset = "negative")
#' }
#'
#' example_fun(1, -9)
#'
#' try(example_fun(1.5, -9))
#'
#' try(example_fun(2, 0))
#'
#' try(example_fun("2", 0))
assert_integer_scalar <- function(arg, subset = "none", optional = FALSE) {
  subsets <- list(
    "positive" = quote(arg > 0L),
    "non-negative" = quote(arg >= 0L),
    "negative" = quote(arg < 0L),
    "none" = quote(TRUE)
  )
  assert_character_scalar(subset, values = names(subsets))
  assert_logical_scalar(optional)

  if (optional && is.null(arg)) {
    return(invisible(arg))
  }

  if (!is_integerish(arg) || length(arg) != 1L || !is.finite(arg) || !eval(subsets[[subset]])) {
    err_msg <- sprintf(
      "`%s` must be %s integer scalar but is %s",
      arg_name(substitute(arg)),
      if (subset == "none") "an" else paste("a", subset),
      what_is_it(arg)
    )
    abort(err_msg)
  }

  invisible(as.integer(arg))
}

#' Is an Argument an Object of a Specific S3 Class?
#'
#' Checks if an argument is an object inheriting from the S3 class specified.
#' @param arg A function argument to be checked
#' @param class The S3 class to check for
#' @param optional Is the checked parameter optional? If set to `FALSE` and `arg`
#'   is `NULL` then an error is thrown
#'
#' @author Thomas Neitmann
#'
#' @return
#' The function throws an error if `arg` is an object which does *not* inherit from `class`
#'
#' @export
#'
#' @keywords assertion
#'
#' @examples
#' example_fun <- function(obj) {
#'   assert_s3_class(obj, "factor")
#' }
#'
#' example_fun(as.factor(letters))
#'
#' try(example_fun(letters))
#'
#' try(example_fun(1:10))
assert_s3_class <- function(arg, class, optional = TRUE) {
  assert_character_scalar(class)
  assert_logical_scalar(optional)

  if (is.null(arg) && optional) {
    return(invisible(arg))
  }

  if (!inherits(arg, class)) {
    err_msg <- sprintf(
      "`%s` must be an object of class '%s' but is %s",
      arg_name(substitute(arg)),
      class,
      what_is_it(arg)
    )
    abort(err_msg)
  }

  invisible(arg)
}

#' Is an Argument a List of Objects of a Specific S3 Class?
#'
#' Checks if an argument is a `list` of objects inheriting from the S3 class specified.
#'
#' @param arg A function argument to be checked
#' @param class The S3 class to check for
#' @param optional Is the checked parameter optional? If set to `FALSE` and `arg`
#'   is `NULL` then an error is thrown
#'
#' @author Thomas Neitmann
#'
#' @return
#' The function throws an error if `arg` is not a list or if `arg` is a list but its
#' elements are not objects inheriting from `class`
#'
#' @export
#'
#' @keywords assertion
#'
#' @examples
#' example_fun <- function(list) {
#'   assert_list_of(list, "data.frame")
#' }
#'
#' example_fun(list(mtcars, iris))
#'
#' try(example_fun(list(letters, 1:10)))
#'
#' try(example_fun(c(TRUE, FALSE)))
assert_list_of <- function(arg, class, optional = TRUE) {
  assert_character_scalar(class)
  assert_logical_scalar(optional)

  if (is.null(arg) && optional) {
    return(invisible(arg))
  }

  assert_s3_class(arg, "list")

  is_class <- map_lgl(arg, inherits, class)
  if (!all(is_class)) {
    info_msg <- paste(
      sprintf("\u2716 Element %s is %s", which(!is_class), map_chr(arg[!is_class], what_is_it)),
      collapse = "\n"
    )
    err_msg <- sprintf(
      "Each element of `%s` must be an object of class '%s' but the following are not:\n%s",
      arg_name(substitute(arg)),
      class,
      info_msg
    )
    abort(err_msg)
  }

  invisible(arg)
}

assert_named_exprs <- function(arg, optional = FALSE) {
  assert_logical_scalar(optional)

  if (optional && is.null(arg)) {
    return(invisible(arg))
  }

  if (!is.list(arg) || !all(map_lgl(arg, is.language)) || any(names(arg) == "")) {
    err_msg <- sprintf(
      "`%s` must be a named list of expressions created using `exprs()` but is %s",
      arg_name(substitute(arg)),
      what_is_it(arg)
    )
    abort(err_msg)
  }

  invisible(arg)
}

assert_list_of_formulas <- function(arg, optional = FALSE) {
  assert_logical_scalar(optional)

  if (optional && is.null(arg)) {
    return(invisible(arg))
  }

  if (!is.list(arg) || !all(map_lgl(arg, ~is_formula(.x, lhs = TRUE))) || !all(map_lgl(arg, ~is.symbol(.x[[2L]])))) {
    err_msg <- paste(
      backquote(arg_name(substitute(arg))),
      "must be a list of formulas where each formula's left-hand side is a single",
      "variable name and each right-hand side is a function, e.g. `list(AVAL ~ mean)`"
    )
    abort(err_msg)
  }

  invisible(arg)
}

#' Does a Dataset Contain All Required Variables?
#'
#' Checks if a dataset contains all required variables
#'
#' @param dataset A `data.frame`
#' @param required_vars A `character` vector of variable names
#'
#' @author Thomas Neitmann
#'
#' @return The function throws an error if any of the required variables are
#' missing in the input dataset
#'
#' @export
#'
#' @keywords assertion
#'
#' @examples
#' data(dm)
#' assert_has_variables(dm, "STUDYID")
#' \dontrun{
#' assert_has_variables(dm, "AVAL")
#' }
assert_has_variables <- function(dataset, required_vars) {
  is_missing <- !required_vars %in% colnames(dataset)
  if (any(is_missing)) {
    missing_vars <- required_vars[is_missing]
    if (length(missing_vars) == 1L) {
      err_msg <- paste0("Required variable `", missing_vars, "` is missing.")
    } else {
      err_msg <- paste0(
        "Required variables ",
        enumerate(missing_vars),
        " are missing."
      )
    }
    abort(err_msg)
  }
}

#' Asserts That a Parameter is Provided in the Expected Unit
#'
#' Checks if a parameter (`PARAMCD`) in a dataset is provided in the expected
#' unit.
#'
#' @param dataset A `data.frame`
#' @param param
#'   Parameter code of the parameter to check
#' @param unit
#'   Expected unit
#'
#' @param unit_var
#'   Variable providing the unit
#'
#' @author Stefan Bundfuss
#'
#' @return The function throws an error if the unit variable differs from the
#'   unit for any observation of the parameter in the input dataset
#'
#' @export
#'
#' @keywords assertion
#'
#' @examples
#' data(advs)
#' assert_unit(advs, param = "WEIGHT", unit = "kg", unit_var = AVALU)
#' \dontrun{
#' assert_unit(advs, param = "WEIGHT", unit = "g", unit_var = AVALU)
#' }
assert_unit <- function(dataset, param, unit_var, unit) {
  unit_var <- assert_symbol(enquo(unit_var))
  assert_data_frame(dataset, required_vars = vars(PARAMCD, !!unit_var))
  units <-
    unique(filter(dataset, PARAMCD == param &
                    !is.na(!!unit_var))[[as_string(quo_get_expr(unit_var))]])
  if (length(units) != 1 || units != unit) {
    abort(
      paste0(
        "It is expected that ",
        param,
        " is measured in ",
        unit,
        ".\n",
        "In the input dataset it is measured in ",
        enumerate(units),
        "."
      )
    )
  }
}

#' Asserts That a Parameter Does not Exist in the Dataset
#'
#' Checks if a parameter (`PARAMCD`) does not exist in a dataset.
#'
#' @param dataset A `data.frame`
#' @param param
#'   Parameter code to check
#'
#' @author Stefan Bundfuss
#'
#' @return The function throws an error if the parameter exists in the input
#'   dataset
#'
#' @export
#'
#' @keywords assertion
#'
#' @examples
#' data(advs)
#' assert_param_does_not_exist(advs, param = "BSA")
#' \dontrun{
#' assert_param_does_not_exist(advs, param = "WEIGHT")
#' }
assert_param_does_not_exist <- function(dataset, param) {
  assert_data_frame(dataset, required_vars = vars(PARAMCD))
  if (param %in% unique(dataset$PARAMCD)) {
    abort(
      paste0(
        "The parameter code ",
        param,
        " does already exist in `",
        arg_name(substitute(dataset)),
        "`."
      )
    )
  }
}


#' Is Date/Date-time?
#'
#' Checks if a date or date-time vector was specified
#'
#' @param arg The argument to check
#'
#' @author Stefan Bundfuss
#'
#' @return `TRUE` if the argument is a date or date-time, `FALSE` otherwise
#'
#' @keywords check
#'
#' @export
#'
#' @examples
#' refdate <- lubridate::ymd("2020-01-02")
#' date <- lubridate::ymd("2020-02-03")
#' assertthat::assert_that(is_date(refdate), is_date(date))
is_date <- function(arg) {
  is.instant(arg)
}
on_failure(is_date) <- function(call, env) {
  evld <- eval(call$arg, envir = env)
  len <- length(evld)
  msg <- if (len == 0) {
    deparse(evld)
  } else if (len == 1) {
    evld
  } else {
    paste0("c(", paste(head(evld, 5), collapse = ", "), `if`(len > 5, ", ..."), ")")
  }
  paste0(
    "Argument ",
    deparse(call$arg),
    " = ",
    msg,
    " is not a lubridate date."
  )
}

#' Is Time Unit?
#'
#' Checks if a string is a time unit, i.e., 'years', 'months', 'days', 'hours',
#' 'minutes', or 'seconds'.
#'
#' @param arg The argument to check
#'
#' @author Stefan Bundfuss
#'
#' @return `TRUE` if the argument is a time unit, `FALSE` otherwise
#'
#' @keywords check
#'
#' @export
#'
#' @examples
#' unit <- "days"
#' assertthat::assert_that(is_timeunit(unit))
is_timeunit <- function(arg) {
  arg %in% c("years", "months", "days", "hours", "minutes", "seconds")
}
on_failure(is_timeunit) <- function(call, env) {
  paste0(
    "Argument ",
    deparse(call$arg),
    " = ",
    eval(call$arg, envir = env),
    " is not a valid time unit.",
    " Valid time units are 'years', 'months', 'days', 'hours', 'minutes', and 'seconds'."
  )
}

#' Check Validity of the Date Imputation Input
#'
#' Date_imputation format should be specified as "dd-mm" (e.g. "01-01")
#' or as a keyword: "FISRT", "MID", "LAST"
#'
#' @param arg The argument to check
#'
#' @author Samia Kabi
#'
#' @return `TRUE` if the argument is a valid date_imputation input, `FALSE` otherwise
#'
#' @keywords check
#'
#' @export
#'
#' @examples
#' assertthat::assert_that(is_valid_date_entry("01-02"))
#' assertthat::assert_that(is_valid_date_entry("FIRST"))
is_valid_date_entry <- function(arg) {
  pattern <- "^(01|02|03|04|05|06|07|08|09|10|11|12)-([0-9]{2})$"
  grepl(pattern, arg) | str_to_upper(arg) %in% c("FIRST", "MID", "LAST")
}
on_failure(is_valid_date_entry) <- function(call, env) {
  paste0(
    "Argument ",
    deparse(call$arg),
    " = ",
    eval(call$arg, envir = env),
    " is not a valid date entry.\n",
    "date_imputation should be specified as 'mm-dd' (e.g. '01-21') or ",
    "'FIRST', 'MID', 'LAST' to get the first/mid/last day/month"
  )
}

#' Check Validity of the Time Imputation Input
#'
#' Time_imputation format should be specified as "hh:mm:ss" (e.g. "00:00:00")
#' or as a keyword: "FISRT", "LAST"
#'
#' @param arg The argument to check
#'
#' @author Samia Kabi
#'
#' @return `TRUE` if the argument is a valid time_imputation input, `FALSE` otherwise
#'
#' @keywords check
#'
#' @export
#'
#' @examples
#' assertthat::assert_that(is_valid_time_entry("23:59:59"))
#' assertthat::assert_that(is_valid_time_entry("FIRST"))
is_valid_time_entry <- function(arg) {
  pattern <- "^([0-9]{2}):([0-9]{2}):([0-9]{2})$"
  grepl(pattern, arg) | str_to_upper(arg) %in% c("FIRST", "LAST")
}
on_failure(is_valid_time_entry) <- function(call, env) {
  paste0(
    "Argument ",
    deparse(call$arg),
    " = ",
    eval(call$arg, envir = env),
    " is not a valid time entry.\n",
    "time_imputation should be specified as 'hh:mm:ss' (e.g. '00:00:00') or ",
    "'FIRST', 'LAST' to get the first/last time of the day"
  )
}

#' Check Validity of the Minute/Second Portion of the Time Input
#'
#' Minutes and seconds are expected to range from 0 to 59
#'
#' @param arg The argument to check
#'
#' @author Samia Kabi
#'
#' @return `TRUE` if the argument is a valid min/sec input, `FALSE` otherwise
#'
#' @keywords check
#'
#' @export
#'
#' @examples
#' assertthat::assert_that(is_valid_sec_min(59))
is_valid_sec_min <- function(arg) {
  arg %in% 0:59
}
on_failure(is_valid_sec_min) <- function(call, env) {
  paste0(
    "Argument ",
    deparse(call$arg),
    " = ",
    eval(call$arg, envir = env),
    " is not a valid min/sec.\n",
    "Values must be between between 0-59"
  )
}

#' Check Validity of the Hour Portion in the Time Input
#'
#' Hours are expected to range from 0 to 23
#'
#' @param arg The argument to check
#'
#' @author Samia Kabi
#'
#' @return `TRUE` if the argument is a valid hour input, `FALSE` otherwise
#'
#' @keywords check
#'
#' @export
#'
#' @examples
#' assertthat::assert_that(is_valid_hour(20))
is_valid_hour <- function(arg) {
  arg %in% 0:23
}
on_failure(is_valid_hour) <- function(call, env) {
  paste0(
    "Argument ",
    deparse(call$arg),
    "=",
    eval(call$arg, envir = env),
    " is not a valid hour.\n",
    "Values must be between 0-23"
  )
}

#' Check Validity of the Day Portion in the Date Input
#'
#' Days are expected to range from 1 to 31
#'
#' @param arg The argument to check
#'
#' @author Samia Kabi
#'
#' @return `TRUE` if the argument is a day input, `FALSE` otherwise
#'
#' @keywords check
#'
#' @export
#'
#' @examples
#' assertthat::assert_that(is_valid_day(20))
is_valid_day <- function(arg) {
  arg %in% 1:31
}
on_failure(is_valid_day) <- function(call, env) {
  paste0(
    "Argument ",
    deparse(call$arg),
    " = ",
    eval(call$arg, envir = env),
    " is not a valid day.\n",
    "Values must be between 1-31"
  )
}

#' Check Validity of the Month Portion in the Date Input
#'
#' Days are expected to range from 1 to 12
#'
#' @param arg The argument to check
#'
#' @author Samia Kabi
#'
#' @return `TRUE` if the argument is a month input, `FALSE` otherwise
#'
#' @keywords check
#'
#' @export
#'
#' @examples
#' assertthat::assert_that(is_valid_month(12))
is_valid_month <- function(arg) {
  arg %in% 1:12
}
on_failure(is_valid_month) <- function(call, env) {
  paste0(
    "Argument ",
    deparse(call$arg),
    " = ",
    eval(call$arg, envir = env),
    " is not a valid month.\n",
    "Values for month must be between 1-12. ",
    "Please check the date_imputation input: it should be sepcified as 'dd-mm'"
  )
}

#' Is an Argument a Variable-Value List?
#'
#' Checks if the argument is a list of `quosures` where the expressions are
#' variable-value pairs. The value can be a symbol, a string, a numeric, or
#' `NA`. More general expression are not allowed.
#'
#' @param arg A function argument to be checked
#' @param optional Is the checked parameter optional? If set to `FALSE` and `arg`
#' is `NULL` then an error is thrown.
#'
#' @author Stefan Bundfuss, Thomas Neitmann
#'
#' @keywords assertion
#'
#' @export
#'
#' @examples
#' example_fun <- function(vars) {
#'   assert_varval_list(vars)
#' }
#' example_fun(vars(DTHDOM = "AE", DTHSEQ = AESEQ))
#'
#' try(example_fun(vars("AE", DTSEQ = AESEQ)))
assert_varval_list <- function(arg, optional = FALSE) {
  assert_logical_scalar(optional)

  if (optional && is.null(arg)) {
    return(invisible(arg))
  }

  err_msg <- sprintf(
    paste0(
      "`%s` must be a named list of quosures where each element is a symbol, ",
      "character scalar, numeric scalar, or `NA` but it is %s\n",
      "\u2139 To create a list of quosures use `vars()`"
    ),
    arg_name(substitute(arg)),
    what_is_it(arg)
  )

  if (!is_quosures(arg) || !is_named(arg)) {
    abort(err_msg)
  } else {
    expr_list <- map(arg, quo_get_expr)
    invalids <- expr_list[!map_lgl(
      expr_list,
      ~ is.symbol(.x) ||
        is.character(.x) ||
        is.numeric(.x) || is.atomic(.x) && is.na(.x)
    )]
    if (length(invalids) > 0) {
      abort(
        paste(
          "The elements of the list",
          arg_name(substitute(arg)),
          "must be a symbol, a character scalar, a numeric, or `NA`.\n",
          paste(
            names(invalids),
            "=",
            map_chr(invalids, expr_label),
            "is of type",
            map_chr(invalids, typeof),
            collapse = "\n"
          )
        )
      )
    }
  }

  invisible(arg)
}

is_vars <- function(arg) {
  inherits(arg, "quosures") && all(map_lgl(arg, quo_is_symbol))
}
on_failure(is_vars) <- function(call, env) {
  paste0(
    "Argument `",
    deparse(call$arg),
    "` is not a list of variables created using `vars()`"
  )
}

is_order_vars <- function(arg) {
  quo_is_desc_call <- function(quo) {
    expr <- quo_get_expr(quo)
    is_call(expr) &&
      length(expr) == 2L &&
      deparse(expr[[1L]]) == "desc" &&
      is_symbol(expr[[2L]])
  }

  inherits(arg, "quosures") &&
    all(map_lgl(arg, ~quo_is_symbol(.x) || quo_is_desc_call(.x)))
}
on_failure(is_order_vars) <- function(call, env) {
  paste0(
    backquote(deparse(call$arg)),
    " is not a valid input for `order_vars`.",
    " Valid inputs are created using `vars()` and may only contain symbols or calls involving `desc()`.\n\n", # nolint
    "  # Bad:\n",
    "  vars(ADT = impute_dtc(LBDTC), is.na(AVAL))\n\n",
    "  # Good:\n",
    "  vars(AVAL, desc(ADT))"
  )
}

is_unnamed_exprs <- function(arg) {
  is.list(arg) &&
    all(map_lgl(arg, is.language)) &&
    all(names(arg) == "")
}
on_failure(is_unnamed_exprs) <- function(call, env) {
  paste0(
    "Argument `",
    deparse(call$arg),
    "` is not a unnamed list of expressions created using `exprs()`"
  )
}

is_expr <- function(arg) {
  # Note: is.language allows both symbol and language
  !is.list(arg) & is.language(arg)
}
on_failure(is_expr) <- function(call, env) {
  paste0(
    "Argument `",
    deparse(call$arg),
    "` is not an expression created using `expr()`"
  )
}

#' Checks the length of derived records and new values of are equal
#'
#' @param x an R object
#' @param y an R object to compare the length with `x`.
#' @param x_arg Argument name of x.
#' @param y_arg Argument name of y.
#'
#' @return Logical value.
#'
#' @noRd
#'
#' @examples
#' \dontrun{
#' x <- list("x", "y")
#' y <- list("y", "z")
#' assertthat::assert_that(are_records_same(x, y, "x", "y"))
#' }
are_records_same <- function(x, y, x_arg, y_arg) {
  stopifnot(is.vector(x), is.vector(y))
  length(x) == length(y)
}

on_failure(are_records_same) <- function(call, env) {
  str_glue("`{call$x_arg}` must have consistent length to the new derived records
           of `{call$y_arg}` within `by_vars`.")
}

#' Check whether an argument is not a quosure of a missing argument
#'
#' @param x Test object
#'
#' @return TRUE or error.
#'
#' @author Thomas Neitmann, Ondrej Slama
#'
#' @noRd
#'
#' @examples
#' test_fun <- function(x) {x <- rlang::enquo(x); assertthat::assert_that(quo_not_missing(x))}
#' test_fun(my_variable) # no missing argument -> returns TRUE
#' \dontrun{
#' test_fun() # missing argument -> throws error
#' }
quo_not_missing <- function(x) {
  !rlang::quo_is_missing(x)
}
on_failure(quo_not_missing) <- function(call, env) {
  paste0("Argument `", deparse(call$x), "` is missing, with no default")
}
