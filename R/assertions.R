#' Is an Argument a Data Frame?
#'
#' Checks if an argument is a data frame and (optionally) whether is contains
#' a set of required variables
#'
#' @param arg A function argument to be checked
#' @param required_vars A list of variables created using `exprs()`
#' @param check_is_grouped Throw an error is `dataset` is grouped? Defaults to `TRUE`.
#' @param optional Is the checked argument optional? If set to `FALSE` and `arg`
#' is `NULL` then an error is thrown
#'
#' @inheritParams assert_logical_scalar
#'
#' @return
#' The function throws an error if `arg` is not a data frame or if `arg` is a
#' data frame but misses any variable specified in `required_vars`. Otherwise,
#' the input is returned invisibly.
#'
#' @export
#'
#' @keywords assertion
#' @family assertion
#'
#' @examples
#' library(dplyr)
#' library(rlang)
#' dm <- tribble(
#'   ~STUDYID, ~USUBJID,
#'   "XYZ",    "1",
#'   "XYZ",    "2"
#' )
#'
#' example_fun <- function(dataset) {
#'   assert_data_frame(dataset, required_vars = exprs(STUDYID, USUBJID))
#' }
#'
#' example_fun(dm)
#'
#' try(example_fun(select(dm, -STUDYID)))
#'
#' try(example_fun("Not a dataset"))
#'
#' try(example_fun(group_by(dm, USUBJID)))
assert_data_frame <- function(arg,
                              required_vars = NULL,
                              check_is_grouped = TRUE,
                              optional = FALSE,
                              arg_name = rlang::caller_arg(arg),
                              message = NULL,
                              class = "assert_data_frame",
                              call = parent.frame()) {
  assert_vars(required_vars, optional = TRUE)
  assert_logical_scalar(check_is_grouped)
  assert_logical_scalar(optional)

  if (optional && is.null(arg)) {
    return(invisible(arg))
  }

  assert_s3_class(
    arg,
    cls = "data.frame",
    optional = optional,
    arg_name = arg_name,
    message = message,
    class = class,
    call = call
  )

  if (check_is_grouped && dplyr::is_grouped_df(arg)) {
    cli_abort(
      message = message %||%
        "Argument {.arg {arg_name}} must not be a grouped dataset, please `ungroup()` it.",
      class = c(class, "assert-admiraldev"),
      call = call
    )
  }

  if (!is.null(required_vars)) {
    required_vars <- vars2chr(required_vars)
    is_missing <- !required_vars %in% colnames(arg)
    if (any(is_missing)) {
      missing_vars <- required_vars[is_missing]
      if (length(missing_vars) == 1L) {
        err_msg <- "Required variable {.var {missing_vars}} is missing in {.arg {arg_name}}"
      } else {
        err_msg <- "Required variables {.var {missing_vars}} are missing in {.arg {arg_name}}"
      }
      cli_abort(
        message = message %||%
          err_msg,
        class = c(class, "assert-admiraldev"),
        call = call
      )
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
#' @param values A `character` vector of valid values for `arg`.
#' Values is converted to a lower case vector if case_sensitive = FALSE is used.
#' @param case_sensitive Should the argument be handled case-sensitive?
#' If set to `FALSE`, the argument is converted to lower case for checking the
#' permitted values and returning the argument.
#' @param optional Is the checked argument optional? If set to `FALSE` and `arg`
#' is `NULL` then an error is thrown
#' @inheritParams assert_logical_scalar
#'
#'
#' @return
#' The function throws an error if `arg` is not a character vector or if `arg`
#' is a character vector but of length > 1 or if its value is not one of the `values`
#' specified. Otherwise, the input is returned invisibly.
#'
#' @export
#'
#' @keywords assertion
#' @family assertion
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
#'
#' # handling arguments case-insensitive
#' example_fun2 <- function(msg_type) {
#'   msg_type <- assert_character_scalar(
#'     msg_type,
#'     values = c("warning", "error"),
#'     case_sensitive = FALSE
#'   )
#'   if (msg_type == "warning") {
#'     print("A warning was requested.")
#'   }
#' }
#'
#' example_fun2("Warning")
assert_character_scalar <- function(arg,
                                    values = NULL,
                                    case_sensitive = TRUE,
                                    optional = FALSE,
                                    arg_name = rlang::caller_arg(arg),
                                    message = NULL,
                                    class = "assert_character_scalar",
                                    call = parent.frame()) {
  assert_character_vector(values, optional = TRUE)
  assert_logical_scalar(optional)

  if (optional && is.null(arg)) {
    return(invisible(arg))
  }

  # change cli `.val` to end with OR instead of AND
  divid <- cli_div(theme = list(.val = list("vec-last" = ", or ", "vec_sep2" = " or ")))

  # check class and length of `arg`
  if (!is.character(arg)) {
    cli_abort(
      message = message %||%
        "Argument {.arg {arg_name}} must be a scalar of class {.cls character},
         but is {.obj_type_friendly {arg}}.",
      call = call,
      class = c(class, "assert-admiraldev")
    )
  }

  if (length(arg) != 1L) {
    cli_abort(
      message = message %||%
        "Argument {.arg {arg_name}} must be a scalar of class {.cls character},
         but is length {.val {length(arg)}}",
      call = call,
      class = c(class, "assert-admiraldev")
    )
  }

  # Create case_adjusted_arg and case_adjusted_values for the following purpose:
  #
  #   1. To simplify the comparison of arg and values; i.e. the "case_adjusted_"
  #      variables take into consideration whether `case_sensitive = TRUE`, or
  #      `case_sensitive = FALSE`.
  #
  #   2. To avoid overwriting the original "arg" and "values", so that subsequent
  #      code can refer directly to the initial function arguments: this is
  #      required whilst generating an error message if "arg" is not one of the
  #      user-specified valid values.

  if (case_sensitive) {
    case_adjusted_arg <- arg
    if (!is.null(values)) {
      case_adjusted_values <- values
    }
  } else {
    case_adjusted_arg <- tolower(arg)
    if (!is.null(values)) {
      case_adjusted_values <- tolower(values)
    }
  }

  if (!is.null(values) && case_adjusted_arg %notin% case_adjusted_values) {
    cli_abort(
      message = message %||%
        "Argument {.arg {arg_name}} must be equal to one of {.val {values}}.",
      call = call,
      class = c(class, "assert-admiraldev")
    )
  }

  invisible(case_adjusted_arg)
}

#' Is an Argument a Character Vector?
#'
#' Checks if an argument is a character vector
#'
#' @param arg A function argument to be checked
#' @param values A `character` vector of valid values for `arg`
#' @param named If set to `TRUE`, an error is issued if not all elements of the
#'   vector are named.
#' @param optional Is the checked argument optional? If set to `FALSE` and `arg`
#' is `NULL` then an error is thrown
#' @inheritParams assert_logical_scalar
#'
#'
#' @return
#' The function throws an error if `arg` is not a character vector or if
#' any element is not included in the list of valid values. Otherwise, the input
#' is returned invisibly.
#'
#' @export
#'
#' @keywords assertion
#' @family assertion
#' @examples
#' example_fun <- function(chr) {
#'   assert_character_vector(chr)
#' }
#'
#' example_fun(letters)
#'
#' try(example_fun(1:10))
#'
#' example_fun2 <- function(chr) {
#'   assert_character_vector(chr, named = TRUE)
#' }
#'
#' try(example_fun2(c(alpha = "a", "b", gamma = "c")))
assert_character_vector <- function(arg, values = NULL, named = FALSE,
                                    optional = FALSE,
                                    arg_name = rlang::caller_arg(arg),
                                    message = NULL,
                                    class = "assert_character_vector",
                                    call = parent.frame()) {
  assert_logical_scalar(named)
  assert_logical_scalar(optional)

  if (optional && is.null(arg)) {
    return(invisible(arg))
  }

  message <-
    message %||%
    ifelse(
      is.null(values),
      "Argument {.arg {arg_name}} must be {.cls character}, but is {.obj_type_friendly {arg}}.",
      "Argument {.arg {arg_name}} must be {.cls character} with values {.val {values}}."
    )

  if (!is.character(arg) ||
    (!is.null(values) && length(unique(arg[!map_lgl(arg, `%in%`, values)])) > 0L)) {
    cli_abort(
      message = message,
      call = call,
      class = c(class, "assert-admiraldev")
    )
  }

  if (named) {
    assert_named(arg, call = call, class = class, arg_name = arg_name)
  }

  invisible(arg)
}

#' Is an Argument a Logical Scalar (Boolean)?
#'
#' Checks if an argument is a logical scalar
#'
#' @param arg A function argument to be checked
#' @param optional Is the checked argument optional?\cr
#' If set to `FALSE` and `arg` is `NULL` then an error is thrown. Otherwise,
#' `NULL` is considered as valid value.
#' @param arg_name string indicating the label/symbol of the object being checked.
#' @param message string passed to `cli::cli_abort(message)`.
#' When `NULL`, default messaging is used (see examples for default messages).
#' `"{arg_name}"` can be used in messaging.
#' @inheritParams cli::cli_abort
#' @inheritParams rlang::abort
#'
#' @return
#' The function throws an error if `arg` is neither `TRUE` or `FALSE`. Otherwise,
#' the input is returned invisibly.
#'
#' @export
#'
#' @keywords assertion
#' @family assertion
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
assert_logical_scalar <- function(arg, optional = FALSE,
                                  arg_name = rlang::caller_arg(arg),
                                  message = NULL,
                                  class = "assert_logical_scalar",
                                  call = parent.frame()) {
  if (optional && is.null(arg)) {
    return(invisible(arg))
  }

  # set default error message, if not specified
  message <-
    message %||%
    "Argument {.arg {arg_name}} must be either {.val {TRUE}} or
     {.val {FALSE}}, but is {.obj_type_friendly {arg}}."

  if (!is.logical(arg) || length(arg) != 1L || is.na(arg)) {
    cli_abort(
      message = message,
      call = call,
      class = c(class, "assert-admiraldev")
    )
  }

  invisible(arg)
}

#' Is an Argument a Symbol?
#'
#' Checks if an argument is a symbol
#'
#' @param arg A function argument to be checked. Must be a `symbol`. See examples.
#' @param optional Is the checked argument optional? If set to `FALSE` and `arg`
#' is `NULL` then an error is thrown
#' @inheritParams assert_logical_scalar
#'
#'
#' @return
#' The function throws an error if `arg` is not a symbol and returns the input
#' invisibly otherwise.
#'
#' @export
#'
#' @keywords assertion
#' @family assertion
#' @examples
#' library(pharmaversesdtm)
#' library(dplyr, warn.conflicts = FALSE)
#' library(rlang)
#' data(dm)
#'
#' example_fun <- function(dat, var) {
#'   var <- assert_symbol(enexpr(var))
#'   select(dat, !!var)
#' }
#'
#' example_fun(dm, USUBJID)
#'
#' try(example_fun(dm))
#'
#' try(example_fun(dm, "USUBJID"))
#'
#' try(example_fun(dm, toupper(PARAMCD)))
assert_symbol <- function(arg,
                          optional = FALSE,
                          arg_name = rlang::caller_arg(arg),
                          message = NULL,
                          class = "assert_symbol",
                          call = parent.frame()) {
  assert_logical_scalar(optional)

  if (optional && is.null(arg)) {
    return(invisible(arg))
  }

  # set default error message, if not specified
  message <-
    message %||%
    ifelse(
      is_missing(arg),
      "Argument {.arg {arg_name}} must be a {.cls symbol}, but is missing.",
      "Argument {.arg {arg_name}} must be a {.cls symbol}, but is {.obj_type_friendly {arg}}."
    )

  if (is_missing(arg) || !is.symbol(arg)) {
    cli_abort(
      message = message,
      call = call,
      class = c(class, "assert-admiraldev")
    )
  }

  invisible(arg)
}

#' Assert Argument is an Expression
#'
#' @inheritParams assert_data_frame
#' @inheritParams assert_character_scalar
#'
#' @keywords assertion
#' @family assertion
#'
#' @return
#' The function throws an error if `arg` is not an expression, i.e. either
#' a symbol or a call, or returns the input invisibly otherwise
#'
#' @export
assert_expr <- function(arg,
                        optional = FALSE,
                        arg_name = rlang::caller_arg(arg),
                        message = NULL,
                        class = "assert_expr",
                        call = parent.frame()) {
  assert_logical_scalar(optional)

  if (optional && is.null(arg)) {
    return(invisible(arg))
  }

  arg_name <- tryCatch(force(arg_name), error = function(e) "arg")
  if (is_missing(arg)) {
    cli_abort(
      message = message %||% "Argument {.arg {arg_name}} cannot be missing.",
      call = call,
      class = c(class, "assert-admiraldev")
    )
  }

  if (!(is_call(arg) || is_expression(arg))) {
    cli_abort(
      message = message %||%
        "Argument {.arg {arg_name}} must be an expression, but is {.obj_type_friendly {arg}}",
      call = call,
      class = c(class, "assert-admiraldev")
    )
  }

  invisible(arg)
}

#' Is an Argument a Filter Condition?
#'
#' @param arg Quosure - filtering condition.
#' @param optional Logical - is the argument optional? Defaults to `FALSE`.
#' @inheritParams assert_logical_scalar
#'
#' @details Check if `arg` is a suitable filtering condition to be used in
#' functions like `subset` or `dplyr::filter`.
#'
#' @return Performs necessary checks and returns `arg` if all pass.
#' Otherwise throws an informative error.
#'
#' @export
#' @keywords assertion
#' @family assertion
#'
#' @examples
#' library(pharmaversesdtm)
#' library(dplyr, warn.conflicts = FALSE)
#' library(rlang)
#' data(dm)
#'
#' # typical usage in a function as an argument check
#' example_fun <- function(dat, x) {
#'   x <- assert_filter_cond(enexpr(x), arg_name = "x")
#'   filter(dat, !!x)
#' }
#'
#' example_fun(dm, AGE == 64)
#'
#' try(assert_filter_cond(mtcars))
assert_filter_cond <- function(arg,
                               optional = FALSE,
                               arg_name = rlang::caller_arg(arg),
                               message = NULL,
                               class = "assert_filter_cond",
                               call = parent.frame()) {
  assert_logical_scalar(optional)

  if (optional && is.null(arg)) {
    return(invisible(arg))
  }

  provided <- !is_missing(arg)
  if (provided && !(is_call(arg) || is_logical(arg))) {
    cli_abort(
      message = message %||%
        "Argument {.arg {arg_name}} must be a filter condition, but is {.obj_type_friendly {arg}}",
      class = c(class, "assert-admiraldev"),
      call = call
    )
  }

  invisible(arg)
}

#' Is an Argument a List of Variables?
#'
#' Checks if an argument is a valid list of symbols (e.g., created by `exprs()`)
#'
#' @param arg A function argument to be checked
#'
#' @param expect_names If the argument is set to `TRUE`, it is checked if all
#'   variables are named, e.g., `exprs(APERSDT = APxxSDT, APEREDT = APxxEDT)`.
#'
#' @param optional Is the checked argument optional? If set to `FALSE` and `arg`
#' is `NULL` then an error is thrown
#'
#' @inheritParams assert_logical_scalar
#'
#' @return
#' The function throws an error if `arg` is not a list of symbols (e.g., created
#' by `exprs()` and returns the input invisibly otherwise.
#'
#' @export
#'
#' @keywords assertion
#' @family assertion
#' @examples
#' library(dplyr, warn.conflicts = FALSE)
#' library(rlang)
#'
#' example_fun <- function(by_vars) {
#'   assert_vars(by_vars)
#' }
#'
#' example_fun(exprs(USUBJID, PARAMCD))
#'
#' try(example_fun(quos(USUBJID, PARAMCD)))
#'
#' try(example_fun(c("USUBJID", "PARAMCD", "VISIT")))
#'
#' try(example_fun(exprs(USUBJID, toupper(PARAMCD), desc(AVAL))))
#'
#' example_fun_name <- function(by_vars) {
#'   assert_vars(by_vars, expect_names = TRUE)
#' }
#'
#' example_fun_name(exprs(APERSDT = APxxSDT, APEREDT = APxxEDT))
#'
#' try(example_fun_name(exprs(APERSDT = APxxSDT, APxxEDT)))
assert_vars <- function(arg,
                        expect_names = FALSE,
                        optional = FALSE,
                        arg_name = rlang::caller_arg(arg),
                        message = NULL,
                        class = "assert_vars",
                        call = parent.frame()) {
  assert_logical_scalar(expect_names)
  assert_logical_scalar(optional)

  if (isTRUE(tryCatch(force(arg), error = function(e) TRUE))) {
    cli_abort(
      message = message %||%
        "Argument {.arg {arg_name}} must be a list of {.cls symbol},
         e.g., {.code exprs(USUBJID, VISIT)}.",
      class = c(class, "assert-admiraldev"),
      call = call
    )
  }

  assert_list_of(
    arg,
    "symbol",
    named = expect_names,
    optional = optional,
    arg_name = arg_name,
    message = message,
    class = class,
    call = call
  )
}

#' Is an Argument an Integer Scalar?
#'
#' Checks if an argument is an integer scalar
#'
#' @param arg A function argument to be checked
#' @param subset A subset of integers that `arg` should be part of. Should be one
#'   of `"none"` (the default), `"positive"`, `"non-negative"` or `"negative"`.
#' @param optional Is the checked argument optional? If set to `FALSE` and `arg`
#'   is `NULL` then an error is thrown
#' @inheritParams assert_logical_scalar
#'
#' @return
#' The function throws an error if `arg` is not an integer belonging to the
#' specified `subset`. Otherwise, the input is returned invisibly.
#'
#' @export
#'
#' @keywords assertion
#' @family assertion
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
assert_integer_scalar <- function(arg,
                                  subset = "none",
                                  optional = FALSE,
                                  arg_name = rlang::caller_arg(arg),
                                  message = NULL,
                                  class = "assert_integer_scalar",
                                  call = parent.frame()) {
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
    cli_abort(
      message = message %||%
        "Argument {.arg {arg_name}} must be
         {ifelse(subset == 'none', 'an', paste('a', subset))}
         integer scalar.",
      class = c(class, "assert-admiraldev"),
      call = call
    )
  }

  invisible(as.integer(arg))
}

#' Is an Argument a Numeric Vector?
#'
#' Checks if an argument is a numeric vector
#'
#' @param arg A function argument to be checked
#' @param optional Is the checked argument optional? If set to `FALSE` and `arg`
#' is `NULL` then an error is thrown
#' @inheritParams assert_logical_scalar
#'
#'
#' @return
#' The function throws an error if `arg` is not a numeric vector.
#' Otherwise, the input is returned invisibly.
#'
#' @export
#'
#' @keywords assertion
#' @family assertion
#' @examples
#' example_fun <- function(num) {
#'   assert_numeric_vector(num)
#' }
#'
#' example_fun(1:10)
#'
#' try(example_fun(letters))
assert_numeric_vector <- function(arg,
                                  optional = FALSE,
                                  arg_name = rlang::caller_arg(arg),
                                  message = NULL,
                                  class = "assert_numeric_vector",
                                  call = parent.frame()) {
  assert_logical_scalar(optional)

  if (optional && is.null(arg)) {
    return(invisible(arg))
  }

  if (!is.numeric(arg)) {
    cli_abort(
      message = message %||%
        "Argument {.arg {arg_name}} must be a numeric vector, but it {.obj_type_friendly {arg}}.",
      class = c(class, "assert-admiraldev"),
      call = call
    )
  }

  invisible(arg)
}

#' Is an Argument an Atomic Vector?
#'
#' Checks if an argument is an atomic vector
#'
#' @param arg A function argument to be checked
#' @param optional Is the checked argument optional? If set to `FALSE` and `arg`
#' is `NULL` then an error is thrown
#' @inheritParams assert_logical_scalar
#'
#'
#' @return
#' The function throws an error if `arg` is not an atomic vector.
#' Otherwise, the input is returned invisibly.
#'
#' @export
#'
#' @keywords assertion
#' @family assertion
#' @examples
#' example_fun <- function(x) {
#'   assert_atomic_vector(x)
#' }
#'
#' example_fun(1:10)
#'
#' try(example_fun(list(1, 2)))
assert_atomic_vector <- function(arg,
                                 optional = FALSE,
                                 arg_name = rlang::caller_arg(arg),
                                 message = NULL,
                                 class = "assert_atomic_vector",
                                 call = parent.frame()) {
  assert_logical_scalar(optional)

  if (optional && is.null(arg)) {
    return(invisible(arg))
  }

  if (!is.atomic(arg)) {
    cli_abort(
      message = message %||%
        "Argument {.arg {arg_name}} must be an atomic vector, but is {.obj_type_friendly {arg}}.",
      call = call,
      class = c(class, "assert-admiraldev")
    )
  }

  invisible(arg)
}

#' Is an Argument an Object of a Specific S3 Class?
#'
#' Checks if an argument is an object inheriting from the S3 class specified.
#' @param arg A function argument to be checked
#' @param cls The S3 class to check for
#' @param optional Is the checked argument optional? If set to `FALSE` and `arg`
#'   is `NULL` then an error is thrown
#' @inheritParams assert_logical_scalar
#'
#'
#' @return
#' The function throws an error if `arg` is an object which does *not* inherit from `class`.
#' Otherwise, the input is returned invisibly.
#'
#' @export
#'
#' @keywords assertion
#' @family assertion
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
assert_s3_class <- function(arg, cls,
                            optional = FALSE,
                            arg_name = rlang::caller_arg(arg),
                            message = NULL,
                            class = "assert_s3_class",
                            call = parent.frame()) {
  assert_character_scalar(cls)
  assert_logical_scalar(optional)

  if (is.null(arg) && optional) {
    return(invisible(arg))
  }

  messagge <-
    message %||%
    "Argument {.arg {arg_name}} must be class {.cls {cls}}, but is {.obj_type_friendly {arg}}."

  if (!inherits(arg, cls)) {
    cli_abort(
      message = messagge,
      class = c(class, "assert-admiraldev"),
      call = call
    )
  }

  invisible(arg)
}

#' Is an Argument a List of Objects of a Specific S3 Class or Type?
#'
#' Checks if an argument is a `list` of objects inheriting from the S3 class or type specified.
#'
#' @param arg A function argument to be checked
#' @param cls The S3 class or type to check for
#' @param named If set to `TRUE`, an error is issued if not all elements of the
#'   list are named.
#' @param optional Is the checked argument optional? If set to `FALSE` and `arg`
#'   is `NULL` then an error is thrown
#' @inheritParams assert_logical_scalar
#'
#'
#' @return
#' The function throws an error if `arg` is not a list or if `arg` is a list but
#' its elements are not objects inheriting from `class` or of type `class`.
#' Otherwise, the input is returned invisibly.
#'
#' @export
#'
#' @keywords assertion
#' @family assertion
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
#'
#' example_fun2 <- function(list) {
#'   assert_list_of(list, "numeric", named = TRUE)
#' }
#' try(example_fun2(list(1, 2, 3, d = 4)))
assert_list_of <- function(arg, cls,
                           named = FALSE,
                           optional = TRUE,
                           arg_name = rlang::caller_arg(arg),
                           message = NULL,
                           class = "assert_list_of",
                           call = parent.frame()) {
  assert_character_scalar(cls)
  assert_logical_scalar(named)
  assert_logical_scalar(optional)

  if (is.null(arg) && optional) {
    return(invisible(arg))
  }

  # check the passed argument is a list
  assert_s3_class(
    arg,
    cls = "list",
    arg_name = arg_name,
    message = message,
    class = class,
    call = call
  )
  # if list must be named, check this
  if (named) {
    assert_named(
      arg,
      arg_name = arg_name,
      message = message,
      class = class,
      call = call
    )
  }

  # check each list element is the expected class/type
  is_class <- map_lgl(arg, inherits, cls) | map_chr(arg, typeof) == cls
  if (!all(is_class)) {
    # construct supplementary message listing elements that are not correct type
    if (is.null(message)) {
      info_msg <- glue_collapse(
        glue(
          "element {{.val {{{which(!is_class)}}}}} is ",
          "{{.obj_type_friendly {{arg[[{which(!is_class)}]]}}}}"
        ),
        sep = ", ", last = ", and "
      )
      message <- c(
        "Each element of the list in argument {.arg {arg_name}}
         must be class/type {.cls {cls}}.",
        i = paste("But,", info_msg)
      )
    }

    cli_abort(
      message = message,
      class = c(class, "assert-admiraldev"),
      call = call
    )
  }

  invisible(arg)
}

#' Assert Argument is a Named List or Vector
#'
#' Assert that all elements of the argument are named.
#'
#' @inheritParams assert_data_frame
#' @inheritParams assert_logical_scalar
#'
#' @keywords assertion
#' @family assertion
#'
#' @return
#' The function throws an error if `arg` is not a named list or vector or
#' returns the input invisibly otherwise
#'
#' @export
#'
#' @examples
#' example_fun <- function(varval_list) {
#'   assert_named(varval_list)
#' }
#'
#' example_fun(list(var1 = 1, var2 = "x"))
#'
#' try(example_fun(list(1, "x")))
#'
#' try(example_fun(list(var = 1, "x")))
assert_named <- function(arg, optional = FALSE,
                         arg_name = rlang::caller_arg(arg),
                         message = NULL,
                         class = "assert_named",
                         call = parent.frame()) {
  if (optional && is.null(arg)) {
    return(invisible(arg))
  }

  # if argument is greater than length 0 and all element named, return arg invisibly
  any_unnamed <- length(arg) > 0L && !is_named(arg)
  if (isFALSE(any_unnamed)) {
    return(invisible(arg))
  }

  # get the indices of the unnamed elements for using in the error message.
  if (is.null(names(arg))) {
    indices <- seq_along(arg)
  } else {
    indices <- which(names(arg) == "")
  }

  message <- message %||%
    c("All elements of {.arg {arg_name}} argument must be named.",
      i = "The indices of the unnamed elements are {.val {indices}}"
    )

  cli_abort(
    message = message,
    call = call,
    class = c(class, "assert-admiraldev")
  )
}

#' Assert Argument is a Named List of Expressions
#'
#' @description
#' `r lifecycle::badge("deprecated")`
#'
#' This function is *deprecated*, please use `assert_expr_list()` instead.
#'
#' @inheritParams assert_data_frame
#'
#' @keywords deprecated
#' @family deprecated
#'
#' @return
#' The function throws an error if `arg` is not a named `list` of expression or
#' returns the input invisibly otherwise
#'
#' @export
assert_named_exprs <- function(arg, optional = FALSE) {
  deprecate_stop("0.5.0", "assert_named_exprs()", "assert_expr_list()")
}

#' Does a Dataset Contain All Required Variables?
#'
#' Checks if a dataset contains all required variables
#'
#' @description
#' `r lifecycle::badge("deprecated")`
#'
#' This function is *deprecated*, please use `assert_data_frame()` instead.
#'
#' @param dataset A `data.frame`
#' @param required_vars A `character` vector of variable names
#'
#'
#' @return The function throws an error if any of the required variables are
#' missing in the input dataset. Otherwise, the dataset is returned invisibly.
#'
#' @export
#'
#' @keywords deprecated
#' @family deprecated
assert_has_variables <- function(dataset, required_vars) {
  deprecate_stop("0.5.0", "assert_has_variables()", "assert_data_frame()")
}

#' Is Argument a Function?
#'
#' Checks if the argument is a function and if all expected arguments are
#' provided by the function.
#'
#' @param arg A function
#'
#' The function to be checked
#'
#' @param params A character vector
#'
#' A character vector of expected argument names for the aforementioned function in `arg`.
#' If ellipsis, `...`, is included in the function formals of the function in `arg`,
#' this argument, `params` will be ignored, accepting all values of the character vector.
#'
#' @param optional Is the checked argument optional?
#'
#' If set to `FALSE` and `arg` is `NULL` then an error is thrown.
#' @inheritParams assert_logical_scalar
#'
#' @return The function throws an error
#'
#'  - if the argument is not a function or
#'
#'  - if the function does not provide all arguments as specified for the
#'  `params` argument (assuming ellipsis is not in function formals)
#'
#' @export
#'
#' @keywords assertion
#' @family assertion
#'
#' @examples
#' example_fun <- function(fun) {
#'   assert_function(fun, params = c("x"))
#' }
#'
#' example_fun(mean)
#'
#' try(example_fun(1))
#'
#' try(example_fun(sum))
assert_function <- function(arg,
                            params = NULL,
                            optional = FALSE,
                            arg_name = rlang::caller_arg(arg),
                            message = NULL,
                            class = "assert_function",
                            call = parent.frame()) {
  assert_character_vector(params, optional = TRUE)
  assert_logical_scalar(optional)

  if (optional && is.null(arg)) {
    return(invisible(arg))
  }

  arg_name <- tryCatch(force(arg_name), error = function(e) "arg")
  if (missing(arg)) {
    cli_abort(
      message = message %||%
        "Argument {.arg {arg_name}} cannot be missing.",
      class = c(class, "assert-admiraldev"),
      call = call
    )
  }

  if (!is.function(arg)) {
    cli_abort(
      message = message %||%
        "Argument {.arg {arg_name}} must be a function, but is {.obj_type_friendly {arg}}.",
      class = c(class, "assert-admiraldev"),
      call = call
    )
  }
  if (!is.null(params)) {
    if ("..." %in% names(formals(arg))) {
      return(invisible(arg))
    }
    is_param <- params %in% names(formals(arg))

    if (!all(is_param)) {
      cli_abort(
        message = message %||%
          "{.val {params[!is_param]}} {?is/are} not {?an argument/arguments}
        of the function specified for {.arg {arg_name}}.",
        call = call,
        class = c(class, "assert-admiraldev")
      )
    }
  }

  invisible(arg)
}

#' Assert Argument is a Parameter of a Function
#'
#' @description
#' `r lifecycle::badge("deprecated")`
#'
#' This function is *deprecated*, please use `assert_function()` instead.
#'
#' @param arg The name of a function passed as a string
#' @param params A character vector of function parameters
#'
#' @keywords deprecated
#' @family deprecated
#'
#' @return
#' The function throws an error if any elements of `params` is not an argument of
#' the function given by `arg`
#'
#' @export
assert_function_param <- function(arg, params) {
  deprecate_stop("0.5.0", "assert_function_param()", "assert_function()")
}

#' Asserts That a Parameter is Provided in the Expected Unit
#'
#' Checks if a parameter (`PARAMCD`) in a dataset is provided in the expected
#' unit.
#'
#' @param dataset A `data.frame`
#' @param param Parameter code of the parameter to check
#' @param required_unit Expected unit
#' @param get_unit_expr Expression used to provide the unit of `param`
#'
#' @inheritParams assert_logical_scalar
#'
#'
#' @keywords assertion
#' @family assertion
#'
#' @return
#' The function throws an error if the unit variable differs from the
#' unit for any observation of the parameter in the input dataset. Otherwise, the
#' dataset is returned invisibly.
#'
#' @export
#'
#' @examples
#' library(dplyr)
#'
#' advs <- tribble(
#'   ~USUBJID, ~VSTESTCD, ~VSTRESN, ~VSSTRESU, ~PARAMCD, ~AVAL,
#'   "P01",    "WEIGHT",      80.1, "kg",      "WEIGHT",  80.1,
#'   "P02",    "WEIGHT",      85.7, "kg",      "WEIGHT",  85.7
#' )
#'
#' assert_unit(advs, param = "WEIGHT", required_unit = "kg", get_unit_expr = VSSTRESU)
assert_unit <- function(dataset,
                        param,
                        required_unit,
                        get_unit_expr,
                        arg_name = rlang::caller_arg(required_unit),
                        message = NULL,
                        class = "assert_unit",
                        call = parent.frame()) {
  assert_data_frame(dataset, required_vars = exprs(PARAMCD))
  assert_character_scalar(param)
  assert_character_scalar(required_unit)
  get_unit_expr <- enexpr(get_unit_expr)

  units <- dataset %>%
    mutate(`_unit` = !!get_unit_expr) %>%
    filter(PARAMCD == param & !is.na(`_unit`)) %>%
    pull(`_unit`) %>%
    unique()

  if (length(units) != 1L) {
    message <-
      message %||%
      "Multiple units {.val {units}} found for {.val {param}}. Please review and update the units."

    cli_abort(
      message = message,
      call = call,
      class = c(class, "assert-admiraldev")
    )
  }
  if (tolower(units) != tolower(required_unit)) {
    message <-
      message %||%
      "It is expected that {.val {param}} has unit of {.val {required_unit}}.
       In the input dataset the unit is {.val {units}}."

    cli_abort(
      message = message,
      call = call,
      class = c(class, "assert-admiraldev")
    )
  }

  invisible(dataset)
}

#' Asserts That a Parameter Does Not Exist in the Dataset
#'
#' Checks if a parameter (`PARAMCD`) does not exist in a dataset.
#'
#' @param dataset A `data.frame`
#' @param param Parameter code to check
#' @inheritParams assert_logical_scalar
#'
#'
#' @return
#' The function throws an error if the parameter exists in the input
#' dataset. Otherwise, the dataset is returned invisibly.
#'
#' @export
#'
#' @keywords assertion
#' @family assertion
#'
#' @examples
#' library(dplyr)
#'
#' advs <- tribble(
#'   ~USUBJID, ~VSTESTCD, ~VSTRESN, ~VSSTRESU, ~PARAMCD, ~AVAL,
#'   "P01",    "WEIGHT",      80.1, "kg",      "WEIGHT",  80.1,
#'   "P02",    "WEIGHT",      85.7, "kg",      "WEIGHT",  85.7
#' )
#' assert_param_does_not_exist(advs, param = "HR")
#' try(assert_param_does_not_exist(advs, param = "WEIGHT"))
assert_param_does_not_exist <- function(dataset,
                                        param,
                                        arg_name = rlang::caller_arg(dataset),
                                        message = NULL,
                                        class = "assert_param_does_not_exist",
                                        call = parent.frame()) {
  assert_data_frame(dataset, required_vars = exprs(PARAMCD))
  if (param %in% unique(dataset$PARAMCD)) {
    message <-
      message %||%
      "The parameter code {.val {param}} already exists in dataset {.arg {arg_name}}."

    cli_abort(
      message = message,
      class = c(class, "assert-admiraldev"),
      call = call
    )
  }
  invisible(dataset)
}

#' Is an Argument a Variable-Value List?
#'
#' Checks if the argument is a list of expressions where the expressions are
#' variable-value pairs. The value can be a symbol, a string, a numeric, an
#' expression, or `NA`.
#'
#' @param arg A function argument to be checked
#' @param required_elements A `character` vector of names that must be present in `arg`
#' @param accept_expr Should expressions on the right hand side be accepted?
#' @param accept_var Should unnamed variable names (e.g. `exprs(USUBJID)`) on the
#'   right hand side be accepted?
#' @param optional Is the checked argument optional? If set to `FALSE` and `arg`
#' is `NULL` then an error is thrown.
#' @inheritParams assert_logical_scalar
#'
#'
#' @return
#' The function throws an error if `arg` is not a list of variable-value expressions.
#' Otherwise, the input it returned invisibly.
#'
#' @keywords assertion
#' @family assertion
#' @export
#'
#' @examples
#' library(dplyr, warn.conflicts = FALSE)
#' library(rlang)
#'
#' example_fun <- function(vars) {
#'   assert_varval_list(vars)
#' }
#' example_fun(exprs(DTHDOM = "AE", DTHSEQ = AESEQ))
#'
#' try(example_fun(exprs("AE", DTSEQ = AESEQ)))
assert_varval_list <- function(arg, # nolint
                               required_elements = NULL,
                               accept_expr = TRUE,
                               accept_var = FALSE,
                               optional = FALSE,
                               arg_name = rlang::caller_arg(arg),
                               message = NULL,
                               class = "assert_varval_list",
                               call = parent.frame()) {
  assert_logical_scalar(accept_expr)
  assert_logical_scalar(accept_var)
  assert_logical_scalar(optional)
  assert_character_vector(required_elements, optional = TRUE)

  if (optional && is.null(arg)) {
    return(invisible(arg))
  }

  if (accept_expr) {
    valid_vals <- "a symbol, character scalar, numeric scalar, an expression, or {.val {NA}}"
  } else if (accept_var) {
    valid_vals <- "a symbol, character scalar, numeric scalar, variable names or {.val {NA}}"
  } else {
    valid_vals <- "a symbol, character scalar, numeric scalar, or {.val {NA}}"
  }

  if (!accept_var && (!inherits(arg, "list") || !is_named(arg))) {
    cli_abort(
      message = message %||%
        c(paste0("Argument {.arg {arg_name}} must be a named list of expressions
                 where each element is ", valid_vals, ", but is {.obj_type_friendly {arg}}."),
          i = "To create a list of expressions use {.fun exprs}."
        ),
      class = c(class, "assert-admiraldev"),
      call = call
    )
  }

  if (accept_var && (!contains_vars(arg))) {
    cli_abort(
      message = message %||%
        c(paste0("Argument {.arg {arg_name}} must be a list of expressions where
                  each element is ", valid_vals, ", but is {.obj_type_friendly {arg}}."),
          i = "To create a list of expressions use {.fun exprs}."
        ),
      class = c(class, "assert-admiraldev"),
      call = call
    )
  }

  if (!is.null(required_elements)) {
    missing_elements <- setdiff(required_elements, names(arg))
    if (length(missing_elements) >= 1L) {
      cli_abort(
        message = message %||%
          "The following required elements are missing from
           argument {.arg {arg_name}}: {.val {missing_elements}}.",
        class = c(class, "assert-admiraldev"),
        call = call
      )
    }
  }

  if (accept_expr) {
    invalids <- arg[!map_lgl(
      arg,
      ~ is.symbol(.x) ||
        is.character(.x) ||
        is.numeric(.x) ||
        is.language(.x) ||
        is.atomic(.x) && is.na(.x)
    )]
  } else {
    invalids <- arg[!map_lgl(
      arg,
      ~ is.symbol(.x) ||
        is.character(.x) ||
        is.numeric(.x) ||
        is.atomic(.x) && is.na(.x)
    )]
  }
  if (length(invalids) > 0) {
    cli_abort(
      message = message %||%
        c(
          paste0(
            "The elements of the list in argument {.arg {arg_name}} must be ",
            valid_vals, "."
          ),
          i = glue_collapse(
            glue("{{.val {names(invalids)}}} = {{.code {invalids}}} is of type
                  {{.cls {map_chr(invalids, typeof)}}}"),
            sep = ", ",
            last = ", and "
          )
        ),
      class = c(class, "assert-admiraldev"),
      call = call
    )
  }

  invisible(arg)
}

#' Is an Argument a List of Expressions?
#'
#' Checks if the argument is a list of expressions.
#'
#' @param arg A function argument to be checked
#' @param required_elements A `character` vector of names that must be present in `arg`
#' @param named If set to `TRUE`, an error is issued if not all elements of the
#'   list are named.
#' @param optional Is the checked argument optional? If set to `FALSE` and `arg`
#' is `NULL` then an error is thrown.
#' @inheritParams assert_logical_scalar
#' @return
#' The function throws an error if `arg` is not a list of expressions.
#' Otherwise, the input it returned invisibly.
#'
#' @keywords assertion
#' @family assertion
#' @export
#'
#' @examples
#' library(rlang)
#'
#' example_fun <- function(vars) {
#'   assert_expr_list(vars)
#' }
#' example_fun(exprs(DTHDOM = "AE", DTHSEQ = AESEQ))
#'
#' try(example_fun(exprs("AE", DTSEQ = AESEQ, !!list("a"), !!list("a"))))
assert_expr_list <- function(arg, # nolint
                             required_elements = NULL,
                             named = FALSE,
                             optional = FALSE,
                             arg_name = rlang::caller_arg(arg),
                             message = NULL,
                             class = "assert_expr_list",
                             call = parent.frame()) {
  assert_logical_scalar(named)
  assert_logical_scalar(optional)
  assert_character_vector(required_elements, optional = TRUE)

  if (optional && is.null(arg)) {
    return(invisible(arg))
  }

  if (!inherits(arg, "list")) {
    cli_abort(
      message = message %||%
        c("Argument {.arg {arg_name}} must be a list of expressions
          but is {.obj_type_friendly {arg}}.",
          i = "To create a list of expressions use {.fun exprs}."
        ),
      class = c(class, "assert-admiraldev"),
      call = call
    )
  }

  if (named) {
    assert_named(arg)
  }

  if (!is.null(required_elements)) {
    missing_elements <- setdiff(required_elements, names(arg))
    if (length(missing_elements) >= 1L) {
      cli_abort(
        message = message %||%
          "The following required elements are missing from
           argument {.arg {arg_name}}: {.val {missing_elements}}.",
        class = c(class, "assert-admiraldev"),
        call = call
      )
    }
  }

  invalids <- !map_lgl(
    arg,
    ~ is_call(.x) || is_expression(.x)
  )
  invalidargs <- arg[invalids]
  index <- if_else(names(invalidargs) == "", as.character(which(invalids)),
    paste0('"', names(invalidargs), '"')
  )

  if (any(invalids)) {
    cli_abort(
      message = message %||%
        c("All elements of {.arg {arg_name}} must be an expression.",
          i = glue_collapse(
            glue("{{.arg {arg_name}[[{index}]]}} =
            {{.code {invalidargs}}} is of type
                 {{.cls {map_chr(invalidargs, typeof)}}}"),
            sep = ", ",
            last = ", and "
          )
        ),
      class = c(class, "assert-admiraldev"),
      call = call
    )
  }

  invisible(arg)
}

#' Is an Element of a List of Lists/Classes Fulfilling a Condition?
#'
#' Checks if the elements of a list of named lists/classes fulfill a certain
#' condition. If not, an error is issued and all elements of the list not
#' fulfilling the condition are listed.
#'
#' @param list A list to be checked
#'   A list of named lists or classes is expected.
#' @param element The name of an element of the lists/classes
#'   A character scalar is expected.
#' @param condition Condition to be fulfilled
#'   The condition is evaluated for each element of the list. The element of the
#'   lists/classes can be referred to by its name, e.g., `censor == 0` to check
#'   the `censor` field of a class.
#' @param message_text Text to be displayed in the error message above
#'   the listing of values that do not meet the condition.
#'   The text should describe the condition to be fulfilled,
#'   e.g., `"Error in {arg_name}: the censor values must be zero."`.
#'   If `message` argument is specified, that text will be displayed and `message_text`
#'   is ignored.
#' @param ... Objects required to evaluate the condition
#'   If the condition contains objects apart from the element, they have to be
#'   passed to the function. See the second example below.
#' @inheritParams assert_logical_scalar
#'
#' @return
#' An error if the condition is not met. The input otherwise.
#'
#' @keywords assertion
#' @family assertion
#' @export
#'
assert_list_element <- function(list,
                                element,
                                condition,
                                message_text,
                                arg_name = rlang::caller_arg(list),
                                message = NULL,
                                class = "assert_list_element",
                                call = parent.frame(), ...) {
  assert_s3_class(list, "list")
  assert_character_scalar(element)
  condition <- assert_filter_cond(enexpr(condition))

  # store elements of the lists/classes in a vector named as the element #
  rlang::env_poke(current_env(), eval(element), lapply(list, `[[`, element))
  invalids <- !eval(
    condition,
    envir = list(...),
    enclos = current_env()
  )
  if (any(invalids)) {
    invalids_idx <- which(invalids)

    # construct supplementary message listing elements that are not correct type
    if (is.null(message)) {
      info_msg <- glue_collapse(
        glue(
          "{{.code {arg_name}[[{invalids_idx}]]${element} =
           {lapply(list[invalids_idx], `[[`, element)}}}"
        ),
        sep = ", ", last = ", and "
      )
      message <- c(
        message_text,
        i = paste(" But,", info_msg)
      )
    }

    cli::cli_abort(
      message = message,
      class = c(class, "assert-admiraldev"),
      call = call
    )
  }

  invisible(list)
}


#' Is There a One to One Mapping between Variables?
#'
#' Checks if there is a one to one mapping between two lists of variables.
#'
#' @param dataset Dataset to be checked
#'
#'   The variables specified for `vars1` and `vars2` are expected.
#'
#' @param vars1 First list of variables
#'
#' @param vars2 Second list of variables
#'
#' @param message string passed to `cli::cli_abort(message)`. When `NULL`, default messaging
#'   is used (see examples for default messages). `"dataset_name"` can be used in messaging.
#'
#' @param dataset_name string indicating the label/symbol of the object being checked.
#'   Default is `rlang::caller_arg(dataset)`.
#' @inheritParams assert_logical_scalar
#'
#' @return
#' An error if the condition is not meet. The input otherwise.
#'
#' @keywords assertion
#' @family assertion
#' @export
#'
#' @examples
#' library(dplyr)
#' library(rlang)
#'
#' df <- tribble(
#'   ~SPECIES, ~SPECIESN,
#'   "DOG",           1L,
#'   "CAT",           2L,
#'   "DOG",           1L
#' )
#'
#' assert_one_to_one(df, vars1 = exprs(SPECIES), vars2 = exprs(SPECIESN))
#'
#' df_many <- tribble(
#'   ~SPECIES, ~SPECIESN,
#'   "DOG",           1L,
#'   "CAT",           2L,
#'   "DOG",           3L
#' )
#'
#' try(
#'   assert_one_to_one(df_many, vars1 = exprs(SPECIES), vars2 = exprs(SPECIESN))
#' )
#'
#' try(
#'   assert_one_to_one(df_many, vars1 = exprs(SPECIESN), vars2 = exprs(SPECIES))
#' )
assert_one_to_one <- function(dataset,
                              vars1,
                              vars2,
                              dataset_name = rlang::caller_arg(dataset),
                              message = NULL,
                              class = "assert_one_to_one",
                              call = parent.frame()) {
  assert_vars(vars1)
  assert_vars(vars2)
  assert_data_frame(dataset, required_vars = expr_c(vars1, vars2))

  uniques <- unique(select(dataset, !!!vars1, !!!vars2))
  one_to_many <- uniques %>%
    group_by(!!!vars1) %>%
    filter(n() > 1) %>%
    arrange(!!!vars1)

  if (nrow(one_to_many) > 0) {
    admiraldev_environment$one_to_many <- one_to_many

    message <- message %||%
      c("For some values of {.val {vars2chr(vars1)}} there is more than one
           value of {.val {vars2chr(vars2)}}",
        "i" = "Call {.fun get_one_to_many_dataset} to get all one-to-many values."
      )

    cli::cli_abort(
      message = message,
      call = call,
      class = c(class, "assert-admiraldev")
    )
  }

  many_to_one <- uniques %>%
    group_by(!!!vars2) %>%
    filter(n() > 1) %>%
    arrange(!!!vars2)

  if (nrow(many_to_one) > 0) {
    admiraldev_environment$many_to_one <- many_to_one

    message <- message %||%
      c("There is more than one value of {.val {vars2chr(vars1)}} for some
         values of {.val {vars2chr(vars2)}}",
        "i" = "Call {.fun get_many_to_one_dataset} to get all many-to-one values."
      )

    cli::cli_abort(
      message = message,
      call = call,
      class = c(class, "assert-admiraldev")
    )
  }

  invisible(dataset)
}

#' Is a Variable in a Dataset a Date or Datetime Variable?
#'
#' Checks if a variable in a dataset is a date or datetime variable
#'
#' @param dataset The dataset where the variable is expected
#'
#' @param var The variable to check
#'
#' @param dataset_name The name of the dataset. If the argument is specified, the
#'   specified name is displayed in the error message.
#'
#' @param var_name The name of the variable. If the argument is specified, the
#'   specified name is displayed in the error message.
#' @param message (`string`)\cr
#'   string passed to `cli::cli_abort(message)`. When `NULL`, default messaging
#'   is used (see examples for default messages). `"var_name"` and `"dataset_name"`,
#'   can be used in messaging.
#' @inheritParams assert_logical_scalar
#'
#' @return
#' The function throws an error if `var` is not a date or datetime variable in
#' `dataset` and returns the input invisibly otherwise.
#'
#' @export
#'
#'
#' @keywords assertion
#'
#' @examples
#' library(lubridate)
#' library(dplyr)
#' library(rlang)
#'
#' example_fun <- function(dataset, var) {
#'   var <- assert_symbol(enexpr(var))
#'   assert_date_var(dataset = dataset, var = !!var)
#' }
#'
#' my_data <- tribble(
#'   ~USUBJID, ~ADT,
#'   "1",      ymd("2020-12-06"),
#'   "2",      ymd("")
#' )
#'
#' example_fun(
#'   dataset = my_data,
#'   var = ADT
#' )
#'
#' try(example_fun(
#'   dataset = my_data,
#'   var = USUBJID
#' ))
#'
#' example_fun2 <- function(dataset, var) {
#'   var <- assert_symbol(enexpr(var))
#'   assert_date_var(
#'     dataset = dataset,
#'     var = !!var,
#'     dataset_name = "your_data",
#'     var_name = "your_var"
#'   )
#' }
#'
#' try(example_fun2(
#'   dataset = my_data,
#'   var = USUBJID
#' ))
assert_date_var <- function(dataset,
                            var,
                            dataset_name = rlang::caller_arg(dataset),
                            var_name = rlang::caller_arg(var),
                            message = NULL,
                            class = "assert_date_var",
                            call = parent.frame()) {
  var <- assert_symbol(enexpr(var))
  assert_data_frame(dataset, required_vars = exprs(!!var))
  assert_character_scalar(dataset_name)
  assert_character_scalar(var_name)
  column <- pull(dataset, !!var)

  if (!is.instant(column)) {
    message <- message %||%
      "Column {.val {var_name}} in dataset {.code {dataset_name}} must be
       a date or datetime, but is {.obj_type_friendly {column}}."

    cli::cli_abort(
      message = message,
      call = call,
      class = c(class, "assert-admiraldev")
    )
  }

  invisible(dataset)
}

#' Is an object a date or datetime vector?
#'
#' Check if an object/vector is a date or datetime variable without needing a dataset as input
#'
#' @param arg The function argument to be checked
#'
#' @param optional Is the checked argument optional? If set to `FALSE`
#' and `arg` is `NULL` then the function `assert_date_vector` exits early and throw and error.
#' @inheritParams assert_logical_scalar
#'
#' @return
#' The function returns an error if `arg` is missing, or not a date or datetime variable
#' but otherwise returns an invisible output.
#'
#' @export
#'
#'
#' @keywords assertion
#'
#' @family assertion
#'
#' @examples
#' example_fun <- function(arg) {
#'   assert_date_vector(arg)
#' }
#'
#' example_fun(
#'   as.Date("2022-01-30", tz = "UTC")
#' )
#' try(example_fun("1993-07-14"))
assert_date_vector <- function(arg,
                               optional = FALSE,
                               arg_name = rlang::caller_arg(arg),
                               message = NULL,
                               class = "assert_date_vector",
                               call = parent.frame()) {
  assert_logical_scalar(optional)

  if (optional && is.null(arg)) {
    return(invisible(arg))
  }

  if (!is.instant(arg)) {
    cli_abort(
      message = message %||%
        "Argument {.arg {arg_name}} must be a date or datetime, but is {.obj_type_friendly {arg}}.",
      class = c(class, "assert-admiraldev"),
      call = call
    )
  }

  invisible(arg)
}

#' Are All Argument of the Same Type?
#'
#'
#' Checks if all arguments are of the same type.
#'
#' @param ... Arguments to be checked
#' @param .message character vector passed to `cli_abort(message)` when assertion fails.
#' @param .class character vector passed to `cli_abort(class)` when assertion fails.
#' @param .call environment passed to `cli_abort(call)` when assertion fails.
#'
#'
#' @return The function throws an error if not all arguments are of the same type.
#'
#' @export
#'
#' @keywords assertion
#' @family assertion
#'
#' @examples
#' example_fun <- function(true_value, false_value, missing_value) {
#'   assert_same_type(true_value, false_value, missing_value)
#' }
#'
#' example_fun(
#'   true_value = "Y",
#'   false_value = "N",
#'   missing_value = NA_character_
#' )
#'
#' try(example_fun(
#'   true_value = 1,
#'   false_value = 0,
#'   missing_value = "missing"
#' ))
assert_same_type <- function(...,
                             .message = c(
                               "Arguments {.arg {arg_names}} must be the same type.",
                               i = paste(
                                 "Argument types are",
                                 paste0("{.arg ", arg_names, "} {.cls ", types, "}",
                                   collapse = ", "
                                 )
                               )
                             ),
                             .class = "assert_same_type",
                             .call = parent.frame()) {
  args <- rlang::dots_list(..., .named = TRUE)
  arg_names <- names(args)
  types <- lapply(args, typeof)

  # if more than one type resent, return error
  if (length(unique(types)) > 1) {
    cli_abort(
      message = .message,
      class = c(.class, "assert-admiraldev"),
      call = .call
    )
  }

  invisible()
}
