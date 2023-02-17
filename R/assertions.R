#' Is an Argument a Data Frame?
#'
#' Checks if an argument is a data frame and (optionally) whether is contains
#' a set of required variables
#'
#' @param arg A function argument to be checked
#' @param required_vars A list of variables created using `exprs()`
#' @param check_is_grouped Throw an error is `dataset` is grouped? Defaults to `TRUE`.
#' @param optional Is the checked parameter optional? If set to `FALSE` and `arg`
#' is `NULL` then an error is thrown
#'
#'
#' @return
#' The function throws an error if `arg` is not a data frame or if `arg`
#' is a data frame but misses any variable specified in `required_vars`. Otherwise,
#' the input is returned invisibly.
#'
#' @export
#'
#' @keywords assertion
#' @family assertion
#'
#' @examples
#' library(admiral.test)
#' library(dplyr, warn.conflicts = FALSE)
#' library(rlang)
#' data(admiral_dm)
#'
#' example_fun <- function(dataset) {
#'   assert_data_frame(dataset, required_vars = exprs(STUDYID, USUBJID))
#' }
#'
#' example_fun(admiral_dm)
#'
#' try(example_fun(select(admiral_dm, -STUDYID)))
#'
#' try(example_fun("Not a dataset"))
assert_data_frame <- function(arg,
                              required_vars = NULL,
                              check_is_grouped = TRUE,
                              optional = FALSE) {
  assert_vars(required_vars, optional = TRUE)
  assert_logical_scalar(check_is_grouped)
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

  if (check_is_grouped && dplyr::is_grouped_df(arg)) {
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
#' @param values A `character` vector of valid values for `arg`.
#' Values is converted to a lower case vector if case_sensitive = FALSE is used.
#' @param case_sensitive Should the argument be handled case-sensitive?
#' If set to `FALSE`, the argument is converted to lower case for checking the
#' permitted values and returning the argument.
#' @param optional Is the checked parameter optional? If set to `FALSE` and `arg`
#' is `NULL` then an error is thrown
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
#' # handling parameters case-insensitive
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
                                    optional = FALSE) {
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
    err_msg <- sprintf(
      "`%s` must be one of %s but is '%s'",
      arg_name(substitute(arg)),
      enumerate(values, quote_fun = squote, conjunction = "or"),
      arg
    )
    abort(err_msg)
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
#' @param optional Is the checked parameter optional? If set to `FALSE` and `arg`
#' is `NULL` then an error is thrown
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
assert_character_vector <- function(arg, values = NULL, named = FALSE, optional = FALSE) {
  assert_logical_scalar(named)
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
      abort(paste0(
        "`", arg_name(substitute(arg)),
        "` contains invalid values:\n",
        enumerate(mismatches), "\n",
        "Valid values:\n",
        enumerate(values)
      ))
    }
  }

  if (named && length(arg) > 0) {
    if (is.null(names(arg))) {
      abort(paste0(
        "All elements of ", arg_name(substitute(arg)), " must be named.\n",
        "No element is named."
      ))
    }
    unnamed <- which(names(arg) == "")
    if (length(unnamed) > 0) {
      abort(paste0(
        "All elements of ", arg_name(substitute(arg)), " must be named.\n",
        "The following elements are not named: ", enumerate(unnamed, quote_fun = NULL)
      ))
    }
  }
}

#' Is an Argument a Logical Scalar (Boolean)?
#'
#' Checks if an argument is a logical scalar
#'
#' @param arg A function argument to be checked
#'
#' @param optional Is the checked parameter optional?
#'
#' If set to `FALSE` and `arg` is `NULL` then an error is thrown. Otherwise,
#' `NULL` is considered as valid value.
#'
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
assert_logical_scalar <- function(arg, optional = FALSE) {
  if (optional && is.null(arg)) {
    return(invisible(arg))
  }

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
#' @param arg A function argument to be checked. Must be a `symbol`. See examples.
#' @param optional Is the checked parameter optional? If set to `FALSE` and `arg`
#' is `NULL` then an error is thrown
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
#' library(admiral.test)
#' library(dplyr, warn.conflicts = FALSE)
#' library(rlang)
#' data(admiral_dm)
#'
#' example_fun <- function(dat, var) {
#'   var <- assert_symbol(enexpr(var))
#'   select(dat, !!var)
#' }
#'
#' example_fun(admiral_dm, USUBJID)
#'
#' try(example_fun(admiral_dm))
#'
#' try(example_fun(admiral_dm, "USUBJID"))
#'
#' try(example_fun(admiral_dm, toupper(PARAMCD)))
assert_symbol <- function(arg, optional = FALSE) {
  assert_logical_scalar(optional)

  if (optional && is.null(arg)) {
    return(invisible(arg))
  }

  if (is_missing(arg)) {
    err_msg <- sprintf("Argument `%s` missing, with no default", arg_name(substitute(arg)))
    abort(err_msg)
  }

  if (!is.symbol(arg)) {
    err_msg <- sprintf(
      "`%s` must be a symbol but is %s",
      arg_name(substitute(arg)),
      what_is_it(arg)
    )
    abort(err_msg)
  }

  invisible(arg)
}

#' Assert Argument is an Expression
#'
#' @inheritParams assert_data_frame
#'
#' @keywords assertion
#' @family assertion
#'
#' @return
#' The function throws an error if `arg` is not an expression, i.e. either
#' a symbol or a call, or returns the input invisibly otherwise
#'
#' @export
assert_expr <- function(arg, optional = FALSE) {
  assert_logical_scalar(optional)

  if (optional && is.null(arg)) {
    return(invisible(arg))
  }

  if (is_missing(arg)) {
    err_msg <- sprintf("Argument `%s` missing, with no default", arg_name(substitute(arg)))
    abort(err_msg)
  }

  if (!is_symbolic(arg)) {
    err_msg <- sprintf(
      "`%s` must be an expression but is %s",
      arg_name(substitute(arg)),
      what_is_it(arg)
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
#' @family assertion
#'
#' @examples
#' library(admiral.test)
#' library(dplyr, warn.conflicts = FALSE)
#' library(rlang)
#' data(admiral_dm)
#'
#' # typical usage in a function as a parameter check
#' example_fun <- function(dat, x) {
#'   x <- assert_filter_cond(enquo(x))
#'   filter(dat, !!x)
#' }
#'
#' example_fun(admiral_dm, AGE == 64)
#'
#' try(example_fun(admiral_dm, USUBJID))
assert_filter_cond <- function(arg, optional = FALSE) {
  assert_logical_scalar(optional)

  if (optional && is.null(arg)) {
    return(invisible(arg))
  }

  provided <- !is_missing(arg)
  if (provided & !(is_call(arg) | is_logical(arg))) {
    err_msg <- sprintf(
      "`%s` must be a filter condition but is %s",
      arg_name(substitute(arg)),
      what_is_it(arg)
    )
    abort(err_msg)
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
#' @param optional Is the checked parameter optional? If set to `FALSE` and `arg`
#' is `NULL` then an error is thrown
#'
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
assert_vars <- function(arg, expect_names = FALSE, optional = FALSE) {
  assert_logical_scalar(expect_names)
  assert_logical_scalar(optional)

  default_err_msg <- sprintf(
    "`%s` must be a list of symbols, e.g. `exprs(USUBJID, VISIT)`",
    arg_name(substitute(arg))
  )

  if (isTRUE(tryCatch(force(arg), error = function(e) TRUE))) {
    abort(default_err_msg)
  }

  assert_list_of(arg, "symbol", named = expect_names, optional = optional)
}

#' Is an Argument a List of Order Variables?
#'
#' Checks if an argument is a valid list of order variables created using `exprs()`
#'
#' @param arg A function argument to be checked
#' @param optional Is the checked parameter optional? If set to `FALSE` and `arg`
#' is `NULL` then an error is thrown
#'
#'
#' @return
#' The function throws an error if `arg` is not a list of variables or `desc()`
#' calls created using `exprs()` and returns the input invisibly otherwise.
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
#'   assert_order_vars(by_vars)
#' }
#'
#' example_fun(exprs(USUBJID, PARAMCD, desc(AVISITN)))
#'
#' try(example_fun(quos(USUBJID, PARAMCD)))
#'
#' try(example_fun(c("USUBJID", "PARAMCD", "VISIT")))
#'
#' try(example_fun(exprs(USUBJID, toupper(PARAMCD), -AVAL)))
assert_order_vars <- function(arg, optional = FALSE) {
  assert_logical_scalar(optional)

  default_err_msg <- paste(
    backquote(arg_name(substitute(arg))),
    "must be a list of unquoted variable names or `desc()` calls,",
    "e.g. `exprs(USUBJID, desc(VISITNUM))`"
  )

  if (isTRUE(tryCatch(force(arg), error = function(e) TRUE))) {
    abort(default_err_msg)
  }

  if (optional && is.null(arg)) {
    return(invisible(arg))
  }

  if (!inherits(arg, "list")) {
    abort(default_err_msg)
  }

  if (isFALSE(is_order_vars(arg))) {
    abort(default_err_msg)
  }

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
      ifelse(subset == "none", "an", paste("a", subset)),
      what_is_it(arg)
    )
    abort(err_msg)
  }

  invisible(as.integer(arg))
}

#' Is an Argument a Numeric Vector?
#'
#' Checks if an argument is a numeric vector
#'
#' @param arg A function argument to be checked
#' @param optional Is the checked parameter optional? If set to `FALSE` and `arg`
#' is `NULL` then an error is thrown
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
assert_numeric_vector <- function(arg, optional = FALSE) {
  assert_logical_scalar(optional)

  if (optional && is.null(arg)) {
    return(invisible(arg))
  }

  if (!is.numeric(arg)) {
    err_msg <- sprintf(
      "`%s` must be a numeric vector but is %s",
      arg_name(substitute(arg)),
      what_is_it(arg)
    )
    abort(err_msg)
  }
}

#' Is an Argument an Atomic Vector?
#'
#' Checks if an argument is an atomic vector
#'
#' @param arg A function argument to be checked
#' @param optional Is the checked parameter optional? If set to `FALSE` and `arg`
#' is `NULL` then an error is thrown
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
assert_atomic_vector <- function(arg, optional = FALSE) {
  assert_logical_scalar(optional)

  if (optional && is.null(arg)) {
    return(invisible(arg))
  }

  if (!is.atomic(arg)) {
    err_msg <- sprintf(
      "`%s` must be an atomic vector but is %s",
      arg_name(substitute(arg)),
      what_is_it(arg)
    )
    abort(err_msg)
  }
}

#' Is an Argument an Object of a Specific S3 Class?
#'
#' Checks if an argument is an object inheriting from the S3 class specified.
#' @param arg A function argument to be checked
#' @param class The S3 class to check for
#' @param optional Is the checked parameter optional? If set to `FALSE` and `arg`
#'   is `NULL` then an error is thrown
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
assert_s3_class <- function(arg, class, optional = FALSE) {
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

#' Is an Argument a List of Objects of a Specific S3 Class or Type?
#'
#' Checks if an argument is a `list` of objects inheriting from the S3 class or type specified.
#'
#' @param arg A function argument to be checked
#' @param class The S3 class or type to check for
#' @param named If set to `TRUE`, an error is issued if not all elements of the
#'   list are named.
#' @param optional Is the checked parameter optional? If set to `FALSE` and `arg`
#'   is `NULL` then an error is thrown
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
assert_list_of <- function(arg, class, named = FALSE, optional = TRUE) {
  assert_character_scalar(class)
  assert_logical_scalar(named)
  assert_logical_scalar(optional)

  if (is.null(arg) && optional) {
    return(invisible(arg))
  }

  assert_s3_class(arg, "list")

  is_class <- map_lgl(arg, inherits, class) | map_chr(arg, typeof) == class
  if (!all(is_class)) {
    info_msg <- paste(
      sprintf("\u2716 Element %s is %s", which(!is_class), map_chr(arg[!is_class], what_is_it)),
      collapse = "\n"
    )
    err_msg <- sprintf(
      "Each element of `%s` must be an object of class/type '%s' but the following are not:\n%s",
      arg_name(substitute(arg)),
      class,
      info_msg
    )
    abort(err_msg)
  }

  if (named && length(arg) > 0) {
    if (is.null(names(arg))) {
      abort(paste0(
        "All elements of ", arg_name(substitute(arg)), " must be named.\n",
        "No element is named."
      ))
    }
    unnamed <- which(names(arg) == "")
    if (length(unnamed) > 0) {
      abort(paste0(
        "All elements of ", arg_name(substitute(arg)), " must be named.\n",
        "The following elements are not named: ", enumerate(unnamed, quote_fun = NULL)
      ))
    }
  }

  invisible(arg)
}

#' Assert Argument is a Named List of Expressions
#'
#' @inheritParams assert_data_frame
#'
#' @keywords assertion
#' @family assertion
#'
#' @return
#' The function throws an error if `arg` is not a named `list` of expression or
#' returns the input invisibly otherwise
#'
#' @export
assert_named_exprs <- function(arg, optional = FALSE) {
  assert_logical_scalar(optional)

  if (optional && is.null(arg)) {
    return(invisible(arg))
  }

  if (!is.list(arg) ||
    !all(map_lgl(arg, ~ is.language(.x) | is.logical(.x))) ||
    any(names(arg) == "")) {
    err_msg <- sprintf(
      "`%s` must be a named list of expressions created using `rlang::exprs()` but is %s",
      arg_name(substitute(arg)),
      what_is_it(arg)
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
#'
#' @return The function throws an error if any of the required variables are
#' missing in the input dataset. Otherwise, the dataset is returned invisibly.
#'
#' @export
#'
#' @keywords assertion
#' @family assertion
#' @examples
#' library(admiral.test)
#' data(admiral_dm)
#'
#' assert_has_variables(admiral_dm, "STUDYID")
#'
#' try(assert_has_variables(admiral_dm, "AVAL"))
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
  invisible(dataset)
}

#' Is Argument a Function?
#'
#' Checks if the argument is a function and if all expected parameters are
#' provided by the function.
#'
#' @param arg A function argument to be checked
#'
#' @param params A character vector of expected parameter names
#'
#' @param optional Is the checked parameter optional?
#'
#' If set to `FALSE` and `arg` is `NULL` then an error is thrown.
#'
#'
#' @return The function throws an error
#'
#'  - if the argument is not a function or
#'
#'  - if the function does not provide all parameters as specified for the
#'  `params` parameter.
#'
#' @export
#'
#' @keywords assertion
#' @family assertion
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
assert_function <- function(arg, params = NULL, optional = FALSE) {
  assert_character_vector(params, optional = TRUE)
  assert_logical_scalar(optional)

  if (optional && is.null(arg)) {
    return(invisible(arg))
  }

  if (missing(arg)) {
    err_msg <- sprintf(
      "Argument `%s` missing, with no default",
      arg_name(substitute(arg))
    )
    abort(err_msg)
  }

  if (!is.function(arg)) {
    err_msg <- sprintf(
      "`%s` must be a function but is %s",
      arg_name(substitute(arg)),
      what_is_it(arg)
    )
    abort(err_msg)
  }
  if (!is.null(params)) {
    is_param <- params %in% names(formals(arg))
    if (!all(is_param)) {
      txt <- if (sum(!is_param) == 1L) {
        "%s is not a parameter of the function specified for `%s`"
      } else {
        "%s are not parameters of the function specified for `%s`"
      }
      err_msg <- sprintf(txt, enumerate(params[!is_param]), arg_name(substitute(arg)))
      abort(err_msg)
    }
  }
  invisible(arg)
}

#' Assert Argument is a Parameter of a Function
#'
#' @param arg The name of a function passed as a string
#' @param params A character vector of function parameters
#'
#' @keywords assertion
#' @family assertion
#'
#' @return
#' The function throws an error if any elements of `params` is not a parameter of
#' the function given by `arg`
#'
#' @export
#'
#' @examples
#' hello <- function(name) {
#'   print(sprintf("Hello %s", name))
#' }
#'
#' assert_function_param("hello", "name")
#'
#' try(assert_function_param("hello", "surname"))
assert_function_param <- function(arg, params) {
  assert_character_scalar(arg)
  assert_character_vector(params)
  fun <- match.fun(arg)

  is_param <- params %in% names(formals(fun))
  if (!all(is_param)) {
    txt <- if (sum(!is_param) == 1L) {
      "%s is not a parameter of `%s()`"
    } else {
      "%s are not parameters of `%s()`"
    }
    err_msg <- sprintf(txt, enumerate(params[!is_param]), arg)
    abort(err_msg)
  }

  invisible(arg)
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
#' library(tibble)
#' advs <- tribble(
#'   ~USUBJID, ~VSTESTCD, ~VSTRESN, ~VSSTRESU, ~PARAMCD, ~AVAL,
#'   "P01",    "WEIGHT",      80.1, "kg",      "WEIGHT",  80.1,
#'   "P02",    "WEIGHT",      85.7, "kg",      "WEIGHT",  85.7
#' )
#'
#' assert_unit(advs, param = "WEIGHT", required_unit = "kg", get_unit_expr = VSSTRESU)
assert_unit <- function(dataset, param, required_unit, get_unit_expr) {
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
    abort(
      paste0(
        "Multiple units ",
        enumerate(units, quote_fun = squote),
        " found for ",
        squote(param),
        ".\n",
        "Please review and update the units."
      )
    )
  }
  if (tolower(units) != tolower(required_unit)) {
    abort(
      paste0(
        "It is expected that ",
        squote(param),
        " is measured in ",
        squote(required_unit),
        ".\n",
        "In the input dataset it is measured in ",
        enumerate(units, quote_fun = squote),
        "."
      )
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
#' library(tibble)
#' advs <- tribble(
#'   ~USUBJID, ~VSTESTCD, ~VSTRESN, ~VSSTRESU, ~PARAMCD, ~AVAL,
#'   "P01",    "WEIGHT",      80.1, "kg",      "WEIGHT",  80.1,
#'   "P02",    "WEIGHT",      85.7, "kg",      "WEIGHT",  85.7
#' )
#' assert_param_does_not_exist(advs, param = "HR")
#' try(assert_param_does_not_exist(advs, param = "WEIGHT"))
assert_param_does_not_exist <- function(dataset, param) {
  assert_data_frame(dataset, required_vars = exprs(PARAMCD))
  if (param %in% unique(dataset$PARAMCD)) {
    abort(
      paste0(
        "The parameter code ",
        squote(param),
        " does already exist in `",
        arg_name(substitute(dataset)),
        "`."
      )
    )
  }
  invisible(dataset)
}

#' Is an Argument a Variable-Value List?
#'
#' Checks if the argument is a list of expressions where the expressions are
#' variable-value pairs. The value can be a symbol, a string, a numeric, or
#' `NA`. More general expression are not allowed.
#'
#' @param arg A function argument to be checked
#' @param required_elements A `character` vector of names that must be present in `arg`
#' @param accept_expr Should expressions on the right hand side be accepted?
#' @param accept_var Should unnamed variable names (e.g. `exprs(USUBJID)`) on the
#'   right hand side be accepted?
#' @param optional Is the checked parameter optional? If set to `FALSE` and `arg`
#' is `NULL` then an error is thrown.
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
                               accept_expr = FALSE,
                               accept_var = FALSE,
                               optional = FALSE) {
  assert_logical_scalar(accept_expr)
  assert_logical_scalar(accept_var)
  assert_logical_scalar(optional)
  assert_character_vector(required_elements, optional = TRUE)

  if (optional && is.null(arg)) {
    return(invisible(arg))
  }

  if (accept_expr) {
    valid_vals <- "a symbol, character scalar, numeric scalar, an expression, or `NA`"
  } else if (accept_var) {
    valid_vals <- "a symbol, character scalar, numeric scalar, variable names or `NA`"
  } else {
    valid_vals <- "a symbol, character scalar, numeric scalar, or `NA`"
  }

  if (!accept_var & (!inherits(arg, "list") || !is_named(arg))) {
    err_msg <- sprintf(
      paste0(
        "`%s` must be a named list of expressions where each element is ",
        valid_vals,
        " but it is %s\n",
        "\u2139 To create a list of expressions use `exprs()`"
      ),
      arg_name(substitute(arg)),
      what_is_it(arg)
    )
    abort(err_msg)
  }

  if (accept_var & (!contains_vars(arg))) {
    err_msg <- sprintf(
      paste0(
        "`%s` must be a list of expressions where each element is ",
        valid_vals,
        " but it is %s\n",
        "\u2139 To create a list of expressions use `exprs()`"
      ),
      arg_name(substitute(arg)),
      what_is_it(arg)
    )
    abort(err_msg)
  }

  if (!is.null(required_elements)) {
    missing_elements <- setdiff(required_elements, names(arg))
    if (length(missing_elements) >= 1L) {
      err_msg <- sprintf(
        "The following required elements are missing in `%s`: %s",
        arg_name(substitute(arg)),
        enumerate(missing_elements, quote_fun = squote)
      )
      abort(err_msg)
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
    abort(
      paste0(
        "The elements of the list ",
        arg_name(substitute(arg)),
        " must be ",
        valid_vals,
        ".\n",
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

  invisible(arg)
}

#' Is an Element of a List of Lists/Classes Fulfilling a Condition?
#'
#' Checks if the elements of a list of named lists/classes fulfill a certain
#' condition. If not, an error is issued and all elements of the list not
#' fulfilling the condition are listed.
#'
#' @param list A list to be checked
#'
#'   A list of named lists or classes is expected.
#'
#' @param element The name of an element of the lists/classes
#'
#'   A character scalar is expected.
#'
#' @param condition Condition to be fulfilled
#'
#'   The condition is evaluated for each element of the list. The element of the
#'   lists/classes can be referred to by its name, e.g., `censor == 0` to check
#'   the `censor` field of a class.
#'
#' @param message_text Text to be displayed in the message
#'
#'   The text should describe the condition to be fulfilled, e.g., "For events
#'   the censor values must be zero.".
#'
#' @param ... Objects required to evaluate the condition
#'
#'   If the condition contains objects apart from the element, they have to be
#'   passed to the function. See the second example below.
#'
#'
#' @return
#' An error if the condition is not meet. The input otherwise.
#'
#' @keywords assertion
#' @family assertion
#' @export
#'
assert_list_element <- function(list, element, condition, message_text, ...) {
  assert_s3_class(list, "list")
  assert_character_scalar(element)
  condition <- assert_filter_cond(enexpr(condition))
  assert_character_scalar(message_text)
  # store elements of the lists/classes in a vector named as the element #
  rlang::env_poke(current_env(), eval(element), lapply(list, `[[`, element))
  invalids <- !eval(
    condition,
    envir = list(...),
    enclos = current_env()
  )
  if (any(invalids)) {
    invalids_idx <- which(invalids)
    abort(
      paste0(
        message_text,
        "\n",
        paste0(
          arg_name(substitute(list)),
          "[[", invalids_idx, "]]$", element,
          " = ",
          lapply(list[invalids_idx], `[[`, element),
          collapse = "\n"
        )
      )
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
#'
#' @return
#' An error if the condition is not meet. The input otherwise.
#'
#' @keywords assertion
#' @family assertion
#' @export
#'
assert_one_to_one <- function(dataset, vars1, vars2) {
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
    abort(
      paste0(
        "For some values of ",
        vars2chr(vars1),
        " there is more than one value of ",
        vars2chr(vars2),
        ".\nCall `get_one_to_many_dataset()` to get all one to many values."
      )
    )
  }
  many_to_one <- uniques %>%
    group_by(!!!vars2) %>%
    filter(n() > 1) %>%
    arrange(!!!vars2)
  if (nrow(many_to_one) > 0) {
    admiraldev_environment$many_to_one <- many_to_one
    abort(
      paste0(
        "There is more than one value of ",
        vars2chr(vars1),
        " for some values of ",
        vars2chr(vars2),
        ".\nCall `get_many_to_one_dataset()` to get all many to one values."
      )
    )
  }
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
#' library(tibble)
#' library(lubridate)
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
assert_date_var <- function(dataset, var, dataset_name = NULL, var_name = NULL) {
  var <- assert_symbol(enexpr(var))
  assert_data_frame(dataset, required_vars = exprs(!!var))
  assert_character_scalar(dataset_name, optional = TRUE)
  assert_character_scalar(var_name, optional = TRUE)
  column <- pull(dataset, !!var)
  if (is.null(dataset_name)) {
    dataset_name <- arg_name(substitute(dataset))
  }
  if (is.null(var_name)) {
    var_name <- as_label(var)
  }
  if (!is.instant(column)) {
    abort(paste0(
      "`",
      var_name,
      "` in dataset `",
      dataset_name,
      "` is not a date or datetime variable but is ",
      friendly_type_of(column)
    ))
  }
}

#' Is an object a date or datetime vector?
#'
#' Check if an object/vector is a date or datetime variable without needing a dataset as input
#'
#' @param arg The function argument to be checked
#'
#' @param optional Is the checked parameter optional? If set to `FALSE`
#' and `arg` is `NULL` then the function `assert_date_vector` exits early and throw and error.
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
assert_date_vector <- function(arg, optional = FALSE) {
  assert_logical_scalar(optional)

  if (optional && is.null(arg)) {
    return(invisible(arg))
  }

  if (!is.instant(arg)) {
    abort(paste0(
      "`",
      deparse(substitute(arg)),
      "` must be a date or datetime variable but it's `",
      friendly_type_of(arg),
      "`"
    ))
  }
}

#' Are All Argument of the Same Type?
#'
#'
#' Checks if all arguments are of the same type.
#'
#' @param ... Arguments to be checked
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
assert_same_type <- function(...) {
  args <- rlang::dots_list(..., .named = TRUE)
  arg_names <- lapply(args, function(x) deparse(substitute(x)))
  types <- lapply(args, typeof)

  if (length(unique(types)) > 1) {
    abort(
      paste(
        "All arguments must be of the same type.",
        "Argument: Type",
        "--------------",
        paste0(names(args), ": ", types, collapse = "\n"),
        sep = "\n"
      )
    )
  }
}
