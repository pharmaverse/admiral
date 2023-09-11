#' Extract Unit From Parameter Description
#'
#' Extract the unit of a parameter from a description like "Param (unit)".
#'
#' @param x A parameter description
#'
#' @return A string
#'
#' @export
#'
#' @keywords utils_help
#' @family utils_help
#'
#' @examples
#' extract_unit("Height (cm)")
#'
#' extract_unit("Diastolic Blood Pressure (mmHg)")
extract_unit <- function(x) {
  assert_character_vector(x)

  x %>%
    str_extract("\\(.+\\)") %>%
    str_remove_all("\\(|\\)")
}

#' Convert Blank Strings Into NAs
#'
#' Turn SAS blank strings into proper R `NA`s.
#'
#' @param x Any R object
#'
#' @details
#' The default methods simply returns its input unchanged. The `character` method
#' turns every instance of `""` into `NA_character_` while preserving *all* attributes.
#' When given a data frame as input the function keeps all non-character columns
#' as is and applies the just described logic to `character` columns. Once again
#' all attributes such as labels are preserved.
#'
#' @return An object of the same class as the input
#'
#'
#' @family utils_fmt
#' @keywords utils_fmt
#'
#' @export
#'
#' @examples
#' library(tibble)
#'
#' convert_blanks_to_na(c("a", "b", "", "d", ""))
#'
#' df <- tribble(
#'   ~USUBJID,   ~RFICDTC,
#'   "1001", "2000-01-01",
#'   "1002", "2001-01-01",
#'   "1003",           ""
#' )
#' print(df)
#' convert_blanks_to_na(df)
convert_blanks_to_na <- function(x) {
  UseMethod("convert_blanks_to_na")
}

#' @export
#' @rdname convert_blanks_to_na
convert_blanks_to_na.default <- function(x) {
  x
}

#' @export
#' @rdname convert_blanks_to_na
convert_blanks_to_na.character <- function(x) {
  do.call(structure, c(list(if_else(x == "", NA_character_, x)), attributes(x)))
}

#' @export
#' @rdname convert_blanks_to_na
convert_blanks_to_na.list <- function(x) {
  lapply(x, convert_blanks_to_na)
}

#' @export
#' @rdname convert_blanks_to_na
convert_blanks_to_na.data.frame <- function(x) { # nolint
  x[] <- lapply(x, convert_blanks_to_na)
  x
}


#' Convert NAs Into Blank Strings
#'
#' Turn `NA`s to blank strings .
#'
#' @param x Any R object
#'
#' @details
#' The default methods simply returns its input unchanged. The `character` method
#' turns every instance of `NA_character_` or `NA` into `""` while preserving *all* attributes.
#' When given a data frame as input the function keeps all non-character columns
#' as is and applies the just described logic to `character`
#' all attributes such as labels are preserved.
#'
#' @return An object of the same class as the input
#'
#'
#' @family utils_fmt
#' @keywords utils_fmt
#'
#' @export
#'
#' @examples
#' library(tibble)
#'
#' convert_na_to_blanks(c("a", "b", NA, "d", NA))
#'
#' df <- tribble(
#'   ~USUBJID,   ~RFICDTC,
#'   "1001", "2000-01-01",
#'   "1002", "2001-01-01",
#'   "1003",           NA
#' )
#' print(df)
#' convert_na_to_blanks(df)
convert_na_to_blanks <- function(x) {
  UseMethod("convert_na_to_blanks")
}

#' @export
#' @rdname convert_na_to_blanks
convert_na_to_blanks.default <- function(x) {
  x
}

#' @export
#' @rdname convert_na_to_blanks
convert_na_to_blanks.character <- function(x) {
  do.call(structure, c(list(if_else(is.na(x), "", x)), attributes(x)))
}

#' @export
#' @rdname convert_na_to_blanks
convert_na_to_blanks.list <- function(x) {
  lapply(x, convert_na_to_blanks)
}

#' @export
#' @rdname convert_na_to_blanks
convert_na_to_blanks.data.frame <- function(x) { # nolint
  x_out <- x %>%
    mutate(across(everything(), convert_na_to_blanks))
  x_out
}


#' Turn a Character Vector into a List of Expressions
#'
#' Turn a character vector into a list of expressions
#'
#' @param chr A character vector
#'
#' @return A `list` of expressions as returned by [`exprs()`]
#'
#'
#' @export
#'
#' @keywords utils_quo
#' @family utils_quo
#'
#' @examples
#' chr2vars(c("USUBJID", "AVAL"))
chr2vars <- function(chr) {
  assert_character_vector(chr, optional = TRUE)
  set_names(
    exprs(!!!syms(chr)),
    names(chr)
  )
}

#' Get One to Many Values that Led to a Prior Error
#'
#' @export
#'
#'
#' @details
#' If `assert_one_to_one()` detects an issue, the one to many values are stored
#' in a dataset. This dataset can be retrieved by `get_one_to_many_dataset()`.
#'
#' Note that the function always returns the one to many values from the last
#' error that has been thrown in the current R session. Thus, after restarting
#' the R sessions `get_one_to_many_dataset()` will return `NULL` and after a
#' second error has been thrown, the dataset of the first error can no longer be
#' accessed (unless it has been saved in a variable).
#'
#' @return A `data.frame` or `NULL`
#'
#' @family utils_ds_chk
#' @keywords utils_ds_chk
#'
#' @examples
#' library(admiraldev, warn.conflicts = FALSE)
#' data(admiral_adsl)
#'
#' try(
#'   assert_one_to_one(admiral_adsl, exprs(STUDYID), exprs(SITEID))
#' )
#'
#' get_one_to_many_dataset()
get_one_to_many_dataset <- function() {
  get_dataset("one_to_many")
}

#' Get Many to One Values that Led to a Prior Error
#'
#' @export
#'
#'
#' @details
#' If `assert_one_to_one()` detects an issue, the many to one values are stored
#' in a dataset. This dataset can be retrieved by `get_many_to_one_dataset()`.
#'
#' Note that the function always returns the many to one values from the last
#' error that has been thrown in the current R session. Thus, after restarting
#' the R sessions `get_many_to_one_dataset()` will return `NULL` and after a
#' second error has been thrown, the dataset of the first error can no longer be
#' accessed (unless it has been saved in a variable).
#'
#' @return A `data.frame` or `NULL`
#'
#' @family utils_ds_chk
#' @keywords utils_ds_chk
#'
#' @examples
#' library(admiraldev, warn.conflicts = FALSE)
#' data(admiral_adsl)
#'
#' try(
#'   assert_one_to_one(admiral_adsl, exprs(SITEID), exprs(STUDYID))
#' )
#'
#' get_many_to_one_dataset()
get_many_to_one_dataset <- function() {
  get_dataset("many_to_one")
}

#' Map `"Y"` and `"N"` to Numeric Values
#'
#' Map `"Y"` and `"N"` to numeric values.
#'
#' @param arg Character vector
#'
#'
#' @keywords utils_fmt
#' @family utils_fmt
#'
#' @export
#'
#' @return `1` if `arg` equals `"Y"`, `0` if `arg` equals `"N"`, `NA_real_` otherwise
#'
#' @examples
#'
#' yn_to_numeric(c("Y", "N", NA_character_))
yn_to_numeric <- function(arg) {
  assert_character_vector(arg)
  case_when(
    arg == "Y" ~ 1,
    arg == "N" ~ 0,
    TRUE ~ NA_real_
  )
}

#' Print `source` Objects
#'
#' @param x An `source` object
#' @param ... If `indent = <numeric value>` is specified the output is indented
#'   by the specified number of characters.
#'
#' @return No return value, called for side effects
#'
#'
#' @keywords internal
#' @family utils_print
#'
#' @export
#'
#' @examples
#' print(death_event)
print.source <- function(x, ...) {
  args <- list(...)
  if ("indent" %in% names(args)) {
    indent <- args[["indent"]]
  } else {
    indent <- 0
  }
  cat(strrep(" ", indent), "<", attr(x, "class")[1], "> object\n", sep = "")
  print_named_list(x, indent = indent)
}


#' Print Named List
#'
#' @param list A named list
#' @param indent Indent
#'
#'   The output is indented by the specified number of characters.
#'
#' @return No return value, called for side effects
#'
#'
#' @keywords internal
#' @family utils_print
#'
#' @export
#'
#' @examples
#' print_named_list(death_event)
print_named_list <- function(list, indent = 0) {
  names <- names(list)
  if (is.null(names)) {
    names <- seq_len(length.out = length(list))
  }
  for (name in names) {
    if (inherits(list[[name]], "source")) {
      cat(strrep(" ", indent), name, ":\n", sep = "")
      print(list[[name]], indent = indent + 2)
    } else if (is.data.frame(list[[name]])) {
      cat(strrep(" ", indent), name, ":\n", sep = "")
      print(list[[name]])
    } else if (is.list(list[[name]])) {
      cat(strrep(" ", indent), name, ":\n", sep = "")
      if (is_named(list[[name]])) {
        print_named_list(list[[name]], indent = indent + 2)
      } else {
        for (item in list[[name]]) {
          if (is.character(item)) {
            chr_val <- dquote(item)
          } else if (is_expression(item)) {
            chr_val <- format(item)
          } else {
            chr_val <- item
          }
          cat(strrep(" ", indent + 2), paste0(chr_val, collapse = "\n"), "\n", sep = "")
        }
      }
    } else {
      if (is.character(list[[name]])) {
        chr_val <- dquote(list[[name]])
      } else if (is_expression(list[[name]])) {
        chr_val <- format(list[[name]])
      } else {
        chr_val <- list[[name]]
      }
      cat(strrep(" ", indent), name, ": ", paste0(chr_val, collapse = "\n"), "\n", sep = "")
    }
  }
}

#' Negate List of Variables
#'
#' The function adds a minus sign as prefix to each variable.
#'
#' This is useful if a list of variables should be removed from a dataset,
#' e.g., `select(!!!negate_vars(by_vars))` removes all by variables.
#'
#' @param vars List of variables created by `exprs()`
#'
#' @return A list of expressions
#'
#'
#' @export
#'
#' @keywords utils_quo
#' @family utils_quo
#'
#' @examples
#' negate_vars(exprs(USUBJID, STUDYID))
negate_vars <- function(vars = NULL) {
  assert_vars(vars, optional = TRUE)
  if (is.null(vars)) {
    NULL
  } else {
    lapply(vars, function(var) expr(-!!var))
  }
}
