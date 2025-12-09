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
  do.call(structure, c(list(if_else(x == "", NA_character_, x)), attributes(x))) # nolint: undesirable_function_linter
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
  do.call(structure, c(list(if_else(is.na(x), "", x)), attributes(x))) # nolint: undesirable_function_linter
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
#' @keywords utils_quo
#' @family utils_quo
#'
#' @export
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

#' Convert XXTPT Strings to Hours
#'
#' Converts timepoint strings (e.g., `PCTPT`, `VSTPT`, `EGTPT`, `ISTPT`, `LBTPT`) into
#' numeric hours. The function handles common dose-centric formats like
#' Pre-dose, Post-dose (hours/minutes), days, and special time markers.
#'
#' @param xxtpt A character vector of timepoint descriptions (e.g., from `PCTPT`, `VSTPT`,
#'   `EGTPT`, `ISTPT`, `LBTPT`, or any other `--TPT` variable)
#'
#' @details
#' The function recognizes the following patterns (all case-insensitive):
#'
#' **Special Cases:**
#' * `"Screening"` → -1
#' * `"Pre-dose"`, `"Predose"`, `"Pre-infusion"`, `"Infusion"`, `"0H"`, `"0 H"` → 0
#' * `"EOI"`, `"End of Infusion"` → 1
#' * `"Morning"`, `"Evening"` → `NA_real_`
#'
#' **Time-based Conversions:**
#' * Days: `"2D"`, `"Day 2"`, `"2 days"` → multiply by 24 (e.g., 48 hours)
#' * Hours + Minutes: `"1H30M"`, `"1 hour 30 min"` → hours + minutes/60 (e.g., 1.5)
#' * Hours only: `"1H"`, `"2 hours"`, `"0.5HR"` → as-is (e.g., 2)
#' * Minutes only: `"30M"`, `"45 min"` → divide by 60 (e.g., 0.5)
#'
#' The function supports various unit formats (H/h, HR/hr, HOUR/hour, M/m, MIN/min,
#' MINUTE/minute, D/d, DAY/day) with optional plurals and whitespace.
#'
#' @return A numeric vector of timepoints in hours
#'
#' @keywords utils_fmt
#' @family utils_fmt
#'
#' @export
#'
#' @examples
#' # Special cases
#' convert_xxtpt_to_hours(c("Screening", "Pre-dose", "EOI", "Morning"))
#'
#' # Days
#' convert_xxtpt_to_hours(c("Day 2", "2D", "3 days"))
#'
#' # Hours and minutes
#' convert_xxtpt_to_hours(c("1H", "30MIN", "1H30M", "2 hours"))
#'
#' # With Post-dose suffix
#' convert_xxtpt_to_hours(c("1h Post-dose", "30 Min Post-dose"))
convert_xxtpt_to_hours <- function(xxtpt) {
  assert_character_vector(xxtpt)

  # Handle empty input
  if (length(xxtpt) == 0) {
    return(numeric(0))
  }

  # Initialize result vector with NA
  result <- rep(NA_real_, length(xxtpt))

  # Handle NA inputs
  na_idx <- is.na(xxtpt)

  # 1. Check special cases first (exact matches, case-insensitive)

  # Screening -> -1
  screening_pattern <- regex("^screening$", ignore_case = TRUE)
  screening_idx <- str_detect(xxtpt, screening_pattern) & !na_idx
  result[screening_idx] <- -1

  # Pre-dose, Predose, Pre-infusion, Infusion, 0H, 0 H -> 0
  zero_pattern <- regex(
    "^(pre-?dose|pre-?infusion|infusion|0\\s*h(?:r|our)?s?)$",
    ignore_case = TRUE
  )
  zero_idx <- str_detect(xxtpt, zero_pattern) & is.na(result) & !na_idx
  result[zero_idx] <- 0

  # EOI, End of Infusion -> 1
  eoi_pattern <- regex("^(eoi|end\\s+of\\s+infusion)$", ignore_case = TRUE)
  eoi_idx <- str_detect(xxtpt, eoi_pattern) & is.na(result) & !na_idx
  result[eoi_idx] <- 1

  # Morning, Evening -> NA (already NA, no action needed)

  # 2. Check days (convert to hours by multiplying by 24)
  # Matches: "2D", "2 days", "Day 2", "2 day", "1.5 days"
  days_pattern <- regex(
    paste0(
      "^(?:day\\s+)?(\\d+(?:\\.\\d+)?)\\s*", # optional "day " prefix, then number
      "(?:d|day|days)",                       # d/day/days
      "(?:\\s+post-?dose)?$"                  # optional post-dose suffix
    ),
    ignore_case = TRUE,
    comments = TRUE
  )
  days_matches <- str_match(xxtpt, days_pattern)
  days_idx <- !is.na(days_matches[, 1]) & is.na(result) & !na_idx
  if (any(days_idx)) {
    result[days_idx] <- as.numeric(days_matches[days_idx, 2]) * 24
  }

  # 3. Check hours+minutes combinations
  # Matches: "1H30M", "1 hour 30 min", "1h 30m", "2HR15MIN"
  hm_pattern <- regex(
    paste0(
      "^(\\d+(?:\\.\\d+)?)\\s*",       # hours number
      "h(?:r|our)?s?",                 # h/hr/hour/hours
      "\\s*",                          # optional space
      "(\\d+(?:\\.\\d+)?)\\s*",        # minutes number
      "m(?:in|inute)?s?",              # m/min/minute/minutes
      "(?:\\s+post-?dose)?$"           # optional post-dose suffix
    ),
    ignore_case = TRUE,
    comments = TRUE
  )
  hm_matches <- str_match(xxtpt, hm_pattern)
  hm_idx <- !is.na(hm_matches[, 1]) & is.na(result) & !na_idx
  if (any(hm_idx)) {
    hours <- as.numeric(hm_matches[hm_idx, 2])
    minutes <- as.numeric(hm_matches[hm_idx, 3])
    result[hm_idx] <- hours + minutes / 60
  }

  # 4. Check time ranges (e.g., "0-6h Post-dose" -> 6, "0.5 - 6.5h" -> 6.5)
  # Process before simple hours to avoid conflicts
  # Fixed: Allow spaces around the dash and decimal numbers
  range_pattern <- regex(
    paste0(
      "^\\d+(?:\\.\\d+)?\\s*-\\s*(\\d+(?:\\.\\d+)?)\\s*", # range with end value captured, spaces allowed
      "h(?:r|our)?s?",                                     # h/hr/hour/hours
      "(?:\\s+post-?dose)?$"                               # optional post-dose suffix
    ),
    ignore_case = TRUE,
    comments = TRUE
  )
  range_matches <- str_match(xxtpt, range_pattern)
  range_idx <- !is.na(range_matches[, 1]) & is.na(result) & !na_idx
  if (any(range_idx)) {
    result[range_idx] <- as.numeric(range_matches[range_idx, 2])
  }

  # 5. Check hours only
  # Matches: "1H", "2 hours", "0.5HR", "1h Post-dose", "8 hour", "12 HOURS"
  hours_pattern <- regex(
    paste0(
      "^(\\d+(?:\\.\\d+)?)\\s*",       # number
      "h(?:r|our)?s?",                 # h/hr/hour/hours
      "(?:\\s+post-?dose)?$"           # optional post-dose suffix
    ),
    ignore_case = TRUE,
    comments = TRUE
  )
  hours_matches <- str_match(xxtpt, hours_pattern)
  hours_idx <- !is.na(hours_matches[, 1]) & is.na(result) & !na_idx
  if (any(hours_idx)) {
    result[hours_idx] <- as.numeric(hours_matches[hours_idx, 2])
  }

  # 6. Check minutes only
  # Matches: "30M", "45 min", "30 Min Post-dose", "2.5 Min"
  minutes_pattern <- regex(
    paste0(
      "^(\\d+(?:\\.\\d+)?)\\s*",       # number
      "m(?:in|inute)?s?",              # m/min/minute/minutes
      "(?:\\s+post-?dose)?$"           # optional post-dose suffix
    ),
    ignore_case = TRUE,
    comments = TRUE
  )
  minutes_matches <- str_match(xxtpt, minutes_pattern)
  minutes_idx <- !is.na(minutes_matches[, 1]) & is.na(result) & !na_idx
  if (any(minutes_idx)) {
    result[minutes_idx] <- as.numeric(minutes_matches[minutes_idx, 2]) / 60
  }

  result
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
    names <- seq_along(list)
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
