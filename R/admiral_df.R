# admiral guidelines loaded
#' Tag a Dataset with the `admiral_df` Class
#'
#' Adds the `admiral_df` class to a data frame (unless it is already present),
#' preserving the existing classes such as `tbl_df`/`data.frame`. This allows a
#' dedicated `summary()` method to provide a quick diagnostic of an admiral
#' dataset (see [summary.admiral_df()]).
#'
#' @param dataset A data frame (or `NULL`)
#'
#' @return
#'   If `dataset` is `NULL`, `NULL` is returned unchanged. Otherwise `dataset`
#'   with `"admiral_df"` prepended to its class attribute.
#'
#' @keywords internal
#' @family internal
as_admiral_df <- function(dataset) {
  if (is.null(dataset)) {
    return(dataset)
  }
  if (!inherits(dataset, "admiral_df")) {
    class(dataset) <- c("admiral_df", class(dataset))
  }
  dataset
}

#' Determine the ADaM Structure of a Dataset
#'
#' Classifies a data frame as one of the common ADaM structures based on the
#' variables it contains. This is used by [summary.admiral_df()] to tailor the
#' diagnostic output.
#'
#' @param dataset A data frame
#'
#' @details
#'   The following precedence is used (the order matters as, e.g., both `"BDS"`
#'   and `"TTE"` datasets contain `PARAMCD`):
#'
#'   1. `"TTE"` -- contains `PARAMCD`, `CNSR`, and `STARTDT`.
#'   2. `"OCCDS"` -- does not contain `PARAMCD` but contains a `--DECOD`/`--TERM`
#'      variable and/or `TRTEMFL`/`AOCCxxFL`.
#'   3. `"BDS"` -- contains `PARAMCD` and `AVAL` or `AVALC`.
#'   4. `"ADSL"` -- does not contain `PARAMCD` and has one record per subject
#'      (with respect to `get_admiral_option("subject_keys")`) or contains
#'      typical subject-level treatment variables (`TRT01P`/`TRTSDT`).
#'   5. `"other"` -- none of the above.
#'
#' @return A character scalar: one of `"ADSL"`, `"BDS"`, `"OCCDS"`, `"TTE"`, or
#'   `"other"`.
#'
#' @keywords internal
#' @family internal
get_admiral_df_type <- function(dataset) {
  cols <- colnames(dataset)
  has <- function(x) all(x %in% cols)
  has_paramcd <- has("PARAMCD")

  is_tte <- has_paramcd && has("CNSR") && has("STARTDT")
  is_occds <- !has_paramcd && any(str_detect(cols, "(DECOD|TERM)$|^AOCC.*FL$|^TRTEMFL$"))
  is_bds <- has_paramcd && (has("AVAL") || has("AVALC"))
  is_adsl <- !has_paramcd && is_adsl_structure(dataset, cols)

  case_when(
    is_tte ~ "TTE",
    is_occds ~ "OCCDS",
    is_bds ~ "BDS",
    is_adsl ~ "ADSL",
    TRUE ~ "other"
  )
}

#' Check Whether a Dataset Has an ADSL (Subject-Level) Structure
#'
#' @param dataset A data frame
#' @param cols The column names of `dataset`
#'
#' @return `TRUE` if the dataset has one record per subject (with respect to
#'   `get_admiral_option("subject_keys")`) or contains typical subject-level
#'   treatment variables (`TRT01P`/`TRTSDT`), `FALSE` otherwise.
#'
#' @keywords internal
#' @family internal
is_adsl_structure <- function(dataset, cols) {
  subject_keys <- intersect(vars2chr(get_admiral_option("subject_keys")), cols)
  one_row_per_subject <-
    length(subject_keys) > 0 &&
      nrow(dataset) == nrow(distinct(dataset, !!!syms(subject_keys)))
  has_adsl_vars <- "TRT01P" %in% cols || "TRTSDT" %in% cols
  one_row_per_subject || has_adsl_vars
}

#' Summarize an `admiral_df` Dataset
#'
#' Provides a quick diagnostic of an admiral dataset, such as the number of
#' subjects and observations, the list of parameters (`PARAMCD`/`PARAM`), and
#' the list of analysis visits (`AVISIT`). The information shown is tailored to
#' the ADaM structure of the dataset (`ADSL`, `BDS`, `OCCDS`, `TTE`, or `other`),
#' as determined by [get_admiral_df_type()].
#'
#' This is useful for checking whether a derivation, e.g.,
#' [derive_param_computed()], added the records that were expected.
#'
#' @param object An `admiral_df` object
#' @param ... Not used
#'
#' @return
#'   An object of class `summary_admiral_df` (a named list) which is printed by
#'   [print.summary_admiral_df()]. It always contains `type`, `n_obs`, and
#'   `n_subjects`, and, depending on the dataset structure, `params`, `avisits`,
#'   `n_events`/`n_censored`, and `n_terms`.
#'
#' @keywords utils_print
#' @family utils_print
#'
#' @export
#'
#' @examplesx
#'
#' @caption Summarizing a computed BDS parameter
#'
#' @code
#' library(tibble)
#'
#' advs <- tribble(
#'   ~USUBJID,      ~PARAMCD, ~AVAL, ~AVISIT,
#'   "01-701-1015", "DIABP",     51, "BASELINE",
#'   "01-701-1015", "SYSBP",    121, "BASELINE",
#'   "01-701-1028", "DIABP",     79, "BASELINE",
#'   "01-701-1028", "SYSBP",    130, "BASELINE"
#' ) %>%
#'   mutate(STUDYID = "PILOT01")
#'
#' map <- derive_param_computed(
#'   advs,
#'   by_vars = exprs(USUBJID, AVISIT),
#'   parameters = c("SYSBP", "DIABP"),
#'   set_values_to = exprs(
#'     AVAL = (AVAL.SYSBP + 2 * AVAL.DIABP) / 3,
#'     PARAMCD = "MAP"
#'   )
#' )
#'
#' summary(map)
summary.admiral_df <- function(object, ...) {
  # work on a plain tibble so derived sub-datasets don't inherit `admiral_df`
  class(object) <- setdiff(class(object), "admiral_df")
  cols <- colnames(object)
  type <- get_admiral_df_type(object)

  # `USUBJID` identifies a subject; prefer it for counting because other
  # subject_keys (e.g. `STUDYID`) may be `NA` on records added by derivations
  # which only populate `by_vars` (see `derive_param_computed()`).
  if ("USUBJID" %in% cols) {
    subject_keys <- "USUBJID"
  } else {
    subject_keys <- intersect(vars2chr(get_admiral_option("subject_keys")), cols)
  }
  n_subjects <- if (length(subject_keys) > 0) {
    nrow(distinct(object, !!!syms(subject_keys)))
  } else {
    NA_integer_
  }

  out <- list(
    type = type,
    n_obs = nrow(object),
    n_subjects = n_subjects,
    subject_keys = subject_keys
  )

  if ("PARAMCD" %in% cols) {
    params <- distinct(object, PARAMCD)
    if ("PARAM" %in% cols) {
      params <- distinct(object, PARAMCD, PARAM)
    }
    out$params <- arrange(params, PARAMCD)
  }

  if ("AVISIT" %in% cols) {
    if ("AVISITN" %in% cols) {
      visits <- object %>%
        distinct(AVISIT, AVISITN) %>%
        arrange(AVISITN)
    } else {
      visits <- object %>%
        distinct(AVISIT) %>%
        arrange(AVISIT)
    }
    out$avisits <- visits$AVISIT
  }

  if (type == "TTE" && "CNSR" %in% cols) {
    out$n_events <- sum(object$CNSR == 0, na.rm = TRUE)
    out$n_censored <- sum(object$CNSR != 0, na.rm = TRUE)
  }

  if (type == "OCCDS") {
    decod_var <- str_subset(cols, "DECOD$")[1]
    if (!is.na(decod_var)) {
      out$n_terms <- n_distinct(object[[decod_var]])
      out$decod_var <- decod_var
    }
  }

  class(out) <- c("summary_admiral_df", "list")
  out
}

#' Print a Summary of an `admiral_df` Dataset
#'
#' @param x A `summary_admiral_df` object created by [summary.admiral_df()]
#' @param ... Not used
#'
#' @return No return value, called for side effects. The `summary_admiral_df`
#'   object is returned invisibly.
#'
#' @keywords utils_print
#' @family utils_print
#'
#' @export
#'
#' @examplesx
#'
#' @caption Printing a summary of an `admiral_df` dataset
#'
#' @code
#' library(tibble)
#'
#' advs <- tribble(
#'   ~STUDYID,  ~USUBJID,      ~PARAMCD, ~AVAL, ~AVISIT,
#'   "PILOT01", "01-701-1015", "DIABP",     51, "BASELINE",
#'   "PILOT01", "01-701-1015", "SYSBP",    121, "BASELINE"
#' ) %>%
#'   as_admiral_df()
#'
#' print(summary(advs))
print.summary_admiral_df <- function(x, ...) {
  cli_text("<admiral_df> summary  --  type: {x$type}")

  subj_label <- if (length(x$subject_keys) > 0) {
    paste0(" (", paste(x$subject_keys, collapse = ", "), ")")
  } else {
    ""
  }
  cli_text("Subjects{subj_label}: {x$n_subjects}")
  cli_text("Observations: {x$n_obs}")

  if (!is.null(x$params)) {
    params <- paste(x$params$PARAMCD, collapse = ", ")
    cli_text("Parameters (PARAMCD): {params}")
  }

  if (!is.null(x$avisits)) {
    avisits <- paste(x$avisits, collapse = ", ")
    cli_text("Analysis visits (AVISIT): {avisits}")
  }

  if (!is.null(x$n_events)) {
    cli_text("Events: {x$n_events} | Censored: {x$n_censored}")
  }

  if (!is.null(x$n_terms)) {
    cli_text("Distinct terms ({x$decod_var}): {x$n_terms}")
  }

  invisible(x)
}
