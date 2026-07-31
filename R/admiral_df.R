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
#' @param keys Optional character vector of the intended key variables of
#'   `dataset`, i.e. the variables it should have one record per, e.g. the
#'   `by_vars` (plus `PARAMCD`) of the derivation that produced it. When
#'   supplied, it is stored in the `"admiral_keys"` attribute and used by
#'   [summary.admiral_df()] to check the record structure instead of guessing
#'   it. `NULL` leaves any existing attribute untouched.
#'
#' @return
#'   If `dataset` is `NULL`, `NULL` is returned unchanged. Otherwise `dataset`
#'   with `"admiral_df"` prepended to its class attribute (and, if `keys` is
#'   supplied, an `"admiral_keys"` attribute).
#'
#' @keywords internal
#' @family internal
as_admiral_df <- function(dataset, keys = NULL) {
  if (is.null(dataset)) {
    return(dataset)
  }
  if (!inherits(dataset, "admiral_df")) {
    class(dataset) <- c("admiral_df", class(dataset))
  }
  if (!is.null(keys)) {
    attr(dataset, "admiral_keys") <- keys
  }
  dataset
}

#' Determine the ADaM Dataset Type
#'
#' Classifies a data frame as one of the common ADaM dataset types based on the
#' variables it contains. This is used by [summary.admiral_df()] to tailor the
#' diagnostic output.
#'
#' Note that this is the *type* of the dataset (`ADSL`, `BDS`, ...), not its
#' record structure: the latter is the set of key variables it has one record
#' per, see [infer_admiral_keys()].
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

#' Extract the Key Variables of a Dataset from a `metacore` Specification
#'
#' Reads the key variables (the variables the dataset should have one record
#' per) for a given dataset from a `{metacore}` specification object. The keys
#' are taken from the `key_seq` column of the `ds_vars` table and returned in
#' key order. These can be used, e.g., to check whether a derivation preserved
#' the expected record structure (see [summary.admiral_df()]).
#'
#' @param metacore A `{metacore}` object (of class `Metacore`), typically created
#'   by `metacore::spec_to_metacore()` or `metacore::metacore()`.
#'
#' @param dataset_name The name of the dataset whose keys should be extracted,
#'   e.g. `"ADVS"`. If the `metacore` object describes a single dataset (for
#'   example after `metacore::select_dataset()`), this may be left `NULL` and the
#'   sole dataset is used.
#'
#' @details
#'   The keys themselves are retrieved with `metacore::get_keys()`, so admiral
#'   tracks `{metacore}`'s own definition of the key variables (this function is a
#'   thin wrapper that returns the key variable names as a character vector and
#'   adds a `dataset_name = NULL` convenience for single-dataset specs).
#'
#'   `{metacore}` is only needed when this function is called; it is a suggested,
#'   not a hard, dependency of admiral. If the package is not installed an
#'   informative error (with an install prompt) is raised.
#'
#' @return
#'   A character vector of key variable names in key order. If the dataset has no
#'   key variables defined in the specification, a zero-length character vector is
#'   returned (with a warning).
#'
#' @keywords utils_help
#' @family utils_help
#'
#' @export
#'
#' @examples
#' # extracting the keys defined in an ADaM specification
#' \dontrun{
#' spec <- metacore::spec_to_metacore(metacore::metacore_example("pilot_ADaM.xlsx"))
#' get_admiral_keys(spec, "ADSL")
#' }
get_admiral_keys <- function(metacore, dataset_name = NULL) {
  rlang::check_installed(
    "metacore",
    reason = "to extract dataset keys from an ADaM specification."
  )
  if (!inherits(metacore, "Metacore")) {
    cli_abort(
      "{.arg metacore} must be a {.cls Metacore} object,
       not {.obj_type_friendly {metacore}}."
    )
  }
  assert_character_scalar(dataset_name, optional = TRUE)

  # allow `dataset_name = NULL` when the spec describes a single dataset
  if (is.null(dataset_name)) {
    available <- unique(metacore$ds_vars$dataset)
    if (length(available) != 1L) {
      cli_abort(c(
        "{.arg dataset_name} must be supplied when the {.cls Metacore} object
         describes more than one dataset.",
        i = "Available datasets: {.val {available}}."
      ))
    }
    dataset_name <- available
  }

  # delegate to metacore so we track its key-variable definition; splice the
  # value in with `!!` because `get_keys()` captures its `dataset` argument via
  # NSE and would otherwise filter for the literal symbol name. `as.character()`
  # drops the `label` attribute carried on the `variable` column.
  keys <- as.character(
    rlang::inject(metacore::get_keys(metacore, !!dataset_name))$variable
  )

  if (length(keys) == 0L) {
    cli_warn(
      "No key variables ({.var key_seq}) are defined for dataset
       {.val {dataset_name}} in the {.cls Metacore} object."
    )
  }

  keys
}

#' Tag a Dataset with its Key Variables from a `metacore` Specification
#'
#' Convenience wrapper that extracts the key variables for `dataset_name` from a
#' `{metacore}` specification (via [get_admiral_keys()]) and stores them on
#' `dataset` in the `"admiral_keys"` attribute, also tagging it as an
#' `admiral_df` (see [as_admiral_df()]). [summary.admiral_df()] uses this
#' attribute to check the record structure, i.e. that the dataset really does
#' have one record per key.
#'
#' @param dataset A data frame
#' @param metacore A `{metacore}` object, see [get_admiral_keys()]
#' @param dataset_name The dataset name, see [get_admiral_keys()]
#'
#' @return `dataset` with an `"admiral_keys"` attribute and the `admiral_df`
#'   class.
#'
#' @details
#'   The core `{dplyr}` verbs (`mutate()`, `filter()`, `arrange()`, `select()`,
#'   `rename()`, `slice()`, `distinct()`, joins where `dataset` is `x`, ...)
#'   restore data-frame attributes, so the keys and the `admiral_df` class
#'   normally survive a derivation pipeline. They are *not* preserved by
#'   operations which build a new data frame, in particular:
#'
#'   * `summarise()`, `group_by()`/`ungroup()` and `rowwise()`,
#'   * `tidyr::pivot_longer()`/`pivot_wider()`, `nest()`/`unnest()`,
#'     `complete()`,
#'   * `bind_rows()` when `dataset` is not the *first* argument,
#'   * joins where `dataset` is the `y` (rather than the `x`) argument.
#'
#'   Call this function afterwards (or re-apply it as needed) if the dataset
#'   passes through any of those.
#'
#'   Because the keys are stored as plain variable names, they are not updated
#'   by `rename()` and not removed by `select()`. [summary.admiral_df()] warns
#'   when keys are no longer present in the dataset.
#'
#' @keywords utils_help
#' @family utils_help
#'
#' @export
set_admiral_keys <- function(dataset, metacore, dataset_name = NULL) {
  assert_data_frame(dataset)
  keys <- get_admiral_keys(metacore, dataset_name)
  attr(dataset, "admiral_keys") <- keys
  as_admiral_df(dataset)
}

#' Find the Minimal Set of Variables that Uniquely Identifies Rows
#'
#' Starting from `must_have`, adds variables from `optional` one at a time (in
#' the given order) until the combination is a unique key of `dataset`, and
#' returns the smallest such combination. If uniqueness is never reached, the
#' full set (`must_have` plus all of `optional`) is returned -- the caller can
#' detect this because the returned key still yields duplicates.
#'
#' This is used by [infer_admiral_keys()] to discover the record structure from
#' the data rather than assuming a fixed key. Only *semantic* key variables
#' should be passed as `optional`; surrogate/sequence keys (e.g. `ASEQ`) must be
#' excluded, otherwise uniqueness is reached trivially and genuine duplicates are
#' masked.
#'
#' @param dataset A data frame
#' @param must_have Character vector of variables always kept in the key
#' @param optional Character vector of candidate variables to add, in priority
#'   order
#'
#' @return A character vector of key variable names.
#'
#' @keywords internal
#' @family internal
minimal_unique_key <- function(dataset, must_have, optional) {
  is_unique <- function(key) {
    length(key) > 0 && nrow(dataset) == nrow(distinct(dataset, !!!syms(key)))
  }
  key <- must_have
  if (is_unique(key)) {
    return(key)
  }
  for (v in setdiff(optional, key)) {
    key <- c(key, v)
    if (is_unique(key)) {
      return(key)
    }
  }
  key
}

#' Infer the Record Structure of a Dataset from its Type and Data
#'
#' Works out which variables a dataset appears to have one record per. This is
#' the fallback used by [summary.admiral_df()] when no spec-defined or
#' derivation keys are available (i.e. the dataset was not tagged via
#' [set_admiral_keys()] or by a `derive_*` function). It starts from the
#' semantic core of the detected ADaM dataset type (e.g. `USUBJID` + `PARAMCD`
#' for BDS) and adds standard ADaM key variables (analysis visit, timepoint,
#' date, period, derivation type) only as far as needed to make the key unique
#' (see [minimal_unique_key()]).
#'
#' For `BDS`/`TTE`, surrogate/sequence keys such as `ASEQ` are deliberately
#' *not* used, so that unintended duplicate records remain detectable. For
#' `OCCDS`, however, `ASEQ` *is* the intended record key (there is no
#' analysis-value structure to fall back on), so it is required: if `ASEQ` is
#' absent the record structure cannot be checked and a warning is issued. This
#' is still a heuristic; spec-defined or `by_vars`-derived keys should be
#' preferred.
#'
#' @param type The dataset type, see [get_admiral_df_type()]
#' @param subject_keys The subject-level key variables present in the dataset
#' @param cols The column names of the dataset
#' @param dataset The dataset itself (used to test uniqueness)
#'
#' @return A character vector of inferred key variable names, or a zero-length
#'   vector when no plausible record structure can be determined.
#'
#' @keywords internal
#' @family internal
infer_admiral_keys <- function(type, subject_keys, cols, dataset) {
  # OCCDS has no analysis-value structure to fall back on; `ASEQ` is its genuine
  # record key, so (unlike BDS/TTE) it is required rather than a surrogate
  if (type == "OCCDS") {
    if (!"ASEQ" %in% cols) {
      cli_warn(c(
        "Cannot check the record structure of an {.val OCCDS} dataset without
         {.var ASEQ}.",
        i = "Attach the intended keys with {.fun set_admiral_keys} to check
             that the dataset has one record per key."
      ))
      return(character(0))
    }
    return(intersect(c(subject_keys, "ASEQ"), cols))
  }

  core <- switch(type,
    ADSL = subject_keys,
    BDS = c(subject_keys, "PARAMCD"),
    TTE = c(subject_keys, "PARAMCD"),
    return(character(0))
  )
  core <- intersect(core, cols)
  if (length(core) == 0) {
    return(character(0))
  }

  # standard "within-core" key variables, in ADaM precedence order.
  # NOTE: no surrogate/sequence keys (ASEQ, SRCSEQ, ...) -- see minimal_unique_key()
  extra <- c(
    "AVISITN", "AVISIT", "ATPTN", "ATPT",
    "ADTM", "ADT", "APERIOD", "APERIODC", "ASPID"
  )
  minimal_unique_key(
    dataset,
    must_have = core,
    optional = intersect(extra, cols)
  )
}

#' Check the `admiral_df` Tag of a Dataset
#'
#' Reports what is left of the `admiral_df` tag on a dataset: whether the class
#' is still present, which key variables are declared (or would be inferred),
#' and whether those keys are still in the dataset. Use it before relying on
#' [summary.admiral_df()], or after an operation which may have dropped the tag.
#'
#' @param dataset A data frame
#'
#' @details
#'   The `admiral_df` class and the `"admiral_keys"` attribute are preserved
#'   independently, so a dataset can carry one without the other and each case
#'   behaves differently:
#'
#'   * **both present** -- the record structure is checked against the declared
#'     keys.
#'   * **class only** -- [summary.admiral_df()] falls back to
#'     [infer_admiral_keys()], so the record structure is still checked, but
#'     against a guessed key (reported as `key_source = "inferred"`).
#'   * **class only, but `"admiral_keys"` set to a zero-length vector** -- the
#'     inference fallback is skipped *and* no key is left to check, so the
#'     record structure is not checked at all. `set_admiral_keys()` produces
#'     this when the specification defines no key variables.
#'   * **keys only** -- `summary()` no longer dispatches to
#'     [summary.admiral_df()] and the diagnostic is silently lost. Re-tagging
#'     with `as_admiral_df(dataset)` restores the class but keeps the existing
#'     (possibly stale) attribute; pass `keys` explicitly to replace it.
#'
#'   Because keys are stored as plain variable names, `dplyr::rename()` and
#'   `dplyr::select()` can leave them pointing at variables which are no longer
#'   there. Those are reported in `stale_keys`, and the record structure is then
#'   checked against fewer variables than declared (or not at all).
#'
#' @return
#'   An object of class `admiral_df_check` (a named list) which is printed by
#'   [print.admiral_df_check()], with elements:
#'
#'   * `is_admiral_df` -- whether the `admiral_df` class is present.
#'   * `keys` -- the key variables which would be used by
#'     [summary.admiral_df()], i.e. the declared keys still present in the
#'     dataset, or the inferred ones. May be zero-length.
#'   * `key_source` -- `"declared"`, `"inferred"`, or `"none"`.
#'   * `declared_keys` -- the `"admiral_keys"` attribute as stored, or `NULL`.
#'   * `stale_keys` -- declared keys which are not in the dataset.
#'   * `structure_checked` -- whether [summary.admiral_df()] would check the
#'     record structure, i.e. report on one record per key.
#'
#' @seealso [as_admiral_df()], [set_admiral_keys()], [summary.admiral_df()]
#'
#' @keywords utils_print
#' @family utils_print
#'
#' @export
#'
#' @examplesx
#'
#' @caption Checking the tag before relying on it
#'
#' @code
#' library(tibble)
#'
#' advs <- as_admiral_df(
#'   tribble(
#'     ~USUBJID, ~PARAMCD, ~AVAL,
#'     "1",      "DIABP",     51,
#'     "1",      "SYSBP",    121
#'   ),
#'   keys = c("USUBJID", "PARAMCD")
#' )
#'
#' check_admiral_df(advs)
#'
#' @caption Keys which no longer match the dataset
#'
#' @info The keys are stored as plain variable names, so they are not updated by
#'   `dplyr::rename()`.
#'
#' @code
#' check_admiral_df(dplyr::rename(advs, SUBJID = USUBJID))
#'
#' @caption Using the result programmatically
#'
#' @code
#' if (!check_admiral_df(advs)$structure_checked) {
#'   print("re-tag before summarizing")
#' }
check_admiral_df <- function(dataset) {
  assert_data_frame(dataset)

  cols <- colnames(dataset)
  declared_keys <- attr(dataset, "admiral_keys")

  # mirror the key resolution of `summary.admiral_df()`: an attribute which is
  # set (even to a zero-length vector) suppresses the inference fallback, so
  # `NULL` and `character(0)` are NOT interchangeable here
  if (is.null(declared_keys)) {
    subject_keys <- intersect(vars2chr(get_admiral_option("subject_keys")), cols)
    if ("USUBJID" %in% cols) {
      subject_keys <- "USUBJID"
    }
    # a pure reporter must not emit the warning `infer_admiral_keys()` raises
    # for an OCCDS dataset without `ASEQ`; the empty result is reported instead
    keys <- suppressWarnings(infer_admiral_keys(
      get_admiral_df_type(dataset), subject_keys, cols, dataset
    ))
    key_source <- if (length(keys) > 0) "inferred" else "none"
    stale_keys <- character(0)
  } else {
    stale_keys <- setdiff(declared_keys, cols)
    keys <- intersect(declared_keys, cols)
    key_source <- if (length(keys) > 0) "declared" else "none"
  }

  out <- list(
    is_admiral_df = inherits(dataset, "admiral_df"),
    keys = keys,
    key_source = key_source,
    declared_keys = declared_keys,
    stale_keys = stale_keys,
    structure_checked = inherits(dataset, "admiral_df") && length(keys) > 0
  )
  class(out) <- c("admiral_df_check", "list")
  out
}

#' Print an `admiral_df` Tag Check
#'
#' @param x An `admiral_df_check` object created by [check_admiral_df()]
#' @param ... Not used
#'
#' @return No return value, called for side effects. The `admiral_df_check`
#'   object is returned invisibly.
#'
#' @keywords utils_print
#' @family utils_print
#'
#' @export
#'
#' @examplesx
#'
#' @caption Printing a tag check
#'
#' @code
#' library(tibble)
#'
#' advs <- as_admiral_df(
#'   tribble(
#'     ~USUBJID, ~PARAMCD, ~AVAL,
#'     "1",      "DIABP",     51
#'   ),
#'   keys = c("USUBJID", "PARAMCD")
#' )
#'
#' print(check_admiral_df(advs))
print.admiral_df_check <- function(x, ...) {
  declared <- if (is.null(x$declared_keys)) {
    if (x$key_source == "inferred") {
      paste0("none -- inferred: ", paste(x$keys, collapse = ", "))
    } else {
      "none -- and none could be inferred"
    }
  } else if (length(x$declared_keys) == 0) {
    "none -- attribute set but empty, so nothing is inferred"
  } else {
    paste(x$declared_keys, collapse = ", ")
  }

  present <- if (length(x$declared_keys) == 0) {
    "n/a"
  } else if (length(x$stale_keys) == 0) {
    cli::col_green("all present")
  } else if (length(x$keys) > 0) {
    cli::col_red(paste0(
      "stale: ", paste(x$stale_keys, collapse = ", "),
      " missing; checked on ", paste(x$keys, collapse = ", "), " only"
    ))
  } else {
    cli::col_red("stale: no declared key left in the dataset")
  }

  # `cli_text()` collapses runs of spaces, which would unalign the padded
  # labels; `cli_verbatim()` passes the line through as it is (and does not
  # interpret inline markup, hence the explicit `cli::col_*()` above)
  field <- function(label, value) {
    cli_verbatim(paste0("  ", sprintf("%-20s", label), value))
  }

  cli_text("<admiral_df> tag check")
  field("is an admiral_df:", if (x$is_admiral_df) {
    cli::col_green("TRUE")
  } else {
    cli::col_red("FALSE")
  })
  field("key variables:", declared)
  field("keys present:", present)
  field("structure checked:", if (x$structure_checked) {
    cli::col_green(paste0(
      "yes -- one record per ", paste(x$keys, collapse = ", "),
      " (", x$key_source, ")"
    ))
  } else {
    cli::col_red("no")
  })

  # the class and the attribute travel independently; an orphaned attribute is
  # the case which fails silently, so it is worth spelling out
  if (!x$is_admiral_df && !is.null(x$declared_keys)) {
    cli_bullets(c("!" = "The {.var admiral_keys} attribute is orphaned:
       {.fun summary} dispatches to {.fun summary.data.frame}, so the
       diagnostic is lost. Re-tag with {.fun as_admiral_df}, passing
       {.arg keys} explicitly to replace the existing attribute."))
  }

  invisible(x)
}

#' Summarize an `admiral_df` Dataset
#'
#' Provides a quick diagnostic of an admiral dataset, such as the number of
#' subjects and observations, the list of parameters (`PARAMCD`/`PARAM`), and
#' the list of analysis visits (`AVISIT`). It also reports the record structure,
#' i.e. whether the dataset has exactly one record per key variable combination.
#' The information shown is tailored to the ADaM type of the dataset (`ADSL`,
#' `BDS`, `OCCDS`, `TTE`, or `other`), as determined by [get_admiral_df_type()].
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

  # record structure check: prefer keys attached by the spec
  # (`set_admiral_keys()`) or by the `derive_*` function that produced the
  # dataset (its `by_vars`); otherwise fall back to keys inferred from the
  # dataset type and validated against the data
  keys <- attr(object, "admiral_keys")
  key_source <- "declared"
  if (is.null(keys)) {
    keys <- infer_admiral_keys(type, subject_keys, cols, object)
    key_source <- "inferred"
  }
  # declared keys are stored as plain names, so they can go stale if the dataset
  # was renamed or the key variables were dropped; the record structure check
  # would then be silently weakened (or skipped), which must not pass unnoticed
  stale_keys <- setdiff(keys, cols)
  if (key_source == "declared" && length(stale_keys) > 0) {
    cli_warn(c(
      "The declared key variable{?s} {.var {stale_keys}} {?is/are} not in the
       dataset, so the record structure was not checked as declared.",
      i = "{.var admiral_keys} is not updated by {.fun dplyr::rename} or
           {.fun dplyr::select}.",
      i = if (length(stale_keys) < length(keys)) {
        "The dataset is only checked for one record per remaining key
         variable{?s} {.var {intersect(keys, cols)}}."
      } else {
        "No key variable is left, the record structure is not checked."
      }
    ))
  }

  keys <- intersect(keys, cols)
  if (length(keys) > 0) {
    n_unique <- nrow(distinct(object, !!!syms(keys)))
    out$keys <- keys
    out$key_source <- key_source
    out$n_duplicate_keys <- nrow(object) - n_unique
  }

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

  if (!is.null(x$keys)) {
    keys <- paste(x$keys, collapse = ", ")
    src <- if (x$key_source == "inferred") " (inferred)" else ""
    if (x$n_duplicate_keys == 0) {
      cli_text(cli::col_green(
        "Structure{src}: one record per {keys} {cli::symbol$tick}"
      ))
    } else {
      cli_text(cli::col_red(
        "Structure{src}: expected one record per {keys}, found
         {x$n_duplicate_keys} extra record{?s} {cli::symbol$cross}"
      ))
    }
  }

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
