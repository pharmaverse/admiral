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
#' if (requireNamespace("metacore", quietly = TRUE)) {
#'   load(metacore::metacore_example("pilot_ADaM.rda"))
#'   get_admiral_keys(metacore, "ADSL")
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
#' @return `dataset` with an `"admiral_keys"` attribute, an `"admiral_ds_name"`
#'   attribute (shown in the heading of [summary.admiral_df()]), and the
#'   `admiral_df` class.
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
  # `get_admiral_keys()` has validated that a `NULL` dataset_name means the
  # specification describes exactly one dataset
  if (is.null(dataset_name)) {
    dataset_name <- unique(metacore$ds_vars$dataset)
  }
  attr(dataset, "admiral_keys") <- keys
  attr(dataset, "admiral_ds_name") <- dataset_name
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

  cli_rule(left = "admiral_df tag check")
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

#' Summarize the ADSL-Specific Content of a Dataset
#'
#' Collects the subject-level information [summary.admiral_df()] reports for a
#' dataset of type `"ADSL"`: the treatment and population breakdowns a
#' programmer checks first, plus a small number of subject-level consistency
#' checks.
#'
#' @param dataset An ADSL-like data frame
#' @param cols The column names of `dataset`
#'
#' @details
#'   Every item is optional: an element is only returned when the variables it
#'   needs are in the dataset, so a partially built ADSL is summarized as far as
#'   it goes. A variable which is present but entirely `NA` is treated as not
#'   yet derived, i.e. the same as absent. The checks (`arm_mismatch`,
#'   `saffl_no_trtsdt`) are diagnostics rather than validations -- they are
#'   reported, never signalled.
#'
#' @return
#'   A named list with any of the following elements, to be spliced into the
#'   result of [summary.admiral_df()]:
#'
#'   * `trt_var`, `trt` -- the treatment variable used and the subject count per
#'     treatment.
#'   * `pops` -- named count of `"Y"` per population flag.
#'   * `funnel` -- named subject counts for the randomized / treated / completed
#'     stages.
#'   * `trtdur` -- median and range of `TRTDURD`.
#'   * `eosstt` -- subject count per end-of-study status.
#'   * `n_dth`, `dthcaus` -- death count, and the cause breakdown among the
#'     deaths.
#'   * `age`, `sex` -- age range and sex breakdown.
#'   * `arm_mismatch` -- subjects whose `ARM` and `ACTARM` differ.
#'   * `saffl_no_trtsdt` -- subjects in the safety population with no treatment
#'     start date.
#'   * `trt_end_before_start` -- subjects with `TRTEDT` before `TRTSDT`.
#'   * `flag_domain`, `flag_domain_vars` -- count of population flag values
#'     outside `Y`/`N`/`NA`, and the flags containing them.
#'
#' @keywords internal
#' @family internal
summarize_adsl <- function(dataset, cols) {
  out <- list()

  # planned treatment is what most outputs are summarized by; fall back to the
  # actual treatment and then to `ARM` so a pre-`TRT01P` ADSL still reports; a
  # candidate which is present but entirely `NA` has not been derived yet, so
  # it is skipped in favour of the next one
  trt_var <- intersect(c("TRT01P", "TRT01A", "ARM"), cols)
  trt_var <- trt_var[vapply(
    trt_var,
    function(v) any(!is.na(dataset[[v]])),
    logical(1)
  )]
  if (length(trt_var) > 0) {
    out$trt_var <- trt_var[1]
    out$trt <- count(dataset, !!sym(trt_var[1]), name = "n")
  }

  # population flags drive the denominators of nearly every downstream table
  pop_vars <- intersect(
    c("RANDFL", "SAFFL", "ITTFL", "FASFL", "EFFFL", "PPROTFL", "COMPLFL"),
    cols
  )
  if (length(pop_vars) > 0) {
    out$pops <- vapply(
      pop_vars,
      function(v) sum(dataset[[v]] == "Y", na.rm = TRUE),
      integer(1)
    )
  }

  # the subject funnel: how many made it to randomization, treatment, and study
  # completion -- the shape of the study, and it makes subjects who fell out
  # before a stage (e.g. screen failures) visible as the difference to the
  # subject count
  funnel <- c(
    if ("RANDDT" %in% cols) {
      c(randomized = sum(!is.na(dataset$RANDDT)))
    } else if ("RANDFL" %in% cols) {
      c(randomized = sum(dataset$RANDFL == "Y", na.rm = TRUE))
    },
    if ("TRTSDT" %in% cols) {
      c(treated = sum(!is.na(dataset$TRTSDT)))
    } else if ("TRTSDTM" %in% cols) {
      c(treated = sum(!is.na(dataset$TRTSDTM)))
    },
    if ("EOSSTT" %in% cols) {
      c(completed = sum(dataset$EOSSTT == "COMPLETED", na.rm = TRUE))
    }
  )
  if (length(funnel) > 0) {
    out$funnel <- funnel
  }

  if ("TRTDURD" %in% cols) {
    dur <- dataset$TRTDURD[!is.na(dataset$TRTDURD)]
    if (length(dur) > 0) {
      out$trtdur <- c(
        median = median(dur), min = min(dur), max = max(dur)
      )
    }
  }

  if ("EOSSTT" %in% cols && any(!is.na(dataset$EOSSTT))) {
    out$eosstt <- count(dataset, EOSSTT, name = "n")
  }

  # deaths are the first safety question; the cause breakdown is restricted to
  # the deaths so an alive majority does not print as `NA: n`
  if ("DTHFL" %in% cols && any(!is.na(dataset$DTHFL))) {
    out$n_dth <- sum(dataset$DTHFL == "Y", na.rm = TRUE)
    if ("DTHCAUS" %in% cols && out$n_dth > 0) {
      out$dthcaus <- count(filter(dataset, DTHFL == "Y"), DTHCAUS, name = "n")
    }
  }

  if ("AGE" %in% cols && any(!is.na(dataset$AGE))) {
    out$age <- c(
      median = median(dataset$AGE, na.rm = TRUE),
      min = min(dataset$AGE, na.rm = TRUE),
      max = max(dataset$AGE, na.rm = TRUE)
    )
  }
  if ("SEX" %in% cols && any(!is.na(dataset$SEX))) {
    out$sex <- count(dataset, SEX, name = "n")
  }

  # planned and actual treatment differing means a subject was mis-dosed, which
  # changes which population they belong to -- always worth surfacing
  if (all(c("ARM", "ACTARM") %in% cols)) {
    out$arm_mismatch <- sum(
      !is.na(dataset$ARM) & !is.na(dataset$ACTARM) &
        dataset$ARM != dataset$ACTARM
    )
  }

  # a subject in the safety population with no treatment start date cannot have
  # treatment-emergent records derived for them downstream
  if (all(c("SAFFL", "TRTSDT") %in% cols)) {
    out$saffl_no_trtsdt <- sum(
      dataset$SAFFL == "Y" & is.na(dataset$TRTSDT),
      na.rm = TRUE
    )
  }

  # treatment ending before it started is always a derivation defect
  if (all(c("TRTSDT", "TRTEDT") %in% cols)) {
    out$trt_end_before_start <- sum(
      !is.na(dataset$TRTSDT) & !is.na(dataset$TRTEDT) &
        dataset$TRTEDT < dataset$TRTSDT
    )
  }

  # values like "Yes", "y", or "" are silently dropped by `== "Y"` subsetting,
  # so a flag outside the Y/N/NA domain poisons every downstream denominator
  if (length(pop_vars) > 0) {
    bad_flags <- vapply(
      pop_vars,
      function(v) sum(!dataset[[v]] %in% c("Y", "N", NA)),
      integer(1)
    )
    out$flag_domain <- sum(bad_flags)
    out$flag_domain_vars <- names(bad_flags)[bad_flags > 0]
  }

  out
}

#' Summarize the BDS-Specific Content of a Dataset
#'
#' Collects the parameter-level information [summary.admiral_df()] reports for a
#' dataset of type `"BDS"`: a per-parameter table (most BDS defects are visible
#' in a per-parameter row), the derived records by `DTYPE`, and a set of
#' record-level consistency checks.
#'
#' @param dataset A BDS-like data frame
#' @param cols The column names of `dataset`
#'
#' @details
#'   As for [summarize_adsl()], every item is optional: an element is only
#'   returned when the variables it needs are in the dataset, and the checks are
#'   diagnostics rather than validations -- reported, never signalled.
#'
#' @return
#'   A named list with any of the following elements:
#'
#'   * `params` -- per-parameter table: records, subjects, visits, missing and
#'     min / median / max of `AVAL`, depending on the variables present.
#'   * `dtype` -- count of derived records per `DTYPE`.
#'   * `multiple_baselines` -- subject/parameter (and `BASETYPE`/`ATPT` where
#'     present) combinations with more than one `ABLFL == "Y"` record.
#'   * `chg_no_base` -- records with `CHG`/`PCHG` but no `BASE`.
#'   * `avisit_mismatch` -- `AVISIT` values mapping to more than one `AVISITN`,
#'     plus the reverse.
#'   * `param_inconsistent` -- `PARAMCD` values with more than one distinct
#'     `PARAM` label or `AVALU` unit.
#'
#' @keywords internal
#' @family internal
summarize_bds <- function(dataset, cols) {
  out <- list()

  # the per-parameter table replaces the flat PARAMCD list: records, subjects,
  # visits, and the AVAL distribution per parameter is where most BDS defects
  # (wrong shape, empty parameter, absurd values) become visible
  if ("PARAMCD" %in% cols) {
    safe_min <- function(v) if (all(is.na(v))) NA_real_ else min(v, na.rm = TRUE)
    safe_max <- function(v) if (all(is.na(v))) NA_real_ else max(v, na.rm = TRUE)
    stats <- list(records = quote(n()))
    if ("USUBJID" %in% cols) {
      stats$subjects <- quote(n_distinct(USUBJID))
    }
    if ("AVISIT" %in% cols) {
      stats$visits <- quote(n_distinct(AVISIT))
    }
    if ("AVAL" %in% cols) {
      stats$missing <- quote(sum(is.na(AVAL)))
      stats$min <- quote(safe_min(AVAL))
      stats$median <- quote(median(AVAL, na.rm = TRUE))
      stats$max <- quote(safe_max(AVAL))
    }
    out$params <- dataset %>%
      group_by(PARAMCD) %>%
      summarise(!!!stats, .groups = "drop") %>%
      arrange(PARAMCD)
  }

  # what the derivations added, by type (LOCF, AVERAGE, MAXIMUM, ...)
  if ("DTYPE" %in% cols && any(!is.na(dataset$DTYPE))) {
    out$dtype <- count(filter(dataset, !is.na(DTYPE)), DTYPE, name = "n")
  }

  # more than one baseline per subject and parameter silently corrupts BASE,
  # CHG, and PCHG downstream, and the record structure check does not see it
  if (all(c("USUBJID", "PARAMCD", "ABLFL") %in% cols)) {
    base_by <- intersect(c("USUBJID", "PARAMCD", "BASETYPE", "ATPT"), cols)
    out$multiple_baselines <- dataset %>%
      filter(ABLFL == "Y") %>%
      count(!!!syms(base_by), name = "n") %>%
      filter(n > 1) %>%
      nrow()
  }

  # a change from baseline without a baseline means CHG was derived before (or
  # despite) BASE -- the values cannot be trusted
  chg_vars <- intersect(c("CHG", "PCHG"), cols)
  if (length(chg_vars) > 0 && "BASE" %in% cols) {
    out$chg_no_base <- sum(
      rowSums(!is.na(dataset[chg_vars])) > 0 & is.na(dataset$BASE)
    )
  }

  # AVISIT and AVISITN must map one-to-one; a mismatch (including one side
  # missing) splits or merges visits in every by-visit output
  if (all(c("AVISIT", "AVISITN") %in% cols)) {
    pairs <- distinct(dataset, AVISIT, AVISITN)
    out$avisit_mismatch <-
      nrow(filter(count(filter(pairs, !is.na(AVISIT)), AVISIT), n > 1)) +
      nrow(filter(count(filter(pairs, !is.na(AVISITN)), AVISITN), n > 1))
  }

  # PARAMCD is the key; two PARAM labels or two units under one code usually
  # means two different measurements were merged into one parameter
  incons_vars <- intersect(c("PARAM", "AVALU"), cols)
  if ("PARAMCD" %in% cols && length(incons_vars) > 0) {
    out$param_inconsistent <- sum(vapply(
      incons_vars,
      function(v) {
        dataset %>%
          filter(!is.na(PARAMCD), !is.na(!!sym(v))) %>%
          distinct(PARAMCD, !!sym(v)) %>%
          count(PARAMCD, name = "n") %>%
          filter(n > 1) %>%
          nrow()
      },
      integer(1)
    ))
  }

  out
}

#' Summarize the OCCDS-Specific Content of a Dataset
#'
#' Collects the occurrence-level information [summary.admiral_df()] reports for
#' a dataset of type `"OCCDS"`: the distinct terms at each coding level, the
#' treatment-emergent and severity/seriousness breakdowns, and the occurrence
#' consistency checks -- in particular the occurrence flag integrity check,
#' which the record structure check cannot catch (an `AOCCFL` duplicated per
#' subject does not change the record count).
#'
#' @param dataset An OCCDS-like data frame
#' @param cols The column names of `dataset`
#'
#' @details
#'   As for [summarize_adsl()], every item is optional: an element is only
#'   returned when the variables it needs are in the dataset, and the checks are
#'   diagnostics rather than validations -- reported, never signalled.
#'
#'   The occurrence flag check derives its grouping from the flag name: a flag
#'   whose infix (the part between `AOCC` and `FL`) contains `"S"` (e.g.
#'   `AOCCSFL`) is checked per subject and body system (`--BODSYS`), one
#'   containing `"P"` (e.g. `AOCCPFL`) per subject and dictionary term
#'   (`--DECOD`), and any other (`AOCCFL`, `AOCC02FL`, `AOCCIFL`, ...) per
#'   subject. Only *more than one* `"Y"` per group is a defect; zero is
#'   legitimate, because occurrence flags are typically restricted to a subset
#'   such as the treatment-emergent records. A flag whose level variable is
#'   absent from the dataset is skipped.
#'
#' @return
#'   A named list with any of the following elements:
#'
#'   * `terms` -- named count of distinct values per coding level (`--BODSYS`,
#'     `--DECOD`, `--TERM`, whichever are present).
#'   * `trtem` -- treatment-emergent (`TRTEMFL == "Y"`) and total record
#'     counts.
#'   * `sev_var`, `sev` -- the severity/toxicity-grade variable used and the
#'     record count per value.
#'   * `n_serious` -- records with `AESER == "Y"`.
#'   * `occ_flag_dups`, `occ_flag_vars` -- occurrence flag groups with more
#'     than one `"Y"`, and the flags containing them.
#'   * `missing_astdt` -- records with no `ASTDT`.
#'   * `pre_trt_emergent` -- treatment-emergent records starting before
#'     `TRTSDT`.
#'
#' @keywords internal
#' @family internal
summarize_occds <- function(dataset, cols) {
  out <- list()

  # the distinct value count per coding level (body system, dictionary term,
  # reported term) shows the coding granularity at a glance; only one --DECOD
  # level was reported before
  term_vars <- c(
    str_subset(cols, "BODSYS$")[1],
    str_subset(cols, "DECOD$")[1],
    str_subset(cols, "TERM$")[1]
  )
  term_vars <- term_vars[!is.na(term_vars)]
  if (length(term_vars) > 0) {
    out$terms <- vapply(
      setNames(term_vars, term_vars),
      function(v) n_distinct(dataset[[v]], na.rm = TRUE),
      integer(1)
    )
  }

  # treatment-emergent vs all records: most occurrence outputs are restricted
  # to TRTEMFL == "Y", so this is the denominator shift a programmer checks
  if ("TRTEMFL" %in% cols && any(!is.na(dataset$TRTEMFL))) {
    out$trtem <- c(
      emergent = sum(dataset$TRTEMFL == "Y", na.rm = TRUE),
      total = nrow(dataset)
    )
  }

  # severity/grade distribution: prefer the analysis version over the SDTM
  # copy, and severity over toxicity grade; an all-NA candidate has not been
  # derived yet and is skipped in favour of the next one
  sev_var <- intersect(c("ASEV", "AESEV", "ATOXGR", "AETOXGR"), cols)
  sev_var <- sev_var[vapply(
    sev_var,
    function(v) any(!is.na(dataset[[v]])),
    logical(1)
  )]
  if (length(sev_var) > 0) {
    out$sev_var <- sev_var[1]
    out$sev <- count(
      filter(dataset, !is.na(!!sym(sev_var[1]))),
      !!sym(sev_var[1]),
      name = "n"
    )
  }

  if ("AESER" %in% cols && any(!is.na(dataset$AESER))) {
    out$n_serious <- sum(dataset$AESER == "Y", na.rm = TRUE)
  }

  # occurrence flag integrity: more than one "Y" per subject (and level)
  # silently duplicates subjects in every "n (%) of subjects" count downstream,
  # and the record structure check does not see it
  occ_flags <- str_subset(cols, "^AOCC.*FL$")
  if (length(occ_flags) > 0 && "USUBJID" %in% cols) {
    bodsys_var <- str_subset(cols, "BODSYS$")[1]
    decod_var <- str_subset(cols, "DECOD$")[1]
    count_dups <- function(flag) {
      # the flag infix encodes the level it is unique per: "S" = body system,
      # "P" = dictionary term, anything else (including numbered sponsor
      # flags) = subject
      infix <- str_remove(str_remove(flag, "^AOCC"), "FL$")
      group <- "USUBJID"
      if (str_detect(infix, "S")) {
        if (is.na(bodsys_var)) {
          return(NA_integer_)
        }
        group <- c(group, bodsys_var)
      } else if (str_detect(infix, "P")) {
        if (is.na(decod_var)) {
          return(NA_integer_)
        }
        group <- c(group, decod_var)
      }
      dataset %>%
        filter(!!sym(flag) == "Y") %>%
        count(!!!syms(group), name = "n") %>%
        filter(n > 1) %>%
        nrow()
    }
    flag_dups <- vapply(occ_flags, count_dups, integer(1))
    checked <- !is.na(flag_dups)
    if (any(checked)) {
      out$occ_flag_dups <- sum(flag_dups[checked])
      out$occ_flag_vars <- names(flag_dups)[checked & flag_dups > 0]
    }
  }

  # an occurrence without a start date cannot be classified as treatment
  # emergent, so it silently drops out of every TRTEMFL-restricted output
  if ("ASTDT" %in% cols) {
    out$missing_astdt <- sum(is.na(dataset$ASTDT))
  }

  # a treatment-emergent record starting before treatment contradicts its own
  # flag -- either the flag or the date imputation is wrong
  if (all(c("TRTEMFL", "ASTDT", "TRTSDT") %in% cols)) {
    out$pre_trt_emergent <- sum(
      dataset$TRTEMFL == "Y" & dataset$ASTDT < dataset$TRTSDT,
      na.rm = TRUE
    )
  }

  out
}

#' Summarize a Dataset Against its ADSL
#'
#' Collects the cross-dataset information [summary.admiral_df()] reports when a
#' subject-level dataset is supplied via its `adsl` argument. ADSL is the
#' subject-level source of truth -- the full subject universe, the
#' authoritative treatment/population values, and the subject-level events
#' (death) a record-level dataset must respect -- so a dataset can be wrong
#' *relative to ADSL* in ways it can never reveal on its own: orphan subjects,
#' stale merged variables, records after death, and it can never supply an
#' "n (%) of subjects" denominator by itself.
#'
#' @param dataset A data frame with a `USUBJID` variable
#' @param cols The column names of `dataset`
#' @param adsl A subject-level data frame with one record per `USUBJID`
#' @param type The dataset type, see [get_admiral_df_type()]
#'
#' @details
#'   As for [summarize_adsl()], every item is optional: an element is only
#'   returned when the variables it needs are present (on both sides, where the
#'   comparison needs both), and the checks are diagnostics rather than
#'   validations -- reported, never signalled.
#'
#'   Subjects are matched by `USUBJID` alone (other subject keys may be `NA` on
#'   derivation-added records). Shared subject-level variables are compared as
#'   character, per record, with `NA` equal to `NA`; a subject not in ADSL is
#'   reported as an orphan, not additionally as a mismatch. Mismatches are
#'   reported as "disagrees with ADSL" without presuming which side is wrong:
#'   the cause may be an outdated dataset *or* a newer ADSL.
#'
#' @return
#'   A named list with any of the following elements:
#'
#'   * `n_adsl`, `n_common` -- subjects in `adsl`, and subjects present in both
#'     datasets.
#'   * `saffl` -- `covered` / `total` subject counts within the safety
#'     population (`SAFFL == "Y"` in `adsl`).
#'   * `arm_var`, `arm_coverage` -- the treatment variable used (first
#'     non-empty of `TRT01A`, `TRT01P`, `ARM` in `adsl`) and, per arm, the
#'     subjects with records vs the arm total.
#'   * `n_orphans`, `orphans` -- subjects in `dataset` but not in `adsl`; the
#'     `USUBJID` values are kept because orphans caused by formatting drift are
#'     only diagnosable from examples.
#'   * `stale_total`, `stale_vars` -- subject/variable combinations where a
#'     shared subject-level variable (`TRTSDT`, `TRTEDT`, `TRT01P`, `TRT01A`,
#'     `SAFFL`, `ITTFL`, `DTHDT`, `EOSDT`, `AGE`, `SEX`, `RACE`) disagrees with
#'     ADSL, and the per-variable breakdown.
#'   * `after_death_var`, `n_after_death` -- the analysis date variable used
#'     (`ADT` or `ASTDT`) and the records dated after the subject's `DTHDT` in
#'     `adsl`.
#'   * `incidence`, `incidence_denom` -- (OCCDS) subjects with at least one
#'     `TRTEMFL == "Y"` record over the safety-population (or, failing that,
#'     ADSL) denominator.
#'   * `n_emergent_untreated` -- (OCCDS) subjects with treatment-emergent
#'     records whose `TRTSDT` in `adsl` is missing.
#'
#' @keywords internal
#' @family internal
summarize_vs_adsl <- function(dataset, cols, adsl, type) {
  out <- list()
  adsl_cols <- colnames(adsl)

  ds_subj <- unique(dataset$USUBJID)
  ds_subj <- ds_subj[!is.na(ds_subj)]
  adsl_subj <- adsl$USUBJID
  # row index of each record's subject in `adsl`; NA for orphans
  adsl_idx <- match(dataset$USUBJID, adsl_subj)

  out$n_adsl <- length(adsl_subj)
  out$n_common <- length(intersect(ds_subj, adsl_subj))

  # subjects with records but no ADSL entry: wrong ADSL version, screen
  # failures or test subjects left in, or USUBJID formatting drift -- the
  # values are kept because formatting drift is only diagnosable from examples
  orphans <- setdiff(ds_subj, adsl_subj)
  out$n_orphans <- length(orphans)
  out$orphans <- orphans

  # the safety population is the denominator most record-level outputs use
  saf_subj <- NULL
  if ("SAFFL" %in% adsl_cols && any(!is.na(adsl$SAFFL))) {
    saf_subj <- adsl_subj[!is.na(adsl$SAFFL) & adsl$SAFFL == "Y"]
    out$saffl <- c(
      covered = length(intersect(ds_subj, saf_subj)),
      total = length(saf_subj)
    )
  }

  # subjects with records per arm vs the arm total: an empty or thin arm is
  # the shape of every downstream incidence table; prefer actual treatment,
  # and skip a candidate which is present but entirely NA
  arm_var <- intersect(c("TRT01A", "TRT01P", "ARM"), adsl_cols)
  arm_var <- arm_var[vapply(
    arm_var,
    function(v) any(!is.na(adsl[[v]])),
    logical(1)
  )]
  if (length(arm_var) > 0) {
    out$arm_var <- arm_var[1]
    out$arm_coverage <- adsl %>%
      group_by(!!sym(arm_var[1])) %>%
      summarise(
        subjects = sum(USUBJID %in% ds_subj),
        total = n(),
        .groups = "drop"
      )
  }

  # shared subject-level variables are copies of ADSL, so any disagreement
  # means a stale merge (ADSL refreshed after this dataset was built) or a
  # merge by the wrong keys; compared as character so dates, factors, and
  # numerics all compare without type gymnastics, with NA equal to NA
  shared_vars <- intersect(
    c(
      "TRTSDT", "TRTEDT", "TRT01P", "TRT01A", "SAFFL", "ITTFL",
      "DTHDT", "EOSDT", "AGE", "SEX", "RACE"
    ),
    intersect(cols, adsl_cols)
  )
  if (length(shared_vars) > 0) {
    mismatched_subjects <- function(v) {
      child <- as.character(dataset[[v]])
      ref <- as.character(adsl[[v]])[adsl_idx]
      differ <- (is.na(child) != is.na(ref)) |
        (!is.na(child) & !is.na(ref) & child != ref)
      # orphans are reported by their own check, not per variable
      differ[is.na(adsl_idx)] <- FALSE
      n_distinct(dataset$USUBJID[differ])
    }
    stale <- vapply(shared_vars, mismatched_subjects, integer(1))
    out$stale_total <- sum(stale)
    out$stale_vars <- tibble(
      variable = names(stale),
      subjects = unname(stale)
    ) %>%
      filter(subjects > 0)
  }

  # records dated after the subject's death are always a finding
  date_var <- intersect(c("ADT", "ASTDT"), cols)
  if (length(date_var) > 0 && "DTHDT" %in% adsl_cols) {
    dth <- adsl$DTHDT[adsl_idx]
    dv <- dataset[[date_var[1]]]
    out$after_death_var <- date_var[1]
    out$n_after_death <- sum(!is.na(dv) & !is.na(dth) & dv > dth)
  }

  if (type == "OCCDS" && "TRTEMFL" %in% cols) {
    em_subj <- unique(dataset$USUBJID[
      !is.na(dataset$USUBJID) &
        !is.na(dataset$TRTEMFL) & dataset$TRTEMFL == "Y"
    ])

    # the incidence line every AE table leads with; the dataset alone cannot
    # compute it because subjects with zero occurrences are not in it
    if (!is.null(saf_subj)) {
      out$incidence <- c(
        subjects = length(intersect(em_subj, saf_subj)),
        total = length(saf_subj)
      )
      out$incidence_denom <- "safety population"
    } else {
      out$incidence <- c(
        subjects = length(intersect(em_subj, adsl_subj)),
        total = out$n_adsl
      )
      out$incidence_denom <- "ADSL"
    }

    # a treatment-emergent record for a subject ADSL says was never treated
    # contradicts the flag; checked against the authoritative TRTSDT, so it
    # works even when (and especially when) the merged copy is absent or stale
    if ("TRTSDT" %in% adsl_cols) {
      untreated <- adsl_subj[is.na(adsl$TRTSDT)]
      out$n_emergent_untreated <- length(intersect(em_subj, untreated))
    }
  }

  out
}

#' Print the ADSL-Comparison Part of a Summary
#'
#' @param x The `vs_adsl` element of a `summary_admiral_df` object, see
#'   [summarize_vs_adsl()]
#'
#' @details
#'   Facts (coverage, per-arm coverage, treatment-emergent incidence) are
#'   printed first, then the checks via [print_summary_checks()]: a failing
#'   check is printed as its own bullet with the offending count, the passing
#'   checks are collapsed into a single confirmation line naming each of them,
#'   and a check whose variables are absent is not mentioned at all.
#'
#' @return No return value, called for side effects.
#'
#' @keywords internal
#' @family internal
print_vs_adsl_summary <- function(x) {
  pct_str <- function(n, d) {
    if (d > 0) {
      sprintf(" (%s%%)", formatC(100 * n / d, format = "f", digits = 1))
    } else {
      ""
    }
  }

  saf <- if (!is.null(x$saffl)) {
    paste0(
      " | safety population ", x$saffl[["covered"]], " of ", x$saffl[["total"]]
    )
  } else {
    ""
  }
  cli_text(
    "Compared with ADSL: {x$n_common} of {x$n_adsl} subjects have
     records{pct_str(x$n_common, x$n_adsl)}{saf}"
  )

  if (!is.null(x$arm_coverage)) {
    cli_text(
      "By arm ({x$arm_var}): {paste0(x$arm_coverage[[x$arm_var]], ' ',
       x$arm_coverage$subjects, '/', x$arm_coverage$total, collapse = ' | ')}"
    )
  }

  if (!is.null(x$incidence)) {
    n_em <- x$incidence[["subjects"]]
    d_em <- x$incidence[["total"]]
    em_pct <- if (d_em > 0) {
      sprintf("%s%% of ", formatC(100 * n_em / d_em, format = "f", digits = 1))
    } else {
      ""
    }
    cli_text(
      "Subjects with a treatment-emergent record: {n_em} of {d_em}
       ({em_pct}{x$incidence_denom})"
    )
  }

  print_summary_checks(x, checks = list(
    n_orphans = list(
      failed = "{x$n_orphans} subject{?s} {?is/are} not in ADSL
                (e.g. {.val {head(x$orphans, 3)}})",
      passed = "all subjects are in ADSL"
    ),
    stale_total = list(
      failed = "{x$stale_total} subject/variable combination{?s}
                disagree{?s/} with ADSL (in {.var {x$stale_vars$variable}})",
      passed = "shared subject-level variables match ADSL"
    ),
    n_after_death = list(
      failed = "{x$n_after_death} record{?s} ({.var {x$after_death_var}})
                dated after death ({.var DTHDT} in ADSL)",
      passed = "no records after death"
    ),
    n_emergent_untreated = list(
      failed = "{x$n_emergent_untreated} subject{?s} with treatment-emergent
                records {?has/have} no {.var TRTSDT} in ADSL",
      passed = "no treatment-emergent records for untreated subjects"
    )
  ))

  invisible(NULL)
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
#'
#' @param keys Optional intended key variables of `object`, i.e. the variables
#'   it should have one record per: either a character vector, or a `{metacore}`
#'   specification object from which they are read with [get_admiral_keys()]
#'   (the specification must describe a single dataset -- select one with
#'   `metacore::select_dataset()` first, or pass
#'   `get_admiral_keys(spec, "ADVS")` instead). When supplied, `keys` overrides
#'   both the `"admiral_keys"` attribute and key inference, and the record
#'   structure is reported with `key_source = "supplied"`. A zero-length vector
#'   skips the record structure check.
#'
#' @param adsl Optional subject-level dataset to compare `object` against: a
#'   data frame with one record per `USUBJID`. When supplied (for a non-ADSL
#'   type), the summary gains a comparison section -- subject coverage overall,
#'   within the safety population, and per arm, plus checks for subjects not in
#'   ADSL, shared subject-level variables disagreeing with ADSL, records dated
#'   after death, and (OCCDS) the treatment-emergent incidence and
#'   treatment-emergent records for subjects ADSL says were never treated. See
#'   [summarize_vs_adsl()] for the details of each item. For a dataset of type
#'   `"ADSL"` the argument is ignored with a warning.
#'
#' @param ... Not used; must be empty
#'
#' @return
#'   An object of class `summary_admiral_df` (a named list) which is printed by
#'   [print.summary_admiral_df()]. It always contains `type`, `n_obs`, and
#'   `n_subjects`, and, depending on the dataset structure,
#'   `keys`/`key_source`/`n_duplicate_keys`, `params`, `avisits`,
#'   `n_events`/`n_censored`, the per-type elements `adsl`
#'   ([summarize_adsl()]), `bds` ([summarize_bds()]), and `occds`
#'   ([summarize_occds()]), and -- when the `adsl` argument is supplied --
#'   `vs_adsl` ([summarize_vs_adsl()]). `ds_name` -- the dataset name shown
#'   in the printed heading -- is the name stored by [set_admiral_keys()] or,
#'   failing that, the variable name `summary()` was called on (it is absent
#'   when neither is available, e.g. under the `magrittr` pipe).
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
#'
#' @caption Supplying the intended keys
#'
#' @info The record structure is then checked against the supplied keys instead
#'   of the tagged or inferred ones. The keys can also be read from a
#'   `{metacore}` specification, see [get_admiral_keys()].
#'
#' @code
#' summary(map, keys = c("USUBJID", "AVISIT", "PARAMCD"))
#'
#' @caption Comparing against ADSL
#'
#' @info Supplying `adsl` adds the subject-coverage facts and the referential
#'   checks only a subject-level source of truth makes possible.
#'
#' @code
#' adsl <- tribble(
#'   ~USUBJID,      ~SAFFL, ~TRT01A,
#'   "01-701-1015", "Y",    "Placebo",
#'   "01-701-1028", "Y",    "Xanomeline",
#'   "01-701-1034", "N",    NA_character_
#' )
#'
#' summary(map, adsl = adsl)
summary.admiral_df <- function(object, keys = NULL, adsl = NULL, ...) {
  rlang::check_dots_empty()
  # a malformed `adsl` argument is a caller error, unlike anything found in the
  # data: comparisons against a subject-level dataset which is not unique by
  # subject would be ambiguous, so this aborts rather than reports
  if (!is.null(adsl)) {
    assert_data_frame(adsl, required_vars = exprs(USUBJID))
    n_dup_adsl <- nrow(adsl) - n_distinct(adsl$USUBJID)
    if (n_dup_adsl > 0) {
      cli_abort(c(
        "{.arg adsl} must have one record per {.var USUBJID}.",
        i = "Found {n_dup_adsl} duplicate subject record{?s}."
      ))
    }
  }
  # the dataset name for the heading: prefer the name stored by
  # `set_admiral_keys()`, else the variable the summary was called on -- unless
  # that is a pipe placeholder or a more complex expression than a name
  ds_name <- attr(object, "admiral_ds_name")
  if (is.null(ds_name)) {
    obj_expr <- substitute(object)
    if (is.name(obj_expr) && !as.character(obj_expr) %in% c(".", "_")) {
      ds_name <- as.character(obj_expr)
    }
  }
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
  out$ds_name <- ds_name

  # record structure check: keys supplied to this call win, then keys attached
  # by the spec (`set_admiral_keys()`) or by the `derive_*` function that
  # produced the dataset (its `by_vars`); otherwise fall back to keys inferred
  # from the dataset type and validated against the data
  if (!is.null(keys)) {
    if (inherits(keys, "Metacore")) {
      keys <- get_admiral_keys(keys)
    }
    assert_character_vector(keys)
    key_source <- "supplied"
  } else {
    keys <- attr(object, "admiral_keys")
    key_source <- "declared"
    if (is.null(keys)) {
      keys <- infer_admiral_keys(type, subject_keys, cols, object)
      key_source <- "inferred"
    }
  }
  # declared/supplied keys are plain names, so they can be absent from the
  # dataset (e.g. after a rename, or a spec entry the derivation has not reached
  # yet); the record structure check would then be silently weakened (or
  # skipped), which must not pass unnoticed
  stale_keys <- setdiff(keys, cols)
  if (key_source != "inferred" && length(stale_keys) > 0) {
    cli_warn(c(
      "The {key_source} key variable{?s} {.var {stale_keys}} {?is/are} not in
       the dataset, so the record structure was not checked as {key_source}.",
      i = if (key_source == "declared") {
        "{.var admiral_keys} is not updated by {.fun dplyr::rename} or
         {.fun dplyr::select}."
      },
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

  if (type == "ADSL") {
    out$adsl <- summarize_adsl(object, cols)
  }

  if (type == "BDS") {
    out$bds <- summarize_bds(object, cols)
  }

  if (type == "TTE" && "CNSR" %in% cols) {
    out$n_events <- sum(object$CNSR == 0, na.rm = TRUE)
    out$n_censored <- sum(object$CNSR != 0, na.rm = TRUE)
  }

  if (type == "OCCDS") {
    out$occds <- summarize_occds(object, cols)
  }

  if (!is.null(adsl)) {
    if (type == "ADSL") {
      cli_warn(
        "{.arg adsl} is ignored for a dataset of type {.val ADSL}; it is meant
         for comparing a record-level dataset against its ADSL."
      )
    } else if (!"USUBJID" %in% cols) {
      cli_warn(
        "{.arg adsl} is ignored because the dataset has no {.var USUBJID} to
         match subjects by."
      )
    } else {
      out$vs_adsl <- summarize_vs_adsl(object, cols, adsl, type)
    }
  }

  class(out) <- c("summary_admiral_df", "list")
  out
}

#' Print the Check Section of a Summary
#'
#' Shared by the per-type printers ([print_adsl_summary()],
#' [print_bds_summary()]): partitions the checks defined in `checks` into
#' failed / passed / not run, prints a detailed red bullet per failing check and
#' one green confirmation line naming each passed check. A check whose element
#' is absent from `x` did not run and is not mentioned at all.
#'
#' @param x The per-type element of a `summary_admiral_df` object; each check is
#'   read from `x[[name]]` and fails when its value is greater than zero
#' @param checks A named list; each element is a list with a `failed` cli
#'   message (may use inline markup and reference `x`) and a `passed` plain-text
#'   label (no markup -- cli does not interpret markup in interpolated values)
#'
#' @return No return value, called for side effects.
#'
#' @keywords internal
#' @family internal
print_summary_checks <- function(x, checks) {
  ran <- names(checks)[!vapply(x[names(checks)], is.null, logical(1))]
  failed <- ran[vapply(ran, function(chk) isTRUE(x[[chk]] > 0), logical(1))]
  passed <- setdiff(ran, failed)

  if (length(failed) > 0) {
    msgs <- vapply(checks[failed], function(chk) chk$failed, character(1))
    cli_bullets(setNames(msgs, rep("x", length(msgs))))
  }
  if (length(passed) > 0) {
    passed_labels <- paste(
      vapply(checks[passed], function(chk) chk$passed, character(1)),
      collapse = "; "
    )
    cli_bullets(c(v = "Checks passed: {passed_labels}"))
  }

  invisible(NULL)
}

#' Print the ADSL-Specific Part of a Summary
#'
#' @param x The `adsl` element of a `summary_admiral_df` object, see
#'   [summarize_adsl()]
#'
#' @details
#'   Facts are printed first, then the checks. A failing check is printed as its
#'   own bullet with the offending count; the checks which passed are collapsed
#'   into a single confirmation line which names each of them, so a clean ADSL
#'   stays one line without hiding what was checked. A check whose variables are
#'   absent is not mentioned at all -- not run is different from passed.
#'
#' @return No return value, called for side effects.
#'
#' @keywords internal
#' @family internal
print_adsl_summary <- function(x) {
  collapse_counts <- function(counts, var) {
    paste0(counts[[var]], ": ", counts$n, collapse = " | ")
  }

  if (!is.null(x$trt)) {
    cli_text("Treatment ({x$trt_var}): {collapse_counts(x$trt, x$trt_var)}")
  }
  if (!is.null(x$pops)) {
    cli_text("Populations: {paste0(names(x$pops), ' ', x$pops, collapse = ' | ')}")
  }
  if (!is.null(x$funnel)) {
    cli_text(
      "Subject flow: {paste0(names(x$funnel), ' ', x$funnel, collapse = ' | ')}"
    )
  }
  if (!is.null(x$trtdur)) {
    cli_text(
      "Treatment duration (days): median {x$trtdur[['median']]}
       (range {x$trtdur[['min']]}-{x$trtdur[['max']]})"
    )
  }
  if (!is.null(x$eosstt)) {
    cli_text("End of study (EOSSTT): {collapse_counts(x$eosstt, 'EOSSTT')}")
  }
  if (!is.null(x$n_dth)) {
    caus <- if (!is.null(x$dthcaus)) {
      paste0(" (", collapse_counts(x$dthcaus, "DTHCAUS"), ")")
    } else {
      ""
    }
    cli_text("Deaths (DTHFL): {x$n_dth}{caus}")
  }
  if (!is.null(x$age) || !is.null(x$sex)) {
    age <- if (!is.null(x$age)) {
      paste0(
        "age median ", x$age[["median"]],
        " (range ", x$age[["min"]], "-", x$age[["max"]], ")"
      )
    }
    sex <- if (!is.null(x$sex)) collapse_counts(x$sex, "SEX")
    cli_text("Demographics: {paste(c(age, sex), collapse = ' | ')}")
  }

  # see print_summary_checks() for how failed/passed/not-run are reported
  print_summary_checks(x, checks = list(
    arm_mismatch = list(
      failed = "{x$arm_mismatch} subject{?s} {?has/have} different {.var ARM}
                and {.var ACTARM}",
      passed = "ARM matches ACTARM"
    ),
    saffl_no_trtsdt = list(
      failed = "{x$saffl_no_trtsdt} subject{?s} in the safety population
                {?has/have} no {.var TRTSDT}",
      passed = "no missing TRTSDT in the safety population"
    ),
    trt_end_before_start = list(
      failed = "{x$trt_end_before_start} subject{?s} {?has/have} {.var TRTEDT}
                before {.var TRTSDT}",
      passed = "no TRTEDT before TRTSDT"
    ),
    flag_domain = list(
      failed = "{x$flag_domain} population flag value{?s} outside Y/N/NA
                (in {.var {x$flag_domain_vars}})",
      passed = "population flags only contain Y/N/NA"
    )
  ))

  invisible(NULL)
}

#' Print the BDS-Specific Part of a Summary
#'
#' @param x The `bds` element of a `summary_admiral_df` object, see
#'   [summarize_bds()]
#'
#' @details
#'   The per-parameter table is truncated to the first `max_params` rows on
#'   printing (the full table is always in the object); numeric statistics are
#'   shown to 4 significant digits with `.` for a missing value. The checks are
#'   printed by [print_summary_checks()].
#'
#' @param max_params Maximum number of parameter rows to print
#'
#' @return No return value, called for side effects.
#'
#' @keywords internal
#' @family internal
print_bds_summary <- function(x, max_params = 10) {
  if (!is.null(x$params)) {
    cli_text("Parameters ({nrow(x$params)}):")
    render_summary_table(head(x$params, max_params))
    if (nrow(x$params) > max_params) {
      cli_text(
        "... and {nrow(x$params) - max_params} more parameter{?s}
         (in {.code $bds$params})"
      )
    }
  }

  if (!is.null(x$dtype)) {
    cli_text(
      "Derived records (DTYPE):
       {paste0(x$dtype$DTYPE, ': ', x$dtype$n, collapse = ' | ')}"
    )
  }

  print_summary_checks(x, checks = list(
    multiple_baselines = list(
      failed = "{x$multiple_baselines} subject-parameter combination{?s}
                {?has/have} more than one baseline record ({.var ABLFL})",
      passed = "at most one baseline per subject and parameter"
    ),
    chg_no_base = list(
      failed = "{x$chg_no_base} record{?s} {?has/have} {.var CHG}/{.var PCHG}
                but no {.var BASE}",
      passed = "no CHG/PCHG without BASE"
    ),
    avisit_mismatch = list(
      failed = "{x$avisit_mismatch} {.var AVISIT}/{.var AVISITN} value{?s}
                map{?s/} to more than one counterpart",
      passed = "AVISIT and AVISITN are consistent"
    ),
    param_inconsistent = list(
      failed = "{x$param_inconsistent} {.var PARAMCD}{?s} {?has/have} more than
                one {.var PARAM} or {.var AVALU} value",
      passed = "one PARAM/AVALU per PARAMCD"
    )
  ))

  invisible(NULL)
}

#' Print the OCCDS-Specific Part of a Summary
#'
#' @param x The `occds` element of a `summary_admiral_df` object, see
#'   [summarize_occds()]
#'
#' @details
#'   Facts are printed first, then the checks via [print_summary_checks()]: a
#'   failing check is printed as its own bullet with the offending count, the
#'   passing checks are collapsed into a single confirmation line naming each of
#'   them, and a check whose variables are absent is not mentioned at all.
#'
#' @return No return value, called for side effects.
#'
#' @keywords internal
#' @family internal
print_occds_summary <- function(x) {
  if (!is.null(x$terms)) {
    cli_text(
      "Distinct terms: {paste0(names(x$terms), ' ', x$terms, collapse = ' | ')}"
    )
  }
  if (!is.null(x$trtem)) {
    cli_text(
      "Treatment emergent (TRTEMFL): {x$trtem[['emergent']]} of
       {x$trtem[['total']]} records"
    )
  }
  if (!is.null(x$sev)) {
    cli_text(
      "Severity/grade ({x$sev_var}):
       {paste0(x$sev[[x$sev_var]], ': ', x$sev$n, collapse = ' | ')}"
    )
  }
  if (!is.null(x$n_serious)) {
    cli_text("Serious (AESER): {x$n_serious}")
  }

  print_summary_checks(x, checks = list(
    occ_flag_dups = list(
      failed = "{x$occ_flag_dups} occurrence flag group{?s} {?has/have} more
                than one {.val Y} (in {.var {x$occ_flag_vars}})",
      passed = "occurrence flags are unique per subject and level"
    ),
    missing_astdt = list(
      failed = "{x$missing_astdt} record{?s} {?has/have} no {.var ASTDT}",
      passed = "no missing ASTDT"
    ),
    pre_trt_emergent = list(
      failed = "{x$pre_trt_emergent} treatment-emergent record{?s} start{?s/}
                before {.var TRTSDT}",
      passed = "no treatment-emergent record before TRTSDT"
    )
  ))

  invisible(NULL)
}

#' Print a Data Frame as an Aligned Text Table
#'
#' Minimal fixed-width rendering for the tables inside a summary: header row
#' plus one line per row, character columns left-aligned, numeric columns
#' right-aligned and shown to 4 significant digits, `NA` as `.`.
#'
#' @param tbl A data frame
#' @param indent Prefix prepended to every line
#'
#' @return No return value, called for side effects.
#'
#' @keywords internal
#' @family internal
render_summary_table <- function(tbl, indent = "  ") {
  fmt_cell <- function(v) {
    if (is.na(v)) {
      "."
    } else if (is.numeric(v)) {
      format(signif(v, 4), trim = TRUE, scientific = FALSE)
    } else {
      as.character(v)
    }
  }
  cells <- lapply(tbl, function(col) vapply(col, fmt_cell, character(1)))
  right <- vapply(tbl, is.numeric, logical(1))
  widths <- mapply(
    function(nm, cell) max(nchar(c(nm, cell), type = "width")),
    names(cells), cells
  )
  pad <- function(txt, w, r) formatC(txt, width = if (r) w else -w)
  line <- function(txt_by_col) {
    paste0(indent, paste(mapply(pad, txt_by_col, widths, right), collapse = "  "))
  }

  # `cli_verbatim()` because `cli_text()` collapses the alignment spaces
  cli_verbatim(line(names(cells)))
  for (i in seq_len(nrow(tbl))) {
    cli_verbatim(line(vapply(cells, `[[`, character(1), i)))
  }

  invisible(NULL)
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
  # the type is what tailors everything below, so it earns the headline; "other"
  # would read oddly there ("other summary"), hence the generic fallback
  title <- if (x$type == "other") "Dataset summary" else paste(x$type, "summary")
  if (!is.null(x$ds_name)) {
    title <- paste0(x$ds_name, ": ", title)
  }
  cli_rule(left = title)

  subj_label <- if (length(x$subject_keys) > 0) {
    paste0(" (", paste(x$subject_keys, collapse = ", "), ")")
  } else {
    ""
  }
  cli_text("Subjects{subj_label}: {x$n_subjects}")
  cli_text("Observations: {x$n_obs}")

  if (!is.null(x$keys)) {
    keys <- paste(x$keys, collapse = ", ")
    src <- switch(x$key_source,
      inferred = " (inferred)",
      supplied = " (supplied)",
      ""
    )
    if (x$n_duplicate_keys == 0) {
      cli_bullets(c(v = "Structure{src}: one record per {keys}"))
    } else {
      cli_bullets(c(x = "Structure{src}: expected one record per {keys}, found
         {x$n_duplicate_keys} extra record{?s}"))
    }
  }

  # for BDS the flat parameter list is superseded by the per-parameter table
  if (!is.null(x$params) && is.null(x$bds)) {
    params <- paste(x$params$PARAMCD, collapse = ", ")
    cli_text("Parameters (PARAMCD): {params}")
  }

  if (!is.null(x$avisits)) {
    avisits <- paste(x$avisits, collapse = ", ")
    cli_text("Analysis visits (AVISIT): {avisits}")
  }

  if (!is.null(x$adsl)) {
    print_adsl_summary(x$adsl)
  }

  if (!is.null(x$bds)) {
    print_bds_summary(x$bds)
  }

  if (!is.null(x$occds)) {
    print_occds_summary(x$occds)
  }

  if (!is.null(x$n_events)) {
    cli_text("Events: {x$n_events} | Censored: {x$n_censored}")
  }

  if (!is.null(x$vs_adsl)) {
    print_vs_adsl_summary(x$vs_adsl)
  }

  invisible(x)
}
