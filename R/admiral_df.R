#' Tag a Dataset with the `admiral_df` Class
#'
#' Adds the `admiral_df` class to a data frame (unless it is already
#' present), preserving the existing classes such as `tbl_df`/`data.frame`.
#' This is the tagging primitive that admiral's `derive_*()` and related
#' dataset-producing functions apply to their output, so that a dataset
#' built with admiral can be recognized as such downstream (e.g. by a
#' future dedicated `summary()` method).
#'
#' @param dataset A data frame (or `NULL`)
#'
#' @return
#'   If `dataset` is `NULL`, `NULL` is returned unchanged. Otherwise
#'   `dataset` with `"admiral_df"` prepended to its class attribute.
#'
#' @keywords utils_help
#' @family utils_help
#'
#' @export
#'
#' @examples
#' library(tibble)
#'
#' adsl <- tribble(
#'   ~USUBJID, ~AGE,
#'   "1",      63,
#'   "2",      71
#' )
#'
#' class(as_admiral_df(adsl))
as_admiral_df <- function(dataset) {
  if (is.null(dataset)) {
    return(dataset)
  }
  if (!inherits(dataset, "admiral_df")) {
    class(dataset) <- c("admiral_df", class(dataset))
  }
  dataset
}

#' Tag a Dataset with its Key Variables
#'
#' Stores the key variables of a dataset -- the variables it is intended to have
#' exactly one record per -- in its `"admiral_keys"` attribute, and tags it with
#' the `admiral_df` class (see [as_admiral_df()]). Downstream admiral tooling
#' reads the attribute to check that the dataset really does have one record per
#' key, e.g. after a derivation which was expected to preserve the record
#' structure.
#'
#' @param dataset A data frame
#'
#' @permitted a data frame
#'
#' @param keys The key variables of `dataset`, in key order, e.g.
#'   `exprs(USUBJID, PARAMCD, AVISIT)` or
#'   `c("USUBJID", "PARAMCD", "AVISIT")`. A zero-length vector records that the
#'   dataset has no key variables defined, which suppresses the record structure
#'   check rather than leaving the keys to be guessed.
#'
#' @permitted list of variables created by `exprs()`, or a character vector
#'
#' @param dataset_name The name of the dataset, e.g. `"ADVS"`, stored in the
#'   `"admiral_ds_name"` attribute and used as a heading by the admiral tooling
#'   which reports on the dataset. `NULL` leaves a name stored by a previous
#'   call untouched, so re-keying a dataset does not lose its name.
#'
#' @permitted character scalar or `NULL`
#'
#' @return `dataset` with an `"admiral_keys"` attribute, an `"admiral_ds_name"`
#'   attribute if `dataset_name` was supplied or had been stored by a previous
#'   call, and the `admiral_df` class.
#'
#' @details
#'   The keys can come from anywhere: the `by_vars` of the derivation which
#'   produced the dataset, a dataset specification, or a hand-written vector.
#'   They are accepted both as `exprs()` -- matching `by_vars` and the
#'   `subject_keys` admiral option -- and as a character vector, which is what a
#'   specification provides. admiral does not read specification objects itself;
#'   to use the keys defined in a `{metacore}` specification, extract them in the
#'   calling code, e.g. `metacore::get_keys(spec, "ADVS")$variable`.
#'
#'   The core `{dplyr}` verbs (`mutate()`, `filter()`, `arrange()`, `select()`,
#'   `rename()`, `slice()`, `distinct()`, joins where `dataset` is `x`, ...)
#'   restore data frame attributes, so the keys and the `admiral_df` class
#'   normally survive a derivation pipeline. They are *not* preserved by
#'   operations which build a new data frame, in particular:
#'
#'   * `summarise()`,
#'   * `tidyr::pivot_longer()`/`pivot_wider()`, `nest()`/`unnest()`,
#'     `complete()`,
#'   * `bind_rows()` when `dataset` is not the *first* argument,
#'   * joins where `dataset` is the `y` (rather than the `x`) argument.
#'
#'   `group_by()` and `rowwise()` are a case of their own: they keep the keys but
#'   drop the `admiral_df` class, leaving the attribute behind on a dataset the
#'   admiral tooling no longer recognizes.
#'
#'   Call this function afterwards (or re-apply it as needed) if the dataset
#'   passes through any of those.
#'
#'   Because the keys are stored as plain variable names, they are not updated by
#'   `rename()` and not removed by `select()`. Keys which are not (or not yet) in
#'   the dataset are accepted without complaint, so that a partially derived
#'   dataset can declare the structure it is being built towards; it is the
#'   admiral tooling which reads the attribute that reports keys no longer
#'   matching the dataset, once the derivation is finished and a missing key is
#'   unambiguously wrong.
#'
#' @keywords utils_help
#' @family utils_help
#'
#' @export
#'
#' @examplesx
#'
#' @caption Tagging a dataset with its keys
#'
#' @code
#' library(tibble)
#'
#' advs <- tribble(
#'   ~USUBJID, ~PARAMCD, ~AVISIT,    ~AVAL,
#'   "1",      "DIABP",  "BASELINE",    51,
#'   "1",      "SYSBP",  "BASELINE",   121,
#'   "2",      "DIABP",  "BASELINE",    79
#' )
#'
#' advs <- set_admiral_keys(
#'   advs,
#'   keys = exprs(USUBJID, PARAMCD, AVISIT),
#'   dataset_name = "ADVS"
#' )
#'
#' attr(advs, "admiral_keys")
#' class(advs)
#'
#' @caption Keys as a character vector
#'
#' @info Keys are equally accepted as a character vector, which is the form a
#'   dataset specification provides them in. admiral does not read specification
#'   objects itself, so they are extracted in the calling code; with a
#'   `{metacore}` specification that is `metacore::get_keys()`:
#'
#'   ```r
#'   set_admiral_keys(advs, metacore::get_keys(spec, "ADVS")$variable, "ADVS")
#'   ```
#'
#' @code
#' advs <- set_admiral_keys(advs, keys = c("USUBJID", "PARAMCD", "AVISIT"))
#'
#' # `dataset_name` was not repeated, so the name from the previous call stands
#' attr(advs, "admiral_ds_name")
#'
#' @caption Declaring that a dataset has no key variables
#'
#' @info A zero-length `keys` is stored all the same. It is not the same as
#'   leaving the attribute unset: it records that no key variables are defined,
#'   so the record structure is not checked against a guess.
#'
#' @code
#' advs_no_keys <- set_admiral_keys(advs, keys = character(0))
#'
#' attr(advs_no_keys, "admiral_keys")
set_admiral_keys <- function(dataset, keys, dataset_name = NULL) {
  assert_data_frame(dataset)
  # keys are accepted as `exprs()` (matching `by_vars` and the `subject_keys`
  # option) as well as a character vector (the form a specification provides),
  # but are always stored as character: they are metadata about the dataset, not
  # expressions to evaluate against it
  if (is.list(keys)) {
    keys <- vars2chr(assert_vars(keys))
  } else {
    assert_character_vector(keys)
  }
  assert_character_scalar(dataset_name, optional = TRUE)

  # `as.character()` drops attributes a key vector carries when it was pulled
  # out of a specification table (e.g. a variable label), so that the stored
  # keys compare equal to a plain character vector
  attr(dataset, "admiral_keys") <- as.character(keys)
  if (!is.null(dataset_name)) {
    attr(dataset, "admiral_ds_name") <- dataset_name
  }
  as_admiral_df(dataset)
}

#' Determine the ADaM Dataset Type
#'
#' Classifies a data frame as one of the common ADaM dataset types based on the
#' variables it contains. This lets admiral tooling tailor what it reports to the
#' kind of dataset it has been given.
#'
#' Note that this is the *type* of the dataset (`ADSL`, `BDS`, ...), not its
#' record structure: the latter is the set of key variables it has one record
#' per, see [set_admiral_keys()].
#'
#' @param dataset A data frame
#'
#' @details
#'   The following precedence is used. The order matters, because the types are
#'   not mutually exclusive by variable content -- a `TTE` dataset contains
#'   `PARAMCD` just as a `BDS` one does, and an occurrence dataset which has had
#'   a parameter added is structurally both.
#'
#'   1. `"TTE"` -- contains `PARAMCD`, `CNSR`, and `STARTDT`. Tested first
#'      because a time-to-event dataset is a `BDS` dataset by variable content;
#'      only the censoring and origin-date variables distinguish it.
#'   2. `"OCCDS"` -- does not contain `PARAMCD`, identifies a subject, and
#'      contains a `--DECOD`/`--TERM` variable, an `AOCCxxFL` occurrence flag, or
#'      `TRTEMFL`.
#'   3. `"BDS"` -- contains `PARAMCD` and `AVAL` or `AVALC`.
#'   4. `"ADSL"` -- does not contain `PARAMCD` and has an ADSL structure, see
#'      [is_adsl_structure()].
#'   5. `"other"` -- none of the above.
#'
#'   The `--DECOD`/`--TERM` test requires the two-letter SDTM domain prefix the
#'   convention gives those variables (`AEDECOD`, `CMTERM`, ...) and requires the
#'   dataset to identify a subject as well. Matching a bare `TERM$` suffix with
#'   no structural corroboration classified any data frame with a variable whose
#'   name happens to end in `TERM` as an occurrence dataset.
#'
#'   A dataset which has not yet been derived far enough to carry any of these
#'   signals is `"other"`, which is a deliberate under-claim: reporting on an
#'   unrecognized dataset is confined to what holds for any data frame, whereas
#'   claiming the wrong type would report the wrong things about it.
#'
#' @return A character scalar: one of `"ADSL"`, `"BDS"`, `"OCCDS"`, `"TTE"`, or
#'   `"other"`.
#'
#' @keywords internal
#' @family internal
get_admiral_df_type <- function(dataset) {
  # a grouped or rowwise dataset is classified rather than refused: this is a
  # diagnostic, and `is_adsl_structure()` ungroups before testing the structure
  assert_data_frame(dataset, check_is_grouped = FALSE, check_is_rowwise = FALSE)
  cols <- colnames(dataset)
  has <- function(x) all(x %in% cols)
  has_paramcd <- has("PARAMCD")
  # occurrence datasets are record-level, so a subject identifier is what
  # corroborates the naming signals below
  has_subject <-
    length(intersect(vars2chr(get_admiral_option("subject_keys")), cols)) > 0

  is_tte <- has_paramcd && has("CNSR") && has("STARTDT")
  is_occds <- !has_paramcd && has_subject &&
    any(str_detect(cols, "^[A-Z]{2}(DECOD|TERM)$|^AOCC[0-9A-Z]*FL$|^TRTEMFL$"))
  is_bds <- has_paramcd && (has("AVAL") || has("AVALC"))
  is_adsl <- !has_paramcd && is_adsl_structure(dataset)

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
#'
#' @details
#'   A dataset is subject-level if it has one record per subject, with respect to
#'   `get_admiral_option("subject_keys")`. The records decide this, not the
#'   variables: treatment variables are no evidence of a subject-level structure,
#'   because the first step of a typical `BDS`/`OCCDS` derivation is to merge
#'   them on from ADSL (see the `admiral` template scripts), which would make
#'   every findings dataset look subject-level between that merge and the point
#'   at which it acquires `PARAMCD`.
#'
#'   The exception is a dataset with no records, which cannot demonstrate its
#'   structure either way -- record and subject counts are trivially equal. There
#'   the declared variables are all there is to go on, so a zero-row dataset is
#'   subject-level if it carries a period-numbered treatment variable (`TRT01P`,
#'   `TRT02A`, ...) or a treatment start date.
#'
#'   A partially derived ADSL which has not yet reached one record per subject is
#'   therefore not recognized here, and is typed `"other"` by
#'   [get_admiral_df_type()]. This is a known limitation rather than an
#'   oversight: at that point nothing in the dataset distinguishes it from any
#'   other record-level data frame, and guessing would misreport a dataset which
#'   is genuinely still being built.
#'
#' @return `TRUE` if the dataset has a subject-level structure, `FALSE`
#'   otherwise.
#'
#' @keywords internal
#' @family internal
is_adsl_structure <- function(dataset) {
  cols <- colnames(dataset)

  # with no records there is nothing to test the structure against, so the
  # variables the dataset declares are all there is to go on
  if (nrow(dataset) == 0) {
    return(
      any(str_detect(cols, "^TRT[0-9]{2}[PA]$")) ||
        any(c("TRTSDT", "TRTSDTM") %in% cols)
    )
  }

  subject_keys <- intersect(vars2chr(get_admiral_option("subject_keys")), cols)
  # `ungroup()` because `distinct()` on a grouped dataset silently adds the
  # grouping variables, which would test uniqueness within group rather than
  # overall and report a grouped findings dataset as subject-level
  length(subject_keys) > 0 &&
    nrow(dataset) == nrow(distinct(ungroup(dataset), !!!syms(subject_keys)))
}
