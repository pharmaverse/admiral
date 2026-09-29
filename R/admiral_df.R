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

#' Find the Minimal Set of Variables that Uniquely Identifies Rows
#'
#' Starting from `must_have`, adds variables from `optional` one at a time (in
#' the given order) until the combination is a unique key of `dataset`, then
#' drops any added variable which is not needed for uniqueness (lowest priority
#' first) and returns what is left. Without that second pass the result would
#' depend on where a genuinely needed variable sits in `optional`: everything
#' ahead of it would be carried along whether it discriminates or not, and the
#' reported record structure would overstate the key. If uniqueness is never
#' reached, the full set (`must_have` plus all of `optional`) is returned -- the
#' caller can detect this because the returned key still yields duplicates.
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
  # `ungroup()` because `distinct()` on a grouped dataset silently adds the
  # grouping variables, which would make any key look unique within group
  dataset <- ungroup(dataset)
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
      # the walk stops at the first unique *prefix* of `optional`, so the key
      # can carry variables passed over on the way which contribute nothing to
      # uniqueness. Drop those, lowest priority first, so the reported record
      # structure is the one the dataset actually has rather than an artefact
      # of the candidate ordering.
      for (redundant in rev(setdiff(key, must_have))) {
        if (is_unique(setdiff(key, redundant))) {
          key <- setdiff(key, redundant)
        }
      }
      return(key)
    }
  }
  key
}

#' Infer the Record Structure of a Dataset from its Type and Data
#'
#' Works out which variables a dataset appears to have one record per. This is
#' the fallback for when no keys have been declared with [set_admiral_keys()],
#' and is the least trustworthy of the ways the record structure can be
#' established -- it reports the structure the data happens to have, which is not
#' necessarily the structure the dataset was meant to have.
#'
#' It starts from the semantic core of the detected ADaM dataset type (e.g.
#' `USUBJID` + `PARAMCD` for `BDS`) and adds standard ADaM key variables
#' (analysis visit, timepoint, relative time, date, interval start/end, period,
#' derivation type) only as far as needed to make the key unique, then drops
#' those which turn out not to discriminate (see [minimal_unique_key()]). The
#' candidate list has to span several dataset shapes: findings data is keyed by
#' visit and timepoint, exposure by interval start, population-PK by relative
#' time and not by visit at all, and any of them may carry derived records
#' separated only by `DTYPE`.
#'
#' @param dataset A data frame
#' @param type The dataset type, see [get_admiral_df_type()]
#'
#' @details
#'   For `BDS`/`TTE`, surrogate/sequence keys such as `ASEQ` are deliberately
#'   *not* used, so that unintended duplicate records remain detectable. For
#'   `OCCDS`, however, the sequence number *is* the intended record key (there is
#'   no analysis-value structure to fall back on: two adverse events for the same
#'   subject need not differ in any analysis variable), so it is required rather
#'   than excluded. `ASEQ` is preferred, falling back to the SDTM domain sequence
#'   (`AESEQ`, `CMSEQ`, `MHSEQ`, ...) which is what occurrence datasets most
#'   often carry; the two-letter domain prefix is what distinguishes those from
#'   provenance variables such as `SRCSEQ`. A sequence variable which is present
#'   but entirely `NA` has not been derived yet and is passed over, as an empty
#'   sequence would otherwise make the key unique without meaning anything. If no
#'   sequence is left -- or several domain sequences are, making the record key
#'   ambiguous -- the record structure cannot be checked and a warning is issued.
#'
#'   Note that a well-formed sequence makes the key unique by construction, so
#'   the `OCCDS` structure check is not looking for semantic duplicates but for
#'   whole records duplicated by a fanned-out merge.
#'
#'   A dataset with no records cannot show what it has one record per -- every
#'   candidate key is trivially unique -- so nothing is inferred for it.
#'
#' @return A character vector of inferred key variable names, or a zero-length
#'   vector when no plausible record structure can be determined.
#'
#' @keywords internal
#' @family internal
infer_admiral_keys <- function(dataset, type = get_admiral_df_type(dataset)) {
  cols <- colnames(dataset)

  # with no records every candidate key is trivially unique, so inference would
  # report a structure the dataset has not demonstrated
  if (nrow(dataset) == 0) {
    return(character(0))
  }

  # `USUBJID` identifies a subject on its own, and is preferred over the full
  # set of subject keys because the others (e.g. `STUDYID`) may be `NA` on
  # records added by a derivation which only populates its `by_vars`
  subject_keys <- if ("USUBJID" %in% cols) {
    "USUBJID"
  } else {
    intersect(vars2chr(get_admiral_option("subject_keys")), cols)
  }

  # OCCDS has no analysis-value structure to fall back on; the sequence number
  # is its genuine record key, so (unlike BDS/TTE) it is required rather than
  # excluded as a surrogate
  if (type == "OCCDS") {
    seq_var <- occds_seq_var(dataset, cols)
    # without the sequence there is no record key at all: the subject keys alone
    # would claim one record per subject, which an occurrence dataset never has
    if (length(seq_var) == 0) {
      return(character(0))
    }
    return(intersect(c(subject_keys, seq_var), cols))
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

  # standard "within-core" key variables, in ADaM precedence order: visit,
  # timepoint, relative time (PK), analysis date, interval start/end, period,
  # then `DTYPE` last -- a derived record (LOCF, AVERAGE, ...) shares every
  # analysis variable with the record it was derived from, so `DTYPE` is what
  # separates them, but only after the semantic variables have had their turn.
  # NOTE: no surrogate/sequence keys (ASEQ, SRCSEQ, ...) -- see minimal_unique_key()
  extra <- c(
    "AVISITN", "AVISIT", "ATPTN", "ATPT",
    "NFRLT", "AFRLT",
    "ADTM", "ADT", "ASTDTM", "ASTDT", "AENDT",
    "APERIOD", "APERIODC", "ASPID", "DTYPE"
  )
  minimal_unique_key(
    dataset,
    must_have = core,
    optional = intersect(extra, cols)
  )
}

#' Find the Record Key of an Occurrence Dataset
#'
#' @param dataset A data frame
#' @param cols The column names of `dataset`
#'
#' @return The name of the sequence variable which keys `dataset`, or a
#'   zero-length vector (with a warning) when none can be determined.
#'
#' @seealso [infer_admiral_keys()], which documents the choice this makes
#'
#' @keywords internal
#' @family internal
occds_seq_var <- function(dataset, cols) {
  # a sequence which is present but entirely `NA` has not been derived yet;
  # treating it as the record key would make the key unique without meaning
  # anything, which is the opposite of what the structure check is for
  is_populated <- function(v) any(!is.na(dataset[[v]]))

  # `ASEQ` is the ADaM analysis sequence, but many occurrence datasets carry
  # only the SDTM domain sequence (`AESEQ`, `CMSEQ`, `MHSEQ`, ...). The
  # two-letter domain prefix is what distinguishes those from provenance
  # variables such as `SRCSEQ`, which is not a record key.
  if ("ASEQ" %in% cols && is_populated("ASEQ")) {
    return("ASEQ")
  }
  domain_seq <- str_subset(cols, "^[A-Z]{2}SEQ$")
  domain_seq <- domain_seq[vapply(domain_seq, is_populated, logical(1))]
  if (length(domain_seq) == 1) {
    return(domain_seq)
  }

  cli_warn(c(
    "Cannot determine the record structure of an {.val OCCDS} dataset without
     {.var ASEQ} or a single populated domain sequence variable.",
    i = if (length(domain_seq) > 1) {
      "{.var {domain_seq}} are all candidates, so the record key is ambiguous."
    },
    i = "Declare the intended keys with {.fun set_admiral_keys} to check that
         the dataset has one record per key."
  ))
  character(0)
}
