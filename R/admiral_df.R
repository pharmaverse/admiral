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
#'   The keys are stored as plain variable names, so they are not updated by
#'   `rename()` or removed by `select()`. Keys which are not (or not yet) in the
#'   dataset are accepted, so that a partially derived dataset can declare the
#'   structure it is being built towards; the admiral tooling reading the
#'   attribute is what reports keys which no longer match.
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
  # keys are always stored as character: they are metadata about the dataset,
  # not expressions to evaluate against it
  if (is.list(keys)) {
    keys <- vars2chr(assert_vars(keys))
  } else {
    assert_character_vector(keys)
  }
  assert_character_scalar(dataset_name, optional = TRUE)

  # `as.character()` drops attributes a key vector pulled out of a specification
  # table carries (e.g. a variable label), which would break comparison against
  # a plain character vector
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
#'   convention gives those variables (`AEDECOD`, `CMTERM`, ...), so that a
#'   variable whose name merely ends in `TERM` is not taken as a signal.
#'
#'   A dataset which has not yet been derived far enough to carry any of these
#'   signals is `"other"`; reporting on it is then confined to what holds for
#'   any data frame.
#'
#' @return A character scalar: one of `"ADSL"`, `"BDS"`, `"OCCDS"`, `"TTE"`, or
#'   `"other"`.
#'
#' @keywords internal
#' @family internal
get_admiral_df_type <- function(dataset) {
  # a grouped or rowwise dataset is classified rather than refused: this is a
  # diagnostic
  assert_data_frame(dataset, check_is_grouped = FALSE, check_is_rowwise = FALSE)
  cols <- colnames(dataset)
  has <- function(x) x %in% cols
  has_paramcd <- has("PARAMCD")
  # occurrence datasets are record-level, so a subject identifier is what
  # corroborates the naming signals below
  has_subject <-
    length(intersect(vars2chr(get_admiral_option("subject_keys")), cols)) > 0

  is_tte <- has_paramcd && has("CNSR") && has("STARTDT")
  is_occds <- !has_paramcd && has_subject &&
    any(str_detect(cols, "^[A-Z]{2}(DECOD|TERM)$|^AOCC[0-9A-Z]*FL$|^TRTEMFL$"))
  is_bds <- has_paramcd && (has("AVAL") || has("AVALC"))

  # Check dataset type
  if (is_tte) return("TTE")
  if (is_occds) return("OCCDS")
  if (is_bds) return("BDS")
  if (!has_paramcd && is_adsl_structure(dataset)) return("ADSL")
  return("other")
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
#'   [get_admiral_df_type()]: nothing in it yet distinguishes it from any other
#'   record-level data frame.
#'
#' @return `TRUE` if the dataset has a subject-level structure, `FALSE`
#'   otherwise.
#'
#' @keywords internal
#' @family internal
is_adsl_structure <- function(dataset) {
  cols <- colnames(dataset)

  # with no records there is nothing to test the structure against
  if (nrow(dataset) == 0) {
    return(
      any(str_detect(cols, "^TRT[0-9]{2}[PA]$")) ||
        any(c("TRTSDT", "TRTSDTM") %in% cols)
    )
  }

  subject_keys <- intersect(vars2chr(get_admiral_option("subject_keys")), cols)
  # `ungroup()` because `distinct()` silently adds the grouping variables, which
  # would test uniqueness within group rather than overall
  length(subject_keys) > 0 &&
    nrow(dataset) == nrow(distinct(ungroup(dataset), !!!syms(subject_keys)))
}

#' Find the Minimal Set of Variables that Uniquely Identifies Rows
#'
#' Starting from `must_have`, adds variables from `optional` one at a time (in
#' the given order) until the combination is a unique key of `dataset`, then
#' drops any added variable which is not needed for uniqueness. If uniqueness is
#' never reached, the full set (`must_have` plus all of `optional`) is returned
#' -- the caller can detect this because the returned key still yields
#' duplicates.
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
  # `ungroup()` for the reason given in `is_adsl_structure()`
  dataset <- ungroup(dataset)
  # `n_distinct()` rather than `nrow(distinct())`: the count is all that is
  # wanted, and this runs often enough (up to twice per candidate) that
  # materializing the deduplicated data frame each time is worth avoiding
  is_unique <- function(key) {
    length(key) > 0 && nrow(dataset) == n_distinct(dataset[key])
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
      # uniqueness. Drop those, lowest priority first -- but not `v` itself,
      # which is what made the key unique, so the key without it is the one
      # just tested and found not to be.
      for (redundant in rev(setdiff(key, c(must_have, v)))) {
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
#'   than excluded; see [occds_seq_var()] for how it is chosen. A well-formed
#'   sequence makes the key unique by construction, so the `OCCDS` structure
#'   check is not looking for semantic duplicates but for whole records
#'   duplicated by a fanned-out merge.
#'
#'   The subject part of the key follows the `subject_keys` admiral option, as
#'   [is_adsl_structure()] does, reduced to those of its variables which
#'   actually distinguish subjects in the dataset at hand. With the default
#'   option that is `USUBJID` alone for a single-study dataset, because
#'   `STUDYID` is constant there and would only overstate the record structure;
#'   a pooled dataset whose subject identifiers repeat across studies keeps
#'   `STUDYID` as well.
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
  # `type` is forced before `dataset` is ungrouped below, so that the default
  # argument classifies the dataset it was called with
  force(type)
  # `n_distinct()` on a grouped dataset would count within group
  dataset <- ungroup(dataset)

  # with no records every candidate key is trivially unique
  if (nrow(dataset) == 0) {
    return(character(0))
  }

  # the subject keys as configured, minus any which do not distinguish subjects
  # in this dataset: `STUDYID` is constant within a single study, and carrying
  # it would overstate the record structure in the way `minimal_unique_key()`
  # exists to avoid. So the default keys reduce to `USUBJID`, while a pooled
  # dataset whose subject identifiers repeat across studies keeps both.
  subject_keys <- intersect(vars2chr(get_admiral_option("subject_keys")), cols)
  if (length(subject_keys) > 1) {
    n_subjects <- n_distinct(dataset[subject_keys])
    for (v in subject_keys) {
      rest <- setdiff(subject_keys, v)
      if (length(rest) > 0 && n_distinct(dataset[rest]) == n_subjects) {
        subject_keys <- rest
      }
    }
  }

  if (type == "OCCDS") {
    seq_var <- occds_seq_var(dataset)
    # without the sequence there is no record key at all: the subject keys alone
    # would claim one record per subject, which an occurrence dataset never has
    if (length(seq_var) == 0) {
      return(character(0))
    }
    return(c(subject_keys, seq_var))
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
  # analysis variable with the record it was derived from, so it is only
  # `DTYPE` which separates them.
  # NOTE: no surrogate/sequence/identifier keys (ASEQ, SRCSEQ, ASPID, ...) --
  # see minimal_unique_key()
  extra <- c(
    "AVISITN", "AVISIT", "ATPTN", "ATPT",
    "NFRLT", "AFRLT",
    "ADTM", "ADT", "ASTDTM", "ASTDT", "AENDTM", "AENDT",
    "APERIOD", "APERIODC", "DTYPE"
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
#'
#' @details
#'   `ASEQ` is preferred, falling back to a populated SDTM domain sequence. If
#'   no sequence is left -- or several domain sequences are, making the record
#'   key ambiguous -- a warning is issued and nothing is returned.
#'
#' @return The name of the sequence variable which keys `dataset`, or a
#'   zero-length vector (with a warning) when none can be determined.
#'
#' @seealso [infer_admiral_keys()], which uses this for `OCCDS` datasets
#'
#' @keywords internal
#' @family internal
occds_seq_var <- function(dataset) {
  cols <- colnames(dataset)
  # a sequence which is present but entirely `NA` has not been derived yet;
  # treating it as the record key would make the key unique without meaning
  # anything
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
