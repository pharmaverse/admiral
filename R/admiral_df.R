#' Tag a Dataset with the `admiral_df` Class
#'
#' Adds the `admiral_df` class to a data frame (unless it is already
#' present), preserving the existing classes such as `tbl_df`/`data.frame`.
#' This is the tagging primitive that admiral's `derive_*()` and related
#' dataset-producing functions apply to their output, so that a dataset
#' built with admiral can be recognized as such downstream, e.g. by
#' [summary.admiral_df()].
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
#' Stores the key variables of a dataset -- the variables which uniquely identify
#' its records -- in its `"admiral_keys"` attribute, and tags it with
#' the `admiral_df` class (see [as_admiral_df()]). The keys describe the record
#' structure the dataset is intended to have; they are not checked against its
#' records when they are stored.
#'
#' @param dataset A data frame. It does not need to have the `admiral_df` class
#'   already; the class is added if it is missing.
#'
#' @permitted a data frame
#'
#' @param keys The key variables of `dataset`, in key order, e.g.
#'   `exprs(USUBJID, PARAMCD, AVISIT)`. They must be variables of `dataset`:
#'   if any is not, a warning is issued and the keys are not set, so any keys
#'   set by a previous call are kept. An empty list, `exprs()`, records that
#'   the dataset has no key variables defined, which is distinct from its keys
#'   never having been set.
#'
#' @permitted list of variables created by `exprs()`
#'
#' @param dataset_name The name of the dataset, e.g. `"ADVS"`, stored in the
#'   `"admiral_ds_name"` attribute. `NULL` leaves a name stored by a previous
#'   call untouched, so re-keying a dataset does not lose its name.
#'
#' @permitted character scalar or `NULL`
#'
#' @return `dataset` with an `"admiral_keys"` attribute (unless a key is not in
#'   `dataset`), an `"admiral_ds_name"` attribute if `dataset_name` was supplied
#'   or had been stored by a previous call, and the `admiral_df` class.
#'
#' @details
#'   The keys can come from anywhere: the `by_vars` of the derivation which
#'   produced the dataset, a dataset specification, or a hand-written list.
#'   A specification provides them as a character vector, which [chr2vars()]
#'   converts. `{admiral}` does not read specification objects itself; to use the
#'   keys defined in a `{metacore}` specification, extract them in the calling
#'   code, e.g. `chr2vars(metacore::get_keys(spec, "ADVS")$variable)`.
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
#'   `group_by()` and `rowwise()` are a case of their own: they keep the
#'   `"admiral_keys"` and `"admiral_ds_name"` attributes but drop the
#'   `admiral_df` class.
#'
#'   Call this function afterwards (or re-apply it as needed) if the dataset
#'   passes through any of those.
#'
#'   The keys are stored as plain variable names, so they are not updated by
#'   `rename()` or removed by `select()`: a key which is no longer in the
#'   dataset means it has changed since its keys were set.
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
#' attr(advs, "admiral_ds_name")
#' class(advs)
#'
#' @caption Keys from a dataset specification
#'
#' @info A dataset specification provides the keys as a character vector, which
#'   [chr2vars()] converts. `{admiral}` does not read specification objects
#'   itself, so they are extracted in the calling code; e.g. with a
#'   `{metacore}` specification that is `metacore::get_keys()`:
#'
#'   ```r
#'   set_admiral_keys(advs, chr2vars(metacore::get_keys(spec, "ADVS")$variable), "ADVS")
#'   ```
#'
#' @code
#' advs <- set_admiral_keys(advs, keys = chr2vars(c("USUBJID", "PARAMCD", "AVISIT")))
#'
#' # `dataset_name` was not repeated, so the name from the previous call stands
#' attr(advs, "admiral_ds_name")
#'
#' @caption Declaring that a dataset has no key variables
#'
#' @info An empty `keys` is stored all the same. It is not the same as leaving
#'   the attribute unset: it records that no key variables are defined.
#'
#' @code
#' advs_no_keys <- set_admiral_keys(advs, keys = exprs())
#'
#' attr(advs_no_keys, "admiral_keys")
set_admiral_keys <- function(dataset, keys, dataset_name = NULL) {
  # keys taken from a specification are a character vector, so that case gets
  # a pointer to the conversion rather than the generic "must be a list"
  if (is.character(keys)) {
    cli_abort(c(
      paste(
        "Argument {.arg keys} must be a list of variables created by {.fun exprs},",
        "but is {.obj_type_friendly {keys}}."
      ),
      i = paste(
        "Use {.fun chr2vars} to convert a character vector, e.g.",
        "{.code chr2vars(metacore::get_keys(spec, \"ADVS\")$variable)}."
      )
    ))
  }
  assert_vars(keys)
  assert_data_frame(dataset)
  assert_character_scalar(dataset_name, optional = TRUE)

  # keys are stored as character: they are metadata about the dataset, not
  # expressions to evaluate against it. `unname()` because `vars2chr()` returns
  # a named vector, which would not compare `identical()` to a plain one
  keys <- unname(vars2chr(keys))

  # keys which are not in the dataset are not stored, so that a typo is caught
  # here rather than when the attribute is read. A warning rather than an error,
  # as the rest of the call -- the name and the class -- is still valid
  missing_keys <- setdiff(keys, colnames(dataset))
  if (length(missing_keys) > 0) {
    previous_keys <- attr(dataset, "admiral_keys")
    cli_warn(c(
      paste(
        "{qty(missing_keys)}Key variable{?s} {.var {missing_keys}}",
        "{qty(missing_keys)}{?is/are} not in {.arg dataset}, so no keys were set."
      ),
      i = if (length(previous_keys) > 0) {
        "The keys set previously are kept: {.var {previous_keys}}."
      }
    ))
  } else {
    attr(dataset, "admiral_keys") <- keys
  }
  if (!is.null(dataset_name)) {
    attr(dataset, "admiral_ds_name") <- dataset_name
  }
  as_admiral_df(dataset)
}

#' Print `admiral_df` Objects
#'
#' Prints the dataset name and key variables stored by [set_admiral_keys()] in a
#' header line, followed by the dataset as it would print without the
#' `admiral_df` class. The header is only printed if a name or keys have been
#' stored, so that the output of a derivation which merely carries the class
#' (see [as_admiral_df()]) prints as before.
#'
#' The header also shows the ADaM dataset type (`ADSL`, `BDS`, `OCCDS`, or
#' `TTE`) if it can be determined from the dataset, see
#' [get_admiral_df_type()].
#'
#' @param x An `admiral_df` object
#' @param ... Passed on to the next print method
#'
#' @return `x`, invisibly
#'
#' @keywords internal
#' @family utils_print
#'
#' @export
#'
#' @examples
#' library(tibble)
#'
#' advs <- tribble(
#'   ~STUDYID, ~USUBJID, ~PARAMCD, ~AVISIT,    ~AVAL,
#'   "P",      "1",      "DIABP",  "BASELINE",    51,
#'   "P",      "1",      "SYSBP",  "BASELINE",   121,
#'   "P",      "2",      "DIABP",  "BASELINE",    79
#' )
#'
#' set_admiral_keys(advs, keys = exprs(USUBJID, PARAMCD, AVISIT), dataset_name = "ADVS")
print.admiral_df <- function(x, ...) {
  name <- attr(x, "admiral_ds_name")
  keys <- attr(x, "admiral_keys")
  # the type is only worked out when there is a header to add it to, so that
  # printing the output of a derivation which merely carries the class costs
  # nothing extra
  if (!is.null(name) || !is.null(keys)) {
    type <- get_admiral_df_type(x)
    details <- c(
      if (!is.null(name)) paste0("name: ", name),
      # only a recognized type is shown: "other" and `NA` say nothing about
      # the dataset worth a place in the header
      if (type %in% c("ADSL", "BDS", "OCCDS", "TTE")) paste0("type: ", type),
      # an empty attribute records that no keys are defined, so it is shown
      # rather than omitted as if the keys had never been set
      if (!is.null(keys)) {
        paste0("keys: ", if (length(keys) > 0) paste(keys, collapse = ", ") else "none")
      }
    )
    cat("# An admiral dataset (", paste(details, collapse = ", "), ")\n", sep = "")
  }
  NextMethod()
}

#' Determine the ADaM Dataset Type
#'
#' Attempts to classify a data frame as one of the common ADaM dataset types based
#' on the variables it contains.
#'
#' Note that this is the *type* of the dataset (`ADSL`, `BDS`, ...), not its
#' record structure: the latter is the set of key variables which uniquely
#' identify its records, see [set_admiral_keys()].
#'
#' @param dataset A data frame
#'
#' @details
#'   Every ADaM dataset type is subject data, so a dataset which does not contain
#'   all the subject keys (`get_admiral_option("subject_keys")`, by default
#'   `STUDYID` and `USUBJID`) is `"other"`. This keeps a data frame which merely
#'   shares variable names with an ADaM dataset -- a coding dictionary with
#'   `AETERM`/`AEDECOD` columns, a summary table with `PARAMCD` and `AVAL` --
#'   from being taken for one.
#'
#'   Otherwise the following precedence is used. The order matters, because the
#'   types are not mutually exclusive by variable content -- a `TTE` dataset
#'   contains `PARAMCD` just as a `BDS` one does, and an occurrence dataset which
#'   has had a parameter added is structurally both.
#'
#'   1. `"TTE"` -- contains `PARAMCD`, `CNSR`, and `STARTDT`. Tested first
#'      because a time-to-event dataset is a `BDS` dataset by variable content;
#'      only the censoring and origin-date variables distinguish it.
#'   2. `"OCCDS"` -- does not contain `PARAMCD`, and contains a `--DECOD`/`--TERM`
#'      variable, an `AOCCxxFL` occurrence flag, or `TRTEMFL`.
#'   3. `"BDS"` -- contains `PARAMCD` and `AVAL` or `AVALC`.
#'   4. `"ADSL"` -- does not contain `PARAMCD` and has one record per subject,
#'      i.e. per combination of the subject keys.
#'   5. `"other"` -- none of the above.
#'
#'   The `--DECOD`/`--TERM` test requires the two-letter SDTM domain prefix the
#'   convention gives those variables (`AEDECOD`, `CMTERM`, ...), so that a
#'   variable whose name merely ends in `TERM` is not taken as a signal.
#'
#'   `ADSL` is the one type decided by the records rather than the variables:
#'   treatment variables are no evidence of it, because the first step of a
#'   typical `BDS`/`OCCDS` derivation is to merge them on from ADSL (see the
#'   `admiral` template scripts). A dataset with no records cannot show whether
#'   it has one record per subject, so its type is unknown (`NA`) unless its
#'   variables already decide it -- e.g. an empty dataset with `PARAMCD` and
#'   `AVAL` is still `"BDS"`, and one without the subject keys is still
#'   `"other"`.
#'
#'   A dataset which has not yet been derived far enough to carry any of these
#'   signals is `"other"`.
#'
#' @return A character scalar: one of `"ADSL"`, `"BDS"`, `"OCCDS"`, `"TTE"`, or
#'   `"other"`; or `NA` if the dataset has no records and its variables do not
#'   decide the type.
#'
#' @keywords internal
#' @family internal
get_admiral_df_type <- function(dataset) {
  # a grouped or rowwise dataset is classified rather than refused: this is a
  # diagnostic
  assert_data_frame(dataset, check_is_grouped = FALSE, check_is_rowwise = FALSE)
  cols <- colnames(dataset)
  subject_keys <- vars2chr(get_admiral_option("subject_keys"))

  # every ADaM dataset type is subject data, so a data frame without the subject
  # keys -- a lookup table, a summary table -- is none of them, whatever
  # variables it shares with them. All the keys are required, as only together
  # do they identify a subject (`STUDYID` alone identifies nothing)
  if (!all(subject_keys %in% cols)) {
    return("other")
  }

  # `PARAMCD` decides between the parameter-based types and the rest; TTE is
  # tested before BDS because a TTE dataset is a BDS dataset by variable content
  if ("PARAMCD" %in% cols) {
    if (all(c("CNSR", "STARTDT") %in% cols)) {
      return("TTE")
    }
    if (any(c("AVAL", "AVALC") %in% cols)) {
      return("BDS")
    }
    return("other")
  }

  if (any(str_detect(cols, "^[A-Z]{2}(DECOD|TERM)$|^AOCC[0-9A-Z]*FL$|^TRTEMFL$"))) {
    return("OCCDS")
  }

  # ADSL is decided by the records: one record per subject. With no records
  # that cannot be shown either way -- the record and subject counts are
  # trivially equal -- so the type is unknown rather than "other"
  if (nrow(dataset) == 0) {
    return(NA_character_)
  }
  # `ungroup()` because `distinct()` silently adds the grouping variables, which
  # would test uniqueness within group rather than overall
  if (nrow(dataset) == nrow(distinct(ungroup(dataset), !!!syms(subject_keys)))) {
    return("ADSL")
  }

  "other"
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
#' @examples
#' library(tibble)
#'
#' advs <- tribble(
#'   ~USUBJID, ~PARAMCD, ~AVISIT,    ~ATPT,     ~AVAL,
#'   "1",      "SYSBP",  "BASELINE", "MORNING",   121,
#'   "1",      "SYSBP",  "BASELINE", "EVENING",   125,
#'   "1",      "SYSBP",  "WEEK 2",   "MORNING",   118,
#'   "1",      "SYSBP",  "WEEK 2",   "EVENING",   120
#' )
#'
#' # candidates are added in order until the key is unique
#' admiral:::minimal_unique_key(
#'   advs,
#'   must_have = c("USUBJID", "PARAMCD"),
#'   optional = c("AVISIT", "ATPT")
#' )
#'
#' # an added candidate which turns out not to be needed is dropped again:
#' # with one timepoint only, `ATPT` distinguishes nothing
#' admiral:::minimal_unique_key(
#'   advs[advs$ATPT == "MORNING", ],
#'   must_have = c("USUBJID", "PARAMCD"),
#'   optional = c("ATPT", "AVISIT")
#' )
#'
#' # if uniqueness is never reached, the full set is returned; the duplicated
#' # record shows up as duplicates of that key
#' admiral:::minimal_unique_key(
#'   rbind(advs, advs[1, ]),
#'   must_have = c("USUBJID", "PARAMCD"),
#'   optional = c("AVISIT", "ATPT")
#' )
#'
#' @keywords internal
#' @family internal
minimal_unique_key <- function(dataset, must_have, optional) {
  # `ungroup()` for the reason given in `get_admiral_df_type()`
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
#' Works out which variables uniquely identify the records of a dataset. This is
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
#'   [get_admiral_df_type()] does, reduced to those of its variables which
#'   actually distinguish subjects in the dataset at hand. With the default
#'   option that is `USUBJID` alone for a single-study dataset, because
#'   `STUDYID` is constant there and would only overstate the record structure;
#'   a pooled dataset whose subject identifiers repeat across studies keeps
#'   `STUDYID` as well.
#'
#'   A dataset with no records cannot show which variables uniquely identify its
#'   records -- every candidate key is trivially unique -- so nothing is
#'   inferred for it.
#'
#' @return A character vector of inferred key variable names, or a zero-length
#'   vector when no plausible record structure can be determined.
#'
#' @examples
#' library(tibble)
#'
#' # BDS: subject and parameter, plus the visit needed for uniqueness. `STUDYID`
#' # is constant in a single study, so it is left out
#' advs <- tribble(
#'   ~STUDYID, ~USUBJID, ~PARAMCD, ~AVISIT,    ~AVAL,
#'   "P",      "P-1",    "SYSBP",  "BASELINE",   121,
#'   "P",      "P-1",    "SYSBP",  "WEEK 2",     118,
#'   "P",      "P-1",    "DIABP",  "BASELINE",    79,
#'   "P",      "P-2",    "SYSBP",  "BASELINE",   130
#' )
#' admiral:::infer_admiral_keys(advs)
#'
#' # OCCDS: keyed by its sequence variable
#' adae <- tribble(
#'   ~STUDYID, ~USUBJID, ~AEDECOD,   ~AESEQ,
#'   "P",      "P-1",    "HEADACHE",      1,
#'   "P",      "P-1",    "NAUSEA",        2,
#'   "P",      "P-2",    "HEADACHE",      1
#' )
#' admiral:::infer_admiral_keys(adae)
#'
#' # ADSL: one record per subject
#' adsl <- tribble(
#'   ~STUDYID, ~USUBJID, ~AGE,
#'   "P",      "P-1",      63,
#'   "P",      "P-2",      71
#' )
#' admiral:::infer_admiral_keys(adsl)
#'
#' # nothing is inferred for a dataset whose type is not recognized
#' admiral:::infer_admiral_keys(tibble(FOO = 1:2))
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

#' Summarize an `admiral_df` Dataset
#'
#' Reports what a dataset contains and whether it has one record per key
#' variable combination. This is useful after a derivation, e.g. to check that
#' [derive_param_computed()] added the records that were expected.
#'
#' The summary always shows the dataset name and type, the number of subjects,
#' observations and variables, and the record structure. It then adds the
#' information for the dataset type (see [get_admiral_df_type()]):
#'
#' * `ADSL`: subjects per treatment arm and the population flag counts.
#' * `BDS`: records and subjects per parameter, and the analysis visits.
#' * `OCCDS`: the number of distinct terms, and the treatment-emergent vs total
#'   records.
#' * `TTE`: events and censored per parameter.
#'
#' @param object An `admiral_df` object
#' @param ... Not used
#'
#' @details
#'   The summary reports facts only, apart from the record structure check.
#'   Each item is only reported when the variables it needs are in the dataset,
#'   so a dataset which is still being derived gives a shorter summary rather
#'   than an error.
#'
#'   The record structure is checked against the keys stored by
#'   [set_admiral_keys()]. If no keys were stored, they are inferred from the
#'   dataset type and its data (see [infer_admiral_keys()]); inferred keys
#'   describe the structure the data happens to have, which is not necessarily
#'   the structure it was meant to have.
#'
#'   Subjects are counted by `USUBJID` if it is present, as the other subject
#'   keys (e.g. `STUDYID`) may not be populated on records added by a
#'   derivation. Otherwise the subject keys
#'   (`get_admiral_option("subject_keys")`) are used.
#'
#'   The dataset name is the one stored by [set_admiral_keys()] or, failing
#'   that, the name of the variable `summary()` was called on.
#'
#' @return An object of class `summary_admiral_df`, a named list which is
#'   printed by [print.summary_admiral_df()].
#'
#' @keywords utils_ds_chk
#' @family utils_ds_chk
#'
#' @export
#'
#' @examplesx
#'
#' @caption Summarizing a BDS dataset after a derivation
#'
#' @code
#' library(tibble)
#' library(dplyr, warn.conflicts = FALSE)
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
#' advs <- derive_param_computed(
#'   advs,
#'   by_vars = exprs(STUDYID, USUBJID, AVISIT),
#'   parameters = c("SYSBP", "DIABP"),
#'   set_values_to = exprs(
#'     AVAL = (AVAL.SYSBP + 2 * AVAL.DIABP) / 3,
#'     PARAMCD = "MAP"
#'   )
#' )
#'
#' summary(advs)
#'
#' @caption Checking against declared keys
#'
#' @info The record structure is checked against the keys stored by
#'   [set_admiral_keys()] rather than inferred ones.
#'
#' @code
#' advs <- set_admiral_keys(
#'   advs,
#'   keys = exprs(USUBJID, PARAMCD, AVISIT),
#'   dataset_name = "ADVS"
#' )
#'
#' summary(advs)
summary.admiral_df <- function(object, ...) {
  name <- attr(object, "admiral_ds_name")
  if (is.null(name)) {
    # the variable `summary()` was called on, unless it is a pipe placeholder
    # or a more complex expression than a name
    obj_expr <- substitute(object)
    if (is.name(obj_expr) && !as.character(obj_expr) %in% c(".", "_")) {
      name <- as.character(obj_expr)
    }
  }
  declared_keys <- attr(object, "admiral_keys")

  # a plain, ungrouped data frame, so that the tables derived from it do not
  # carry the class
  class(object) <- setdiff(class(object), "admiral_df")
  dataset <- ungroup(object)
  cols <- colnames(dataset)
  type <- get_admiral_df_type(dataset)

  subject_vars <- if ("USUBJID" %in% cols) {
    "USUBJID"
  } else {
    intersect(vars2chr(get_admiral_option("subject_keys")), cols)
  }

  out <- list(
    name = name,
    type = type,
    n_subjects = if (length(subject_vars) > 0) {
      n_distinct(dataset[subject_vars])
    } else {
      NA_integer_
    },
    n_obs = nrow(dataset),
    n_vars = ncol(dataset),
    subject_vars = subject_vars,
    structure = summarize_structure(dataset, declared_keys, type)
  )

  # the type is `NA` for a dataset with no records whose variables do not
  # decide it, which has nothing type-specific to report
  out$details <- switch(coalesce(type, "other"),
    ADSL = summarize_adsl(dataset, subject_vars),
    BDS = summarize_bds(dataset, subject_vars),
    OCCDS = summarize_occds(dataset),
    TTE = summarize_tte(dataset),
    list()
  )

  class(out) <- c("summary_admiral_df", "list")
  out
}

#' Check the Record Structure of a Dataset
#'
#' @param dataset A data frame
#' @param declared_keys The keys stored by [set_admiral_keys()], or `NULL`
#' @param type The dataset type, see [get_admiral_df_type()]
#'
#' @return A list with the elements `keys`, `source` (`"declared"` or
#'   `"inferred"`), `missing_keys` (declared keys which are not in `dataset`),
#'   and `n_duplicates` (the number of records beyond one per key, `NA` if
#'   the structure was not checked).
#'
#' @keywords internal
#' @family internal
summarize_structure <- function(dataset, declared_keys, type) {
  if (!is.null(declared_keys)) {
    keys <- declared_keys
    source <- "declared"
  } else {
    keys <- infer_admiral_keys(dataset, type = coalesce(type, "other"))
    source <- "inferred"
  }

  # declared keys are plain names, so a `rename()` or `select()` since they
  # were set can leave them pointing at nothing; the check is skipped rather
  # than run on part of the key
  missing_keys <- setdiff(keys, colnames(dataset))
  n_duplicates <- if (length(keys) > 0 && length(missing_keys) == 0) {
    nrow(dataset) - n_distinct(dataset[keys])
  } else {
    NA_integer_
  }

  list(
    keys = keys,
    source = source,
    missing_keys = missing_keys,
    n_duplicates = n_duplicates
  )
}

#' Summarize the Type-Specific Content of a Dataset
#'
#' Collect the information [summary.admiral_df()] reports for a dataset of
#' type `ADSL`, `BDS`, `OCCDS`, or `TTE`. An element is only returned when
#' the variables it needs are in the dataset.
#'
#' @param dataset A data frame
#' @param subject_vars The variables subjects are counted by
#'
#' @return A named list:
#'
#'   * `summarize_adsl()`: `trt_var` and `trt` (subjects per treatment arm),
#'     `pops` (named count of `"Y"` per population flag).
#'   * `summarize_bds()`: `params` (records and subjects per `PARAMCD`),
#'     `visits` (the analysis visits, in `AVISITN` order if available).
#'   * `summarize_occds()`: `terms` (named count of distinct values per coding
#'     variable), `trtem` (treatment-emergent and total records).
#'   * `summarize_tte()`: `events` (events and censored per `PARAMCD`).
#'
#' @keywords internal
#' @family internal
#' @rdname summarize_admiral_df_type
summarize_adsl <- function(dataset, subject_vars) {
  out <- list()

  # planned treatment first, falling back to actual treatment and `ARM`; a
  # variable which is entirely `NA` has not been derived yet
  trt_var <- intersect(c("TRT01P", "TRT01A", "ARM"), colnames(dataset))
  trt_var <- keep(trt_var, ~ any(!is.na(dataset[[.x]])))
  if (length(trt_var) > 0) {
    out$trt_var <- trt_var[1]
    out$trt <- dataset %>%
      distinct(!!sym(trt_var[1]), !!!syms(subject_vars)) %>%
      count(!!sym(trt_var[1]), name = "subjects")
  }

  pop_vars <- intersect(
    c("RANDFL", "ENRLFL", "SAFFL", "ITTFL", "FASFL", "EFFFL", "PPROTFL", "COMPLFL"),
    colnames(dataset)
  )
  if (length(pop_vars) > 0) {
    out$pops <- vapply(
      set_names(pop_vars),
      function(v) sum(dataset[[v]] == "Y", na.rm = TRUE),
      integer(1)
    )
  }

  out
}

#' @rdname summarize_admiral_df_type
summarize_bds <- function(dataset, subject_vars) {
  out <- list()
  cols <- colnames(dataset)

  if (nrow(dataset) > 0) {
    params <- count(dataset, PARAMCD, name = "records")
    if (length(subject_vars) > 0) {
      params <- left_join(
        params,
        dataset %>%
          distinct(PARAMCD, !!!syms(subject_vars)) %>%
          count(PARAMCD, name = "subjects"),
        by = "PARAMCD"
      )
    }
    out$params <- params
  }

  if ("AVISIT" %in% cols && any(!is.na(dataset$AVISIT))) {
    visits <- filter(dataset, !is.na(AVISIT))
    visits <- if ("AVISITN" %in% cols) {
      arrange(distinct(visits, AVISITN, AVISIT), AVISITN, AVISIT)
    } else {
      arrange(distinct(visits, AVISIT), AVISIT)
    }
    out$visits <- unique(visits$AVISIT)
  }

  out
}

#' @rdname summarize_admiral_df_type
summarize_occds <- function(dataset) {
  out <- list()
  cols <- colnames(dataset)

  # the same `--BODSYS`/`--DECOD`/`--TERM` naming `get_admiral_df_type()` relies
  # on, from the highest coding level to the reported term
  term_vars <- c(
    str_subset(cols, "^[A-Z]{2}BODSYS$"),
    str_subset(cols, "^[A-Z]{2}DECOD$"),
    str_subset(cols, "^[A-Z]{2}TERM$")
  )
  if (length(term_vars) > 0) {
    out$terms <- vapply(
      set_names(term_vars),
      function(v) n_distinct(dataset[[v]], na.rm = TRUE),
      integer(1)
    )
  }

  if ("TRTEMFL" %in% cols) {
    out$trtem <- c(
      emergent = sum(dataset$TRTEMFL == "Y", na.rm = TRUE),
      total = nrow(dataset)
    )
  }

  out
}

#' @rdname summarize_admiral_df_type
summarize_tte <- function(dataset) {
  out <- list()

  # each parameter is a separate endpoint, so a single total would mislead
  if (nrow(dataset) > 0) {
    out$events <- dataset %>%
      group_by(PARAMCD) %>%
      summarise(
        events = sum(CNSR == 0, na.rm = TRUE),
        censored = sum(CNSR != 0, na.rm = TRUE)
      )
  }

  out
}

#' Print a Summary of an `admiral_df` Dataset
#'
#' @param x A `summary_admiral_df` object created by [summary.admiral_df()]
#' @param ... Not used
#'
#' @return `x`, invisibly
#'
#' @keywords internal
#' @family utils_print
#'
#' @export
#'
#' @examples
#' library(tibble)
#'
#' adsl <- tribble(
#'   ~STUDYID, ~USUBJID, ~TRT01P,   ~SAFFL,
#'   "P",      "P-1",    "Placebo", "Y",
#'   "P",      "P-2",    "Drug A",  "Y",
#'   "P",      "P-3",    "Drug A",  "N"
#' )
#'
#' print(summary(set_admiral_keys(adsl, keys = exprs(USUBJID), dataset_name = "ADSL")))
print.summary_admiral_df <- function(x, ...) {
  type <- if (is.na(x$type)) "unknown type" else x$type
  title <- if (is.null(x$name)) type else paste0(x$name, " (", type, ")")
  cat("# Summary of an admiral dataset: ", title, "\n", sep = "")

  dims <- c(
    if (!is.na(x$n_subjects)) paste0("Subjects: ", x$n_subjects),
    paste0("Observations: ", x$n_obs),
    paste0("Variables: ", x$n_vars)
  )
  cat(paste(dims, collapse = " | "), "\n", sep = "")
  cat("Record structure: ", format_structure(x$structure), "\n", sep = "")

  d <- x$details
  if (!is.null(d[["trt"]])) {
    cat(
      "Subjects per treatment arm (", d[["trt_var"]], "): ",
      paste0(d[["trt"]][[d[["trt_var"]]]], ": ", d[["trt"]]$subjects, collapse = " | "), "\n",
      sep = ""
    )
  }
  if (!is.null(d[["pops"]])) {
    cat("Population flags (\"Y\"): ", format_counts(d[["pops"]]), "\n", sep = "")
  }
  if (!is.null(d[["params"]])) {
    cat("Parameters (", nrow(d[["params"]]), "):\n", sep = "")
    print_summary_table(d[["params"]])
  }
  if (!is.null(d[["visits"]])) {
    cat(
      "Analysis visits (", length(d[["visits"]]), "): ",
      paste(d[["visits"]], collapse = ", "), "\n",
      sep = ""
    )
  }
  if (!is.null(d[["terms"]])) {
    cat("Distinct terms: ", format_counts(d[["terms"]]), "\n", sep = "")
  }
  if (!is.null(d[["trtem"]])) {
    cat(
      "Treatment-emergent records (TRTEMFL): ", d[["trtem"]][["emergent"]], " of ",
      d[["trtem"]][["total"]], "\n",
      sep = ""
    )
  }
  if (!is.null(d[["events"]])) {
    cat("Events and censored per parameter:\n")
    print_summary_table(d[["events"]])
  }

  invisible(x)
}

#' Format the Record Structure of a Summary
#'
#' @param structure The `structure` element of a `summary_admiral_df` object,
#'   see [summarize_structure()]
#'
#' @return A character scalar
#'
#' @keywords internal
#' @family internal
format_structure <- function(structure) {
  keys <- paste(structure$keys, collapse = ", ")
  source <- paste0(" (", structure$source, " keys)")
  if (length(structure$missing_keys) > 0) {
    paste0(
      "not checked, key variable(s) ", paste(structure$missing_keys, collapse = ", "),
      " not in the dataset", source
    )
  } else if (length(structure$keys) == 0) {
    # empty declared keys record that no keys are defined, which is not the
    # same as being unable to infer any
    if (structure$source == "declared") "no key variables defined" else "could not be determined"
  } else if (structure$n_duplicates == 0) {
    paste0("one record per ", keys, source)
  } else {
    paste0(
      "not one record per ", keys, ", ", structure$n_duplicates,
      " duplicate record", if (structure$n_duplicates > 1) "s", source
    )
  }
}

#' Format Named Counts as `name n | name n`
#'
#' @param counts A named numeric vector
#'
#' @return A character scalar
#'
#' @keywords internal
#' @family internal
format_counts <- function(counts) {
  paste(names(counts), counts, collapse = " | ")
}

#' Print a Table of a Summary Indented and Without Row Names
#'
#' @param tbl A data frame
#'
#' @return No return value, called for side effects
#'
#' @keywords internal
#' @family internal
print_summary_table <- function(tbl) {
  lines <- capture.output(print(as.data.frame(tbl), row.names = FALSE))
  cat(paste0("  ", lines), sep = "\n")
}
