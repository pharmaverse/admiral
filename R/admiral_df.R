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
#' @param keys The declared key variables of `dataset`, see
#'   [set_admiral_keys()], or `NULL` if none were declared
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
#'   4. `"ADSL"` -- does not contain `PARAMCD` and is one record per subject,
#'      i.e. per combination of the subject keys. If `keys` are declared, they
#'      decide this: `"ADSL"` if they are all subject keys. Otherwise the
#'      records decide it.
#'   5. `"other"` -- none of the above.
#'
#'   The `--DECOD`/`--TERM` test requires the two-letter SDTM domain prefix the
#'   convention gives those variables (`AEDECOD`, `CMTERM`, ...), so that a
#'   variable whose name merely ends in `TERM` is not taken as a signal.
#'
#'   `ADSL` is the one type decided by the records rather than the variables:
#'   treatment variables are no evidence of it, because the first step of a
#'   typical `BDS`/`OCCDS` derivation is to merge them on from ADSL (see the
#'   `admiral` template scripts). Declared keys take precedence over the
#'   records because they state the intended structure: an ADSL with duplicated
#'   records is still an ADSL, and the duplicates are what
#'   [summary.admiral_df()] reports. Without declared keys, a dataset with no
#'   records cannot show whether it has one record per subject, so its type is
#'   unknown (`NA`) unless its variables already decide it -- e.g. an empty
#'   dataset with `PARAMCD` and `AVAL` is still `"BDS"`, and one without the
#'   subject keys is still `"other"`.
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
get_admiral_df_type <- function(dataset, keys = attr(dataset, "admiral_keys")) {
  # a grouped or rowwise dataset is classified rather than refused: this is a
  # diagnostic
  assert_data_frame(dataset, check_is_grouped = FALSE, check_is_rowwise = FALSE)
  assert_character_vector(keys, optional = TRUE)
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

  # declared keys state the intended record structure, so they decide ADSL
  # rather than the records: duplicated subject records are a defect of an
  # ADSL, not evidence that it is something else. Empty keys declare no
  # structure, so the records decide as if none were declared
  if (length(keys) > 0) {
    return(if (all(keys %in% subject_keys)) "ADSL" else "other")
  }

  # otherwise ADSL is decided by the records: one record per subject. With no
  # records that cannot be shown either way -- the record and subject counts are
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
#' * `ADSL`: subjects per treatment arm, the population flag counts, the subject
#'   flow (randomized, treated, completed), the treatment duration, the end of
#'   study status, the deaths and their causes, and the age and sex of the
#'   subjects.
#' * `BDS`: records, subjects, visits, and the distribution of `AVAL` per
#'   parameter, the analysis visits, and the derived records per `DTYPE`.
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
  # the keys are passed explicitly as `ungroup()` need not keep the attribute
  type <- get_admiral_df_type(dataset, keys = declared_keys)

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
#' @details
#'   The ADSL items other than the treatment arms and population flags are
#'   taken from the first record of each subject, so that duplicated subject
#'   records -- which the record structure check reports -- do not inflate them.
#'
#' @return A named list:
#'
#'   * `summarize_adsl()`: `trt_var` and `trt` (subjects per treatment arm),
#'     `pops` (named count of subjects with `"Y"` per population flag), `flow`
#'     (subjects randomized, treated, and completed), `trtdur` (median, minimum,
#'     and maximum of `TRTDURD`), `eosstt` (subjects per end of study status),
#'     `n_deaths` and `dthcaus` (subjects who died, and per cause of death),
#'     `age` (median, minimum, and maximum of `AGE`), `sex` (subjects per
#'     `SEX`).
#'   * `summarize_bds()`: `params` (per `PARAMCD`: the `PARAM` values, the
#'     records, subjects, and analysis visits, and the missing values, minimum,
#'     median, and maximum of a numeric `AVAL`), `visits` (the analysis visits,
#'     in `AVISITN` order if available), `dtype` (records per `DTYPE`).
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
  # subjects rather than records are counted, as the treatment arms are, so
  # that duplicated subject records do not inflate the counts
  if (length(pop_vars) > 0) {
    out$pops <- vapply(
      set_names(pop_vars),
      function(v) n_distinct(dataset[dataset[[v]] %in% "Y", subject_vars]),
      integer(1)
    )
  }

  # one record per subject for the remaining items, so that each counts
  # subjects, and so that a subject's age or duration is only counted once
  subjects <- distinct(dataset, !!!syms(subject_vars), .keep_all = TRUE)
  cols <- colnames(subjects)

  out$flow <- summarize_subject_flow(subjects)

  if ("TRTDURD" %in% cols) {
    out$trtdur <- summarize_distribution(subjects$TRTDURD)
  }
  if ("EOSSTT" %in% cols && any(!is.na(subjects$EOSSTT))) {
    out$eosstt <- count_values(subjects$EOSSTT)
  }

  # the causes only among the deaths, so that the subjects alive do not show up
  # as a missing cause
  if ("DTHFL" %in% cols && any(!is.na(subjects$DTHFL))) {
    deaths <- filter(subjects, DTHFL %in% "Y")
    out$n_deaths <- nrow(deaths)
    if ("DTHCAUS" %in% cols && nrow(deaths) > 0) {
      out$dthcaus <- count_values(deaths$DTHCAUS)
    }
  }

  if ("AGE" %in% cols) {
    out$age <- summarize_distribution(subjects$AGE)
  }
  if ("SEX" %in% cols && any(!is.na(subjects$SEX))) {
    out$sex <- count_values(subjects$SEX)
  }

  out
}

#' @rdname summarize_admiral_df_type
summarize_bds <- function(dataset, subject_vars) {
  out <- list()
  cols <- colnames(dataset)

  # per parameter, as that is where most BDS defects show: a parameter with too
  # few subjects or visits, values which were not derived, implausible values
  if (nrow(dataset) > 0) {
    stats <- exprs(records = n())
    if ("PARAM" %in% cols) {
      # all values rather than the first, so that a `PARAMCD` with more than one
      # `PARAM` is visible
      stats$PARAM <- expr(paste(unique(PARAM[!is.na(PARAM)]), collapse = " / "))
    }
    if (length(subject_vars) > 0) {
      stats$subjects <- expr(n_distinct(pick(all_of(!!subject_vars))))
    }
    if ("AVISIT" %in% cols) {
      stats$visits <- expr(n_distinct(AVISIT, na.rm = TRUE))
    }
    if ("AVAL" %in% cols && is.numeric(dataset$AVAL)) {
      stats$missing <- expr(sum(is.na(AVAL)))
      stats$min <- expr(if (all(is.na(AVAL))) NA_real_ else min(AVAL, na.rm = TRUE))
      stats$median <- expr(median(AVAL, na.rm = TRUE))
      stats$max <- expr(if (all(is.na(AVAL))) NA_real_ else max(AVAL, na.rm = TRUE))
    }
    out$params <- dataset %>%
      group_by(PARAMCD) %>%
      summarise(!!!stats, .groups = "drop") %>%
      select(PARAMCD, any_of("PARAM"), everything())
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

  # what the derivations added, by type (`AVERAGE`, `LOCF`, ...)
  if ("DTYPE" %in% cols && any(!is.na(dataset$DTYPE))) {
    out$dtype <- count_values(dataset$DTYPE[!is.na(dataset$DTYPE)])
  }

  out
}

#' Count the Subjects Randomized, Treated, and Completed
#'
#' @param subjects A data frame with one record per subject
#'
#' @return A named integer vector with any of the elements `randomized` (from
#'   `RANDDT`, or `RANDFL` if not available), `treated` (from `TRTSDT`, or
#'   `TRTSDTM`), and `completed` (from `EOSSTT`), depending on the variables in
#'   `subjects`; `NULL` if there are none.
#'
#' @keywords internal
#' @family internal
summarize_subject_flow <- function(subjects) {
  cols <- colnames(subjects)
  # how far subjects got: the differences to the subject count are the screen
  # failures, the subjects randomized but not treated, and the dropouts
  c(
    if ("RANDDT" %in% cols) {
      c(randomized = sum(!is.na(subjects$RANDDT)))
    } else if ("RANDFL" %in% cols) {
      c(randomized = sum(subjects$RANDFL %in% "Y"))
    },
    if ("TRTSDT" %in% cols) {
      c(treated = sum(!is.na(subjects$TRTSDT)))
    } else if ("TRTSDTM" %in% cols) {
      c(treated = sum(!is.na(subjects$TRTSDTM)))
    },
    if ("EOSSTT" %in% cols) {
      c(completed = sum(subjects$EOSSTT %in% "COMPLETED"))
    }
  )
}

#' Summarize the Distribution of a Numeric Variable
#'
#' @param x A numeric vector
#'
#' @return A named numeric vector with the elements `median`, `min`, and `max`,
#'   or `NULL` if `x` has no non-missing values.
#'
#' @keywords internal
#' @family internal
summarize_distribution <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) == 0) {
    return(NULL)
  }
  c(median = median(x), min = min(x), max = max(x))
}

#' Count the Values of a Vector
#'
#' @param x A vector
#'
#' @return A named integer vector of counts per value, in value order with `NA`
#'   (named `"NA"`) last.
#'
#' @keywords internal
#' @family internal
count_values <- function(x) {
  counts <- table(x, useNA = "ifany")
  set_names(as.integer(counts), coalesce(names(counts), "NA"))
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
  # the type is what tailors the summary, so it heads the title; "other" and an
  # unknown type would read oddly there ("other summary"), hence the fallback
  title <- if (x$type %in% c("ADSL", "BDS", "OCCDS", "TTE")) {
    paste(x$type, "summary")
  } else {
    "Dataset summary"
  }
  # a name which only repeats the type ("ADSL: ADSL summary") is left out
  if (!is.null(x$name) && !identical(x$name, x$type)) {
    title <- paste0(x$name, ": ", title)
  }
  cli_rule(left = title)

  if (!is.na(x$n_subjects)) {
    subject_vars <- paste(x$subject_vars, collapse = ", ")
    cli_text("Subjects ({subject_vars}): {x$n_subjects}")
  }
  cli_text("Observations: {x$n_obs} | Variables: {x$n_vars}")
  print_summary_structure(x$structure)

  d <- x$details
  print_adsl_details(d)
  print_bds_details(d)
  if (!is.null(d[["terms"]])) {
    terms <- format_counts(d[["terms"]])
    cli_text("Distinct terms: {terms}")
  }
  if (!is.null(d[["trtem"]])) {
    cli_text(
      "Treatment emergent (TRTEMFL): {d[['trtem']][['emergent']]} of
       {d[['trtem']][['total']]} records"
    )
  }
  if (!is.null(d[["events"]])) {
    cli_text("Events and censored per parameter (CNSR):")
    print_summary_table(d[["events"]])
  }

  invisible(x)
}

#' Print the Type-Specific Part of a Summary
#'
#' Print the lines of [print.summary_admiral_df()] for the type-specific items
#' collected by [summarize_adsl()] and [summarize_bds()]. Each line is only
#' printed if its item is present.
#'
#' @param d The `details` element of a `summary_admiral_df` object
#'
#' @return No return value, called for side effects
#'
#' @keywords internal
#' @family internal
#' @rdname print_admiral_df_details
print_adsl_details <- function(d) {
  if (!is.null(d[["trt"]])) {
    trt <- format_counts(set_names(d[["trt"]]$subjects, d[["trt"]][[d[["trt_var"]]]]), ": ")
    cli_text("Subjects per treatment arm ({d[['trt_var']]}): {trt}")
  }
  if (!is.null(d[["pops"]])) {
    pops <- format_counts(d[["pops"]])
    cli_text("Populations: {pops}")
  }
  if (!is.null(d[["flow"]])) {
    flow <- format_counts(d[["flow"]])
    cli_text("Subject flow: {flow}")
  }
  if (!is.null(d[["trtdur"]])) {
    trtdur <- format_distribution(d[["trtdur"]])
    cli_text("Treatment duration (TRTDURD, days): {trtdur}")
  }
  if (!is.null(d[["eosstt"]])) {
    eosstt <- format_counts(d[["eosstt"]], ": ")
    cli_text("End of study (EOSSTT): {eosstt}")
  }
  if (!is.null(d[["n_deaths"]])) {
    dthcaus <- if (!is.null(d[["dthcaus"]])) {
      paste0(" (", format_counts(d[["dthcaus"]], ": "), ")")
    } else {
      ""
    }
    cli_text("Deaths (DTHFL): {d[['n_deaths']]}{dthcaus}")
  }
  if (!is.null(d[["age"]]) || !is.null(d[["sex"]])) {
    demographics <- paste(
      c(
        if (!is.null(d[["age"]])) paste("AGE", format_distribution(d[["age"]])),
        if (!is.null(d[["sex"]])) paste("SEX", format_counts(d[["sex"]], ": "))
      ),
      collapse = " | "
    )
    cli_text("Demographics: {demographics}")
  }

  invisible(NULL)
}

#' @rdname print_admiral_df_details
print_bds_details <- function(d) {
  if (!is.null(d[["params"]])) {
    cli_text("Parameters ({nrow(d[['params']])}):")
    # `PARAM` is left out as it would make the table too wide for the console;
    # it is in the summary object
    print_summary_table(
      select(d[["params"]], -any_of("PARAM")),
      max_rows = 10,
      what = "parameter"
    )
  }
  if (!is.null(d[["visits"]])) {
    visits <- paste(d[["visits"]], collapse = ", ")
    cli_text("Analysis visits (AVISIT): {visits}")
  }
  if (!is.null(d[["dtype"]])) {
    dtype <- format_counts(d[["dtype"]], ": ")
    cli_text("Derived records (DTYPE): {dtype}")
  }

  invisible(NULL)
}

#' Print the Record Structure of a Summary
#'
#' Prints the record structure check as a bullet: a tick if the dataset has one
#' record per key, a cross with the number of extra records if it does not, and
#' a warning sign if the check could not be run.
#'
#' @param structure The `structure` element of a `summary_admiral_df` object,
#'   see [summarize_structure()]
#'
#' @return No return value, called for side effects
#'
#' @keywords internal
#' @family internal
print_summary_structure <- function(structure) {
  keys <- paste(structure$keys, collapse = ", ")
  # declared keys are what the dataset is meant to have, so only inferred ones
  # are qualified
  src <- if (structure$source == "inferred") " (inferred)" else ""
  n_extra <- structure$n_duplicates
  if (length(structure$missing_keys) > 0) {
    missing_keys <- paste(structure$missing_keys, collapse = ", ")
    cli_bullets(c("!" = "Structure: not checked,
       {qty(length(structure$missing_keys))}declared key{?s} {missing_keys}
       {?is/are} not in the dataset"))
  } else if (length(structure$keys) == 0) {
    # empty declared keys record that no keys are defined, which is not the
    # same as being unable to infer any
    if (structure$source == "declared") {
      cli_bullets(c(i = "Structure: no key variables defined"))
    } else {
      cli_bullets(c(i = "Structure: could not be determined"))
    }
  } else if (n_extra == 0) {
    cli_bullets(c(v = "Structure{src}: one record per {keys}"))
  } else {
    cli_bullets(c(x = "Structure{src}: expected one record per {keys}, found
       {n_extra} extra record{?s}"))
  }

  invisible(NULL)
}

#' Format a Number for a Summary
#'
#' @param x A numeric scalar
#'
#' @return A character scalar: a whole number in full (rounding a count would
#'   misreport it), any other number to 4 significant digits, and `NA` as `.`.
#'
#' @keywords internal
#' @family internal
format_number <- function(x) {
  if (is.na(x)) {
    "."
  } else {
    format(if (x %% 1 == 0) x else signif(x, 4), scientific = FALSE)
  }
}

#' Format a Distribution as `median m (range min-max)`
#'
#' @param stats A named numeric vector created by [summarize_distribution()]
#'
#' @return A character scalar
#'
#' @keywords internal
#' @family internal
format_distribution <- function(stats) {
  paste0(
    "median ", format_number(stats[["median"]]),
    " (range ", format_number(stats[["min"]]), "-", format_number(stats[["max"]]), ")"
  )
}

#' Format Named Counts as `name n | name n`
#'
#' @param counts A named numeric vector
#' @param sep The separator between a name and its count
#'
#' @return A character scalar
#'
#' @keywords internal
#' @family internal
format_counts <- function(counts, sep = " ") {
  paste0(names(counts), sep, counts, collapse = " | ")
}

#' Print a Table of a Summary as Aligned Text
#'
#' Prints a header line and one line per row, indented, with character columns
#' left-aligned and numeric columns right-aligned. Fractional numbers are
#' shown to 4 significant digits, whole numbers in full, and missing values as
#' `.`.
#'
#' @param tbl A data frame
#' @param max_rows The maximum number of rows to print. The number of rows left
#'   out is reported below the table.
#' @param what What a row is, used to report the rows left out
#'
#' @return No return value, called for side effects
#'
#' @keywords internal
#' @family internal
print_summary_table <- function(tbl, max_rows = Inf, what = "row") {
  n_more <- nrow(tbl) - max_rows
  tbl <- tbl[seq_len(min(nrow(tbl), max_rows)), , drop = FALSE]

  # formatted per value rather than per column, as `format()` on a column pads
  # every value to the decimals of the most precise one
  cells <- lapply(tbl, function(col) {
    out <- if (is.numeric(col)) {
      # whole numbers (the counts) are shown in full: rounding them would
      # misreport them
      vapply(col, format_number, character(1))
    } else {
      as.character(col)
    }
    out[is.na(col)] <- "."
    out
  })
  right <- vapply(tbl, is.numeric, logical(1))
  widths <- vapply(
    names(cells),
    function(nm) max(nchar(c(nm, cells[[nm]]), type = "width")),
    integer(1)
  )
  line <- function(txt_by_col) {
    # `formatC()` takes a single width for character input, so one column at a
    # time
    padded <- vapply(
      seq_along(txt_by_col),
      function(j) formatC(txt_by_col[[j]], width = if (right[[j]]) widths[[j]] else -widths[[j]]),
      character(1)
    )
    paste0("  ", paste(padded, collapse = "  "))
  }

  # `cli_verbatim()` because `cli_text()` collapses the alignment spaces
  cli_verbatim(line(names(cells)))
  for (i in seq_len(nrow(tbl))) {
    cli_verbatim(line(vapply(cells, `[[`, character(1), i)))
  }
  if (n_more > 0) {
    # the plural is built here as cli would take the quantity from `what`
    what <- if (n_more > 1) paste0(what, "s") else what
    cli_text("... and {n_more} more {what}")
  }

  invisible(NULL)
}
