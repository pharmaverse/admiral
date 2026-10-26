#' Construct an `admiral_df` Object
#'
#' Low-level constructor of the `admiral_df` class: adds the class to a data
#' frame (unless it is already present), preserving the existing classes such
#' as `tbl_df`/`data.frame`, and stores any metadata supplied. The arguments are
#' not checked, so this is what admiral's `derive_*()` and related
#' dataset-producing functions apply to their output. Users call
#' [set_admiral_attrs()], which validates its arguments before calling this.
#'
#' @param x A data frame (or `NULL`)
#' @param keys Character vector of key variables, stored in the
#'   `"admiral_keys"` attribute
#' @param dataset_name Character scalar, stored in the `"admiral_ds_name"`
#'   attribute
#' @param dataset_type Character scalar, stored in the `"admiral_ds_type"`
#'   attribute
#'
#' @details
#'   An argument left as `NULL` leaves the corresponding attribute untouched, so
#'   that tagging a dataset which already carries metadata does not remove it.
#'
#' @return
#'   If `x` is `NULL`, `NULL` is returned unchanged. Otherwise `x` with the
#'   supplied attributes set and `"admiral_df"` prepended to its class
#'   attribute.
#'
#' @keywords internal
#' @family internal
new_admiral_df <- function(x, keys = NULL, dataset_name = NULL, dataset_type = NULL) {
  if (is.null(x)) {
    return(x)
  }
  if (!is.null(keys)) {
    attr(x, "admiral_keys") <- keys
  }
  if (!is.null(dataset_name)) {
    attr(x, "admiral_ds_name") <- dataset_name
  }
  if (!is.null(dataset_type)) {
    attr(x, "admiral_ds_type") <- dataset_type
  }
  if (!inherits(x, "admiral_df")) {
    class(x) <- c("admiral_df", class(x))
  }
  x
}

# the ADaM dataset types `set_admiral_attrs()` accepts, and those of them which
# are medical device types
admiral_md_ds_types <- c("ADDL", "MDBDS", "MDOCCDS")
admiral_ds_types <- c("ADSL", "BDS", "OCCDS", "TTE", admiral_md_ds_types)

#' Set the Metadata of a Dataset
#'
#' @description
#' `r lifecycle::badge("experimental")`
#'
#' Tags a data frame with the `admiral_df` class and stores its metadata: the
#' key variables -- the variables which uniquely identify its records -- the
#' dataset name, and the ADaM dataset type. Call it on the finished dataset,
#' once it has been derived, e.g. just before summarising it: the metadata is
#' not guaranteed to survive the derivations which build the dataset (see
#' Details).
#'
#' The keys describe the record structure the dataset is intended to have; they
#' are not checked against its records when they are stored. `{admiral}` never
#' infers the keys of a dataset: if they have not been set with this function,
#' functions which need them ask for them to be set. The type, by contrast, is
#' inferred from the dataset if it has not been set (see
#' [guess_admiral_df_type()]), but a type set here always takes precedence.
#'
#' @param dataset A data frame. It does not need to have the `admiral_df` class
#'   already; the class is added if it is missing.
#'
#' @permitted a data frame
#'
#' @param keys The key variables of `dataset`, in key order, e.g.
#'   `exprs(USUBJID, PARAMCD, AVISIT)`, stored in the `"admiral_keys"`
#'   attribute. They must be variables of `dataset`: if any is not, a warning is
#'   issued and the keys are not set, so any keys set by a previous call are
#'   kept. An empty list, `exprs()`, records that the dataset has no key
#'   variables defined, which is distinct from its keys never having been set.
#'   `NULL` leaves keys stored by a previous call untouched.
#'
#' @permitted list of variables created by `exprs()` or `NULL`
#'
#' @param dataset_name The name of the dataset, e.g. `"ADVS"`, stored in the
#'   `"admiral_ds_name"` attribute. `NULL` leaves a name stored by a previous
#'   call untouched, so re-keying a dataset does not lose its name.
#'
#' @permitted character scalar or `NULL`
#'
#' @param dataset_type The ADaM dataset type, stored in the `"admiral_ds_type"`
#'   attribute. `"ADDL"`, `"MDBDS"`, and `"MDOCCDS"` are the medical device
#'   types (device-level, and the device versions of `"BDS"` and `"OCCDS"`);
#'   a device dataset is identified by `SPDEVID`, so if `dataset` does not
#'   contain it a warning is issued, but the type is stored all the same.
#'   `NULL` leaves a type stored by a previous call untouched.
#'
#' @permitted `"ADSL"`, `"BDS"`, `"OCCDS"`, `"TTE"`, `"ADDL"`, `"MDBDS"`,
#'   `"MDOCCDS"`, or `NULL`
#'
#' @return `dataset` with the `admiral_df` class, and the `"admiral_keys"`,
#'   `"admiral_ds_name"`, and `"admiral_ds_type"` attributes for those of
#'   `keys`, `dataset_name`, and `dataset_type` which were supplied (or had been
#'   stored by a previous call).
#'
#' @details
#'   The keys can come from anywhere: the `by_vars` of the derivation which
#'   produced the dataset, a dataset specification, or a hand-written list.
#'   A specification provides them as a character vector, which [chr2vars()]
#'   converts. `{admiral}` does not read specification objects itself; to use the
#'   keys defined in a `{metacore}` specification, extract them in the calling
#'   code, e.g. `chr2vars(metacore::get_keys(spec, "ADVS")$variable)`.
#'
#'   All `{admiral}` dataset functions tag their output with the `admiral_df`
#'   class, but they do not carry the metadata through. The core `{dplyr}` verbs
#'   (`mutate()`, `filter()`, `arrange()`, `select()`, `rename()`, `slice()`,
#'   `distinct()`, joins where `dataset` is `x`, ...) restore data frame
#'   attributes, so the metadata survives those. It is *not* preserved by
#'   operations which build a new data frame, in particular:
#'
#'   * `summarise()`,
#'   * `tidyr::pivot_longer()`/`pivot_wider()`, `nest()`/`unnest()`,
#'     `complete()`,
#'   * `bind_rows()` when `dataset` is not the *first* argument,
#'   * joins where `dataset` is the `y` (rather than the `x`) argument,
#'   * `ungroup()`, after `group_by()` or `rowwise()` (which themselves keep the
#'     attributes but drop the `admiral_df` class).
#'
#'   Many derivations use these internally, so metadata set at the start of a
#'   derivation pipeline is lost along the way. Set it once the dataset is
#'   finished, or re-apply this function if the dataset changes afterwards.
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
#' @caption Tagging a dataset with its metadata
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
#' advs <- set_admiral_attrs(
#'   advs,
#'   keys = exprs(USUBJID, PARAMCD, AVISIT),
#'   dataset_name = "ADVS",
#'   dataset_type = "BDS"
#' )
#'
#' attr(advs, "admiral_keys")
#' attr(advs, "admiral_ds_name")
#' attr(advs, "admiral_ds_type")
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
#'   set_admiral_attrs(
#'     advs,
#'     keys = chr2vars(metacore::get_keys(spec, "ADVS")$variable),
#'     dataset_name = "ADVS"
#'   )
#'   ```
#'
#' @code
#' advs <- set_admiral_attrs(
#'   advs,
#'   keys = chr2vars(c("USUBJID", "PARAMCD", "AVISIT"))
#' )
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
#' advs_no_keys <- set_admiral_attrs(advs, keys = exprs())
#'
#' attr(advs_no_keys, "admiral_keys")
set_admiral_attrs <- function(dataset,
                              keys = NULL,
                              dataset_name = NULL,
                              dataset_type = NULL) {
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
  assert_vars(keys, optional = TRUE)
  assert_data_frame(dataset)
  assert_character_scalar(dataset_name, optional = TRUE)
  assert_character_scalar(dataset_type, values = admiral_ds_types, optional = TRUE)

  # a declared type is not checked against the inferred one, but a device type
  # declared for a dataset without the device identifier is most likely the
  # wrong type (or the wrong dataset). Only a warning, as the declaration is the
  # user's to make
  if (!is.null(dataset_type) && dataset_type %in% admiral_md_ds_types &&
    !"SPDEVID" %in% colnames(dataset)) {
    cli_warn(c(
      paste(
        "Dataset type {.val {dataset_type}} is a medical device type, but",
        "{.arg dataset} does not contain {.var SPDEVID}."
      ),
      i = "Check that {.arg dataset_type} is correct."
    ))
  }

  if (!is.null(keys)) {
    # keys are stored as character: they are metadata about the dataset, not
    # expressions to evaluate against it. `unname()` because `vars2chr()`
    # returns a named vector, which would not compare `identical()` to a plain
    # one
    keys <- unname(vars2chr(keys))

    # keys which are not in the dataset are not stored, so that a typo is
    # caught here rather than when the attribute is read. A warning rather than
    # an error, as the rest of the call -- the name, type and class -- is still
    # valid
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
      keys <- NULL
    }
  }

  new_admiral_df(
    dataset,
    keys = keys,
    dataset_name = dataset_name,
    dataset_type = dataset_type
  )
}

#' Print `admiral_df` Objects
#'
#' @description
#' `r lifecycle::badge("experimental")`
#'
#' Prints the dataset name, type, and key variables stored by [set_admiral_attrs()]
#' in a header line, followed by the dataset as it would print without the
#' `admiral_df` class. The header is only printed if any of them have been
#' stored, so that the output of a derivation which merely carries the class
#' prints as before.
#'
#' If no type has been stored, the type inferred from the dataset (see
#' [guess_admiral_df_type()]) is shown instead, marked as inferred, provided it
#' is one of the ADaM dataset types.
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
#' set_admiral_attrs(
#'   advs,
#'   keys = exprs(USUBJID, PARAMCD, AVISIT),
#'   dataset_name = "ADVS",
#'   dataset_type = "BDS"
#' )
print.admiral_df <- function(x, ...) {
  name <- attr(x, "admiral_ds_name")
  type <- attr(x, "admiral_ds_type")
  keys <- attr(x, "admiral_keys")
  # the header is only printed if something has been stored, so that the output
  # of a derivation which merely carries the class prints as before -- and the
  # type is only inferred then, so that printing such output costs nothing extra
  if (is.null(name) && is.null(type) && is.null(keys)) {
    return(NextMethod())
  }
  if (is.null(type)) {
    inferred <- guess_admiral_df_type(x)
    # only a recognized type is shown: "other" and `NA` say nothing about the
    # dataset worth a place in the header
    if (inferred %in% admiral_ds_types) {
      type <- paste(inferred, "(inferred)")
    }
  }
  details <- c(
    if (!is.null(name)) paste0("name: ", name),
    if (!is.null(type)) paste0("type: ", type),
    # an empty attribute records that no keys are defined, so it is shown
    # rather than omitted as if the keys had never been set
    if (!is.null(keys)) {
      paste0("keys: ", if (length(keys) > 0) paste(keys, collapse = ", ") else "none")
    }
  )
  cat("# An admiral dataset (", paste(details, collapse = ", "), ")\n", sep = "")
  NextMethod()
}

#' Get the ADaM Dataset Type
#'
#' @description
#' `r lifecycle::badge("experimental")`
#'
#' Returns the type declared for a dataset with [set_admiral_attrs()] if there
#' is one, and otherwise the type inferred from its variables and records, see
#' [guess_admiral_df_type()].
#'
#' @param dataset A data frame
#'
#' @return A character scalar: the declared type, or the inferred one (which may
#'   be `"other"` or `NA`, see [guess_admiral_df_type()]).
#'
#' @keywords internal
#' @family internal
get_admiral_df_type <- function(dataset) {
  declared <- attr(dataset, "admiral_ds_type")
  if (!is.null(declared)) {
    return(declared)
  }
  guess_admiral_df_type(dataset)
}

#' Guess the ADaM Dataset Type
#'
#' @description
#' `r lifecycle::badge("experimental")`
#'
#' Attempts to classify a data frame as one of the common ADaM dataset types based
#' on the variables it contains. This is the fallback for when no type has been
#' declared with [set_admiral_attrs()]; use [get_admiral_df_type()] to get the
#' declared type where there is one.
#'
#' Note that this is the *type* of the dataset (`ADSL`, `BDS`, ...), not its
#' record structure: the latter is the set of key variables which uniquely
#' identify its records. Unlike the type, the keys are never inferred; they are
#' only ever declared with [set_admiral_attrs()].
#'
#' @param dataset A data frame
#'
#' @details
#'   A dataset is first classed as device data or subject data:
#'
#'   * **Device data** contains `SPDEVID`, the device identifier, and is not
#'     keyed by subject: it does not contain all the subject keys
#'     (`get_admiral_option("subject_keys")`, by default `STUDYID` and
#'     `USUBJID`), or they are missing for some of its records. Its type is one
#'     of the medical device types `"ADDL"`, `"MDBDS"`, and `"MDOCCDS"`.
#'   * **Subject data** contains all the subject keys. Its type is one of
#'     `"ADSL"`, `"BDS"`, `"OCCDS"`, and `"TTE"`.
#'
#'   A dataset which is neither is `"other"`. This keeps a data frame which merely
#'   shares variable names with an ADaM dataset -- a coding dictionary with
#'   `AETERM`/`AEDECOD` columns, a summary table with `PARAMCD` and `AVAL` --
#'   from being taken for one.
#'
#'   `SPDEVID` is not enough on its own to make a dataset device data, because
#'   subject datasets record it too, e.g. the device involved in an adverse
#'   event. A dataset whose records are all keyed by subject is therefore subject
#'   data whether or not it contains `SPDEVID`; a device dataset which also
#'   identifies the subject of every record is inferred as the corresponding
#'   subject type, so declare its type with [set_admiral_attrs()] instead.
#'
#'   Within each class the following precedence is used. The order matters,
#'   because the types are not mutually exclusive by variable content -- a `TTE`
#'   dataset contains `PARAMCD` just as a `BDS` one does, and an occurrence
#'   dataset which has had a parameter added is structurally both.
#'
#'   1. `"TTE"` -- contains `PARAMCD`, `CNSR`, and `STARTDT`. Tested first
#'      because a time-to-event dataset is a `BDS` dataset by variable content;
#'      only the censoring and origin-date variables distinguish it. There is no
#'      device time-to-event type, so for device data this is `"MDBDS"`.
#'   2. `"BDS"`/`"MDBDS"` -- contains `PARAMCD` and `AVAL` or `AVALC`.
#'   3. `"OCCDS"`/`"MDOCCDS"` -- does not contain `PARAMCD`, and contains a
#'      `--DECOD`/`--TERM` variable, an `AOCCxxFL` occurrence flag, or
#'      `TRTEMFL`.
#'   4. `"ADSL"`/`"ADDL"` -- does not contain `PARAMCD` and has one record per
#'      subject (per combination of the subject keys) or per device (per
#'      combination of `STUDYID`, where present, and `SPDEVID`).
#'   5. `"other"` -- none of the above.
#'
#'   The `--DECOD`/`--TERM` test requires the two-letter SDTM domain prefix the
#'   convention gives those variables (`AEDECOD`, `DETERM`, ...), so that a
#'   variable whose name merely ends in `TERM` is not taken as a signal.
#'
#'   `ADSL`/`ADDL` are the types decided by the records rather than the
#'   variables: treatment variables are no evidence of ADSL, because the first
#'   step of a typical `BDS`/`OCCDS` derivation is to merge them on from ADSL
#'   (see the `admiral` template scripts). A dataset with no records cannot show
#'   whether it has one record per subject or device, so its type is unknown
#'   (`NA`) unless its variables already decide it -- e.g. an empty dataset with
#'   the subject keys, `PARAMCD`, and `AVAL` is still `"BDS"`, and one with
#'   neither the subject keys nor `SPDEVID` is still `"other"`.
#'
#'   A dataset which has not yet been derived far enough to carry any of these
#'   signals is `"other"`.
#'
#' @return A character scalar: one of `"ADSL"`, `"BDS"`, `"OCCDS"`, `"TTE"`,
#'   `"ADDL"`, `"MDBDS"`, `"MDOCCDS"`, or `"other"`; or `NA` if the dataset has
#'   no records and its variables do not decide the type.
#'
#' @keywords internal
#' @family internal
guess_admiral_df_type <- function(dataset) {
  # a grouped or rowwise dataset is classified rather than refused: this is a
  # diagnostic
  assert_data_frame(dataset, check_is_grouped = FALSE, check_is_rowwise = FALSE)
  cols <- colnames(dataset)
  subject_keys <- vars2chr(get_admiral_option("subject_keys"))

  # All the subject keys are required, as only together do they identify a
  # subject (`STUDYID` alone identifies nothing)
  has_subject_keys <- all(subject_keys %in% cols)

  # `SPDEVID` is recorded in subject datasets too, so it only makes a dataset
  # device data where the records are not all keyed by subject
  if ("SPDEVID" %in% cols &&
    !(has_subject_keys && !any(is.na(dataset[subject_keys])))) {
    return(guess_structure_type(
      dataset,
      id_vars = intersect(c("STUDYID", "SPDEVID"), cols),
      types = c(tte = "MDBDS", bds = "MDBDS", occds = "MDOCCDS", level = "ADDL")
    ))
  }

  # without the subject keys -- a lookup table, a summary table -- a data frame
  # is none of the subject types, whatever variables it shares with them
  if (!has_subject_keys) {
    return("other")
  }

  guess_structure_type(
    dataset,
    id_vars = subject_keys,
    types = c(tte = "TTE", bds = "BDS", occds = "OCCDS", level = "ADSL")
  )
}

#' Guess the ADaM Dataset Structure from the Variables and Records
#'
#' Applies the precedence described in [guess_admiral_df_type()] for one class
#' of data, subject or device.
#'
#' @param dataset A data frame
#' @param id_vars Character vector of the variables identifying the unit of the
#'   class (the subject keys, or `STUDYID` and `SPDEVID`)
#' @param types Named character vector of the type to return for each
#'   structure: `tte`, `bds`, `occds`, and `level` (one record per unit)
#'
#' @return A character scalar: one of `types`, `"other"`, or `NA`.
#'
#' @keywords internal
#' @family internal
guess_structure_type <- function(dataset, id_vars, types) {
  cols <- colnames(dataset)

  # `PARAMCD` decides between the parameter-based types and the rest; TTE is
  # tested before BDS because a TTE dataset is a BDS dataset by variable content
  if ("PARAMCD" %in% cols) {
    if (all(c("CNSR", "STARTDT") %in% cols)) {
      return(types[["tte"]])
    }
    if (any(c("AVAL", "AVALC") %in% cols)) {
      return(types[["bds"]])
    }
    return("other")
  }

  if (any(str_detect(cols, "^[A-Z]{2}(DECOD|TERM)$|^AOCC[0-9A-Z]*FL$|^TRTEMFL$"))) {
    return(types[["occds"]])
  }

  # ADSL/ADDL is decided by the records: one record per subject or device. With
  # no records that cannot be shown either way -- the record and unit counts are
  # trivially equal -- so the type is unknown rather than "other"
  if (nrow(dataset) == 0) {
    return(NA_character_)
  }
  # `ungroup()` because `distinct()` silently adds the grouping variables, which
  # would test uniqueness within group rather than overall
  if (nrow(dataset) == nrow(distinct(ungroup(dataset), !!!syms(id_vars)))) {
    return(types[["level"]])
  }

  "other"
}
