#' Construct an `admiral_df` Object
#'
#' Low-level constructor of the `admiral_df` class: adds the class to a data
#' frame (unless it is already present), preserving the existing classes such
#' as `tbl_df`/`data.frame`, and stores any metadata supplied. The arguments are
#' not checked, so this is what admiral's `derive_*()` and related
#' dataset-producing functions apply to their output. Users call
#' [admiral_df()], which validates its arguments before calling this.
#'
#' @param x A data frame (or `NULL`)
#' @param dataset_keys Character vector of key variables, stored in the
#'   `"admiral_dataset_keys"` attribute
#' @param dataset_name Character scalar, stored in the `"admiral_dataset_name"`
#'   attribute
#' @param dataset_class Character scalar, stored in the `"admiral_dataset_class"`
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
new_admiral_df <- function(x, dataset_keys = NULL, dataset_name = NULL, dataset_class = NULL) {
  if (is.null(x)) {
    return(x)
  }
  if (!is.null(dataset_keys)) {
    attr(x, "admiral_dataset_keys") <- dataset_keys
  }
  if (!is.null(dataset_name)) {
    attr(x, "admiral_dataset_name") <- dataset_name
  }
  if (!is.null(dataset_class)) {
    attr(x, "admiral_dataset_class") <- dataset_class
  }
  if (!inherits(x, "admiral_df")) {
    class(x) <- c("admiral_df", class(x))
  }
  x
}

# the ADaM dataset classes `admiral_df()` accepts
admiral_dataset_classes <- c(
  "ADSL", "BDS", "OCCDS", "TTE", "ADDL", "MDBDS", "MDTTE", "MDOCCDS"
)

#' Create an `admiral_df` Dataset
#'
#' @description
#' `r lifecycle::badge("experimental")`
#'
#' Tags a data frame with the `admiral_df` class and stores its metadata: the
#' key variables -- the variables which uniquely identify its records -- the
#' dataset name, and the ADaM dataset class. Call it on the finished dataset,
#' once it has been derived, e.g. just before summarizing it: the metadata is
#' not guaranteed to survive the derivations which build the dataset (see
#' Details). The metadata is read back with [get_dataset_keys()],
#' [get_dataset_name()], and [get_dataset_class()].
#'
#' Called on a dataset which is already an `admiral_df`, it updates the
#' metadata rather than rebuilding it: an argument left as `NULL` keeps the
#' value stored by the previous call, so each item can be set on its own.
#'
#' The keys describe the record structure the dataset is intended to have; they
#' are not checked against its records when they are stored. `{admiral}` never
#' infers the keys of a dataset: if they have not been set with this function,
#' functions which need them ask for them to be set. The class, by contrast, is
#' inferred from the dataset if it has not been set (see
#' [infer_dataset_class()]), but a class set here always takes precedence.
#'
#' @param dataset A data frame. It does not need to have the `admiral_df` class
#'   already; the class is added if it is missing.
#'
#' @permitted a data frame
#'
#' @param dataset_keys The key variables of `dataset`, in key order, e.g.
#'   `exprs(USUBJID, PARAMCD, AVISIT)`, stored in the `"admiral_dataset_keys"`
#'   attribute. They must be variables of `dataset`: if any is not, a warning is
#'   issued and the keys are not set, so any keys set by a previous call are
#'   kept. An empty list, `exprs()`, records that the dataset has no key
#'   variables defined, which is distinct from its keys never having been set.
#'   `NULL` leaves keys stored by a previous call untouched.
#'
#' @permitted list of variables created by `exprs()` or `NULL`
#'
#' @param dataset_name The name of the dataset, e.g. `"ADVS"`, stored in the
#'   `"admiral_dataset_name"` attribute. There are no restrictions on its
#'   casing: it is stored exactly as given. `NULL` leaves a name stored by a
#'   previous call untouched, so re-keying a dataset does not lose its name.
#'
#' @permitted character scalar or `NULL`
#'
#' @param dataset_class The ADaM dataset class, stored in the
#'   `"admiral_dataset_class"` attribute. This is the ADaM class of the
#'   dataset, not its R class (which is `admiral_df`). It must be given in upper
#'   case, as listed below. `"TTE"` is the time-to-event subclass of `"BDS"`.
#'   `"ADDL"`, `"MDBDS"`, `"MDTTE"`, and `"MDOCCDS"` are the medical device
#'   classes (device-level, and the device versions of `"BDS"`, `"TTE"`, and
#'   `"OCCDS"`). The class is not checked against the variables of `dataset`.
#'   `NULL` leaves a class stored by a previous call untouched.
#'
#' @permitted `"ADSL"`, `"BDS"`, `"OCCDS"`, `"TTE"`, `"ADDL"`, `"MDBDS"`,
#'   `"MDTTE"`, `"MDOCCDS"`, or `NULL`
#'
#' @return `dataset` with the `admiral_df` class, and the `"admiral_dataset_keys"`,
#'   `"admiral_dataset_name"`, and `"admiral_dataset_class"` attributes for those of
#'   `dataset_keys`, `dataset_name`, and `dataset_class` which were supplied (or had been
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
#' advs <- admiral_df(
#'   advs,
#'   dataset_keys = exprs(USUBJID, PARAMCD, AVISIT),
#'   dataset_name = "ADVS",
#'   dataset_class = "BDS"
#' )
#'
#' get_dataset_keys(advs)
#' get_dataset_name(advs)
#' get_dataset_class(advs)
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
#'   admiral_df(
#'     advs,
#'     dataset_keys = chr2vars(metacore::get_keys(spec, "ADVS")$variable),
#'     dataset_name = "ADVS"
#'   )
#'   ```
#'
#' @code
#' advs <- admiral_df(
#'   advs,
#'   dataset_keys = chr2vars(c("USUBJID", "PARAMCD", "AVISIT"))
#' )
#'
#' # `dataset_name` was not repeated, so the name from the previous call stands
#' get_dataset_name(advs)
#'
#' @caption Declaring that a dataset has no key variables
#'
#' @info An empty `dataset_keys` is stored all the same. It is not the same as leaving
#'   the attribute unset: it records that no key variables are defined.
#'
#' @code
#' advs_no_keys <- admiral_df(advs, dataset_keys = exprs())
#'
#' get_dataset_keys(advs_no_keys)
admiral_df <- function(dataset,
                       dataset_keys = NULL,
                       dataset_name = NULL,
                       dataset_class = NULL) {
  assert_vars(dataset_keys, optional = TRUE)
  assert_data_frame(dataset)
  assert_character_scalar(dataset_name, optional = TRUE)
  assert_character_scalar(dataset_class, values = admiral_dataset_classes, optional = TRUE)

  if (!is.null(dataset_keys)) {
    # keys are stored as character: they are metadata about the dataset, not
    # expressions to evaluate against it. `unname()` because `vars2chr()`
    # returns a named vector, which would not compare `identical()` to a plain
    # one
    dataset_keys <- unname(vars2chr(dataset_keys))

    # keys which are not in the dataset are not stored, so that a typo is
    # caught here rather than when the attribute is read. A warning rather than
    # an error, as the rest of the call -- the name, ADaM class and R class -- is still
    # valid
    missing_keys <- setdiff(dataset_keys, colnames(dataset))
    if (length(missing_keys) > 0) {
      previous_keys <- attr(dataset, "admiral_dataset_keys")
      cli_warn(c(
        paste(
          "{qty(missing_keys)}Key variable{?s} {.var {missing_keys}}",
          "{qty(missing_keys)}{?is/are} not in {.arg dataset}, so no keys were set."
        ),
        i = if (length(previous_keys) > 0) {
          "The keys set previously are kept: {.var {previous_keys}}."
        }
      ))
      dataset_keys <- NULL
    }
  }

  new_admiral_df(
    dataset,
    dataset_keys = dataset_keys,
    dataset_name = dataset_name,
    dataset_class = dataset_class
  )
}

#' Print `admiral_df` Objects
#'
#' @description
#' `r lifecycle::badge("experimental")`
#'
#' Prints the dataset name, ADaM class, and key variables stored by [admiral_df()] in
#' a header line, followed by the dataset as it would print without the
#' `admiral_df` class. The header is only printed if any of them have been
#' stored, so that the output of a derivation which merely carries the class
#' prints as before.
#'
#' If no class has been stored, the class inferred from the dataset (see
#' [infer_dataset_class()]) is shown instead, marked as inferred, provided it
#' is one of the ADaM dataset classes.
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
#' admiral_df(
#'   advs,
#'   dataset_keys = exprs(USUBJID, PARAMCD, AVISIT),
#'   dataset_name = "ADVS",
#'   dataset_class = "BDS"
#' )
print.admiral_df <- function(x, ...) {
  dataset_name <- attr(x, "admiral_dataset_name")
  dataset_class <- attr(x, "admiral_dataset_class")
  dataset_keys <- attr(x, "admiral_dataset_keys")
  # the header is only printed if something has been stored, so that the output
  # of a derivation which merely carries the class prints as before -- and the
  # class is only inferred then, so that printing such output costs nothing extra
  if (is.null(dataset_name) && is.null(dataset_class) && is.null(dataset_keys)) {
    return(NextMethod())
  }
  if (is.null(dataset_class)) {
    inferred <- infer_dataset_class(x)
    # only a recognized class is shown: "other" and `NA` say nothing about the
    # dataset worth a place in the header
    if (inferred %in% admiral_dataset_classes) {
      dataset_class <- paste(inferred, "(inferred)")
    }
  }
  details <- c(
    if (!is.null(dataset_name)) paste0("name: ", dataset_name),
    if (!is.null(dataset_class)) paste0("class: ", dataset_class),
    # an empty attribute records that no keys are defined, so it is shown
    # rather than omitted as if the keys had never been set
    if (!is.null(dataset_keys)) {
      paste0(
        "keys: ",
        if (length(dataset_keys) > 0) paste(dataset_keys, collapse = ", ") else "none"
      )
    }
  )
  cat("# An admiral dataset (", paste(details, collapse = ", "), ")\n", sep = "")
  NextMethod()
}

#' Get the Key Variables of a Dataset
#'
#' @description
#' `r lifecycle::badge("experimental")`
#'
#' Returns the key variables stored for a dataset with [admiral_df()].
#'
#' @param dataset A data frame
#'
#' @permitted a data frame
#'
#' @return The key variables as a list of symbols, in key order, so that they
#'   can be passed to a `by_vars` argument or back to [admiral_df()]; an empty
#'   list if the dataset has been declared to have no key variables; or `NULL` if
#'   its keys have not been set. The keys are never inferred.
#'
#' @keywords utils_help
#' @family utils_help
#'
#' @export
#'
#' @examples
#' library(tibble)
#'
#' advs <- tribble(
#'   ~USUBJID, ~PARAMCD, ~AVISIT,    ~AVAL,
#'   "1",      "DIABP",  "BASELINE",    51,
#'   "1",      "SYSBP",  "BASELINE",   121
#' )
#'
#' get_dataset_keys(advs)
#'
#' advs <- admiral_df(advs, dataset_keys = exprs(USUBJID, PARAMCD, AVISIT))
#' get_dataset_keys(advs)
get_dataset_keys <- function(dataset) {
  assert_data_frame(dataset, check_is_grouped = FALSE, check_is_rowwise = FALSE)
  dataset_keys <- attr(dataset, "admiral_dataset_keys")
  if (is.null(dataset_keys)) {
    return(NULL)
  }
  chr2vars(dataset_keys)
}

#' Get the Name of a Dataset
#'
#' @description
#' `r lifecycle::badge("experimental")`
#'
#' Returns the dataset name stored for a dataset with [admiral_df()].
#'
#' @param dataset A data frame
#'
#' @permitted a data frame
#'
#' @return A character scalar, or `NULL` if no name has been set.
#'
#' @keywords utils_help
#' @family utils_help
#'
#' @export
#'
#' @examples
#' library(tibble)
#'
#' advs <- tribble(
#'   ~USUBJID, ~PARAMCD, ~AVAL,
#'   "1",      "DIABP",     51
#' )
#'
#' get_dataset_name(advs)
#'
#' advs <- admiral_df(advs, dataset_name = "ADVS")
#' get_dataset_name(advs)
get_dataset_name <- function(dataset) {
  assert_data_frame(dataset, check_is_grouped = FALSE, check_is_rowwise = FALSE)
  attr(dataset, "admiral_dataset_name")
}

#' Get the ADaM Dataset Class
#'
#' @description
#' `r lifecycle::badge("experimental")`
#'
#' Returns the ADaM class declared for a dataset with [admiral_df()] if there is
#' one, and otherwise the class inferred from its variables and records, see
#' [infer_dataset_class()].
#'
#' @param dataset A data frame
#'
#' @permitted a data frame
#'
#' @return A character scalar: the declared class, or the inferred one (which may
#'   be `"other"` or `NA`, see [infer_dataset_class()]).
#'
#' @keywords utils_help
#' @family utils_help
#'
#' @export
#'
#' @examples
#' library(tibble)
#'
#' advs <- tribble(
#'   ~STUDYID, ~USUBJID, ~PARAMCD, ~AVAL,
#'   "P",      "1",      "DIABP",     51,
#'   "P",      "1",      "SYSBP",    121
#' )
#'
#' # no class has been declared, so it is inferred
#' get_dataset_class(advs)
#'
#' # a declared class takes precedence
#' advs <- admiral_df(advs, dataset_class = "TTE")
#' get_dataset_class(advs)
get_dataset_class <- function(dataset) {
  declared <- attr(dataset, "admiral_dataset_class")
  if (!is.null(declared)) {
    return(declared)
  }
  infer_dataset_class(dataset)
}

#' Infer the ADaM Dataset Class
#'
#' @description
#' Attempts to classify a data frame as one of the common ADaM dataset classes
#' from its variables and records. It ignores any class declared with
#' [admiral_df()]: [get_dataset_class()] returns the declared class if there is
#' one, and only falls back to this function otherwise.
#'
#' Note that this is the *class* of the dataset (`ADSL`, `BDS`, ...), not its
#' record structure: the latter is the set of key variables which uniquely
#' identify its records. Unlike the class, the keys are never inferred; they are
#' only ever declared with [admiral_df()].
#'
#' @param dataset A data frame
#'
#' @details
#'   A dataset is first classed as device data or subject data:
#'
#'   * **Device data** contains `SPDEVID`, the device identifier, populated
#'     (neither `NA` nor `""`) for every record. Its class is one of the medical
#'     device classes `"ADDL"`, `"MDBDS"`, `"MDTTE"`, and `"MDOCCDS"`.
#'   * **Subject data** is otherwise data which contains all the subject keys
#'     (`get_admiral_option("subject_keys")`, by default `STUDYID` and
#'     `USUBJID`). Its class is one of `"ADSL"`, `"BDS"`, `"OCCDS"`, and `"TTE"`.
#'
#'   A dataset which is neither is `"other"`. This keeps a data frame which merely
#'   shares variable names with an ADaM dataset -- a coding dictionary with
#'   `AETERM`/`AEDECOD` columns, a summary table with `PARAMCD` and `AVAL` --
#'   from being taken for one.
#'
#'   The subject keys do not decide the class: in the ADaMIG-MD `STUDYID` is
#'   required in a device dataset and `USUBJID` is conditional, so a device
#'   dataset may identify the subject of every record. Nor does the presence of
#'   `SPDEVID`, as subject datasets record it too, e.g. the device involved in
#'   an adverse event; but there it is populated only for the records which
#'   involve a device. A subject dataset in which every record involves a device
#'   (e.g. an `ADAE` restricted to device-related events) is therefore inferred
#'   as device data, so declare its class with [admiral_df()] instead.
#'
#'   Within each class the following precedence is used. The order matters,
#'   because the classes are not mutually exclusive by variable content -- a `TTE`
#'   dataset contains `PARAMCD` just as a `BDS` one does, and an occurrence
#'   dataset which has had a parameter added is structurally both.
#'
#'   1. `"TTE"`/`"MDTTE"` -- contains `PARAMCD`, `CNSR`, and `STARTDT`. Tested
#'      first because a time-to-event dataset is a `BDS` dataset by variable
#'      content; only the censoring and origin-date variables distinguish it.
#'   2. `"BDS"`/`"MDBDS"` -- contains `PARAMCD` and `AVAL` or `AVALC`.
#'   3. `"OCCDS"`/`"MDOCCDS"` -- does not contain `PARAMCD`, and contains a
#'      `--DECOD`/`--TERM` variable with a two-letter SDTM domain prefix
#'      (`AEDECOD`, `DETERM`, ...; a variable merely ending in `TERM` does not
#'      count), an `AOCCxxFL` occurrence flag, or `TRTEMFL`.
#'   4. `"ADSL"`/`"ADDL"` -- does not contain `PARAMCD` and has one record per
#'      subject (per combination of the subject keys) or per device (per
#'      combination of `STUDYID` and `SPDEVID`). This is decided by the records:
#'      treatment variables are no evidence of ADSL, as `BDS`/`OCCDS`
#'      derivations merge them on from ADSL.
#'   5. `"other"` -- none of the above, if there is at least one record.
#'   6. `NA` -- none of the above.
#'
#' @return A character scalar: one of `"ADSL"`, `"BDS"`, `"OCCDS"`, `"TTE"`,
#'   `"ADDL"`, `"MDBDS"`, `"MDTTE"`, `"MDOCCDS"`, or `"other"`; or `NA` if the
#'   dataset has no records and its variables do not decide the class.
#'
#' @keywords internal
#' @family internal
infer_dataset_class <- function(dataset) {
  # a grouped or rowwise dataset is classified rather than refused: this is a
  # diagnostic
  assert_data_frame(dataset, check_is_grouped = FALSE, check_is_rowwise = FALSE)
  cols <- colnames(dataset)
  subject_keys <- vars2chr(get_admiral_option("subject_keys"))

  # All the subject keys are required, as only together do they identify a
  # subject (`STUDYID` alone identifies nothing)
  has_subject_keys <- all(subject_keys %in% cols)

  # `SPDEVID` is recorded in subject datasets too, but there only for the
  # records involving a device, so it makes a dataset device data only where it
  # is populated for every record. The subject keys do not decide it, as
  # `USUBJID` is conditional in a device dataset. A blank counts as missing, as
  # data read from SAS has `""` rather than `NA`
  if ("SPDEVID" %in% cols &&
    !any(is.na(dataset$SPDEVID) | dataset$SPDEVID == "")) {
    return(infer_class_by_structure(
      dataset,
      id_vars = intersect(c("STUDYID", "SPDEVID"), cols),
      classes = c(tte = "MDTTE", bds = "MDBDS", occds = "MDOCCDS", level = "ADDL")
    ))
  }

  # without the subject keys -- a lookup table, a summary table -- a data frame
  # is none of the subject classes, whatever variables it shares with them
  if (!has_subject_keys) {
    return("other")
  }

  infer_class_by_structure(
    dataset,
    id_vars = subject_keys,
    classes = c(tte = "TTE", bds = "BDS", occds = "OCCDS", level = "ADSL")
  )
}

#' Infer the ADaM Dataset Class of Subject or Device Data
#'
#' Applies the precedence described in [infer_dataset_class()] for one class
#' of data, subject or device.
#'
#' @param dataset A data frame
#' @param id_vars Character vector of the variables identifying the unit of the
#'   class (the subject keys, or `STUDYID` and `SPDEVID`)
#' @param classes Named character vector of the class to return for each
#'   structure: `tte`, `bds`, `occds`, and `level` (one record per unit)
#'
#' @return A character scalar: one of `classes`, `"other"`, or `NA`.
#'
#' @keywords internal
#' @family internal
infer_class_by_structure <- function(dataset, id_vars, classes) {
  cols <- colnames(dataset)

  # `PARAMCD` decides between the parameter-based classes and the rest; TTE is
  # tested before BDS because a TTE dataset is a BDS dataset by variable content
  if ("PARAMCD" %in% cols) {
    if (all(c("CNSR", "STARTDT") %in% cols)) {
      return(classes[["tte"]])
    }
    if (any(c("AVAL", "AVALC") %in% cols)) {
      return(classes[["bds"]])
    }
    # "other" is a conclusion drawn from the records, so a dataset without any
    # is unknown instead
    return(if (nrow(dataset) == 0) NA_character_ else "other")
  }

  if (any(str_detect(cols, "^[A-Z]{2}(DECOD|TERM)$|^AOCC[0-9A-Z]*FL$|^TRTEMFL$"))) {
    return(classes[["occds"]])
  }

  # ADSL/ADDL is decided by the records: one record per subject or device. With
  # no records that cannot be shown either way -- the record and unit counts are
  # trivially equal -- so the class is unknown rather than "other"
  if (nrow(dataset) == 0) {
    return(NA_character_)
  }
  # `ungroup()` because `distinct()` silently adds the grouping variables, which
  # would test uniqueness within group rather than overall
  if (nrow(dataset) == nrow(distinct(ungroup(dataset), !!!syms(id_vars)))) {
    return(classes[["level"]])
  }

  "other"
}
