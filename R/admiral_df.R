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

#' Summarize an `admiral_df` Dataset
#'
#' Reports what a dataset contains and whether it has one record per key
#' variable combination. This is useful after a derivation, e.g. to check that
#' [derive_param_computed()] added the records that were expected.
#'
#' The summary always shows the dataset name and ADaM class, the number of
#' subjects, observations and variables, and the record structure. It then adds
#' the information for the dataset class (see [get_dataset_class()]):
#'
#' * `ADSL`: subjects per treatment arm, the population flag counts, the subject
#'   flow (randomized, treated, completed), the treatment duration, the end of
#'   study status, the deaths and their causes, and the age and sex of the
#'   subjects.
#' * `BDS`/`MDBDS`: records, subjects, visits, and the distribution of `AVAL`
#'   per parameter, the analysis visits, and the derived records per `DTYPE`.
#' * `OCCDS`/`MDOCCDS`: the number of distinct terms, and the
#'   treatment-emergent vs total records.
#' * `TTE`/`MDTTE`: events and censored per parameter.
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
#'   The record structure is checked against the keys stored by [admiral_df()].
#'   Keys are never inferred, so if none were stored the structure is not
#'   checked.
#'
#'   Subjects are counted by `USUBJID` if it is present, as the other subject
#'   keys (e.g. `STUDYID`) may not be populated on records added by a
#'   derivation. Otherwise the subject keys
#'   (`get_admiral_option("subject_keys")`) are used.
#'
#'   The dataset name is the one stored by [admiral_df()] or, failing that, the
#'   name of the variable `summary()` was called on.
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
#' @info The record structure is only checked once the keys have been stored
#'   with [admiral_df()].
#'
#' @code
#' advs <- admiral_df(
#'   advs,
#'   dataset_keys = exprs(USUBJID, PARAMCD, AVISIT),
#'   dataset_name = "ADVS"
#' )
#'
#' summary(advs)
summary.admiral_df <- function(object, ...) {
  name <- attr(object, "admiral_dataset_name")
  if (is.null(name)) {
    # the variable `summary()` was called on, unless it is a pipe placeholder
    # or a more complex expression than a name
    obj_expr <- substitute(object)
    if (is.name(obj_expr) && !as.character(obj_expr) %in% c(".", "_")) {
      name <- as.character(obj_expr)
    }
  }
  dataset_keys <- attr(object, "admiral_dataset_keys")
  # read before `ungroup()`, which need not keep the attributes
  dataset_class <- get_dataset_class(object)

  # a plain, ungrouped data frame, so that the tables derived from it do not
  # carry the class
  class(object) <- setdiff(class(object), "admiral_df")
  dataset <- ungroup(object)
  cols <- colnames(dataset)

  subject_vars <- if ("USUBJID" %in% cols) {
    "USUBJID"
  } else {
    intersect(vars2chr(get_admiral_option("subject_keys")), cols)
  }

  out <- list(
    name = name,
    dataset_class = dataset_class,
    n_subjects = if (length(subject_vars) > 0) {
      n_distinct(dataset[subject_vars])
    } else {
      NA_integer_
    },
    n_obs = nrow(dataset),
    n_vars = ncol(dataset),
    subject_vars = subject_vars,
    structure = summarize_structure(dataset, dataset_keys)
  )

  # the class is `NA` for a dataset with no records whose variables do not
  # decide it, which has nothing class-specific to report. The medical device
  # classes share the content of their subject counterparts
  out$details <- switch(coalesce(dataset_class, "other"),
    ADSL = summarize_adsl(dataset, subject_vars),
    BDS = ,
    MDBDS = summarize_bds(dataset, subject_vars),
    OCCDS = ,
    MDOCCDS = summarize_occds(dataset),
    TTE = ,
    MDTTE = summarize_tte(dataset),
    list()
  )

  class(out) <- c("summary_admiral_df", "list")
  out
}

#' Check the Record Structure of a Dataset
#'
#' @param dataset A data frame
#' @param keys The keys stored by [admiral_df()], as a character vector, or
#'   `NULL` if none were stored
#'
#' @return A list with the elements `keys`, `missing_keys` (keys which are not
#'   in `dataset`), and `n_duplicates` (the number of records beyond one per
#'   key, `NA` if the structure was not checked).
#'
#' @keywords internal
#' @family internal
summarize_structure <- function(dataset, keys) {
  # stored keys are plain names, so a `rename()` or `select()` since they
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
    missing_keys = missing_keys,
    n_duplicates = n_duplicates
  )
}

#' Summarize the Type-Specific Content of a Dataset
#'
#' Collect the information [summary.admiral_df()] reports for a dataset of
#' class `ADSL`, `BDS`, `OCCDS`, or `TTE` (or their medical device
#' counterparts). An element is only returned when the variables it needs are in
#' the dataset. This matters because a class set with [admiral_df()] is not
#' checked against the variables, so e.g. a dataset declared as `"TTE"` need
#' not have `CNSR`.
#'
#' @param dataset A data frame
#' @param subject_vars The variables subjects are counted by
#'
#' @details
#'   The ADSL items other than the treatment arms and population flags are
#'   taken from the first record of each subject, so that duplicated subject
#'   records -- which the record structure check reports -- do not inflate them.
#'   Every ADSL item counts subjects, so none is returned for a dataset without
#'   any of the variables subjects are counted by.
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

  # every item counts subjects, which cannot be identified without these
  if (length(subject_vars) == 0) {
    return(out)
  }

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
  if ("PARAMCD" %in% cols && nrow(dataset) > 0) {
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

  # the same `--BODSYS`/`--DECOD`/`--TERM` naming `infer_dataset_class()` relies
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
  if (all(c("PARAMCD", "CNSR") %in% colnames(dataset)) && nrow(dataset) > 0) {
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
#' print(summary(admiral_df(adsl, dataset_keys = exprs(USUBJID), dataset_name = "ADSL")))
print.summary_admiral_df <- function(x, ...) {
  # the class is what tailors the summary, so it heads the title; "other" and an
  # unknown class would read oddly there ("other summary"), hence the fallback
  title <- if (x$dataset_class %in% admiral_dataset_classes) {
    paste(x$dataset_class, "summary")
  } else {
    "Dataset summary"
  }
  # a name which only repeats the class ("ADSL: ADSL summary") is left out
  if (!is.null(x$name) && !identical(x$name, x$dataset_class)) {
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
#' Print the lines of [print.summary_admiral_df()] for the class-specific items
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
  n_extra <- structure$n_duplicates
  if (length(structure$missing_keys) > 0) {
    missing_keys <- paste(structure$missing_keys, collapse = ", ")
    cli_bullets(c("!" = "Structure: not checked,
       {qty(length(structure$missing_keys))}key{?s} {missing_keys}
       {?is/are} not in the dataset"))
  } else if (is.null(structure$keys)) {
    # keys are never inferred, so without stored keys there is nothing to check
    cli_bullets(c(i = "Structure: not checked, no keys set (see {.fun admiral_df})"))
  } else if (length(structure$keys) == 0) {
    # empty keys record that no keys are defined, which is not the same as
    # never having set them
    cli_bullets(c(i = "Structure: no key variables defined"))
  } else if (n_extra == 0) {
    cli_bullets(c(v = "Structure: one record per {keys}"))
  } else {
    cli_bullets(c(x = "Structure: expected one record per {keys}, found
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
