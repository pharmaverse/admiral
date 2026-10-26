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
#'   that a derivation which merely tags its output does not lose the metadata
#'   set by the user on its input.
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
#' dataset name, and the ADaM dataset type. All `{admiral}` dataset functions
#' tag their output with the class, and keep any metadata set on their input.
#'
#' The keys describe the record structure the dataset is intended to have; they
#' are not checked against its records when they are stored. `{admiral}` never
#' infers the keys or the type of a dataset: if they have not been set with
#' this function, functions which need them ask for them to be set.
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
#'   The core `{dplyr}` verbs (`mutate()`, `filter()`, `arrange()`, `select()`,
#'   `rename()`, `slice()`, `distinct()`, joins where `dataset` is `x`, ...)
#'   restore data frame attributes, so the metadata and the `admiral_df` class
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
#'   attributes but drop the `admiral_df` class.
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

  # the type is not inferred, but a device type declared for a dataset without
  # the device identifier is most likely the wrong type (or the wrong dataset).
  # Only a warning, as the declaration is the user's to make
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
  details <- c(
    if (!is.null(name)) paste0("name: ", name),
    if (!is.null(type)) paste0("type: ", type),
    # an empty attribute records that no keys are defined, so it is shown
    # rather than omitted as if the keys had never been set
    if (!is.null(keys)) {
      paste0("keys: ", if (length(keys) > 0) paste(keys, collapse = ", ") else "none")
    }
  )
  if (length(details) > 0) {
    cat("# An admiral dataset (", paste(details, collapse = ", "), ")\n", sep = "")
  }
  NextMethod()
}
