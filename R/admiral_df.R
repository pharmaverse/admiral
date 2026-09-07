#' Tag a Dataset with the `admiral_df` Class
#'
#' Adds the `admiral_df` class to a data frame (unless it is already present),
#' preserving the existing classes such as `tbl_df`/`data.frame`. This allows
#' downstream admiral tooling to recognize a dataset as one produced by an
#' admiral derivation.
#'
#' @param dataset A data frame (or `NULL`)
#'
#' @param keys Optional character vector of the intended key variables of
#'   `dataset`, i.e. the variables it should have one record per, e.g. the
#'   `by_vars` (plus `PARAMCD`) of the derivation that produced it. When
#'   supplied, it is stored in the `"admiral_keys"` attribute. `NULL` leaves
#'   any existing attribute untouched.
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

#' Extract the Key Variables of a Dataset from a `metacore` Specification
#'
#' Reads the key variables (the variables the dataset should have one record
#' per) for a given dataset from a `{metacore}` specification object. The keys
#' are taken from the `key_seq` column of the `ds_vars` table and returned in
#' key order. These can be used, e.g., to check whether a derivation preserved
#' the expected record structure.
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
#' `admiral_df` (see [as_admiral_df()]).
#'
#' @param dataset A data frame
#' @param metacore A `{metacore}` object, see [get_admiral_keys()]
#' @param dataset_name The dataset name, see [get_admiral_keys()]
#'
#' @return `dataset` with an `"admiral_keys"` attribute, an `"admiral_ds_name"`
#'   attribute, and the `admiral_df` class.
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
#'   by `rename()` and not removed by `select()`.
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
