#' Default Format for the Disposition Reason
#'
#' @description
#' `r lifecycle::badge("deprecated")`
#'
#' This function is *deprecated*. This function is a default for `derive_vars_disposition_reason()`
#' for the `format_new_vars` argument. Please use `derive_vars_merged()` and
#' specify the `filter_add` argument to derive the respective variables.
#'
#' Define a function to map the disposition reason, to be used as a parameter in
#' `derive_vars_disposition_reason()`.
#'
#' @param reason the disposition variable used for the mapping (e.g. `DSDECOD`).
#' @param reason_spe the disposition variable used for the mapping of the details
#' if required (e.g. `DSTERM`).
#'
#' @details
#' `format_reason_default(DSDECOD)` returns `DSDECOD` when `DSDECOD` is not `'COMPLETED'` nor `NA`.
#' \cr`format_reason_default(DSDECOD, DSTERM)` returns `DSTERM` when `DSDECOD` is
#' equal to `'OTHER'`.
#' \cr Usually this function can not be used with `%>%`.
#'
#' @return A `character` vector
#'
#' @export
#' @family deprecated
#' @keywords deprecated
#' @seealso [derive_vars_disposition_reason()]
format_reason_default <- function(reason, reason_spe = NULL) {
  ### DEPRECATION
  deprecate_warn("0.10.0",
    "format_reason_default()",
    details = paste(
      "This function is a default for `derive_vars_disposition_reason() and is being deprecated`",
      "Please use `derive_vars_merged()` and",
      "specify the `filter_add` argument to derive the respective variables."
    )
  )

  if (is.null(reason_spe)) {
    if_else(reason != "COMPLETED" & !is.na(reason), reason, NA_character_)
  } else {
    if_else(reason == "OTHER", reason_spe, NA_character_)
  }
}

#' Derive a Disposition Reason at a Specific Timepoint
#'
#' @description
#' `r lifecycle::badge("deprecated")`
#'
#' This function is *deprecated*. Please use `derive_vars_merged()` and
#' specify the `filter_add` argument to derive the respective variables.
#'
#' Derive a disposition reason from the the relevant records in the disposition domain.
#'
#' @param dataset Input dataset
#'
#' @param dataset_ds Dataset containing the disposition information (e.g. `ds`)
#'
#' The dataset must contain:
#' - `STUDYID`, `USUBJID`,
#' - The variable(s) specified in the `reason_var` (and `reason_var_spe`, if required)
#' - The variables used in `filter_ds`.
#'
#' @param new_var Name of the disposition reason variable
#'
#' A variable name is expected (e.g. `DCSREAS`).
#'
#' @param reason_var The variable used to derive the disposition reason
#'
#' A variable name is expected (e.g. `DSDECOD`).
#'
#' @param new_var_spe Name of the disposition reason detail variable
#'
#' A variable name is expected (e.g. `DCSREASP`).
#' If `new_var_spe` is specified, it is expected that `reason_var_spe` is also specified,
#' otherwise an error is issued.
#'
#' Default: NULL
#'
#' @param reason_var_spe The variable used to derive the disposition reason detail
#'
#' A variable name is expected (e.g. `DSTERM`).
#' If `new_var_spe` is specified, it is expected that `reason_var_spe` is also specified,
#' otherwise an error is issued.
#'
#' Default: NULL
#'
#' @param format_new_vars The function used to derive the reason(s)
#'
#' This function is used to derive the disposition reason(s) and must follow the below conventions
#'
#' - If only the main reason for discontinuation needs to be derived (i.e. `new_var_spe` is NULL),
#' the function must have at least one character vector argument, e.g.
#' `format_reason <- function(reason)`
#' and `new_var` will be derived as `new_var = format_reason(reason_var)`.
#' Typically, the content of the function would return `reason_var` or `NA` depending on the
#' value (e.g. `if_else ( reason != "COMPLETED" & !is.na(reason), reason, NA_character_)`).
#' `DCSREAS = format_reason(DSDECOD)` returns `DCSREAS = DSDECOD`
#' when `DSDECOD` is not `'COMPLETED'` nor `NA`, `NA` otherwise.
#'
#' - If both the main reason and the details needs to be derived (`new_var_spe` is specified)
#' the function must have two character vectors argument, e.g.
#' `format_reason2 <- function(reason, reason_spe)` and
#' `new_var` will be derived as `new_var = format_reason(reason_var)`,
#' `new_var_spe` will be derived as `new_var_spe = format_reason(reason_var, reason_var_spe)`.
#' Typically, the content of the function would return `reason_var_spe` or `NA` depending on the
#' `reason_var` value (e.g. `if_else ( reason == "OTHER", reason_spe, NA_character_)`).
#' `DCSREASP = format_reason(DSDECOD, DSTERM)` returns `DCSREASP = DSTERM` when
#' `DSDECOD` is equal to `'OTHER'`.
#'
#' Default: `format_reason_default`, see [`format_reason_default()`] for details.
#'
#' @param filter_ds Filter condition for the disposition data.
#'
#' Filter used to select the relevant disposition data.
#' It is expected that the filter restricts `dataset_ds` such that there is at most
#' one observation per patient. An error is issued otherwise.
#'
#' Permitted Values: logical expression.
#'
#' @param subject_keys Variables to uniquely identify a subject
#'
#' A list of expressions where the expressions are symbols as returned by
#' `exprs()` is expected.
#'
#' @return the input dataset with the disposition reason(s) (`new_var` and
#' if required `new_var_spe`) added.
#'
#' @details
#' This functions returns the main reason for discontinuation (e.g. `DCSREAS` or `DCTREAS`).
#' The reason for discontinuation is derived based on `reason_var` (e.g. `DSDECOD`) and
#' `format_new_vars`.
#' If `new_var_spe` is not NULL, then the function will also return the details associated
#' with the reason for discontinuation (e.g. `DCSREASP`).
#' The details associated with the reason for discontinuation are derived based on
#' `reason_var_spe` (e.g. `DSTERM`), `reason_var` and `format_new_vars`.
#'
#' @family deprecated
#' @seealso [format_reason_default()]
#' @keywords deprecated
#'
#'
#' @export
#'
derive_vars_disposition_reason <- function(dataset,
                                           dataset_ds,
                                           new_var,
                                           reason_var,
                                           new_var_spe = NULL,
                                           reason_var_spe = NULL,
                                           format_new_vars = format_reason_default,
                                           filter_ds,
                                           subject_keys = get_admiral_option("subject_keys")) {
  ### DEPRECATION
  deprecate_warn("0.10.0",
    "derive_vars_disposition_reason()",
    details = paste(
      "Please use `derive_vars_merged()`",
      "and specify the `filter_add` argument to derive the respective variables"
    )
  )

  new_var <- assert_symbol(enexpr(new_var))
  reason_var <- assert_symbol(enexpr(reason_var))
  new_var_spe <- assert_symbol(enexpr(new_var_spe), optional = T)
  reason_var_spe <- assert_symbol(enexpr(reason_var_spe), optional = T)
  assert_s3_class(format_new_vars, "function")
  filter_ds <- assert_filter_cond(enexpr(filter_ds))
  assert_vars(subject_keys)
  assert_data_frame(dataset, required_vars = subject_keys)
  assert_data_frame(
    dataset_ds,
    required_vars = expr_c(subject_keys, reason_var, reason_var_spe)
  )
  warn_if_vars_exist(dataset, as_name(new_var))

  # Additional checks
  if (!is.null(new_var_spe)) {
    if (!is.null(reason_var_spe)) {
      statusvar <- c(as_name(reason_var), as_name(reason_var_spe))
    } else {
      err_msg <- paste(
        "`new_var_spe` is specified as ", as_name(new_var_spe),
        "but `reason_var_spe` is NULL.",
        "Please specify `reason_var_spe` together with `new_var_spe`."
      )
      abort(err_msg)
    }
  } else {
    statusvar <- as_name(reason_var)
  }

  dataset <- dataset %>%
    derive_vars_merged(
      dataset_add = dataset_ds,
      filter_add = !!filter_ds,
      new_vars = expr_c(reason_var, reason_var_spe),
      by_vars = subject_keys
    ) %>%
    mutate(!!new_var := format_new_vars(!!reason_var))

  if (!is.null(new_var_spe)) {
    dataset <- mutate(
      dataset,
      !!new_var_spe := format_new_vars(!!reason_var, !!reason_var_spe)
    )
  }
  select(dataset, -statusvar)
}
