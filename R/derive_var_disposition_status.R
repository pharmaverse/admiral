#' Default Format for Disposition Status
#'
#' @description
#' `r lifecycle::badge("deprecated")`
#'
#' This function is *deprecated*. This function is a default for `derive_var_disposition_status()`
#' for the `format_new_var` argument. Please define your own function and use that as input for the
#' `cat_fun` argument in `derive_var_merged_cat()` instead.
#'
#' Define a function to map the disposition status. To be used as an input for
#' `derive_var_disposition_status()`.
#'
#' @param status the disposition variable used for the mapping (e.g. `DSDECOD`).
#'
#' @return A `character` vector derived based on the values given in `status`:
#'  "NOT STARTED" if `status` is "SCREEN FAILURE" or "SCREENING NOT COMPLETED",
#'  "COMPLETED" if `status` is "COMPLETED",
#'  "DISCONTINUED" if `status` is not in ("COMPLETED","SCREEN FAILURE",
#'  "SCREENING NOT COMPLETED") nor NA,
#'  "ONGOING" otherwise.
#'
#' @details Usually this function can not be used with `%>%`.
#' @export
#' @family deprecated
#' @keywords deprecated
format_eoxxstt_default <- function(status) {
  ### DEPRECATION
  deprecate_stop("0.10.0",
    "format_eoxxstt_default()",
    details = paste(
      "This function is deprecated",
      "Please define your own function and use that as input for the
                    `cat_fun` argument in `derive_var_merged_cat()` instead"
    )
  )
}

#' Derive a Disposition Status at a Specific Timepoint
#'
#' @description
#' `r lifecycle::badge("deprecated")`
#'
#' This function is *deprecated*, Please define your own function and use that as input for the
#' `cat_fun` argument in `derive_var_merged_cat()` instead.
#'
#' Derive a disposition status from the the relevant records in the disposition domain.
#'
#' @param dataset Input dataset.
#'
#' @param dataset_ds Dataset containing the disposition information (e.g.: ds).
#'
#' It must contain:
#' - `STUDYID`, `USUBJID`,
#' - The variable(s) specified in the `status_var`
#' - The variables used in `filter_ds`.
#'
#' @param new_var Name of the disposition status variable.
#'
#' A variable name is expected (e.g. `EOSSTT`).
#'
#' @param status_var The variable used to derive the disposition status.
#'
#'   A variable name is expected (e.g. `DSDECOD`).
#'
#' @param format_new_var The format used to derive the status.
#'
#' Default: `format_eoxxstt_default()` defined as:
#' ``` {r echo=TRUE, eval=FALSE}
#' format_eoxxstt_default <- function(status) {
#'   case_when(
#'     status %in% c("SCREEN FAILURE", "SCREENING NOT COMPLETED") ~ "NOT STARTED",
#'     status == "COMPLETED" ~ "COMPLETED",
#'     !status %in% c("COMPLETED", "SCREEN FAILURE", "SCREENING NOT COMPLETED")
#'     & !is.na(status) ~ "DISCONTINUED",
#'     TRUE ~ "ONGOING"
#'   )
#' }
#' ```
#' where `status` is the `status_var.`
#'
#' @param filter_ds Filter condition for the disposition data.
#'
# 'It is expected that the filter restricts `dataset_ds` such that there is at most
#' one observation per patient. An error is issued otherwise.
#'
#' Permitted Values: logical expression.
#'
#' @param subject_keys Variables to uniquely identify a subject
#'
#'   A list of expressions where the expressions are symbols as returned by
#'   `exprs()` is expected.
#'
#' @return The input dataset with the disposition status (`new_var`) added.
#' `new_var` is derived based on the values given in `status_var` and according to the format
#'  defined by `format_new_var` (e.g. when the default format is used, the function will derive
#'  `new_var` as:
#'  "NOT STARTED" if `status` is "SCREEN FAILURE" or "SCREENING NOT COMPLETED",
#'  "COMPLETED" if `status_var` == "COMPLETED",
#'  "DISCONTINUED" if `status` is not in ("COMPLETED","SCREEN FAILURE",
#'  "SCREENING NOT COMPLETED") nor NA,
#'  "ONGOING" otherwise).
#'
#' @family deprecated
#' @keywords deprecated
#'
#'
#' @export
#'
derive_var_disposition_status <- function(dataset,
                                          dataset_ds,
                                          new_var,
                                          status_var,
                                          format_new_var = format_eoxxstt_default,
                                          filter_ds,
                                          subject_keys = get_admiral_option("subject_keys")) {
  ### DEPRECATION
  deprecate_stop("0.10.0",
    "derive_var_disposition_status()",
    details = "Please use `derive_var_merged_cat()` instead"
  )
}
