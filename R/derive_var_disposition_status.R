#' Default Format for Disposition Status
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
#' @author Samia Kabi
#' @details Usually this function can not be used with `%>%`.
#' @export
#' @family utils_fmt
#' @keywords utils_fmt
#' @seealso [derive_var_disposition_status()]
#' @examples
#' library(dplyr, warn.conflicts = FALSE)
#' library(admiral.test)
#' data("admiral_dm")
#' data("admiral_ds")
#'
#' admiral_dm %>%
#'   derive_var_disposition_status(
#'     dataset_ds = admiral_ds,
#'     new_var = EOSSTT,
#'     status_var = DSDECOD,
#'     format_new_var = format_eoxxstt_default,
#'     filter_ds = DSCAT == "DISPOSITION EVENT"
#'   ) %>%
#'   select(STUDYID, USUBJID, EOSSTT)
format_eoxxstt_default <- function(status) {
  case_when(
    status %in% c("SCREEN FAILURE", "SCREENING NOT COMPLETED") ~ "NOT STARTED",
    status == "COMPLETED" ~ "COMPLETED",
    !status %in% c("COMPLETED", "SCREEN FAILURE", "SCREENING NOT COMPLETED") &
      !is.na(status) ~ "DISCONTINUED",
    TRUE ~ "ONGOING"
  )
}

#' Derive a Disposition Status at a Specific Timepoint
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
#'   A list of quosures where the expressions are symbols as returned by
#'   `vars()` is expected.
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
#' @family der_adsl
#' @keywords der_adsl
#'
#' @author Samia Kabi
#'
#' @export
#'
#' @examples
#' library(dplyr, warn.conflicts = FALSE)
#' library(admiral.test)
#' data("admiral_dm")
#' data("admiral_ds")
#'
#' # Default derivation: EOSSTT =
#' #- NOT STARTED when status_var is SCREEN FAILURE or SCREENING NOT COMPLETED
#' #- COMPLETED when status_var is COMPLETED
#' #- DISCONTINUED when status_var is not COMPLETED nor SCREEN FAILURE nor
#' #  SCREENING NOT COMPLETED nor NA
#' #- ONGOING otherwise
#'
#' admiral_dm %>%
#'   derive_var_disposition_status(
#'     dataset_ds = admiral_ds,
#'     new_var = EOSSTT,
#'     status_var = DSDECOD,
#'     filter_ds = DSCAT == "DISPOSITION EVENT"
#'   ) %>%
#'   select(STUDYID, USUBJID, EOSSTT)
#'
#' # Specific derivation: EOSSTT =
#' #- NOT STARTED when status_var = SCREEN FAILURE
#' #- COMPLETED when status_var = COMPLETED
#' #- DISCONTINUED DUE TO AE when status_var = ADVERSE EVENT
#' #- DISCONTINUED NOT DUE TO AE when status_var != ADVERSE EVENT nor COMPLETED
#' #  nor SCREEN FAILURE nor missing
#' #- ONGOING otherwise
#'
#' format_eoxxstt1 <- function(x) {
#'   case_when(
#'     x == "SCREEN FAILURE" ~ "NOT STARTED",
#'     x == "COMPLETED" ~ "COMPLETED",
#'     x == "ADVERSE EVENT" ~ "DISCONTINUED DUE TO AE",
#'     !(x %in% c("ADVERSE EVENT", "COMPLETED", "SCREEN FAILURE")) & !is.na(x) ~
#'       "DISCONTINUED NOT DUE TO AE",
#'     TRUE ~ "ONGOING"
#'   )
#' }
#'
#' admiral_dm %>%
#'   derive_var_disposition_status(
#'     dataset_ds = admiral_ds,
#'     new_var = EOSSTT,
#'     status_var = DSDECOD,
#'     format_new_var = format_eoxxstt1,
#'     filter_ds = DSCAT == "DISPOSITION EVENT"
#'   ) %>%
#'   select(STUDYID, USUBJID, EOSSTT)
derive_var_disposition_status <- function(dataset,
                                          dataset_ds,
                                          new_var,
                                          status_var,
                                          format_new_var = format_eoxxstt_default,
                                          filter_ds,
                                          subject_keys = get_admiral_option("subject_keys")) {
  new_var <- assert_symbol(enquo(new_var))
  status_var <- assert_symbol(enquo(status_var))
  filter_ds <- assert_filter_cond(enquo(filter_ds))
  assert_s3_class(format_new_var, "function")
  assert_data_frame(dataset)
  assert_data_frame(dataset_ds, quo_c(status_var))
  warn_if_vars_exist(dataset, quo_text(new_var))
  assert_vars(subject_keys)

  # Add the status variable and derive the new dispo status in the input dataset
  dataset %>%
    derive_vars_merged(
      dataset_add = dataset_ds,
      filter_add = !!filter_ds,
      new_vars = vars(!!status_var),
      by_vars = subject_keys
    ) %>%
    mutate(!!new_var := format_new_var(!!status_var)) %>%
    select(-!!status_var)
}
