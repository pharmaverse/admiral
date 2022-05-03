#' Default Format for Disposition Status
#'
#' Define a function to map the disposition status. To be used as an input for
#' `derive_var_disposition_status()`.
#'
#' @param x the disposition variable used for the mapping (e.g. `DSDECOD`).
#'
#' @return A `character` vector derived based on the values given in `x`:
#'  "COMPLETED" if `x` is "COMPLETED",
#'  "DISCONTINUED" if `x` is not "COMPLETED" nor NA,
#'  "ONGOING" otherwise.
#'
#' @author Samia Kabi
#' @export
#' @keywords user_utility adsl computation
#' @seealso [derive_var_disposition_status()]
#' @examples
#' library(dplyr, warn.conflicts = FALSE)
#' library(admiraltest)
#' data("dm")
#' data("ds")
#'
#' dm %>%
#'   derive_var_disposition_status(
#'     dataset_ds = ds,
#'     new_var = EOSSTT,
#'     status_var = DSDECOD,
#'     format_new_var = format_eoxxstt_default,
#'     filter_ds = DSCAT == "DISPOSITION EVENT"
#'   ) %>%
#'   select(STUDYID, USUBJID, EOSSTT)
format_eoxxstt_default <- function(x) {
  case_when(
    x == "COMPLETED" ~ "COMPLETED",
    x != "COMPLETED" & !is.na(x) ~ "DISCONTINUED",
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
#' ```
#' format_eoxxstt_default <- function(x) {
#'   case_when(
#'     x == "COMPLETED" ~ "COMPLETED",
#'     x != "COMPLETED" & !is.na(x) ~ "DISCONTINUED",
#'     TRUE ~ "ONGOING"
#'   )
#' }
#' ```
#' where `x` is the `status_var.`
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
#'  "COMPLETED" if `status_var` == "COMPLETED",
#'  "DISCONTINUED" if `status_var` is not "COMPLETED" nor NA,
#'  "ONGOING" otherwise).
#'
#' @keywords adsl
#'
#' @author Samia Kabi
#'
#' @export
#'
#' @examples
#' library(dplyr, warn.conflicts = FALSE)
#' library(admiraltest)
#' data("dm")
#' data("ds")
#'
#' # Default derivation: EOSSTT =
#' #- COMPLETED when status_var = COMPLETED
#' #- DISCONTINUED when status_var is not COMPLETED nor NA
#' #- ONGOING otherwise
#'
#' dm %>%
#'   derive_var_disposition_status(
#'     dataset_ds = ds,
#'     new_var = EOSSTT,
#'     status_var = DSDECOD,
#'     filter_ds = DSCAT == "DISPOSITION EVENT"
#'   ) %>%
#'   select(STUDYID, USUBJID, EOSSTT)
#'
#' # Specific derivation: EOSSTT =
#' #- COMPLETED when status_var = COMPLETED
#' #- DISCONTINUED DUE TO AE when status_var = ADVERSE EVENT
#' #- DISCONTINUED NOT DUE TO AE when status_var != ADVERSE EVENT nor COMPLETED nor missing
#' #- ONGOING otherwise
#'
#' format_eoxxstt1 <- function(x) {
#'   case_when(
#'     x == "COMPLETED" ~ "COMPLETED",
#'     x == "ADVERSE EVENT" ~ "DISCONTINUED DUE TO AE",
#'     !(x %in% c("ADVERSE EVENT", "COMPLETED")) & !is.na(x) ~ "DISCONTINUED NOT DUE TO AE",
#'     TRUE ~ "ONGOING"
#'   )
#' }
#'
#' dm %>%
#'   derive_var_disposition_status(
#'     dataset_ds = ds,
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
                                      subject_keys = vars(STUDYID, USUBJID)) {
  new_var <- assert_symbol(enquo(new_var))
  status_var <- assert_symbol(enquo(status_var))
  filter_ds <- assert_filter_cond(enquo(filter_ds))
  assert_that(is.function(format_new_var))
  assert_data_frame(dataset)
  assert_data_frame(dataset_ds, quo_c(status_var))
  warn_if_vars_exist(dataset, quo_text(new_var))
  assert_vars(subject_keys)

  # Add the status variable and derive the new dispo status in the input dataset
  dataset %>%
    derive_vars_merged(dataset_add = dataset_ds,
                       filter_add = !!filter_ds,
                       new_vars = vars(!!status_var),
                       by_vars = subject_keys) %>%
    mutate(!!new_var := format_new_var(!!status_var)) %>%
    select(-!!status_var)
}
