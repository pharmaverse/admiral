#' Create an Existence Flag
#'
#' @description Create a flag variable for the input dataset which indicates if
#' there exists at least one observation in another dataset fulfilling a certain
#' condition.
#'
#' **Note:** This is a helper function for `derive_vars_merged_exist_flag` which inputs
#' this result into `derive_vars_merged()`.
#'
#'
#' @param new_var New variable
#'
#'   The specified variable is added to the input dataset.
#'
#' @param condition Condition
#'
#'   The condition is evaluated at the additional dataset (`dataset_add`). For
#'   all by groups where it evaluates as `TRUE` at least once the new variable
#'   is set to the true value (`true_value`). For all by groups where it
#'   evaluates as `FALSE` or `NA` for all observations the new variable is set
#'   to the false value (`false_value`). The new variable is set to the missing
#'   value (`missing_value`) for by groups not present in the additional
#'   dataset.
#'
#' @param true_value True value
#'
#' @param false_value False value
#'
#' @param missing_value Values used for missing information
#'
#'   The new variable is set to the specified value for all by groups without
#'   observations in the additional dataset.
#'
#'   *Permitted Value*: A character scalar
#'
#' @param filter_add Filter for additional data
#'
#'   Only observations fulfilling the specified condition are taken into account
#'   for flagging. If the argument is not specified, all observations are
#'   considered.
#'
#'   *Permitted Values*: a condition
#'
#' @return The output dataset contains all observations and variables of the
#'   input dataset and additionally the variable specified for `new_var` derived
#'   from the additional dataset (`dataset_add`).
#'
#' @details
#'
#'   1. The additional dataset is restricted to the observations matching the
#'   `filter_add` condition.
#'
#'   1. The new variable is created to be added to the input dataset and set to the true value
#'   (`true_value`) if for the by group at least one observation exists in the
#'   (restricted) additional dataset where the condition evaluates to `TRUE`. It
#'   is set to the false value (`false_value`) if for the by group at least one
#'   observation exists and for all observations the condition evaluates to
#'   `FALSE` or `NA`. Otherwise, it is set to the missing value
#'   (`missing_value`).
#'
#'
#' @family der_gen
#' @keywords der_gen
#'
#' @export
#'
#' @examples
#'
#' library(dplyr, warn.conflicts = FALSE)
#'
#'
#' ae <- tribble(
#'   ~STUDYID,  ~DOMAIN,  ~USUBJID,    ~AETERM,     ~AEREL,
#'   "PILOT01",    "AE", "01-1028", "ERYTHEMA", "POSSIBLE",
#'   "PILOT01",    "AE", "01-1028", "PRURITUS", "PROBABLE",
#'   "PILOT01",    "AE", "06-1049",  "SYNCOPE", "POSSIBLE",
#'   "PILOT01",    "AE", "06-1049",  "SYNCOPE", "PROBABLE"
#' )
#'
#'
#' derive_var_merged_exist_flag(
#'   dataset_add = ae,
#'   new_var = exprs(AERELFL),
#'   condition = exprs(AEREL == "PROBABLE")
#' ) %>%
#'   select(STUDYID, USUBJID, AGE, AGEU, AERELFL)
#'
#' vs <- tribble(
#'   ~STUDYID,  ~DOMAIN,  ~USUBJID,      ~VISIT, ~VSTESTCD, ~VSSTRESN, ~VSBLFL,
#'   "PILOT01",    "VS", "01-1028", "SCREENING",  "HEIGHT",     177.8,      NA,
#'   "PILOT01",    "VS", "01-1028", "SCREENING",  "WEIGHT",     98.88,      NA,
#'   "PILOT01",    "VS", "01-1028",  "BASELINE",  "WEIGHT",     99.34,     "Y",
#'   "PILOT01",    "VS", "01-1028",    "WEEK 4",  "WEIGHT",     98.88,      NA,
#'   "PILOT01",    "VS", "04-1127", "SCREENING",  "HEIGHT",     165.1,      NA,
#'   "PILOT01",    "VS", "04-1127", "SCREENING",  "WEIGHT",     42.87,      NA,
#'   "PILOT01",    "VS", "04-1127",  "BASELINE",  "WEIGHT",     41.05,     "Y",
#'   "PILOT01",    "VS", "04-1127",    "WEEK 4",  "WEIGHT",     41.73,      NA,
#'   "PILOT01",    "VS", "06-1049", "SCREENING",  "HEIGHT",    167.64,      NA,
#'   "PILOT01",    "VS", "06-1049", "SCREENING",  "WEIGHT",     57.61,      NA,
#'   "PILOT01",    "VS", "06-1049",  "BASELINE",  "WEIGHT",     57.83,     "Y",
#'   "PILOT01",    "VS", "06-1049",    "WEEK 4",  "WEIGHT",     58.97,      NA
#' )
#' derive_var_merged_exist_flag(
#'   dataset_add = vs,
#'   filter_add = VSTESTCD == "WEIGHT" & VSBLFL == "Y",
#'   new_var = exprs(WTBLHIFL),
#'   condition = exprs(VSSTRESN > 90),
#'   false_value = "N",
#'   missing_value = "M"
#' ) %>%
#'   select(STUDYID, USUBJID, AGE, AGEU, WTBLHIFL)
#'
derive_var_exist_flag <- function(dataset_add,
                                  new_var,
                                  condition,
                                  true_value = "Y",
                                  false_value = NA_character_,
                                  missing_value = NA_character_,
                                  filter_add = NULL) {
  # Very unclear to me why, but without these variables just hanging out here, it does not work
  condition
  new_var
  new_var <- assert_symbol(enexpr(new_var))
  condition <- assert_filter_cond(enexpr(condition))

  filter_add <-
    assert_filter_cond(filter_add, optional = TRUE)
  if (is.null(filter_add)) {
    add_data <- dataset_add %>%
      mutate(!!new_var := if_else(!!condition, 1, 0, 0))
  } else {
    add_data <- filter_if(dataset_add, filter_add) %>%
      mutate(!!new_var := if_else(!!condition, 1, 0, 0))
  }
  add_data
}
