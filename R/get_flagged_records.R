#' Create an Existence Flag
#'
#' @description Create a flag variable for the input dataset which indicates if
#' there exists at least one observation in another dataset fulfilling a certain
#' condition.
#'
#' **Note:** This is a helper function for `derive_vars_merged_exist_flag` which
#' inputs this result into `derive_vars_merged()`.
#'
#' @param dataset Input dataset
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
#' @param filter Filter for additional data
#'
#'   Only observations fulfilling the specified condition are taken into account
#'   for flagging. If the argument is not specified, all observations are
#'   considered.
#'
#'   *Permitted Values*: a condition
#'
#' @return The output dataset contains a filtered version of the
#'   input dataset with the variable specified for `new_var` representing a flag
#'   for the condition
#'
#' @details
#'
#'   1. The additional dataset is restricted to the observations matching the
#'   `filter` condition.
#'
#'
#'
#' @family utils_help
#' @keywords utils_help
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
#' get_flagged_records(
#'   dataset = ae,
#'   new_var = AERELFL,
#'   condition = AEREL == "PROBABLE"
#' ) %>%
#'   select(STUDYID, USUBJID, AERELFL)
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
#' get_flagged_records(
#'   dataset = vs,
#'   new_var = WTBLHIFL,
#'   condition = VSSTRESN > 90,
#'   filter = VSTESTCD == "WEIGHT" & VSBLFL == "Y"
#' ) %>%
#'   select(STUDYID, USUBJID, WTBLHIFL)
#'
get_flagged_records <- function(dataset,
                                new_var,
                                condition,
                                filter = NULL) {
  new_var <- assert_symbol(enexpr(new_var))
  condition <- assert_filter_cond(enexpr(condition))
  filter <-
    assert_filter_cond(enexpr(filter), optional = TRUE)
  filter_if(dataset, filter) %>%
    mutate(!!new_var := if_else(!!condition, 1, 0, 0))
}
