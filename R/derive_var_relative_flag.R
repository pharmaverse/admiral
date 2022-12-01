#' Flag Observations Before or After a Condition is Fulfilled
#'
#' Flag all observations before or after the observation where a specified
#' condition is fulfilled for each by group. For example, the function could be
#' called to flag for each subject all observations before the first disease
#' progression or to flag all AEs after a specific AE.
#'
#' @param dataset Input dataset
#'
#'   The variables specified by the `order` and the `by_vars` argument are
#'   expected.
#'
#' @param by_vars Grouping variables
#'
#'   *Permitted Values:* list of variables created by `vars()`
#'
#' @param order Sort order
#'
#'   Within each by group the observations are ordered by the specified order.
#'
#'   *Permitted Values:* list of variables or `desc(<variable>)` function calls
#'   created by `vars()`, e.g., `vars(ADT, desc(AVAL))`
#'
#' @param new_var New variable
#'
#'   The variable is added to the input dataset and set to `"Y"` for all
#'   observations before or after the condition is fulfilled. For all other
#'   observations it is set to `NA`.
#'
#' @param condition Condition for Reference Observation
#'
#'   The specified condition determines the reference observation. In the output
#'   dataset all observations before or after (`selection` argument)
#'   the reference observation are flagged.
#'
#' @param mode Selection mode (first or last)
#'
#'   If `"first"` is specified, for each by group the observations before or
#'   after (`selection` argument) the observation where the condition
#'   (`condition` argument) is fulfilled the *first* time is flagged in the
#'   output dataset. If `"last"` is specified, for each by group the
#'   observations before or after (`selection` argument) the observation where
#'   the condition (`condition` argument) is fulfilled the *last* time is
#'   flagged in the output dataset.
#'
#'   *Permitted Values:* `"first"`, `"last"`
#'
#' @param selection Flag observations before or after the reference observation?
#'
#'   *Permitted Values:* `"before"`, `"after"`
#'
#' @param inclusive Flag the reference observation?
#'
#'   *Permitted Values:* `TRUE`, `FALSE`
#'
#' @param flag_no_ref_groups Should by groups without reference observation be flagged?
#'
#'   *Permitted Values:* `TRUE`, `FALSE`
#'
#' @param check_type Check uniqueness?
#'
#'   If `"warning"` or `"error"` is specified, the specified message is issued
#'   if the observations of the input dataset are not unique with respect to the
#'   by variables and the order.
#'
#'   *Permitted Values:* `"none"`, `"warning"`, `"error"`
#'
#' @author Stefan Bundfuss
#'
#' @details For each by group (`by_vars` argument) the observations before or
#'   after (`selection` argument) the observations where the condition
#'   (`condition` argument) is fulfilled the first or last time (`order`
#'   argument and `mode` argument) is flagged in the output dataset.
#'
#' @return The input dataset with the new variable (`new_var`) added
#'
#' @keywords der_gen
#' @family der_gen
#'
#' @export
#'
#' @examples
#' library(tibble)
#' library(dplyr, warn.conflicts = FALSE)
#'
#' # Flag all AEs after the first COVID AE
#' adae <- tribble(
#'   ~USUBJID, ~ASTDY, ~ACOVFL, ~AESEQ,
#'   "1",           2, NA,           1,
#'   "1",           5, "Y",          2,
#'   "1",           5, NA,           3,
#'   "1",          17, NA,           4,
#'   "1",          27, "Y",          5,
#'   "1",          32, NA,           6,
#'   "2",           8, NA,           1,
#'   "2",          11, NA,           2,
#' )
#'
#' derive_var_relative_flag(
#'   adae,
#'   by_vars = vars(USUBJID),
#'   order = vars(ASTDY, AESEQ),
#'   new_var = PSTCOVFL,
#'   condition = ACOVFL == "Y",
#'   mode = "first",
#'   selection = "after",
#'   inclusive = FALSE,
#'   flag_no_ref_groups = FALSE
#' )
#'
#' response <- tribble(
#'   ~USUBJID, ~AVISITN, ~AVALC,
#'   "1",      0,        "PR",
#'   "1",      1,        "CR",
#'   "1",      2,        "CR",
#'   "1",      3,        "SD",
#'   "1",      4,        "NE",
#'   "2",      0,        "SD",
#'   "2",      1,        "PD",
#'   "2",      2,        "PD",
#'   "3",      0,        "SD",
#'   "4",      0,        "SD",
#'   "4",      1,        "PR",
#'   "4",      2,        "PD",
#'   "4",      3,        "SD",
#'   "4",      4,        "PR"
#' )
#'
#' # Flag observations up to first PD for each patient
#' response %>%
#'   derive_var_relative_flag(
#'     by_vars = vars(USUBJID),
#'     order = vars(AVISITN),
#'     new_var = ANL02FL,
#'     condition = AVALC == "PD",
#'     mode = "first",
#'     selection = "before",
#'     inclusive = TRUE
#'   )
#'
#' # Flag observations up to first PD excluding baseline (AVISITN = 0) for each patient
#' response %>%
#'   restrict_derivation(
#'     derivation = derive_var_relative_flag,
#'     args = params(
#'       by_vars = vars(USUBJID),
#'       order = vars(AVISITN),
#'       new_var = ANL02FL,
#'       condition = AVALC == "PD",
#'       mode = "first",
#'       selection = "before",
#'       inclusive = TRUE
#'     ),
#'     filter = AVISITN > 0
#'   ) %>%
#'   arrange(USUBJID, AVISITN)
derive_var_relative_flag <- function(dataset,
                                     by_vars,
                                     order,
                                     new_var,
                                     condition,
                                     mode,
                                     selection,
                                     inclusive,
                                     flag_no_ref_groups = TRUE,
                                     check_type = "warning") {
  new_var <- assert_symbol(enquo(new_var))
  condition <- assert_filter_cond(enquo(condition))
  assert_logical_scalar(flag_no_ref_groups)

  # add obs number for merging
  tmp_obs_nr <- get_new_tmp_var(dataset, prefix = "tmp_obs_nr")
  data <- derive_var_obs_number(
    dataset,
    by_vars = by_vars,
    order = order,
    new_var = !!tmp_obs_nr,
    check_type = check_type
  )

  # select observations before/after the condition
  # set check_type to "none" as uniqueness was already checked by derive_var_obs_number()
  flag_obs <- filter_relative(
    data,
    by_vars = by_vars,
    order = order,
    condition = !!condition,
    mode = mode,
    selection = selection,
    inclusive = inclusive,
    keep_no_ref_groups = flag_no_ref_groups,
    check_type = "none"
  )

  # flag observations based on the selected observations
  derive_var_merged_exist_flag(
    data,
    dataset_add = flag_obs,
    by_vars = quo_c(by_vars, quo(!!tmp_obs_nr)),
    new_var = !!new_var,
    condition = TRUE
  ) %>%
    remove_tmp_vars()
}
