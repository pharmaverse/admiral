#' Filter the Observations Before or After a Condition is Fulfilled
#'
#' Filters the observations before or after the observation where a specified
#' condition is fulfilled for each by group. For example, the function could be
#' called to select for each subject all observations before the first disease
#' progression.
#'
#' @param dataset Input dataset
#'
#'   The variables specified by the `order` and the `by_vars` parameter are
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
#' @param condition Condition for Reference Observation
#'
#'   The specified condition determines the reference observation. The output
#'   dataset contains all observations before or after (`selection` parameter)
#'   the reference observation.
#'
#' @param mode Selection mode (first or last)
#'
#'   If `"first"` is specified, for each by group the observations before or
#'   after (`selection` parameter) the observation where the condition
#'   (`condition` parameter) is fulfilled the *first* time is included in the
#'   output dataset. If `"last"` is specified, for each by group the
#'   observations before or after (`selection` parameter) the observation where
#'   the condition (`condition` parameter) is fulfilled the *last* time is
#'   included in the output dataset.
#'
#'   *Permitted Values:* `"first"`, `"last"`
#'
#' @param selection Select observations before or after the reference observation?
#'
#'   *Permitted Values:* `"before"`, `"after"`
#'
#' @param inclusive Include the reference observation?
#'
#'   *Permitted Values:* `TRUE`, `FALSE`
#'
#' @param keep_no_ref_groups Should by groups without reference observation be kept?
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
#' @details For each by group ( `by_vars` parameter) the observations before or
#'   after (`selection` parameter) the observations where the condition
#'   (`condition` parameter) is fulfilled the first or last time (`order`
#'   parameter and `mode` parameter) is included in the output dataset.
#'
#' @author Stefan Bundfuss
#'
#' @return A dataset containing for each by group the observations before or
#'   after the observation where the condition was fulfilled the first or last
#'   time
#'
#' @keywords utils_fil
#' @family utils_fil
#'
#'
#' @export
#'
#' @examples
#' library(tibble)
#'
#' response <- tribble(
#'   ~USUBJID, ~AVISITN, ~AVALC,
#'   "1",      1,        "PR",
#'   "1",      2,        "CR",
#'   "1",      3,        "CR",
#'   "1",      4,        "SD",
#'   "1",      5,        "NE",
#'   "2",      1,        "SD",
#'   "2",      2,        "PD",
#'   "2",      3,        "PD",
#'   "3",      1,        "SD",
#'   "4",      1,        "SD",
#'   "4",      2,        "PR",
#'   "4",      3,        "PD",
#'   "4",      4,        "SD",
#'   "4",      5,        "PR"
#' )
#'
#' # Select observations up to first PD for each patient
#' response %>%
#'   filter_relative(
#'     by_vars = vars(USUBJID),
#'     order = vars(AVISITN),
#'     condition = AVALC == "PD",
#'     mode = "first",
#'     selection = "before",
#'     inclusive = TRUE
#'   )
#'
#' # Select observations after last CR, PR, or SD for each patient
#' response %>%
#'   filter_relative(
#'     by_vars = vars(USUBJID),
#'     order = vars(AVISITN),
#'     condition = AVALC %in% c("CR", "PR", "SD"),
#'     mode = "last",
#'     selection = "after",
#'     inclusive = FALSE
#'   )
#'
#' # Select observations from first response to first PD
#' response %>%
#'   filter_relative(
#'     by_vars = vars(USUBJID),
#'     order = vars(AVISITN),
#'     condition = AVALC %in% c("CR", "PR"),
#'     mode = "first",
#'     selection = "after",
#'     inclusive = TRUE,
#'     keep_no_ref_groups = FALSE
#'   ) %>%
#'   filter_relative(
#'     by_vars = vars(USUBJID),
#'     order = vars(AVISITN),
#'     condition = AVALC == "PD",
#'     mode = "first",
#'     selection = "before",
#'     inclusive = TRUE
#'   )
filter_relative <- function(dataset,
                            by_vars,
                            order,
                            condition,
                            mode,
                            selection,
                            inclusive,
                            keep_no_ref_groups = TRUE,
                            check_type = "warning") {
  assert_vars(by_vars)
  assert_order_vars(order)
  condition <- assert_filter_cond(enquo(condition))
  mode <-
    assert_character_scalar(
      mode,
      values = c("first", "last"),
      case_sensitive = FALSE
    )
  selection <-
    assert_character_scalar(
      selection,
      values = c("before", "after"),
      case_sensitive = FALSE
    )
  assert_logical_scalar(inclusive)
  assert_data_frame(dataset, required_vars = quo_c(by_vars, extract_vars(order)))
  check_type <-
    assert_character_scalar(
      check_type,
      values = c("none", "warning", "error"),
      case_sensitive = FALSE
    )

  data <- dataset %>%
    derive_var_obs_number(
      new_var = tmp_obs_nr_filter_relative,
      order = order,
      by_vars = by_vars,
      check_type = check_type
    )

  cond_nrs <- data %>%
    filter(!!condition) %>%
    filter_extreme(
      by_vars = by_vars,
      order = order,
      mode = mode,
      check_type = check_type
    ) %>%
    select(!!!by_vars, tmp_obs_nr_filter_relative)

  data <- derive_vars_merged(
    data,
    dataset_add = cond_nrs,
    new_vars = vars(tmp_obs_nr_match_filter_relative = tmp_obs_nr_filter_relative),
    by_vars = by_vars
  )

  # Build condition for selecting observations
  if (selection == "before") {
    if (inclusive) {
      operator <- "<="
    } else {
      operator <- "<"
    }
  } else {
    if (inclusive) {
      operator <- ">="
    } else {
      operator <- ">"
    }
  }
  selection_condition <- paste(
    "tmp_obs_nr_filter_relative",
    operator,
    "tmp_obs_nr_match_filter_relative"
  )
  if (keep_no_ref_groups) {
    selection_condition <- paste(selection_condition, "| is.na(tmp_obs_nr_match_filter_relative)")
  }

  data %>%
    filter(!!parse_expr(selection_condition)) %>%
    select(-tmp_obs_nr_match_filter_relative, -tmp_obs_nr_filter_relative)
}
