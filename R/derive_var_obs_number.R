#' Adds a Variable Numbering the Observations Within Each By Group
#'
#' Adds a variable numbering the observations within each by group
#'
#' @param dataset
#'   `r roxygen_param_dataset(expected_vars = c("by_vars", "order"))`
#'
#' @param by_vars Grouping variables
#'
#'   `r roxygen_param_by_vars()`
#'
#'
#' @param order Sort order
#'
#'   Within each by group the observations are ordered by the specified order.
#'
#'   `r roxygen_order_na_handling()`
#'
#' @permitted list of variables or functions of variables
#'
#' @param new_var Name of variable to create
#'
#'   The new variable is set to the observation number for each by group. The
#'   numbering starts with 1.
#'
#' @permitted [var]
#'
#' @param check_type Check uniqueness?
#'
#'   If `"message"`, `"warning"` or `"error"` is specified, the specified
#'   message is issued if the observations of the input dataset are not unique
#'   with respect to the by variables and the order.
#'
#' @permitted [msg_type]
#'
#' @details For each group (with respect to the variables specified for the
#'   `by_vars` parameter) the first or last observation (with respect to the
#'   order specified for the `order` parameter and the mode specified for the
#'   `mode` parameter) is included in the output dataset.
#'
#'
#' @return A dataset containing all observations and variables of the input
#'   dataset and additionally the variable specified by the `new_var` parameter.
#'
#' @family der_gen
#' @keywords der_gen
#'
#' @export
#'
#' @examples
#' library(dplyr, warn.conflicts = FALSE)
#' vs <- tribble(
#'   ~STUDYID,  ~DOMAIN,      ~USUBJID, ~VSTESTCD, ~VISITNUM, ~VSTPTNUM,
#'   "PILOT01",    "VS", "01-703-1182",   "DIABP",         3,       815,
#'   "PILOT01",    "VS", "01-703-1182",   "DIABP",         3,       816,
#'   "PILOT01",    "VS", "01-703-1182",   "DIABP",         4,       815,
#'   "PILOT01",    "VS", "01-703-1182",   "DIABP",         4,       816,
#'   "PILOT01",    "VS", "01-703-1182",   "PULSE",         3,       815,
#'   "PILOT01",    "VS", "01-703-1182",   "PULSE",         3,       816,
#'   "PILOT01",    "VS", "01-703-1182",   "PULSE",         4,       815,
#'   "PILOT01",    "VS", "01-703-1182",   "PULSE",         4,       816,
#'   "PILOT01",    "VS", "01-703-1182",   "SYSBP",         3,       815,
#'   "PILOT01",    "VS", "01-703-1182",   "SYSBP",         3,       816,
#'   "PILOT01",    "VS", "01-703-1182",   "SYSBP",         4,       815,
#'   "PILOT01",    "VS", "01-703-1182",   "SYSBP",         4,       816,
#'   "PILOT01",    "VS", "01-716-1229",   "DIABP",         3,       815,
#'   "PILOT01",    "VS", "01-716-1229",   "DIABP",         3,       816,
#'   "PILOT01",    "VS", "01-716-1229",   "DIABP",         4,       815,
#'   "PILOT01",    "VS", "01-716-1229",   "DIABP",         4,       816,
#'   "PILOT01",    "VS", "01-716-1229",   "PULSE",         3,       815,
#'   "PILOT01",    "VS", "01-716-1229",   "PULSE",         3,       816,
#'   "PILOT01",    "VS", "01-716-1229",   "PULSE",         4,       815,
#'   "PILOT01",    "VS", "01-716-1229",   "PULSE",         4,       816,
#'   "PILOT01",    "VS", "01-716-1229",   "SYSBP",         3,       815,
#'   "PILOT01",    "VS", "01-716-1229",   "SYSBP",         3,       816,
#'   "PILOT01",    "VS", "01-716-1229",   "SYSBP",         4,       815,
#'   "PILOT01",    "VS", "01-716-1229",   "SYSBP",         4,       816
#' )
#' vs %>%
#'   derive_var_obs_number(
#'     by_vars = exprs(USUBJID, VSTESTCD),
#'     order = exprs(VISITNUM, desc(VSTPTNUM))
#'   )
derive_var_obs_number <- function(dataset,
                                  by_vars = NULL,
                                  order = NULL,
                                  new_var = ASEQ,
                                  check_type = "none") {
  # checks and quoting
  new_var <- assert_symbol(enexpr(new_var))
  assert_vars(by_vars, optional = TRUE)
  assert_expr_list(order, optional = TRUE)
  if (!is.null(by_vars)) {
    required_vars <- by_vars
  } else {
    required_vars <- NULL
  }
  if (!is.null(order)) {
    required_vars <- exprs(!!!required_vars, !!!extract_vars(order))
  }
  assert_data_frame(dataset, required_vars = required_vars)
  check_type <-
    assert_character_scalar(
      check_type,
      values = c("none", "warning", "error", "message"),
      case_sensitive = FALSE
    )

  # derivation
  data <- dataset

  if (!is.null(by_vars) || !is.null(order)) {
    # group and sort input dataset
    if (!is.null(by_vars)) {
      data <- data %>%
        group_by(!!!by_vars) %>%
        arrange(!!!order, .by_group = TRUE)

      if (check_type != "none") {
        signal_duplicate_records(
          data,
          by_vars = expr_c(by_vars, order),
          cnd_type = check_type
        )
      }
    } else {
      data <- data %>%
        arrange(!!!order)

      if (check_type != "none") {
        signal_duplicate_records(
          data,
          by_vars = extract_vars(order),
          cnd_type = check_type
        )
      }
    }
  }

  data %>%
    mutate(!!new_var := row_number()) %>%
    ungroup()
}
