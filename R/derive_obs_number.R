#' Adds a variable numbering the observations within by group
#'
#' Adds a variable numbering the observations within by group
#'
#' @param dataset Input dataset
#'
#'   The variables specified by the `order` and the `by_vars` parameter are
#'   expected.
#'
#' @param new_var Name of variable to create
#'
#'   The new variable is set to the observation number for each by group. The
#'   numbering starts with 1.
#'
#'   Default: `ASEQ`

#' @param order Sort order
#'
#'   Within each by group the observations are ordered by the specified order.
#'
#'   Permitted Values: list of variables or functions of variables
#'
#' @param by_vars Grouping variables
#'
#'   Permitted Values: list of variables
#'
#' @param check_type Check uniqueness?
#'
#'   If `"warning"` or `"error"` is specified, the specified message is issued
#'   if the observations of the input dataset are not unique with respect to the
#'   by variables and the order.
#'
#'   Default: `"none"`
#'
#'   Permitted Values: `"none"`, `"warning"`, `"error"`
#'
#' @details For each group (with respect to the variables specified for the
#'   `by_vars` parameter) the first or last observation (with respect to the
#'   order specified for the `order` parameter and the mode specified for the
#'   `mode` parameter) is included in the output dataset.
#'
#' @author Stefan Bundfuss
#'
#' @return A dataset containing all observations and variables of the input
#'   dataset and additionally the variable specified by the `new_var` parameter.
#'
#' @keywords adam derivation
#'
#' @export
#'
#' @examples
#' library(dplyr)
#' data("vs")
#'
#' vs %>%
#'   select(USUBJID, VSTESTCD, VISITNUM, VSTPTNUM) %>%
#'   filter(VSTESTCD %in% c("HEIGHT", "WEIGHT")) %>%
#'   derive_obs_number(
#'     by_vars = vars(USUBJID, VSTESTCD),
#'     order = vars(VISITNUM, VSTPTNUM)
#'   )
derive_obs_number <- function(dataset,
                              new_var = ASEQ,
                              order = NULL,
                              by_vars = NULL,
                              check_type = "none") {
  arg_match(check_type, c("none", "warning", "error"))
  data <- dataset

  if (!is.null(by_vars) | !is.null(order)) {
    # group and sort input dataset
    if (!is.null(by_vars)) {
      assert_has_variables(data, vars2chr(by_vars))

      data <- data %>%
        group_by(!!!by_vars) %>%
        arrange(!!!order, .by_group = TRUE)

      if (check_type != "none") {
        signal_duplicate_records(
          data,
          by_vars = c(by_vars, extract_vars(order)),
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
    mutate(!!enquo(new_var) := row_number()) %>%
    ungroup()
}
