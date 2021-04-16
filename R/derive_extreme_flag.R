#' Adds a variable flagging the first or last observation within each by group
#'
#' Adds a variable flagging the first or last observation within each by group
#'
#' @param dataset Input dataset
#'
#'   The variables specified by the `by_vars` parameter are expected.
#'
#' @param new_var Variable to add
#'
#'   The specified variable is added to the output dataset. It is set to `"Y"`
#'   for the first or last observation (depending on the mode) of each by group.
#'
#'   Permitted Values: list of name-value pairs
#'
#' @param order Sort order
#'
#'   The first or last observation is determined with respect to the specified
#'   order.
#'
#'   Permitted Values: list of variables or functions of variables
#'
#' @param mode Flag mode
#'
#'   Determines of the first or last observation is flagged.
#'
#'   Default: `"last"`
#'
#'   Permitted Values: `"first"`, `"last"`
#'
#' @param by_vars Grouping variables
#'
#'   Permitted Values: list of variables
#'
#' @details For each group (with respect to the variables specified for the
#'   `by_vars` parameter) the first or last observation (with respect to the
#'   order specified for the `order` parameter and the flag mode) is included in
#'   the output dataset.
#'
#' @author Stefan Bundfuss
#'
#' @return The input dataset with the new flag variable added
#'
#' @keywords derivation adam
#'
#' @export
#'
#' @examples
#' library(dplyr)
#' library(magrittr)
#'
#' data("vs")
#'
#' # flag last value for each patient, test, and visit
#' derive_extreme_flag(vs,
#'                     new_var = LASTFL,
#'                     by_vars = rlang::exprs(USUBJID, VSTESTCD, VISIT),
#'                    order = rlang::exprs(VSTPTNUM)) %>%
#'   select(USUBJID, VSTESTCD, VISIT, VSTPTNUM, VSSTRESN, LASTFL)
#'

derive_extreme_flag <- function(dataset,
                                new_var,
                                by_vars,
                                order,
                                mode = "last",
                                check_type = "warning"){
  arg_match(mode, c("first", "last"))
  arg_match(check_type, c("none", "warning", "error"))
  assert_has_variables(dataset, map_chr(by_vars, as_string))

  data <- dataset %>%
    derive_obs_number(order = order,
                      by_vars = by_vars,
                      check_type = check_type)

  if (mode == "first"){
    data <- data %>%
      mutate(!!enquo(new_var) := if_else(temp_obs_nr == 1, "Y", ""))
  }
  else{
    data <- data %>%
      group_by(!!!by_vars) %>%
      mutate(!!enquo(new_var) := if_else(temp_obs_nr == n(), "Y", "")) %>%
      ungroup()
  }
  data %>% select(-temp_obs_nr)
}
