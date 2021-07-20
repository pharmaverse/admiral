#' Adds a Parameter Computed from the Analysis Value of Other Parameters
#'
#' Adds a parameter computed from the analysis value of other parameters. It is
#' expected that the analysis value of the new parameter is defined by a formula
#' using the analysis values of other parameters. For example mean arterial
#' pressure (MAP) can be derived from systolic (SYSBP) and diastolic blood
#' pressure (DIABP) with the formula \deqn{MAP = (SYSBP + 2 DIABP) / 3}
#'
#' @param dataset Input dataset
#'
#'   The variables specified by the `by_vars` parameter, `PARAMCD`, and `AVAL`
#'   are expected.
#'
#'   The variable specified by `by_vars` and `PARAMCD` must be a unique key of
#'   the input dataset after restricting it by the filter condition (`filter`
#'   parameter) and to the parameters specified by `parameters`.
#'
#' @param filter Filter condition
#'
#'   The specified condition is applied to the input dataset before deriving the
#'   new parameter, i.e., only observations fulfilling the condition are taken
#'   into account.
#'
#'   *Permitted Values:* a condition
#'
#' @param parameters Required parameter codes
#'
#'   It is expected that all parameter codes (`PARAMCD`) which are required to
#'   derive the new parameter are specified.
#'
#'   *Permitted Values:* A character vector of `PARAMCD` values
#'
#' @param by_vars Grouping variables
#'
#'   For each group defined by `by_vars` an observation is added to the output
#'   dataset.
#'
#'   *Permitted Values:* list of variables
#'
#' @param analysis_value Definition of the analysis value
#'
#'   A formula defining the analysis value (`AVAL`) of the new parameter is
#'   expected. The analysis values of the parameters specified by `parameters`
#'   can be accessed using `AVAL.<parameter code>`, e.g., `AVAL.SYSBP`.
#'
#'   *Permitted Values:* An unquoted expression
#'
#' @param set_values_to Variables to be set
#'
#'   The specified variables are set to the specified values for the new
#'   observations. For example `vars(PARAMCD = "MAP")` defines the parameter
#'   code for the new parameter.
#'
#'   *Permitted Values:* List of variable-value pairs
#'
#' @details For each group (with respect to the variables specified for the
#'   `by_vars` parameter) an observation is added to the output dataset if the
#'   filtered input dataset contains exactly one observation for each parameter
#'   code specified for `parameters`.
#'
#'   For the new observations `AVAL` is set to the value specified by
#'   `analysis_value` and the variables specified for `set_values_to` are set to
#'   the provided values. The values of the other variables of the input dataset
#'   are retained if they are constant within each by group. Otherwise they are
#'   set to `NA`.
#'
#' @author Stefan Bundfuss
#'
#' @return The input dataset with the new parameter added
#'
#' @keywords derivation bds
#'
#' @export
#'
#' @examples
#' library(dplyr, warn.conflicts = FALSE)
#' library(magrittr)
#' data("advs")
#' advss <- advs %>%
#'   select(USUBJID, PARAMCD, PARAM, AVAL, VSTESTCD, ANL01FL, DTYPE, AVISIT, AVISITN) %>%
#'   filter(
#'     USUBJID %in% c("01-701-1015", "01-701-1363") &
#'       PARAMCD %in% c("SYSBP", "DIABP") &
#'       AVISIT %in% c("Baseline", "Week 2")
#'   )
#'
#' derive_derived_param(
#'   advss,
#'   filter = ANL01FL == "Y" & DTYPE == "AVERAGE",
#'   parameters = c("SYSBP", "DIABP"),
#'   by_vars = vars(USUBJID, AVISIT),
#'   analysis_value = (AVAL.SYSBP + 2 * AVAL.DIABP) / 3,
#'   set_values_to = vars(
#'     PARAMCD = "MAP",
#'     PARAM = "Mean arterial pressure (mmHg)",
#'     DTYPE = NA
#'   )
#' ) %>%
#'   filter(DTYPE == "AVERAGE" | PARAMCD == "MAP")

derive_derived_param <- function(dataset,
                                 filter = NULL,
                                 parameters,
                                 by_vars,
                                 analysis_value,
                                 set_values_to
) {
  # checking and quoting
  assert_vars(by_vars)
  assert_data_frame(dataset,
                    required_vars = vars(!!!by_vars, AVAL))
  filter <- assert_filter_cond(enquo(filter), optional = TRUE)
  assert_character_vector(parameters)

  # select observations and variables required for new observations
  data <- dataset %>%
    filter_if(filter) %>%
    filter(PARAMCD %in% parameters)

  keep_vars <- get_constant_vars(data, by_vars = by_vars)
  data <- data %>%
    select(c(vars2chr(keep_vars),
             "PARAMCD",
             "AVAL"))

  signal_duplicate_records(
    data,
    by_vars = vars(!!!by_vars, PARAMCD),
    msg = paste("The filtered input dataset contains duplicate records with respect to",
                enumerate(c(vars2chr(by_vars), "PARAMCD")),
                "\nPlease ensure that the variables specified for `by_vars` and `PARAMCD`",
                "are a unique key of the input data set restricted by the condition",
                "specified for `filter` and to the parameters specified for `parameters`.")
  )

  # horizontalize data, AVAL for PARAMCD = "PARAMx" -> AVAL.PARAMx
  hori_data <- data %>%
    spread(key = PARAMCD,
           value = AVAL,
           sep = ".")
  names(hori_data) <- map_chr(names(hori_data), str_replace, "PARAMCD.", "AVAL.")

  # add analysis value (AVAL) and parameter variables, e.g., PARAMCD
  hori_data <- hori_data %>%
    mutate(AVAL = !!enquo(analysis_value),
           !!!set_values_to) %>%
    select(-starts_with("AVAL."))

  bind_rows(dataset, hori_data)
}
