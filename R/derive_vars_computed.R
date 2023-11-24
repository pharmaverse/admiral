#' Adds a Variable Computed from the Analysis Value of one or more Parameters
#'
#' Adds a variable computed from the analysis value of one or more parameters.
#' It is expected that the value of the new variable is defined by an expression
#' using the analysis values of other parameters. For example Body Mass Index at
#' Baseline (BMIBL) can be derived from of HEIGHT and WEIGHT parameters.
#'
#' @param dataset
#'
#'   The variables specified by the `by_vars` parameter are expected.
#'
#'   The variable specified by `by_vars` must be a unique key of the input
#'   dataset after restricting it by the filter condition (`filter` parameter).
#'
#' @param dataset_add Additional dataset
#'
#'   The variables specified by the `by_vars` parameter are expected.
#'
#'   The variable specified by `by_vars` and `PARAMCD` must be a unique key of
#'   the additional dataset after restricting it by the filter condition
#'   (`filter_add` parameter) and to the parameters specified by `parameters`.
#'
#'   If the argument is specified, the observations of the additional dataset
#'   are considered in addition to the observations from the input dataset
#'   (`dataset` restricted by `filter`).
#'
#'
#' @param filter_add Filter condition of additional dataset
#'
#'   The specified condition is applied to the additional dataset before
#'   deriving the new variable, i.e., only observations fulfilling the
#'   condition are taken into account.
#'
#'   *Permitted Values:* a condition
#'
#' @param parameters Required parameter codes
#'
#'   It is expected that all parameter codes (`PARAMCD`) which are required to
#'   derive the new variable are specified for this parameter or the
#'   `constant_parameters` parameter.
#'
#'   If observations should be considered which do not have a parameter code,
#'   e.g., if an SDTM dataset is used, temporary parameter codes can be derived
#'   by specifying a list of expressions. The name of the element defines the
#'   temporary parameter code and the expression defines the condition for
#'   selecting the records. For example,
#'   `parameters = exprs(HGHT = VSTESTCD == "HEIGHT")` selects the observations
#'   with `VSTESTCD == "HEIGHT"` from the input data (`dataset` and
#'   `dataset_add`), sets `PARAMCD = "HGHT"` for these observations, and adds
#'   them to the observations to consider.
#'
#'   Unnamed elements in the list of expressions are considered as parameter
#'   codes. For example, `parameters = exprs(WEIGHT, HGHT = VSTESTCD ==
#'   "HEIGHT")` uses the parameter code `"WEIGHT"` and creates a temporary
#'   parameter code `"HGHT"`.
#'
#'   *Permitted Values:* A character vector of `PARAMCD` values or a list of expressions
#'
#'
#' @param by_vars Grouping variables
#'
#'   Grouping variables uniquely identifying a set of records for which
#'   `new_vars` are to be calculated.
#'
#'   *Permitted Values:* list of variables
#'
#' @param constant_parameters Required constant parameter codes
#'
#'   It is expected that all the parameter codes (`PARAMCD`) which are required
#'   to derive the new variable and are measured only once are specified here.
#'   For example if BMI should be derived and height is measured only once while
#'   weight is measured at each visit. Height could be specified in the
#'   `constant_parameters` parameter. (Refer to the Example)
#'
#'   If observations should be considered which do not have a parameter code,
#'   e.g., if an SDTM dataset is used, temporary parameter codes can be derived
#'   by specifying a list of expressions. The name of the element defines the
#'   temporary parameter code and the expression defines the condition for
#'   selecting the records. For example `constant_parameters =
#'   exprs(HGHT = VSTESTCD == "HEIGHT")` selects the observations with
#'   `VSTESTCD == "HEIGHT"` from the input data (`dataset` and `dataset_add`),
#'   sets `PARAMCD = "HGHT"` for these observations, and adds them to the
#'   observations to consider.
#'
#'   Unnamed elements in the list of expressions are considered as parameter
#'   codes. For example, `constant_parameters = exprs(WEIGHT, HGHT = VSTESTCD ==
#'   "HEIGHT")` uses the parameter code `"WEIGHT"` and creates a temporary
#'   parameter code `"HGHT"`.
#'
#'   *Permitted Values:* A character vector of `PARAMCD` values or a list of expressions
#'
#' @param constant_by_vars By variables for constant parameters
#'
#'   The constant parameters (parameters that are measured only once) are merged
#'   to the other parameters using the specified variables.
#'   (Refer to the Example)
#'
#'   *Permitted Values:* list of variables
#'
#'
#' @param new_vars Name of the newly created variables
#'
#'   The specified variables are set to the specified values. The values of
#'   variables of the parameters specified by `parameters` can be accessed using
#'   `<variable name>.<parameter code>`. For example
#'   ```
#'   exprs(
#'     BMIBL = (AVAL.WEIGHT / (AVAL.HEIGHT/100)^2)
#'   )
#'   ```
#'   defines the analysis value and parameter code for the new variable.
#'
#'   Variable names in the expression must not contain more than one dot.
#'
#'   *Permitted Values:* List of variable-value pairs
#'
#' @param keep_nas Keep observations with `NA`s
#'
#'   If the argument is set to `TRUE`, observations are added even if some of
#'   the values contributing to the computed value are `NA`.
#'
#' @details For each group (with respect to the variables specified for the
#'   `by_vars` parameter), the value of the new variable is computed in the
#'   output dataset if the filtered input dataset (`dataset`) or the additional
#'   dataset (`dataset_add`) contains exactly one observation for each parameter
#'   code specified for `parameters`.
#'
#' @return The input dataset with the new variable added.
#'
#' @family der_prm_bds_findings
#'
#' @keywords der_prm_bds_findings
#'
#' @export
#'
#' @examples
#' library(tibble)
#' library(dplyr)
#'
#' # Example 1: Derive BMIBL
#' adsl <- tribble(
#'   ~STUDYID,   ~USUBJID, ~AGE,   ~AGEU,
#'   "PILOT01", "01-1302",   61, "YEARS",
#'   "PILOT01", "17-1344",   64, "YEARS"
#' )
#'
#' advs <- tribble(
#'   ~STUDYID, ~USUBJID, ~PARAMCD, ~PARAM, ~VISIT, ~AVAL, ~AVALU, ~ABLFL,
#'   "PILOT01", "01-1302", "HEIGHT", "Height (cm)", "SCREENING", 177.8, "cm", "Y",
#'   "PILOT01", "01-1302", "WEIGHT", "Weight (kg)", "SCREENING", 81.19, "kg", "N",
#'   "PILOT01", "01-1302", "WEIGHT", "Weight (kg)", "BASELINE", 82.1, "kg", "Y",
#'   "PILOT01", "01-1302", "WEIGHT", "Weight (kg)", "WEEK 2", 81.19, "kg", "N",
#'   "PILOT01", "01-1302", "WEIGHT", "Weight (kg)", "WEEK 4", 82.56, "kg", "N",
#'   "PILOT01", "01-1302", "WEIGHT", "Weight (kg)", "WEEK 6", 80.74, "kg", "N",
#'   "PILOT01", "17-1344", "HEIGHT", "Height (cm)", "SCREENING", 163.5, "cm", "Y",
#'   "PILOT01", "17-1344", "WEIGHT", "Weight (kg)", "SCREENING", 58.06, "kg", "N",
#'   "PILOT01", "17-1344", "WEIGHT", "Weight (kg)", "BASELINE", 58.06, "kg", "Y",
#'   "PILOT01", "17-1344", "WEIGHT", "Weight (kg)", "WEEK 2", 58.97, "kg", "N",
#'   "PILOT01", "17-1344", "WEIGHT", "Weight (kg)", "WEEK 4", 57.97, "kg", "N",
#'   "PILOT01", "17-1344", "WEIGHT", "Weight (kg)", "WEEK 6", 58.97, "kg", "N"
#' )
#'
#' derive_vars_computed(
#'   dataset = adsl,
#'   dataset_add = advs,
#'   by_vars = exprs(STUDYID, USUBJID),
#'   parameters = c("WEIGHT"),
#'   constant_by_vars = exprs(STUDYID, USUBJID),
#'   constant_parameters = c("HEIGHT"),
#'   new_vars = exprs(BMIBL = compute_bmi(height = AVAL.HEIGHT, weight = AVAL.WEIGHT)),
#'   filter_add = ABLFL == "Y"
#' )
derive_vars_computed <- function(dataset,
                                 dataset_add = NULL,
                                 by_vars,
                                 parameters,
                                 new_vars,
                                 filter_add = NULL,
                                 constant_by_vars = NULL,
                                 constant_parameters = NULL,
                                 keep_nas = FALSE) {
  filter_add <- assert_filter_cond(enexpr(filter_add), optional = TRUE)

  dataset_add <- dataset_add %>%
    filter_if(filter_add)

  derive_param_return <- derive_param_computed(
    dataset = dataset,
    dataset_add = dataset_add,
    by_vars = by_vars,
    parameters = parameters,
    set_values_to = new_vars,
    constant_by_vars = constant_by_vars,
    constant_parameters = constant_parameters,
    keep_nas = keep_nas
  )

  if ((names(new_vars)[1]) %in% names(derive_param_return)) {
    derive_param_return <- derive_param_return %>%
      filter(!is.na(!!!parse_exprs(names(new_vars)[1]))) %>%
      select(where(function(x) any(!is.na(x))))

    derive_vars_merged(
      dataset,
      dataset_add = derive_param_return,
      by_vars = by_vars
    )
  } else {
    dataset
  }
}
