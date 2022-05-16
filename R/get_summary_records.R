#' To Get The Summary Record
#'
#' @description
#' It is not uncommon to have an analysis need whereby one needs to derive an
#' analysis value (`AVAL`) from multiple records. The ADaM basic dataset
#' structure variable `DTYPE` is available to indicate when a new derived
#' records has been added to a dataset.
#'
#' @details
#' When all records have same values within `by_vars` then this function will
#' retain those common values in the newly derived records. Otherwise new value
#' will be set to `NA`.
#'
#' @param dataset A data frame.
#'
#' @param by_vars Variables to consider for generation of groupwise summary
#'   records. Providing the names of variables in [vars()] will create a
#'   groupwise summary and generate summary records for the specified groups.
#'
#' @param filter Filter condition as logical expression to apply during
#'   summary calculation. By default, filtering expressions are computed within
#'   `by_vars` as this will help when an aggregating, lagging, or ranking
#'   function is involved.
#'
#'   For example,
#'
#'   + `filter_rows = (AVAL > mean(AVAL, na.rm = TRUE))` will filter all AVAL
#'   values greater than mean of AVAL with in `by_vars`.
#'   + `filter_rows = (dplyr::n() > 2)` will filter n count of `by_vars` greater
#'   than 2.
#'
#' @param analysis_var Analysis variable.
#'
#' @param summary_fun Function that takes as an input the `analysis_var` and
#'   performs the calculation.
#'   This can include built-in functions as well as user defined functions,
#'   for example `mean` or `function(x) mean(x, na.rm = TRUE)`.
#'
#' @param set_values_to A list of variable name-value pairs. Use this argument
#'   if you need to change the values of any newly derived records.
#'
#'   Set a list of variables to some specified value for the new observation(s)
#'   + LHS refer to a variable.
#'   + RHS refers to the values to set to the variable. This can be a string, a symbol, a numeric
#'   value or NA.
#'   (e.g.  `vars(PARAMCD = "TDOSE",PARCAT1 = "OVERALL")`).
#'   More general expression are not allowed.
#'
#' @author Pavan Kumar, updated by Alana Harris
#'
#' @return A data frame with derived records appended to original dataset.
#'
#' @keywords adam user_utility.
#'
#' @export
#'
#' @examples


get_summary_records <- function(dataset,
                                by_vars,
                                filter = NULL,
                                analysis_var,
                                summary_fun,
                                set_values_to){

  # Quotes the filter argument and ensure that it is a condition
  filter <- assert_filter_cond(enquo(filter), optional = TRUE)

  # Summarise the analysis value
  dataset %>%
    group_by(!!!by_vars) %>%
    filter_if(filter) %>%
    summarise(!!analysis_var := summary_fun(!!analysis_var)) %>%
    mutate(!!!set_values_to) %>%
    ungroup()
}
