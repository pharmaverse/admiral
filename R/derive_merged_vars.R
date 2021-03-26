#' Merge variables from a second dataset
#'
#' Merges variables from a second dataset. The observations to merge can be
#' selected by a condition and/or by selecting the first observation in each by
#' group.
#'
#' @param dataset Input dataset
#'
#'   The variables specified by the `by_vars` parameter are expected.
#'
#' @param dataset_add Dataset to add
#'
#'   The variables specified by the `by_vars`, the `filter_add`, and the
#'   `filter_first_order` parameter are expected.
#'
#' @param filter_add Fiter condition for dataset to add
#'
#'   Only observations of the add dataset which fulfill the specified condition
#'   are used for merging.
#'
#'   Permitted Values: logical expression
#'
#' @param new_vars Variable to add
#'
#'   The specified variables are created in the add dataset and then merge to
#'   the input dataset.
#'
#'   Permitted Values: list of name-value pairs
#'
#' @param filter_first_order Sort order
#'
#'   If the parameter is specified, the add dataset is ordered by the specified
#'   order and only the first observation in each by group is used for merging.
#'
#'   Permitted Values: list of variables or functions of variables
#'
#' @param by_vars Grouping variables
#'
#'   Default: `exprs(USUBJID)`
#'   Permitted Values: list of variables
#'
#' @details For each group (with respect to the variables specified for the
#'   `by_vars` parameter) the first observation (with respect to the order
#'   specified for the `order` parameter) is included in the output dataset.
#'
#' @author Stefan Bundfuss
#'
#' @return The input dataset with variables from the add dataset
#'
#' @family {general functions}
#'
#' @export
#'
#' @examples
#' data("ex")
#' data("dm")
#'
#' # adding first treatment start date for each patient
#' derive_merged_vars(dm,
#'                    dataset_add = ex,
#'                    filter_add = rlang::exprs(EXDOSE > 0 | str_detect(EXTRT, 'PLACEBO')),
#'                    new_vars = rlang::exprs(TRTSDT := ymd(EXSTDTC)),
#'                    filter_first_order = rlang::exprs(EXSTDTC, EXSEQ))
#'

derive_merged_vars <- function(dataset,
                               dataset_add,
                               filter_add,
                               new_vars,
                               by_vars = exprs(USUBJID),
                               filter_first_order){
  assert_has_variables(dataset, by_vars)
  assert_has_variables(dataset_add, by_vars)

  if (!missing(filter_add)){
     add <- dataset_add %>%
       filter(!!!filter_add)
  }
  else{
    add <- dataset_add
  }
  if (!missing(filter_first_order)){
    add <- add %>%
      filter_first(order = filter_first_order,
                   by_vars = by_vars)
  }
  add <- add %>% transmute(!!!by_vars, !!!new_vars)
  left_join(dataset, add, by = map_chr(by_vars, as_string))
}
