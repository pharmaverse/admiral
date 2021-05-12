#' Filters observations from more than one dataset
#'
#' Filters observations from more than one dataset. For each dataset the
#' observations can be selected by a condition and/or by selecting the first or
#' last observation in each by group.
#'
#' @param sources List of sources
#'
#'   For each source the observations are selected as specified. The selected
#'   observations of all sources are collected in one dataset and then for each
#'   by group one observation (with respect to the `order` and `mode` parameter)
#'   is selected.
#'
#'   Each element of the list must be a list with the following named elements.
#'
#'   \itemize{
#'   \item \emph{dataset:}
#'     Dataset to add
#'
#'   The variables specified by the \code{by_vars} parameter, the \code{filter}, and the
#'   \code{order} element are expected.
#'
#'   \item \emph{filter:} Filter condition for dataset
#'
#'   Only observations of the dataset which fulfill the specified condition
#'   are used for merging.
#'
#'   \emph{Permitted Values:} logical expression
#'
#'   \item \emph{new_vars:} Variable to add
#'
#'   The specified variables are created in the source dataset and added to the
#'   output dataset.
#'
#'   \emph{Permitted Values:} list of name-value pairs
#'
#'   \item \emph{order:} Sort order
#'
#'   If the parameter is specified, the source dataset is ordered by the
#'   specified order and only the first or last observation (depending on the
#'   mode) in each by group is used for merging.
#'
#'   \emph{Permitted Values:} list of variables or functions of variables
#'
#'   \item \emph{mode:} mode
#'
#'   If the \code{filter_order} parameter is specified, the mode determines if
#'   the first or last observation of each by group is selected.
#'
#'   \emph{Permitted Values:} \code{"first"}, \code{"last"}
#'   }
#'
#' @param order Sort order
#'
#'   If the parameter is specified, all selected observations from the source
#'   dataset are ordered by the specified order and only the first or last
#'   observation (depending on the filter mode) in each by group is used for
#'   merging.
#'
#'   The temporary variable `temp_source_nr`, which is set to the number of the
#'   source, can be used for ordering. It is not included in the output dataset.
#'
#'   Permitted Values: list of variables or functions of variables
#'
#' @param mode Filter mode
#'
#'   If the `order` parameter is specified, the filter mode determines if the
#'   first or last observation of each by group is selected.
#'
#'   Permitted Values: `"first"`, `"last"`
#'
#' @param by_vars Grouping variables
#'
#'   Default: `exprs(USUBJID)`
#'
#'   Permitted Values: list of variables
#'
#' @details The following steps are performed to create the output dataset:
#'
#' \enumerate{
#' \item For each source dataset the observations as specified by the `filter`
#' element are selected. If the `order` element is specified, for each by group
#' the first or last element (depending on the `mode` element) with respect to
#' the specified order is selected.
#'
#' \item The selected observations of all source datasets are combined into a
#' single dataset. The variable `temp_source_nr` is added. It is set to the
#' source number.
#'
#' \item For each group (with respect to the variables specified for the
#' `by_vars` parameter) the first or last observation (with respect to the order
#' specified for the `order` parameter and the `mode` parameter) from the single
#' dataset is added to the output dataset.
#' }
#'
#' @author Stefan Bundfuss
#'
#' @return A dataset containing the first or last observation of each by group.
#'   Only the by variables and the new variables are included.
#'
#' @keywords derivation adam
#'
#' @export
#'
filter_extreme_multisource <- function(sources,
                                       by_vars,
                                       order,
                                       mode){
  add_data <- vector("list", length(sources))
  for (i in seq_along(sources)) {
    if (!is.null(sources[[i]]$filter)) {
      add_data[[i]] <- sources[[i]]$dataset %>%
        filter(!!(sources[[i]]$filter))
    }
    else {
      add_data[[i]] <- sources[[i]]$dataset
    }
    if (!is.null(sources[[i]]$order)){
      add_data[[i]] <- filter_extreme(add_data[[i]],
                                      order = sources[[i]]$order,
                                      by_vars = by_vars,
                                      mode = sources[[i]]$mode)
    }
    add_data[[i]] <- transmute(add_data[[i]], temp_source_nr = i, !!!by_vars, !!!sources[[i]]$vars)
  }
  bind_rows(add_data) %>%
    filter_extreme(by_vars = by_vars,
                   order = order,
                   mode = mode) %>%
    select(-starts_with("temp_"))
}
