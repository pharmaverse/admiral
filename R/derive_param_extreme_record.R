#' Title
#'
#' @param adevent
#' @param sources
#' @param source_datasets
#' @param by_vars
#' @param order
#' @param mode
#' @param set_values_to
#'
#' @return
#' @export
#' @details The function takes as input an ADaM dataset (adevent), a list of sources with their filter and new variable expressions (sources), a list of source datasets (source_datasets), an order expression (order), a mode (mode) to select the extreme record (either "first" or "last"), and a list of values to set in the new row (set_values_to).
#' The function first loads the dplyr and lubridate packages, which are used for data manipulation and date-time conversion, respectively. Then, it applies the following steps:
#' For each source in the sources list, it filters the corresponding dataset in source_datasets based on the filter expression and creates new variables based on the new variable expressions. It returns the selected record with the new variables based on the order expression and the mode.
#' It concatenates the selected records from all sources into a single data frame, sorts them based on the order expression, and selects the extreme record based on the mode.
#' It creates a new row with the values in set_values_to and the values from the extreme record, and appends it to the input ADaM dataset.
#' Note that this implementation assumes that all data frames in source_datasets have the same column names as the corresponding ADaM datasets. If this is not the case, you may need to adjust the code accordingly.

#' @examples
derive_param_extreme_record <- function(aevent,
                                        sources,
                                        source_datasets,
                                        by_vars = NULL,
                                        order,
                                        mode,
                                        set_values_to) {
  # Create Empty list to contain source datasets
  data_list <- vector("list", length(sources))

  # Evaluate the expressions contained in the sources
  for (i in seq_along(sources)) {
    source_dataset <- source_datasets[[sources[[i]]$dataset_name]]
    data_list[[i]] <- source_dataset %>%
      filter_if(sources[[i]]$filter) %>%
      mutate(!!!sources[[i]]$new_vars) %>%
      select(!!!by_vars, ADT, AVALC)
  }

  # Bind the source datasets together and parse out the extreme value
  param_data <- bind_rows(data_list) %>%
    group_by(!!!by_vars)%>%
    distinct() %>%
    ungroup() %>%
    filter_extreme(.,
                   by_vars = by_vars,
                   order = order,
                   mode = mode) %>%
    mutate(!!!set_values_to)

  # Bind the parameter rows back to original adevent dataset
  data <- bind_rows(aevent, param_data)
  return(data)
}
