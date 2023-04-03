#' Title
#'
#' @param adevent
#' @param sources
#' @param source_datasets
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
derive_param_extreme_record <- function(adevent,
                                        sources,
                                        source_datasets,
                                        order,
                                        mode,
                                        set_values_to) {

  # Create a list of data frames with the selected records from each source
  selected_records <- lapply(sources, function(source) {
    dataset_name <- source$dataset_name
    filter_expr <- source$filter
    new_vars <- source$new_vars

    # Filter the source dataset based on the given expression
    source_dataset <- source_datasets[[dataset_name]] %>%
      filter(!!filter_expr)

    # Create new variables based on the given expressions
    new_vars_values <- source_dataset %>%
      transmute_at(vars(all_of(names(new_vars))), .funs = list(~eval(new_vars[[as.character(substitute(.))]])))

    # Return the selected record with the new variables
    new_vars_values[order(new_vars_values$ADT), ][mode == "first", ]
  })

  # Select the extreme record based on the given mode
  extreme_record <- do.call(rbind, selected_records)[order(do.call(c, order)), ][mode == "first", ]

  # Add the new row to the input ADaM dataset with the given values
  new_row <- set_values_to %>%
    mutate_at(vars(all_of(names(extreme_record))), .funs = list(~eval(set_values_to[[as.character(substitute(.))]])))
  bind_rows(adevent, new_row)
}
