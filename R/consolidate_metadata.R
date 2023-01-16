#' Consolidate Multiple Meta Datasets Into a Single One
#'
#' The purpose of the function is to consolidate multiple meta datasets into a
#' single one. For example, from global and project specific parameter mappings
#' a single lookup table can be created.
#'
#' @param datasets List of datasets to consolidate
#'
#'   *Permitted Values*: A named list of datasets
#'
#' @param key_vars Key variables
#'
#'   The specified variables must be a unique of all input datasets.
#'
#'   *Permitted Values*: A list of variables created by `vars()`
#'
#' @param source_var Source variable
#'
#'   The specified variable is added to the output dataset. It is set the name
#'   of the dataset the observation is originating from.
#'
#'   *Permitted Values*: A symbol
#'
#' @author Stefan Bundfuss
#'
#' @details All observations of the input datasets are put together into a
#'   single dataset. If a by group (defined by `key_vars`) exists in more of one
#'   of the input datasets, the observation from the last dataset is selected.
#'
#' @return A dataset which contains one row for each by group occuring in any of
#'   the input datasets.
#'
#' @keywords create_aux
#' @family create_aux
#'
#' @export
#'
#' @examples
#' library(tibble)
#' glob_ranges <- tribble(
#'   ~PARAMCD, ~ANRLO, ~ANRHI,
#'   "PULSE",      60,    100,
#'   "SYSBP",      90,    130,
#'   "DIABP",      60,     80
#' )
#' proj_ranges <- tribble(
#'   ~PARAMCD, ~ANRLO, ~ANRHI,
#'   "SYSBP",     100,    140,
#'   "DIABP",      70,     90
#' )
#' stud_ranges <- tribble(
#'   ~PARAMCD, ~ANRLO, ~ANRHI,
#'   "BMI",        18,     25
#' )
#'
#' consolidate_metadata(
#'   datasets = list(
#'     global = glob_ranges,
#'     project = proj_ranges,
#'     study = stud_ranges
#'   ),
#'   key_vars = vars(PARAMCD)
#' )
consolidate_metadata <- function(datasets,
                                 key_vars,
                                 source_var = SOURCE) {
  assert_list_of(datasets, class = "data.frame", named = TRUE)
  source_var <- assert_symbol(enquo(source_var))
  assert_vars(key_vars)

  data_order <- 1:length(datasets)
  names(data_order) <- names(datasets)
  all_data <- bind_rows(datasets, .id = as_label(source_var))
  temp_ord_var <- get_new_tmp_var(all_data)
  all_data %>% mutate(!!temp_ord_var := data_order[!!source_var]) %>%
    filter_extreme(
      by_vars = key_vars,
      order = vars(!!temp_ord_var),
      mode = "last",
      check_type = "none"
    ) %>%
    remove_tmp_vars()
}
