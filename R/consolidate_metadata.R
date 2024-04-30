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
#'   *Permitted Values*: A list of variables created by `exprs()`
#'
#' @param source_var Source variable
#'
#'   The specified variable is added to the output dataset. It is set the name
#'   of the dataset the observation is originating from.
#'
#'   *Permitted Values*: A symbol
#'
#' @param check_vars Check variables?
#'
#'   If `"message"`, `"warning"`, or `"error"` is specified, a message is issued
#'   if the variable names differ across the input datasets (`datasets`).
#'
#'   *Permitted Values*: `"none"`, `"message"`, `"warning"`, `"error"`
#'
#' @param check_keys Check keys?
#'
#'  `r lifecycle::badge("deprecated")` Please use `check_type` instead.
#'
#'   If `"warning"` or `"error"` is specified, a message is issued if the key
#'   variables (`key_vars`) are not a unique key in all of the input datasets
#'   (`datasets`).
#'
#'   *Permitted Values*: `"none"`, `"warning"`, `"error"`
#'
#' @param check_type Check uniqueness?
#'
#'   If `"warning"` or `"error"` is specified, a message is issued if the key
#'   variables (`key_vars`) are not a unique key in all of the input datasets
#'   (`datasets`).
#'
#'   *Permitted Values*: `"none"`, `"warning"`, `"error"`
#'
#' @details All observations of the input datasets are put together into a
#'   single dataset. If a by group (defined by `key_vars`) exists in more than
#'   one of the input datasets, the observation from the last dataset is
#'   selected.
#'
#' @return A dataset which contains one row for each by group occurring in any
#'   of the input datasets.
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
#'   key_vars = exprs(PARAMCD)
#' )
consolidate_metadata <- function(datasets,
                                 key_vars,
                                 source_var = SOURCE,
                                 check_vars = "warning",
                                 check_keys,
                                 check_type = "error") {
  assert_list_of(datasets, cls = "data.frame", named = TRUE)
  assert_vars(key_vars)
  source_var <- assert_symbol(enexpr(source_var))
  check_vars <-
    assert_character_scalar(
      check_vars,
      values = c("none", "message", "warning", "error"),
      case_sensitive = FALSE
    )
  if (!is_missing(check_keys)) {
    deprecate_stop(
      "1.1.0",
      "consolidate_metadata(check_keys = )",
      "consolidate_metadata(check_type = )"
    )
    check_type <-
      assert_character_scalar(
        check_keys,
        values = c("none", "warning", "error"),
        case_sensitive = FALSE,
        optional = TRUE
      )
  }
  check_type <-
    assert_character_scalar(
      check_type,
      values = c("none", "warning", "error"),
      case_sensitive = FALSE
    )

  if (check_vars != "none") {
    if (length(unique(map(datasets, function(x) sort(names(x))))) > 1) {
      msg_funs <- list(message = cli_inform, warning = cli_warn, error = cli_abort)
      msg_funs[[check_vars]](c(
        "The variable names differ across the input datasets.",
        i = "This message can be suppressed by setting {.code check_vars = \"none\"}."
      ))
    }
  }

  data_order <- seq_len(length(datasets))
  names(data_order) <- names(datasets)
  all_data <- bind_rows(datasets, .id = as_label(source_var))
  tmp_source_ord <- get_new_tmp_var(all_data, prefix = "tmp_source_ord")
  all_data %>%
    mutate(!!tmp_source_ord := data_order[!!source_var]) %>%
    filter_extreme(
      by_vars = key_vars,
      order = exprs(!!tmp_source_ord),
      mode = "last",
      check_type = check_type
    ) %>%
    remove_tmp_vars()
}
