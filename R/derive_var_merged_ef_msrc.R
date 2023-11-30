#' Merge an Existence Flag From Multiple Sources
#'
#' @description Adds a flag variable to the input dataset which indicates if
#'   there exists at least one observation in one of the source datasets
#'   fulfilling a certain condition. For example, if a dose adjustment flag
#'   should be added to `ADEX` but the dose adjustment information is collected
#'   in different datasets, e.g., `EX`, `EC`, and `FA`.
#'
#' @param flag_events Flag events
#'
#'   A list of `flag_event()` objects is expected. For each event the condition
#'   (`condition` field) is evaluated in the source dataset referenced by the
#'   `dataset_name` field. If it evaluates to `TRUE` at least once, the new
#'   variable is set to `true_value`.
#'
#' @param source_datasets Source datasets
#'
#'   A named list of datasets is expected. The `dataset_name` field of
#'   `flag_event()` refers to the dataset provided in the list.
#'
#' @param new_var New variable
#'
#'   The specified variable is added to the input dataset.
#'
#' @param true_value True value
#'
#'   The new variable (`new_var`) is set to the specified value for all by
#'   groups for which at least one of the source object (`sources`) has the
#'   condition evaluate to `TRUE`.
#'
#'   The values of `true_value`, `false_value`, and `missing_value` must be of
#'   the same type.
#'
#' @param false_value False value
#'
#'   The new variable (`new_var`) is set to the specified value for all by
#'   groups which occur in at least one source (`sources`) but the condition
#'   never evaluates to `TRUE`.
#'
#'   The values of `true_value`, `false_value`, and `missing_value` must be of
#'   the same type.
#'
#' @param missing_value Values used for missing information
#'
#'   The new variable is set to the specified value for all by groups without
#'   observations in any of the sources (`sources`).
#'
#'   The values of `true_value`, `false_value`, and `missing_value` must be of
#'   the same type.
#'
#' @inheritParams derive_var_merged_exist_flag
#'
#' @return The output dataset contains all observations and variables of the
#'   input dataset and additionally the variable specified for `new_var`.
#'
#' @details
#'
#'   1. For each `flag_event()` object specified for `flag_events`: The
#'   condition (`condition`) is evaluated in the dataset referenced by
#'   `dataset_name`. If the `by_vars` field is specified the dataset is grouped
#'   by the specified variables for evaluating the condition. If named elements
#'   are used in `by_vars` like `by_vars = exprs(USUBJID, EXLNKID = ECLNKID)`,
#'   the variables are renamed after the evaluation. If the `by_vars` element is
#'   not specified, the observations are grouped by the variables specified for
#'   the `by_vars` argument.
#'
#'   1. The new variable (`new_var`) is added to the input dataset and set to
#'   the true value (`true_value`) if for the by group at least one condition
#'   evaluates to `TRUE` in one of the sources. It is set to the false value
#'   (`false_value`) if for the by group at least one observation exists and for
#'   all observations the condition evaluates to `FALSE` or `NA`. Otherwise, it
#'   is set to the missing value (`missing_value`).
#'
#' @seealso [flag_event()]
#'
#' @family der_gen
#' @keywords der_gen
#'
#' @export
#'
#' @examples
#' library(dplyr)
#'
#' # Derive a flag indicating anti-cancer treatment based on CM and PR
#' adsl <- tribble(
#'   ~USUBJID,
#'   "1",
#'   "2",
#'   "3",
#'   "4"
#' )
#'
#' cm <- tribble(
#'   ~USUBJID, ~CMCAT,        ~CMSEQ,
#'   "1",      "ANTI-CANCER",      1,
#'   "1",      "GENERAL",          2,
#'   "2",      "GENERAL",          1,
#'   "3",      "ANTI-CANCER",      1
#' )
#'
#' # Assuming all records in PR indicate cancer treatment
#' pr <- tibble::tribble(
#'   ~USUBJID, ~PRSEQ,
#'   "2",      1,
#'   "3",      1
#' )
#'
#' derive_var_merged_ef_msrc(
#'   adsl,
#'   by_vars = exprs(USUBJID),
#'   flag_events = list(
#'     flag_event(
#'       dataset_name = "cm",
#'       condition = CMCAT == "ANTI-CANCER"
#'     ),
#'     flag_event(
#'       dataset_name = "pr"
#'     )
#'   ),
#'   source_datasets = list(cm = cm, pr = pr),
#'   new_var = CANCTRFL
#' )
#'
#' # Using different by variables depending on the source
#' # Add a dose adjustment flag to ADEX based on ADEX, EC, and FA
#' adex <- tribble(
#'   ~USUBJID, ~EXLNKID, ~EXADJ,
#'   "1",      "1",      "AE",
#'   "1",      "2",      NA_character_,
#'   "1",      "3",      NA_character_,
#'   "2",      "1",      NA_character_,
#'   "3",      "1",      NA_character_
#' )
#'
#' ec <- tribble(
#'   ~USUBJID, ~ECLNKID, ~ECADJ,
#'   "1",      "3",      "AE",
#'   "3",      "1",      NA_character_
#' )
#'
#' fa <- tribble(
#'   ~USUBJID, ~FALNKID, ~FATESTCD, ~FAOBJ,            ~FASTRESC,
#'   "3",      "1",      "OCCUR",   "DOSE ADJUSTMENT", "Y"
#' )
#'
#' derive_var_merged_ef_msrc(
#'   adex,
#'   by_vars = exprs(USUBJID, EXLNKID),
#'   flag_events = list(
#'     flag_event(
#'       dataset_name = "ex",
#'       condition = !is.na(EXADJ)
#'     ),
#'     flag_event(
#'       dataset_name = "ec",
#'       condition = !is.na(ECADJ),
#'       by_vars = exprs(USUBJID, EXLNKID = ECLNKID)
#'     ),
#'     flag_event(
#'       dataset_name = "fa",
#'       condition = FATESTCD == "OCCUR" & FAOBJ == "DOSE ADJUSTMENT" & FASTRESC == "Y",
#'       by_vars = exprs(USUBJID, EXLNKID = FALNKID)
#'     )
#'   ),
#'   source_datasets = list(ex = adex, ec = ec, fa = fa),
#'   new_var = DOSADJFL
#' )
derive_var_merged_ef_msrc <- function(dataset,
                                      by_vars,
                                      flag_events,
                                      source_datasets,
                                      new_var,
                                      true_value = "Y",
                                      false_value = NA_character_,
                                      missing_value = NA_character_) {
  new_var <- assert_symbol(enexpr(new_var))
  assert_list_of(source_datasets, class = "data.frame", named = TRUE)
  assert_list_of(flag_events, "flag_event")

  source_names <- names(source_datasets)
  assert_list_element(
    list = flag_events,
    element = "dataset_name",
    condition = dataset_name %in% source_names,
    source_names = source_names,
    message_text = paste0(
      "The dataset names must be included in the list specified for the ",
      "`source_datasets` parameter.\n",
      "Following names were provided by `source_datasets`:\n",
      enumerate(source_names, quote_fun = squote)
    )
  )

  tmp_cond_val <- get_new_tmp_var(dataset, prefix = "tmp_cond_val")

  selected_records <- map(
    flag_events,
    function(event) {
      data_event <- source_datasets[[event$dataset_name]]
      if (is.null(event$by_vars)) {
        event_by_vars <- by_vars
      } else {
        event_by_vars <- event$by_vars
      }
      if (is.null(event$condition)) {
        event_condition <- TRUE
      } else {
        event_condition <- event$condition
      }
      data_selected <- data_event %>%
        group_by(!!!unname(event_by_vars)) %>%
        mutate(!!tmp_cond_val := if_else(!!event_condition, TRUE, FALSE, FALSE)) %>%
        ungroup() %>%
        select(!!tmp_cond_val, !!!event_by_vars)
    }
  )

  derive_var_merged_exist_flag(
    dataset,
    dataset_add = bind_rows(selected_records),
    by_vars = by_vars,
    new_var = !!new_var,
    condition = !!tmp_cond_val == TRUE,
    true_value = true_value,
    false_value = false_value,
    missing_value = missing_value
  )
}

#' Create a `flag_event` Object
#'
#' The `flag_event` object is used to define events as input for the
#' `derive_var_merged_ef_msrc()` function.
#'
#' @param dataset_name Dataset name of the dataset to be used as input for the
#'   event. The name refers to the dataset specified for `source_datasets` in
#'   `derive_var_merged_ef_msrc()`.
#'
#'   *Permitted Values*: a character scalar
#'
#' @param condition Condition
#'
#'   The condition is evaluated at the dataset referenced by `dataset_name`. For
#'   all by groups where it evaluates as `TRUE` at least once the new variable
#'   is set to the true value (`true_value`).
#'
#' @param by_vars Grouping variables
#'
#'   If specified, the dataset is grouped by the specified variables before the
#'   condition is evaluated. If named elements are used in `by_vars` like
#'   `by_vars = exprs(USUBJID, EXLNKID = ECLNKID)`, the variables are renamed
#'   after the evaluation. If the `by_vars` element is not specified, the
#'   observations are grouped by the variables specified for the `by_vars`
#'   argument of `derive_var_merged_ef_msrc()`.
#'
#' @seealso [derive_var_merged_ef_msrc()]
#'
#' @family source_specifications
#' @keywords source_specifications
#'
#' @export
flag_event <- function(dataset_name,
                       condition = NULL,
                       by_vars = NULL) {
  out <- list(
    dataset_name = assert_character_scalar(dataset_name),
    condition = assert_filter_cond(enexpr(condition), optional = TRUE),
    by_vars = assert_expr_list(by_vars, optional = TRUE)
  )
  class(out) <- c("flag_event", "source", "list")
  out
}
