#' Merge an Existence Flag From Muliple Sources
#'
#' @description Adds a flag variable to the input dataset which indicates if
#'   there exists at least one observation in one of the source dataset
#'   fulfilling a certain condition.
#'
#' @param dataset_add Additional dataset
#'
#'   The variables specified by the `by_vars` argument are expected.
#'
#' @param by_vars Grouping variables
#'
#'   *Permitted Values*: list of variables
#'
#' @param new_var New variable
#'
#'   The specified variable is added to the input dataset.
#'
#' @param condition Condition
#'
#'   The condition is evaluated at the additional dataset (`dataset_add`). For
#'   all by groups where it evaluates as `TRUE` at least once the new variable
#'   is set to the true value (`true_value`). For all by groups where it
#'   evaluates as `FALSE` or `NA` for all observations the new variable is set
#'   to the false value (`false_value`). The new variable is set to the missing
#'   value (`missing_value`) for by groups not present in the additional
#'   dataset.
#'
#' @inheritParams derive_var_merged_exist_flag
#'
#' @return The output dataset contains all observations and variables of the
#'   input dataset and additionally the variable specified for `new_var` derived
#'   from the additional dataset (`dataset_add`).
#'
#' @details
#'
#'   1. For each `flag_event()` object specified for `sources`: The condition (`condition`) is evaluated in the
#'   dataset referenced by `dataset_name`. If the `by_vars` field is specified the dataset is grouped by the specified variables.
#'
#'   1. The new variable (`new_var`) is added to the input dataset and set to
#'   the true value (`true_value`) if for the by group at least one condition
#'   evaluates to `TRUE` in one of the sources. It is set to the false value
#'   (`false_value`) if for the by group at least one observation exists and for
#'   all observations the condition evaluates to `FALSE` or `NA`. Otherwise, it
#'   is set to the missing value (`missing_value`).
#'
#'
#' @family der_gen
#' @keywords der_gen
#'
#' @export
#'
#' @examples
#'
#' library(dplyr)
#'
#' adsl <- tibble::tribble(
#'   ~USUBJID,
#'   "1",
#'   "2",
#'   "3",
#'   "4"
#' )
#'
#' cm <- tibble::tribble(
#'   ~USUBJID, ~CMCAT,        ~CMSEQ,
#'   "1",      "ANTI-CANCER",      1,
#'   "1",      "GENERAL",          2,
#'   "2",      "GENERAL",          1,
#'   "3",      "ANTI-CANCER",      1
#' )
#'
#' pr<- tibble::tribble(
#'   ~USUBJID, ~PRSEQ,
#'   "1",      1,
#'   "1",      2,
#'   "2",      1,
#'   "3",      1
#' )
#'
#' actual <- derive_var_merged_exist_flag_msrc(
#'   adsl,
#'   sources = list(
#'     flag_event(
#'       dataset_name = "cm",
#'       condition = CMCAT == "ANTI-CANCER"
#'     ),
#'     flag_event(
#'       dataset_name = "pr"
#'     )
#'   ),
#'   source_datasets = list(cm = cm, pr = pr),
#'   by_vars = exprs(USUBJID),
#'   new_var = CANCTRFL
#' )
#'
#' adex <- tibble::tribble(
#'   ~USUBJID, ~EXLNKID, ~EXADJ,
#'   "1",       "1",      "AE",
#'   "1",       "2",      NA_character_,
#'   "1",       "3",      NA_character_,
#'   "2",       "1",      NA_character_,
#'   "3",       "1",      NA_character_
#' )
#'
#' ec <- tibble::tribble(
#'   ~USUBJID, ~ECLNKID, ~ECADJ,
#'   "1",      "3",      "AE",
#'   "3",      "1",      NA_character_
#' )
#'
#' fa <- tibble::tribble(
#'   ~USUBJID, ~FALNKID, ~FATESTCD, ~FAOBJ,            ~FASTRESC,
#'   "3",      "1",      "OCCUR",   "DOSE ADJUSTMENT", "Y"
#' )
#'
#' actual <- derive_var_merged_exist_flag_msrc(
#'   adex,
#'   sources = list(
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
#'   by_vars = exprs(USUBJID, EXLNKID),
#'   new_var = CANCTRFL
#' )
derive_var_merged_exist_flag_msrc <- function(dataset,
                                              source_datasets,
                                              sources,
                                              by_vars,
                                              new_var,
                                              true_value = "Y",
                                              false_value = NA_character_,
                                              missing_value = NA_character_,
                                              filter_add = NULL) {
  new_var <- assert_symbol(enexpr(new_var))
  assert_list_of(source_datasets, class = "data.frame", named = TRUE)
  assert_list_of(sources, "flag_event")

  source_names <- names(source_datasets)
  assert_list_element(
    list = sources,
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
    sources,
    function(source) {
      data_source <- source_datasets[[source$dataset_name]]
      if (is.null(source$by_vars)) {
        source_by_vars <- by_vars
      } else {
        source_by_vars <- source$by_vars
      }
      if (is.null(source$condition)) {
        source_condition <- TRUE
      } else {
        source_condition <- source$condition
      }
      data_selected <- data_source %>%
        group_by(!!!unname(source_by_vars)) %>%
        mutate(!!tmp_cond_val := if_else(!!source_condition, TRUE, FALSE, FALSE)) %>%
        ungroup() %>%
        select(!!tmp_cond_val, !!!source_by_vars)
    }
  )

  derive_var_merged_exist_flag(
    dataset,
    dataset_add = bind_rows(selected_records) ,
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
#' `derive_var_merged_exist_flag_msrc()` function.
#'
#' @param dataset_name Dataset name of the dataset to be used as input for the
#'   event. The name refers to the dataset specified for `source_datasets` in
#'   `derive_var_merged_exist_flag()`.
#'
#'   *Permitted Values*: a character scalar
#'
#' @param condition Condition
#'
#'   The condition is evaluated at the dataset referenced by `dataset_name`. For
#'   all by groups where it evaluates as `TRUE` at least once the new variable
#'   is set to the true value (`true_value`).
#'

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
