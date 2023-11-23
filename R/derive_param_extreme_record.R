#' Adds a Parameter Based on First or Last Record from Multiple Sources
#'
#' @description
#' `r lifecycle::badge("superseded")` The `derive_param_extreme_record()`
#' function has been superseded in favor of `derive_extreme_event()`.
#'
#' Generates parameter based on the first or last observation from multiple
#' source datasets, based on user-defined filter, order and by group criteria.
#' All variables of the selected observation are kept.
#'
#' @param dataset `r roxygen_param_dataset()`
#'
#' @param sources Sources
#'
#'    A list of `records_source()` objects is expected.
#'
#' @param source_datasets Source datasets
#'
#'   A named list of datasets is expected. The `dataset_name` field of
#'   `records_source()` refers to the dataset provided in the list. The variables
#'   specified by the `order` and the `by_vars` arguments are expected after applying `new_vars`.
#'
#' @param by_vars Grouping variables
#'
#'   If the argument is specified, for each by group the observations are
#'   selected separately.
#'
#'   `r roxygen_param_by_vars()`
#'
#' @param order Sort order
#'
#'   If the argument is set to a non-null value, for each by group the first or
#'   last observation from the source datasets is selected with respect to
#'   the specified order. Variables created via `new_vars` e.g., imputed date variables,
#'   can be specified as well (see examples below).
#'
#'   Please note that `NA` is considered as the last value. I.e., if a order
#'   variable is `NA` and `mode = "last"`, this observation is chosen while for
#'   `mode = "first"` the observation is chosen only if there are no
#'   observations where the variable is not `NA`.
#'
#'   *Permitted Values:* list of expressions created by `exprs()`, e.g.,
#'   `exprs(ADT, desc(AVAL))`
#'
#' @param mode Selection mode (first or last)
#'
#'   If `"first"` is specified, for each by group the first observation with
#'   respect to `order` is included in the output dataset. If `"last"` is
#'   specified, the last observation is included in the output dataset.
#'
#'   Permitted Values:  `"first"`, `"last"`
#'
#' @param set_values_to Variables to be set
#'
#'   The specified variables are set to the specified values for the new
#'   observations.
#'
#'   A list of variable name-value pairs is expected.
#'   + LHS refers to a variable.
#'   + RHS refers to the values to set to the variable. This can be a string, a
#'   symbol, a numeric value or `NA`, e.g., `exprs(PARAMCD = "PD", PARAM =
#'   "First Progressive Disease")`.
#'
#' @details The following steps are performed to create the output dataset:
#'
#'   \enumerate{
#'   \item For each source dataset the observations as specified by
#'   the `filter` element are selected.
#'
#'   \item Variables specified by `new_vars` are created for each source dataset.
#'
#'   \item The first or last observation (with respect to the
#'   `order` variable) for each by group (specified by `by_vars`) from multiple sources
#'   is selected and added to the input dataset. }
#'
#' @return
#' The input dataset with the first or last observation of each by group
#' added as new observations.
#'
#' @family superseded
#' @keywords superseded
#'
#' @export
#' @examples
#' aevent_samp <- tibble::tribble(
#'   ~USUBJID, ~PARAMCD,                       ~PARAM,     ~RSSTDTC,
#'   "1",          "PD",  "First Progressive Disease", "2022-04-01",
#'   "2",          "PD",  "First Progressive Disease", "2021-04-01",
#'   "3",          "PD",  "First Progressive Disease", "2023-04-01"
#' )
#'
#' cm <- tibble::tribble(
#'   ~STUDYID, ~USUBJID, ~CMDECOD,     ~CMSTDTC,
#'   "1001",        "1",    "ACT", "2021-12-25"
#' )
#'
#' pr <- tibble::tribble(
#'   ~STUDYID, ~USUBJID, ~PRDECOD,     ~PRSTDTC,
#'   "1001",        "1",    "ACS", "2021-12-27",
#'   "1001",        "2",    "ACS", "2020-12-25",
#'   "1001",        "3",    "ACS", "2022-12-25",
#' )
#' derive_param_extreme_record(
#'   dataset = aevent_samp,
#'   sources = list(
#'     records_source(
#'       dataset_name = "cm",
#'       filter = CMDECOD == "ACT",
#'       new_vars = exprs(
#'         ADT = convert_dtc_to_dt(CMSTDTC),
#'         AVALC = CMDECOD
#'       )
#'     ),
#'     records_source(
#'       dataset_name = "pr",
#'       filter = PRDECOD == "ACS",
#'       new_vars = exprs(
#'         ADT = convert_dtc_to_dt(PRSTDTC),
#'         AVALC = PRDECOD
#'       )
#'     )
#'   ),
#'   source_datasets = list(cm = cm, pr = pr),
#'   by_vars = exprs(USUBJID),
#'   order = exprs(ADT),
#'   mode = "first",
#'   set_values_to = exprs(
#'     PARAMCD = "FIRSTACT",
#'     PARAM = "First Anti-Cancer Therapy"
#'   )
#' )
derive_param_extreme_record <- function(dataset = NULL,
                                        sources,
                                        source_datasets,
                                        by_vars = NULL,
                                        order,
                                        mode,
                                        set_values_to) {
  # Check arguments assertions
  assert_data_frame(dataset, optional = TRUE)
  assert_list_of(sources, "records_source")
  assert_list_of(source_datasets, "data.frame")
  assert_vars(by_vars, optional = TRUE)
  assert_character_scalar(
    mode,
    values = c("first", "last"),
    case_sensitive = FALSE
  )
  assert_varval_list(set_values_to, accept_expr = TRUE, optional = TRUE)

  source_names <- names(source_datasets)

  # Create Empty list to contain source datasets
  data_list <- vector("list", length(sources))

  # Evaluate the expressions contained in the sources
  for (i in seq_along(sources)) {
    source_dataset <- source_datasets[[sources[[i]]$dataset_name]]
    new_vars_colnames <- replace_values_by_names(sources[[i]]$new_vars)
    data_list[[i]] <- source_dataset %>%
      filter_if(sources[[i]]$filter) %>%
      mutate(!!!sources[[i]]$new_vars) %>%
      select(!!!by_vars, !!!new_vars_colnames)
  }

  # Bind the source datasets together and parse out the extreme value
  param_data <- bind_rows(data_list) %>%
    filter_extreme(.,
      by_vars = by_vars,
      order = order,
      mode = mode
    ) %>%
    process_set_values_to(
      set_values_to
    )

  # Bind the parameter rows back to original dataset
  bind_rows(dataset, param_data)
}

#' Create a `records_source` Object
#'
#' The `records_source` object is used to find extreme records of interest.
#'
#' @param dataset_name The name of the source dataset
#'
#'   The name refers to the dataset provided by the `source_datasets` argument
#'   of `derive_param_extreme_record()`.
#'
#' @param filter An unquoted condition for selecting the observations from
#'   `dataset`.
#'
#' @param new_vars Variables to add
#'
#'   The specified variables from the source datasets are added to the output
#'   dataset. Variables can be renamed by naming the element, i.e., `new_vars =
#'   exprs(<new name> = <old name>)`.
#'
#'   For example `new_vars = exprs(var1, var2)` adds variables `var1` and `var2`
#'   from to the input dataset.
#'
#'   And `new_vars = exprs(var1, new_var2 = old_var2)` takes `var1` and
#'   `old_var2` from the source dataset and adds them to the input dataset renaming
#'   `old_var2` to `new_var2`. Expressions can be used to create new variables
#'   (see for example `new_vars` argument in `derive_vars_merged()`).
#'
#'   *Permitted Values:* list of expressions created by `exprs()`, e.g.,
#'   `exprs(ADT, desc(AVAL))`
#'
#' @keywords source_specifications
#' @family source_specifications
#'
#' @seealso [derive_param_extreme_record()]
#'
#' @return  An object of class `records_source`
#' @export
records_source <- function(dataset_name,
                           filter = NULL,
                           new_vars) {
  out <- list(
    dataset_name = assert_character_scalar(dataset_name),
    filter = assert_filter_cond(enexpr(filter), optional = TRUE),
    new_vars = assert_expr_list(new_vars)
  )
  class(out) <- c("records_source", "source", "list")
  out
}
