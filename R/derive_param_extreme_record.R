#' Adds a Parameter based on first or last record from multiple sources
#'
#' Adds a Parameter based on first or last record from multiple sources
#'
#' @param dataset Input dataset
#'
#'   The variables specified by the `order` and the `by_vars` parameter are
#'   expected.
#'
#' @param sources Sources
#'
#'    A list of `records_source()` objects is expected.
#'
#' @param source_datasets Source datasets
#'
#'   A named list of datasets is expected. The `dataset_name` field of
#'   `records_source()` refers to the dataset provided in the list.
#'
#' @param by_vars By variables
#'
#'   If the parameter is specified, for each by group the observations are
#'   selected separately.
#'
#' @param order Sort order
#'
#'   If the argument is set to a non-null value, for each by group the first or
#'   last observation from the additional dataset is selected with respect to
#'   the specified order. The imputed date variable can be specified as well
#'   (see examples below).
#'
#'   Please note that `NA` is considered as the last value. I.e., if a order
#'   variable is `NA` and `mode = "last"`, this observation is chosen while for
#'   `mode = "first"` the observation is chosen only if there are no
#'   observations where the variable is not 'NA'.
#'
#'   *Default*: `NULL`
#'
#'   *Permitted Values*: list of variables or `desc(<variable>)` function calls
#'   created by `exprs()`, e.g., `exprs(ADT, desc(AVAL)` or `NULL`
#'
#' @param mode Selection mode (first or last)
#'
#'   If `"first"` is specified, for each subject the first observation with
#'   respect to the date is included in the output dataset. If `"last"` is
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
#'   symbol, a numeric value or `NA`, e.g., `exprs(PARAMCD = "TDOSE", PARCAT1 =
#'   "OVERALL")`. More general expression are not allowed.
#'
#' @details The following steps are performed to create the output dataset:
#'
#'   \enumerate{
#'   \item For each source dataset the observations as specified by
#'   the `filter` element are selected. Then for each patient the first or last
#'   observation (with respect to `order`) is selected.
#'
#'   \item The `ADT` variable is set to the variable specified by the
#'   \code{date} element. If the date variable is a datetime variable, only
#'   the datepart is copied. If the source variable is a character variable, it
#'   is converted to a date. If the date is incomplete, it is imputed as
#'   the first possible date.
#'
#'   \item For each patient the first or last observation (with respect to the
#'   `ADT` variable) from multiple sources is selected and added back to the original
#'   dataset. }
#'
#' @return
#' The input dataset with the first or last observation of each by group
#' added as new observations.
#'
#' @family der_prm_bds_findings
#' @keywords der_prm_bds_findings
#'
#' @export
#' @examples
#' aevent_samp <- tibble::tribble(
#'   ~USUBJID, ~AESEQ, ~AETERM,     ~AESTDTC,
#'   "1"     ,      1,     "X", "2022-01-01",
#'   "1"     ,      2,     "Y", "2022-01-02",
#'   "1"     ,      3,     "Z", "2022-01-03",
#'   "2"     ,      1,     "A", "2021-01-01",
#'   "2"     ,      2,     "B", "2021-01-02",
#'   "2"     ,      3,     "C", "2021-01-03",
#'   "3"     ,      1,     "J", "2023-01-01",
#'   "3"     ,      2,     "K", "2023-01-02",
#'   "3"     ,      3,     "L", "2023-01-03"
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
#'         AVALC = CMDECOD)
#'     ),
#'     records_source(
#'       dataset_name = "pr",
#'       filter = PRDECOD == "ACS",
#'       new_vars = exprs(
#'         ADT = convert_dtc_to_dt(PRSTDTC),
#'         AVALC = PRDECOD)
#'     )),
#'   source_datasets = list(cm = cm, pr = pr),
#'   by_vars = exprs(USUBJID),
#'   order = exprs(ADT),
#'   mode = "first",
#'   set_values_to = exprs(
#'     PARAMCD = "FIRSTACT",
#'     PARAM = "First Anti-Cancer Therapy"
#'   )
#' )
derive_param_extreme_record <- function(dataset,
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
    filter_extreme(.,
                   by_vars = by_vars,
                   order = order,
                   mode = mode) %>%
    mutate(!!!set_values_to)

  # Bind the parameter rows back to original adevent dataset
  data <- bind_rows(dataset, param_data)
  return(data)
}
#' Create a `records_source` Object
#'
#' The `records_source` object is used to find extreme records of interest.
#'
#' @param dataset_name The name of the source dataset
#'
#'   The name refers to the dataset provided by the `source_datasets` parameter
#'   of `derive_param_extreme_record()`.
#'
#' @param filter An unquoted condition for selecting the observations from
#'   `dataset` which are events or possible censoring time points.
#'
#' @param new_vars Variables to add
#'
#'   The specified variables from the additional dataset are added to the output
#'   dataset. Variables can be renamed by naming the element, i.e., `new_vars =
#'   exprs(<new name> = <old name>)`.
#'
#'   For example `new_vars = exprs(var1, var2)` adds variables `var1` and `var2`
#'   from `dataset_add` to the input dataset.
#'
#'   And `new_vars = exprs(var1, new_var2 = old_var2)` takes `var1` and
#'   `old_var2` from `dataset_add` and adds them to the input dataset renaming
#'   `old_var2` to `new_var2`.
#'
#'   If the argument is not specified or set to `NULL`, all variables from the
#'   additional dataset (`dataset_add`) are added.
#'
#'   *Permitted Values*: list of variables created by `exprs()`
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
                           new_vars = NULL) {
  out <- list(
    dataset_name = assert_character_scalar(dataset_name),
    filter = assert_filter_cond(enexpr(filter), optional = TRUE),
    new_vars = assert_expr_list(new_vars, optional = TRUE)
  )
  class(out) <- c("records_source", "source", "list")
  out
}
