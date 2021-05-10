#' Derive variables that are based on data from a joined dataset
#'
#' Derives new variables which depend on either input dataset or joined dataset or both.
#' Additionally, pre-/post- filters can be supplied to filder the joined/resulting dataset.
#'
#' @param dataset Input dataset
#'
#'   The variables specified by the `by_vars` parameter are expected.
#'
#' @param dataset_join Dataset to join
#'
#'   The variables specified by the `by_vars` and the `pre_filter` are expected.
#'
#'
#' @param by_vars Grouping variables
#'
#'   Default: `exprs(STUDYID, USUBJID)`
#'
#'   Permitted Values: list of variables
#'
#' @param new_vars Variables to derive
#'
#'   Expressions specifying the variables definitions. See Examples.
#'
#' @param order Sort order
#'
#'   Within each by group the observations are ordered by the specified order.
#'
#'   Permitted Values: list of variables or functions of variables
#'
#' @param pre_filter Filter condition before joining the join dataset
#'
#'   Only observations of the joining dataset which fulfill the specified condition(-s)
#'   are used for joining.
#'
#'   Permitted Values: logical expression
#'
#' @param post_filter Filter condition for resulting dataset
#'
#'   Only observations of the joined dataset which fulfill the specified condition(-s)
#'   are returned.
#'
#'   Permitted Values: logical expression
#'
#' @return The input dataset with additional derived variables
#'
#' @keywords derivation adam
#'
#' @author Ondrej Slama
#'
#' @export
#'
#' @examples
#' data(ae)
#' data(ex)
#' ae %>%
#'   derive_vars_dtm(new_vars_prefix = "AST", dtc = AESTDTC) %>%
#'   derive_joined_vars(
#'     dataset_join = ex,
#'     new_vars = exprs(LDOSEDTM = compute_dtm(dtc = EXENDTC,
#'                                             date_imputation = NULL,
#'                                             time_imputation = NULL)),
#'     pre_filter = exprs(is_valid_dose()),
#'     by_vars = exprs(USUBJID),
#'     order = exprs(EXENDTC),
#'     post_filter = exprs(ASTDTM <= c(dtc = EXENDTC))
#'   )
#'
#' data(time_def)
#' ae %>%
#'   derive_vars_dtm(new_vars_prefix = "AST", dtc = AESTDTC) %>%
#'   derive_joined_vars(
#'     dataset_join = time_def,
#'     new_vars = exprs(APERIOD, APERIODC, TRTP, TRTA),
#'     by_vars = exprs(STUDYID, USUBJID),
#'     post_filter = exprs(STARTDTM <= ASTDTM && ASTDTM <= ENDDTM)
#'   )
derive_joined_vars <- function(dataset,
                               dataset_join,
                               by_vars = exprs(STUDYID, USUBJID),
                               new_vars,
                               order = NULL,
                               pre_filter = NULL,
                               post_filter = NULL) {

  assert_has_variables(dataset, map_chr(by_vars, as_string))
  assert_has_variables(dataset_join, map_chr(by_vars, as_string))
  stopifnot(is.list(new_vars) && vapply(new_vars, rlang::is_expression, logical(1)))

  cols_dataset <- colnames(dataset)

  if (!is.null(pre_filter)) {
    dataset_join <- filter(dataset_join, !!!pre_filter)
  }

  dataset_res <- dataset %>%
    left_join(dataset_join, by = map_chr(by_vars, as_string)) %>%
    mutate(!!!new_vars)

  if (!is.null(post_filter)) {
    dataset_res <- filter(dataset_res, !!!post_filter)
  }

  if (!is.null(order)) {
    dataset_res <- dataset_res %>%
      group_by(!!!by_vars) %>%
      arrange(!!!order, .by_group = TRUE)
  }

  if (is.null(names(new_vars)) || all(names(new_vars) == "")) {
    names(new_vars) <- unname(map_chr(new_vars, as_string))
  }

  return(select(dataset_res, c(cols_dataset, names(new_vars))))
}
