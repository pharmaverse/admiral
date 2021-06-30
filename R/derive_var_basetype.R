#' Derive BASETYPE variable
#'
#' Adds the `BASETYPE` variable to a dataset and duplicates records based upon
#' the provided conditions
#'
#' @param dataset Input dataset
#'
#'   The columns specified in the expressions inside `basetypes` are required.
#'
#' @param basetypes A *named* list of expressions created using `exprs()`
#'
#'   The names corresponds to the values of the newly created `BASETYPE` variables
#'   and the expressions are used to subset the input dataset.
#'
#' @details
#' For each element of `basetypes` the input dataset is subset based upon
#' the provided expression and the `BASETYPE` variable is set to the name of the
#' expression. Then, all subsets are stacked. Records which do not match any
#' condition are kept and `BASETYPE` is set to `NA`.
#'
#' @author Thomas Neitmann
#'
#' @keywords bds derivation
#'
#' @export
#'
#' @examples
#' bds <- tibble::tribble(
#'   ~USUBJID, ~EPOCH,         ~PARAMCD,  ~ASEQ, ~AVAL,
#'   "P01",    "RUN-IN",       "PARAM01", 1,     10,
#'   "P01",    "RUN-IN",       "PARAM01", 2,      9.8,
#'   "P01",    "DOUBLE-BLIND", "PARAM01", 3,      9.2,
#'   "P01",    "DOUBLE-BLIND", "PARAM01", 4,     10.1,
#'   "P01",    "OPEN-LABEL",   "PARAM01", 5,     10.4,
#'   "P01",    "OPEN-LABEL",   "PARAM01", 6,      9.9,
#'   "P02",    "RUN-IN",       "PARAM01", 1,     12.1,
#'   "P02",    "DOUBLE-BLIND", "PARAM01", 2,     10.2,
#'   "P02",    "DOUBLE-BLIND", "PARAM01", 3,     10.8,
#'   "P02",    "OPEN-LABEL",   "PARAM01", 4,     11.4,
#'   "P02",    "OPEN-LABEL",   "PARAM01", 5,     10.8
#' )
#'
#' derive_var_basetype(
#'   dataset = bds,
#'   basetypes = exprs(
#'     "RUN-IN" = EPOCH %in% c("RUN-IN", "STABILIZATION", "DOUBLE-BLIND", "OPEN-LABEL"),
#'     "DOUBLE-BLIND" = EPOCH %in% c("DOUBLE-BLIND", "OPEN-LABEL"),
#'     "OPEN-LABEL" = EPOCH == "OPEN-LABEL"
#'   )
#' )
derive_var_basetype <- function(dataset, basetypes) {
  assert_that(
    is.data.frame(dataset),
    is_named_exprs(basetypes)
  )
  assert_has_variables(
    dataset,
    unique(map_chr(basetypes, all.vars))
  )

  records_with_basetype <- map2(names(basetypes), basetypes, function(label, condition) {
    dataset %>%
      filter(!!condition) %>%
      mutate(BASETYPE = label)
  }) %>%
    bind_rows()

  records_without_basetype <- anti_join(dataset, records_with_basetype, by = colnames(dataset))

  bind_rows(records_without_basetype, records_with_basetype)
}
