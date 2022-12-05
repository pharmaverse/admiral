#' Derive Basetype Variable
#'
#' Baseline Type `BASETYPE` is needed when there is more than one definition of
#' baseline for a given Analysis Parameter `PARAM` in the same dataset.  For a
#' given parameter, if Baseline Value `BASE` is populated, and there is more than
#' one definition of baseline, then `BASETYPE` must be non-null on all records of
#' any type for that parameter. Each value of `BASETYPE` refers to a definition of
#' baseline that characterizes the value of `BASE` on that row.  Please see
#' section 4.2.1.6 of the ADaM Implementation Guide, version 1.3 for further
#' background.
#'
#' Adds the `BASETYPE` variable to a dataset and duplicates records based upon
#' the provided conditions.
#'
#' @param dataset Input dataset
#'
#'   The columns specified in the expressions inside `basetypes` are required.
#'
#' @param basetypes A *named* list of expressions created using the
#' `rlang::exprs()` function
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
#' @return The input dataset with variable `BASETYPE` added
#'
#' @author Thomas Neitmann
#'
#' @family der_bds_findings
#'
#' @keywords der_bds_findings
#'
#' @export
#'
#' @examples
#' library(tibble)
#' library(dplyr, warn.conflicts = FALSE)
#' library(rlang)
#'
#' bds <- tribble(
#'   ~USUBJID, ~EPOCH,         ~PARAMCD,  ~ASEQ, ~AVAL,
#'   "P01",    "RUN-IN",       "PARAM01",     1,  10.0,
#'   "P01",    "RUN-IN",       "PARAM01",     2,   9.8,
#'   "P01",    "DOUBLE-BLIND", "PARAM01",     3,   9.2,
#'   "P01",    "DOUBLE-BLIND", "PARAM01",     4,  10.1,
#'   "P01",    "OPEN-LABEL",   "PARAM01",     5,  10.4,
#'   "P01",    "OPEN-LABEL",   "PARAM01",     6,   9.9,
#'   "P02",    "RUN-IN",       "PARAM01",     1,  12.1,
#'   "P02",    "DOUBLE-BLIND", "PARAM01",     2,  10.2,
#'   "P02",    "DOUBLE-BLIND", "PARAM01",     3,  10.8,
#'   "P02",    "OPEN-LABEL",   "PARAM01",     4,  11.4,
#'   "P02",    "OPEN-LABEL",   "PARAM01",     5,  10.8
#' )
#'
#' bds_with_basetype <- derive_var_basetype(
#'   dataset = bds,
#'   basetypes = exprs(
#'     "RUN-IN" = EPOCH %in% c("RUN-IN", "STABILIZATION", "DOUBLE-BLIND", "OPEN-LABEL"),
#'     "DOUBLE-BLIND" = EPOCH %in% c("DOUBLE-BLIND", "OPEN-LABEL"),
#'     "OPEN-LABEL" = EPOCH == "OPEN-LABEL"
#'   )
#' )
#'
#'
#' # Below print statement will print all 23 records in the data frame
#' # bds_with_basetype
#' print(bds_with_basetype, n = Inf)
#'
#' count(bds_with_basetype, BASETYPE, name = "Number of Records")
#'
#' # An example where all parameter records need to be included for 2 different
#' # baseline type derivations (such as LAST and WORST)
#' bds <- tribble(
#'   ~USUBJID, ~EPOCH,         ~PARAMCD,  ~ASEQ, ~AVAL,
#'   "P01",    "RUN-IN",       "PARAM01",     1,  10.0,
#'   "P01",    "RUN-IN",       "PARAM01",     2,   9.8,
#'   "P01",    "DOUBLE-BLIND", "PARAM01",     3,   9.2,
#'   "P01",    "DOUBLE-BLIND", "PARAM01",     4,  10.1
#' )
#'
#' bds_with_basetype <- derive_var_basetype(
#'   dataset = bds,
#'   basetypes = exprs(
#'     "LAST" = TRUE,
#'     "WORST" = TRUE
#'   )
#' )
#'
#' print(bds_with_basetype, n = Inf)
#'
#' count(bds_with_basetype, BASETYPE, name = "Number of Records")
derive_var_basetype <- function(dataset, basetypes) {
  assert_data_frame(dataset)
  assert_named_exprs(basetypes)

  records_with_basetype <- map2(names(basetypes), basetypes, function(label, condition) {
    dataset %>%
      filter(!!condition) %>%
      mutate(BASETYPE = label)
  }) %>%
    bind_rows()

  complementary_condition <- Reduce(function(x, y) bquote(.(x) | .(y)), basetypes)
  records_without_basetype <- filter(dataset, !(!!complementary_condition))

  bind_rows(records_without_basetype, records_with_basetype)
}
