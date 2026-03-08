#' Derive Basetype Variable
#'
#' Baseline Type `BASETYPE` is needed when there is more than one definition of
#' baseline for a given Analysis Parameter `PARAM` in the same dataset.  For a
#' given parameter, if Baseline Value `BASE` or `BASEC` are derived and there
#' is more than one definition of baseline, then `BASETYPE` must be non-null on
#' all records of any type for that parameter where either `BASE` or `BASEC`
#' are also non-null. Each value of `BASETYPE` refers to a definition of
#' baseline that characterizes the value of `BASE` on that row.  Please see
#' section 4.2.1.6 of the ADaM Implementation Guide, version 1.3 for further
#' background.
#'
#' Adds the `BASETYPE` variable to a dataset and duplicates records based upon
#' the provided conditions.
#'
#' @param dataset
#'   `r roxygen_param_dataset(expected_vars = c("basetypes"))`
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
#'
#' @family der_bds_findings
#'
#' @keywords der_bds_findings
#'
#' @export
#'
#' @examplesx
#'
#' @caption Derive `BASETYPE` based on epoch (`basetypes`)
#' @info The `basetypes` argument is a named list of expressions where each name
#' becomes a value of `BASETYPE` and each expression defines which records
#' receive that value. A record can match multiple expressions and will be
#' duplicated once for each matching `BASETYPE`.
#' @code
#' library(tibble)
#' library(dplyr, warn.conflicts = FALSE)
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
#' derive_basetype_records(
#'   dataset = bds,
#'   basetypes = exprs(
#'     "RUN-IN" = EPOCH %in% c("RUN-IN", "STABILIZATION", "DOUBLE-BLIND", "OPEN-LABEL"),
#'     "DOUBLE-BLIND" = EPOCH %in% c("DOUBLE-BLIND", "OPEN-LABEL"),
#'     "OPEN-LABEL" = EPOCH == "OPEN-LABEL"
#'   )
#' )
#'
#' @caption Records not matching any condition are retained with `BASETYPE = NA`
#' @info Records that do not match any condition in `basetypes` are kept in the
#' output dataset with `BASETYPE` set to `NA`. In this example, `SCREENING`
#' records do not match any of the `basetypes` conditions and are therefore
#' retained with `BASETYPE = NA`.
#' @code
#' bds <- tribble(
#'   ~USUBJID, ~EPOCH,         ~PARAMCD,  ~ASEQ, ~AVAL,
#'   "P01",    "SCREENING",    "PARAM01",     1,  10.2,
#'   "P01",    "RUN-IN",       "PARAM01",     2,  10.0,
#'   "P01",    "RUN-IN",       "PARAM01",     3,   9.8,
#'   "P01",    "DOUBLE-BLIND", "PARAM01",     4,   9.2,
#'   "P01",    "DOUBLE-BLIND", "PARAM01",     5,  10.1,
#'   "P02",    "SCREENING",    "PARAM01",     1,  12.2,
#'   "P02",    "RUN-IN",       "PARAM01",     2,  12.1,
#'   "P02",    "DOUBLE-BLIND", "PARAM01",     3,  10.2
#' )
#'
#' derive_basetype_records(
#'   dataset = bds,
#'   basetypes = exprs(
#'     "RUN-IN" = EPOCH %in% c("RUN-IN", "DOUBLE-BLIND"),
#'     "DOUBLE-BLIND" = EPOCH == "DOUBLE-BLIND"
#'   )
#' )
#'
#' @caption Include all records for multiple baseline type derivations (`basetypes = TRUE`)
#' @info When all parameter records need to be included for multiple baseline
#' type derivations (such as `"LAST"` and `"WORST"`), set each expression in
#' `basetypes` to `TRUE`. This duplicates every record once for each named
#' baseline type.
#' @code
#' bds <- tribble(
#'   ~USUBJID, ~EPOCH,         ~PARAMCD,  ~ASEQ, ~AVAL,
#'   "P01",    "RUN-IN",       "PARAM01",     1,  10.0,
#'   "P01",    "RUN-IN",       "PARAM01",     2,   9.8,
#'   "P01",    "DOUBLE-BLIND", "PARAM01",     3,   9.2,
#'   "P01",    "DOUBLE-BLIND", "PARAM01",     4,  10.1
#' )
#'
#' derive_basetype_records(
#'   dataset = bds,
#'   basetypes = exprs(
#'     "LAST" = TRUE,
#'     "WORST" = TRUE
#'   )
#' )
derive_basetype_records <- function(dataset, basetypes) {
  assert_data_frame(dataset)
  assert_expr_list(basetypes, named = TRUE)

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
