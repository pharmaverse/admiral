#' Derive Basetype Variable
#'
#' @description
#' `r lifecycle::badge("deprecated")`
#'
#' This function is *deprecated*, please use `derive_basetype_records()` instead.
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
#'
#' @family deprecated
#'
#' @keywords deprecated
#'
#' @export
derive_var_basetype <- function(dataset, basetypes) {
  deprecate_stop("0.11.0", "derive_var_basetype()", "derive_basetype_records()")
}
