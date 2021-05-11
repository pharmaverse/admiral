#' Merge variables from more than one dataset
#'
#' Merges variables from more than one dataset. For each dataset the
#' observations to merge can be selected by a condition and/or by selecting the
#' first or last observation in each by group.
#'
#' @param dataset Input dataset
#'
#'   The variables specified by the `by_vars` parameter are expected.
#'
#' @param sources List of sources
#'
#'   For each source the observations are selected as specified. The selected
#'   observations of all sources are collected in one dataset and then for each
#'   by group one observation (with respect to the `order` and `mode` parameter)
#'   is selected and merged to the input dataset.
#'
#'   Each element of the list must be a list with the following named elements.
#'
#'   \itemize{
#'   \item \emph{dataset:}
#'     Dataset to add
#'
#'   The variables specified by the \code{by_vars} parameter, the \code{filter}, and the
#'   \code{order} element are expected.
#'
#'   \item \emph{filter:} Filter condition for dataset
#'
#'   Only observations of the dataset which fulfill the specified condition
#'   are used for merging.
#'
#'   \emph{Permitted Values:} logical expression
#'
#'   \item \emph{new_vars:} Variable to add
#'
#'   The specified variables are created in the add dataset and then merge to
#'   the input dataset.
#'
#'   \emph{Permitted Values:} list of name-value pairs
#'
#'   \item \emph{order:} Sort order
#'
#'   If the parameter is specified, the source dataset is ordered by the
#'   specified order and only the first or last observation (depending on the
#'   mode) in each by group is used for merging.
#'
#'   \emph{Permitted Values:} list of variables or functions of variables
#'
#'   \item \emph{mode:} mode
#'
#'   If the \code{filter_order} parameter is specified, the mode determines if
#'   the first or last observation of each by group is selected.
#'
#'   \emph{Permitted Values:} \code{"first"}, \code{"last"}
#'   }
#'
#' @param order Sort order
#'
#'   If the parameter is specified, all selected observations from the source
#'   dataset are ordered by the specified order and only the first or last
#'   observation (depending on the filter mode) in each by group is used for
#'   merging.
#'
#'   The temporary variable `temp_source_nr`, which is set to the number of the
#'   source, can be used for ordering. It is not included in the output dataset.
#'
#'   Permitted Values: list of variables or functions of variables
#'
#' @param mode Filter mode
#'
#'   If the `order` parameter is specified, the filter mode determines if the
#'   first or last observation of each by group is selected.
#'
#'   Permitted Values: `"first"`, `"last"`
#'
#' @param by_vars Grouping variables
#'
#'   Default: `exprs(USUBJID)`
#'
#'   Permitted Values: list of variables
#'
#' @details For each group (with respect to the variables specified for the
#'   `by_vars` parameter) the first observation (with respect to the order
#'   specified for the `order` parameter) is included in the output dataset.
#'
#' @author Stefan Bundfuss
#'
#' @return The input dataset with variables from the source datasets
#'
#' @keywords derivation adam
#'
#' @export
#'
#' @examples
#' library(dplyr, warn.conflict = FALSE)
#' library(rlang, warn.conflict = FALSE)
#' library(stringr)
#' data("dm")
#' data("ds")
#' data("ae")
#' data("lb")
#'
#' # derive last known alive date as maximum AE start date, AE end date,
#' # and lab test date
#' ae_start <- list(dataset = ae,
#'                  vars = exprs(LSTALVDT = convert_dtc_to_dt(AESTDTC)),
#'                  filter = expr(str_length(AESTDTC) >= 10))
#' ae_end <- list(dataset = ae,
#'                vars = exprs(LSTALVDT = convert_dtc_to_dt(AEENDTC)),
#'                filter = expr(str_length(AEENDTC) >= 10))
#' lb_date <- list(dataset = lb,
#'                 vars = exprs(LSTALVDT = convert_dtc_to_dt(LBDTC)),
#'                 filter = expr(str_length(LBDTC) >= 10))
#'
#' derive_merged_multisource(dm,
#'                           list(ae_start, ae_end, lb_date),
#'                           order = exprs(LSTALVDT),
#'                           mode = "last") %>%
#'   select(USUBJID, LSTALVDT)
#'
#' # derive death cause and domain from AE and DS dataset
#' dthae <- list(dataset = ae,
#'               vars = exprs(DTHDOM = "AE", DTHCAUS = AEDECOD),
#'               filter = expr(AEOUT == "FATAL"))
#' dthds <- list(dataset = ds,
#'               vars = exprs(DTHDOM = "DS", DTHCAUS = DSTERM),
#'               filter = expr(DSDECOD == "DEATH"))
#' derive_merged_multisource(dm,
#'                           list(dthae, dthds),
#'                           order = exprs(temp_source_nr),
#'                           mode = "first") %>%
#'   filter(!is.null(DTHDOM)) %>%
#'   select(USUBJID, DTHCAUS, DTHDOM)
derive_merged_multisource <- function(dataset,
                                      sources,
                                      by_vars = exprs(USUBJID),
                                      order,
                                      mode){
  all_data <- filter_extreme_multisource(sources = sources,
                                         by_vars = by_vars,
                                         order = order,
                                         mode = mode)
  left_join(dataset, all_data, by = map_chr(by_vars, as_string)) %>%
    select(-starts_with("temp_"))
}
