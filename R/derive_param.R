#' Derive PARAM and PARAMCD
#'
#' @param dataset The input dataset
#' @param mapping A dataset containing the mapping of one or more SDTM variables
#'        to `PARAM` and `PARAMCD`
#' @param by_vars The variables used to join `dataset` and `mapping`
#'
#' @author Thomas Neitmann
#'
#' @return
#' @export The input dataset with additional columns `PARAM` and `PARAMCD`
#'
#' @examples
#' data(vs)
#' data(advs_param)
#'
#' derive_param(vs, advs_param, by_vars = c("VSTESTCD", "VSPOS"))
#'
#' ## Disregard the `VSPOS` variable for mapping
#' derive_param(vs, advs_param, by_vars = c("VSTESTCD"))
derive_param <- function(dataset, mapping, by_vars) {
  assert_has_variables(dataset, by_vars)
  assert_has_variables(mapping, c("PARAM", "PARAMCD", by_vars))

  keys <- colnames(mapping)[!colnames(mapping) %in% c("PARAM", "PARAMCD")]

  mapping <- mapping %>%
    arrange(!!!syms(keys)) %>%
    group_by(!!!syms(by_vars)) %>%
    slice(1L) %>%
    select(PARAM, PARAMCD, !!!syms(by_vars))

  result <- left_join(dataset, mapping, by = by_vars)

  warn_if_param_missing(result, by_vars)

  result
}
