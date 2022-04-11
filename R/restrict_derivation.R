#' Execute a Derivation on a Subset of the Input Dataset
#'
#' Execute a derivation on a subset of the input dataset.
#'
#' @param dataset Input dataset
#'
#' @param derivation Derivation
#'
#' @param args Arguments of the derivation
#'
#'   A `params()` object is expected.
#'
#' @param filter Filter condition
#'
#' @keywords user_utility high_order_function
#'
#' @author Stefan Bundfuss
#'
#' @seealso [params()] [slice_derivation()]
#'
#' @export
#'
#' @examples
#'
#' library(magrittr)
#' adlb <- tibble::tribble(
#' ~USUBJID, ~AVISITN, ~AVAL, ~ABLFL,
#' "1",      -1,       113,   NA_character_,
#' "1",       0,       113,   "Y",
#' "1",       3,       117,   NA_character_,
#' "2",       0,        95,   "Y",
#' "3",       0,       111,   "Y",
#' "3",       1,       101,   NA_character_,
#' "3",       2,       123,   NA_character_
#' )
#'
#' # derive BASE for post-baseline records only (derive_var_base() can not be used in this case
#' # as it requires the baseline observation to be in the input dataset)
#' restrict_derivation(
#'   adlb,
#'   derivation = derive_vars_merged,
#'   args = params(
#'     by_vars = vars(USUBJID),
#'     dataset_add = adlb,
#'     filter_add = ABLFL == "Y",
#'     new_vars = vars(BASE = AVAL)
#'   ),
#'   filter = AVISITN > 0
#' )
#'
#' # derive BASE for baseline and post-baseline records only
#' restrict_derivation(
#'   adlb,
#'   derivation = derive_var_base,
#'   args = params(
#'     by_vars = vars(USUBJID)
#'   ),
#'   filter = AVISITN >= 0
#' ) %>%
#'
#' # derive CHG for post-baseline records only
#' restrict_derivation(
#'   derivation = derive_var_chg,
#'   filter = AVISITN > 0
#' )
#'
restrict_derivation <- function(dataset,
                                derivation,
                                args = NULL,
                                filter) {
  # check input
  assert_data_frame(dataset)
  assert_function(derivation, params = c("dataset"))
  assert_s3_class(args, "params", optional = TRUE)
  if (!is.null(args)) {
    assert_function_param(deparse(substitute(derivation)), names(args))
  }
  filter <- assert_filter_cond(enquo(filter))

  # split input dataset
  data_ignore <- dataset %>%
    filter(!(!!filter) | is.na(!!filter))
  data <- dataset %>%
    filter(!!filter)

  # call derivation on subset
  call <- as.call(c(substitute(derivation), c(quote(data), args)))
  data <- eval(call, envir = list(data = data), enclos = parent.frame())

  # put datasets together again
  bind_rows(data, data_ignore)
}
