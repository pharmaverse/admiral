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
#' @family high_order_function
#' @keywords high_order_function
#'
#' @author Stefan Bundfuss
#'
#' @seealso [params()] [slice_derivation()]
#'
#' @export
#'
#' @examples
#'
#' library(tibble)
#'
#' adlb <- tribble(
#'   ~USUBJID, ~AVISITN, ~AVAL, ~ABLFL,
#'   "1",            -1,   113, NA_character_,
#'   "1",             0,   113, "Y",
#'   "1",             3,   117, NA_character_,
#'   "2",             0,    95, "Y",
#'   "3",             0,   111, "Y",
#'   "3",             1,   101, NA_character_,
#'   "3",             2,   123, NA_character_
#' )
#'
#' # Derive BASE for post-baseline records only (derive_var_base() can not be used in this case
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
#' # Derive BASE for baseline and post-baseline records only
#' restrict_derivation(
#'   adlb,
#'   derivation = derive_var_base,
#'   args = params(
#'     by_vars = vars(USUBJID)
#'   ),
#'   filter = AVISITN >= 0
#' ) %>%
#'   # Derive CHG for post-baseline records only
#'   restrict_derivation(
#'     derivation = derive_var_chg,
#'     filter = AVISITN > 0
#'   )
restrict_derivation <- function(dataset,
                                derivation,
                                args = NULL,
                                filter) {
  # Check input
  assert_data_frame(dataset)
  assert_function(derivation, params = c("dataset"))
  assert_s3_class(args, "params", optional = TRUE)
  if (!is.null(args)) {
    assert_function_param(deparse(substitute(derivation)), names(args))
  }
  filter <- assert_filter_cond(enquo(filter))

  # Split input dataset
  data_ignore <- dataset %>%
    filter(!(!!filter) | is.na(!!filter))
  data_derive <- dataset %>%
    filter(!!filter)

  # Call derivation on subset
  call <- as.call(c(substitute(derivation), c(quote(data_derive), args)))
  data_derive <-
    eval(
      call,
      envir = list(data_derive = data_derive),
      enclos = parent.frame()
    )

  # Put datasets together again
  bind_rows(data_derive, data_ignore)
}
