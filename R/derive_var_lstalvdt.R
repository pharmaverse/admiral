#' Merge the first or last date from more than one dataset
#'
#' Merges the first or last date from more than one dataset. For each dataset
#' the observations to merge can be selected by a condition and/or by selecting
#' the first or last observation in each by group.
#'
#' @param dataset Input dataset
#'
#'   The variables specified by the `by_vars` parameter are expected.
#'
#' @param sources List of sources
#'
#'   For each source the observations are selected as specified. The selected
#'   observations of all sources are collected in one dataset and then for each
#'   by group one observation (with respect to the `mode` parameter) is selected
#'   and the date is merged to the input dataset.
#'
#'   Each element of the list must be a list with the following named elements.
#'
#'   \itemize{ \item \emph{dataset:} Dataset to add
#'
#'   The variables specified by the \code{by_vars} parameter, the \code{filter},
#'   and the \code{order} element are expected.
#'
#'   \item \emph{filter:} Filter condition for dataset
#'
#'   Only observations of the dataset which fulfill the specified condition are
#'   used for merging.
#'
#'   \emph{Permitted Values:} logical expression
#'
#'   \item \emph{var:} Variable to add
#'
#'   The specified variable is compared with the variables from the other
#'   datasets. A date, a datetime, or a character variable containing ISO 8601
#'   dates can be specified.
#'
#'   \emph{Permitted Values:} variable name (as symbol)
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
#'   If the \code{order} element is specified, the mode determines if the first
#'   or last observation of each by group is selected.
#'
#'   \emph{Permitted Values:} \code{"first"}, \code{"last"} }
#'
#' @details The following steps are performed to create the output dataset:
#'
#'   \enumerate{ \item For each source dataset the observations as specified by
#'   the `filter` element are selected. If the `order` element is specified, for
#'   each by group the first or last element (depending on the `mode` element)
#'   with respect to the specified order is selected.
#'
#'   \item The new variable is set to the variable specified by the \code{var}
#'   element. If the source variable is a datetime variable, only the datepart
#'   is copied. If the source variable is a character variable, it is converted
#'   to a date. If the date is imcomplete, it is set to NA.
#'
#'   \item The selected observations of all source datasets are combined into a
#'   single dataset.
#'
#'   \item For each group (with respect to the variables specified for the
#'   `by_vars` parameter) the first or last observation (with respect to the new
#'   variable and the `mode` parameter) from the single dataset is selected and
#'   the new variable is merged to the input dataset. }
#'
#' @author Stefan Bundfuss
#'
#' @return The input dataset with the additional variable
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
#' data("ae")
#' data("lb")
#' data("adsl")
#' ae_start <- list(dataset = ae,
#'                  var = expr(AESTDTC),
#'                  filter = expr(str_length(AESTDTC) >= 10))
#' ae_end <- list(dataset = ae,
#'                vars = expr(AEENDTC),
#'                filter = expr(str_length(AEENDTC) >= 10))
#' lb_date <- list(dataset = lb,
#'                 var = expr(LBDTC),
#'                 filter = expr(str_length(LBDTC) >= 10))
#'
#' adsl_date <- list(dataset = adsl,
#'                   var = expr(TRTEDT))
#' derive_var_lstalvdt(dm,
#'                     sources = list(ae_start, ae_end, lb_date, adsl_date)) %>%
#'   select(USUBJID, LSTALVDT)
#'
derive_var_lstalvdt <- function(dataset,
                                ...) {
  sources <- list(...)
  add_data <- vector("list", length(sources))
  for (i in seq_along(sources)) {
    if (!is.null(sources[[i]]$filter)) {
      add_data[[i]] <- sources[[i]]$dataset %>%
        filter(!!(sources[[i]]$filter))
    }
    else {
      add_data[[i]] <- sources[[i]]$dataset
    }
    add_data[[i]] <- filter_extreme(add_data[[i]],
                                    order = sources[[i]]$date_var,
                                    by_vars = exprs(USUBJID),
                                    mode = "last",
                                    check_type = "none")
    date_var <- sources[[i]]$date_var
    if (is.Date(add_data[[i]][[as_string(date_var)]])) {
      add_data[[i]] <- transmute(add_data[[i]],
                                 USUBJID,
                                 !!!sources[[i]]$traceability_vars,
                                 LSTALVDT = !!date_var)
    }
    else if (is.instant(add_data[[i]][[as_string(date_var)]])) {
      add_data[[i]] <- transmute(add_data[[i]],
                                 USUBJID,
                                 !!!sources[[i]]$traceability_vars,
                                 LSTALVDT = date(!!date_var))
    }
    else {
      add_data[[i]] <- transmute(add_data[[i]],
                               USUBJID,
                               !!!sources[[i]]$traceability_vars,
                               LSTALVDT = convert_dtc_to_dt(!!date_var))
    }
  }

  all_data <- bind_rows(add_data) %>%
    filter_extreme(by_vars = exprs(USUBJID),
                   order = exprs(LSTALVDT),
                   mode = "last",
                   check_type = "none")

  # remove label to avoid warning from left_join
  attr(dataset$USUBJID, "label") <- NULL
  left_join(dataset, all_data, by = c("USUBJID"))
}

#' Create an `lstalvdt_source` object
#'
#' @param dataset A data.frame containing a source dataset.
#' @param filter An unquoted condition for filtering `dataset`.
#' @param date_var An unquoted symbol to be used for sorting `dataset`.
#' @param traceabilty_vars A named list returned by `vars()` defining the traceability variables, e.g.
#'  `vars(LALVDOM = "AE", LALVSEQ = AESEQ, LALVVAR = "AESTDTC")`.
#'
#' @author Stefan Bundfuss
#'
#' @export
#'
#' @return An object of class "lstalvdt_source".
lstalvdt_source <- function(dataset,
                            filter = NULL,
                            date_var,
                            traceability_vars = NULL) {
  if (is.null(filter)) {
    filter_value <- NULL
  }
  else {
    filter_value <- enquo(filter)
  }
  if (!(is.character(enexpr(date_var)) | is.symbol(enexpr(date_var)))) {
    abort(paste0("Invalid value of date_var parameter:\n",
                capture_output(print(enexpr(date_var))),"\n",
                "The value has to be a symbol or a character.",
                collapse = "\n"))
  }
  out <- list(
    dataset = dataset,
    filter = filter_value,
    date_var = ensym(date_var),
    traceability_vars = traceability_vars
  )
  class(out) <- c("lstalvdt_source", "list")
  validate_lstalvdt_source(out)
}

#' Validate an object is indeed a `lstalvdt_source` object
#'
#' @param obj An object to be validated.
#'
#' @author Stefan Bundfuss
#'
#' @export
#'
#' @return The original object.
validate_lstalvdt_source <- function(obj) {
  assert_that(inherits(obj, "lstalvdt_source"))
  values <- unclass(obj)
  assert_that(is.data.frame(values$dataset))
  assert_that(is.symbol(values$date_var))
  if (!is.null(values$traceability_vars)) {
    assert_that(is_varval_list(values$traceability_vars))
  }
  obj
}
