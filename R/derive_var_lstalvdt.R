#' Derive Last Known Alive Date
#'
#' Add the last known alive date (`LSTALVDT`) to the dataset.
#'
#' @param dataset Input dataset
#'
#'   The `USUBJID` variable is expected.
#'
#' @param ... Source of known alive dates. A `lstalvdt_source()` object is
#'   expected.
#'
#' @param subject_keys Variables to uniquely identify a subject
#'   A list of quosures where the expressions are symbols as returned by
#'   `vars()` is expected.
#'
#' @details The following steps are performed to create the output dataset:
#'
#'   \enumerate{ \item For each source dataset the observations as specified by
#'   the `filter` element are selected. Then for each patient the last
#'   oservation (with respect to `date_var`) is selected.
#'
#'   \item The `LSTALVDT` variable is set to the variable specified by the
#'   \code{date_var} element. If the date variable is a datetime variable, only
#'   the datepart is copied. If the source variable is a character variable, it
#'   is converted to a date. If the date is incomplete, it is imputed as
#'   specified by the \code{date_imputation} element.
#'
#'   \item The variables specified by the \code{traceability_vars} element are
#'   added.
#'
#'   \item The selected observations of all source datasets are combined into a
#'   single dataset.
#'
#'   \item For each patient the last observation (with respect to the `LSTALVDT`
#'   variable) from the single dataset is selected and the new variable is
#'   merged to the input dataset. }
#'
#' @author Stefan Bundfuss
#'
#' @return The input dataset with the `LSTALVDT` variable added.
#'
#' @keywords derivation adsl
#'
#' @export
#'
#' @examples
#' library(dplyr, warn.conflict = FALSE)
#' library(stringr)
#' data("dm")
#' data("ae")
#' data("lb")
#' data("adsl")
#' ae_start <- lstalvdt_source(dataset = ae,
#'                             date_var = AESTDTC,
#'                             date_imputation = "first")
#' ae_end <- lstalvdt_source(dataset = ae,
#'                           date_var = AEENDTC,
#'                           date_imputation = "first")
#' lb_date <- lstalvdt_source(dataset = lb,
#'                            date_var = LBDTC,
#'                            filter = str_length(LBDTC) >= 10)
#' adsl_date <- lstalvdt_source(dataset = adsl,
#'                              date_var = TRTEDT)
#'
#' derive_var_lstalvdt(dm,
#'                     ae_start,
#'                     ae_end,
#'                     lb_date,
#'                     adsl_date) %>%
#'   select(USUBJID, LSTALVDT)
#'
#' # derive last alive date and traceability variables
#' ae_start <- lstalvdt_source(
#'   dataset = ae,
#'   date_var = AESTDTC,
#'   date_imputation = "first",
#'   traceability_vars = vars(
#'     LALVDOM = "AE",
#'     LALVSEQ = AESEQ,
#'     LALVVAR = "AESTDTC"
#'   )
#' )
#'
#' ae_end <- lstalvdt_source(
#'   dataset = ae,
#'   date_var = AEENDTC,
#'   date_imputation = "first",
#'   traceability_vars = vars(
#'     LALVDOM = "AE",
#'     LALVSEQ = AESEQ,
#'     LALVVAR = "AEENDTC"
#'   )
#' )
#' lb_date <- lstalvdt_source(
#'   dataset = lb,
#'   date_var = LBDTC,
#'   filter = str_length(LBDTC) >= 10,
#'   traceability_vars = vars(
#'     LALVDOM = "LB",
#'     LALVSEQ = LBSEQ,
#'     LALVVAR = "LBDTC"
#'   )
#' )
#'
#' adsl_date <- lstalvdt_source(
#'   dataset = adsl,
#'   date_var = TRTEDT,
#'   traceability_vars = vars(
#'     LALVDOM = "ADSL",
#'     LALVSEQ = NA_integer_,
#'     LALVVAR = "TRTEDTM"
#'   )
#' )
#'
#' derive_var_lstalvdt(dm,
#'                     ae_start, ae_end, lb_date, adsl_date) %>%
#'   select(USUBJID, LSTALVDT, LALVDOM, LALVSEQ, LALVVAR)
derive_var_lstalvdt <- function(dataset,
                                ...,
                                subject_keys = vars(STUDYID, USUBJID)) {
  sources <- list(...)
  add_data <- vector("list", length(sources))
  for (i in seq_along(sources)) {
    if (i > 1) {
      warn_if_inconsistent_list(
        base = sources[[i - 1]]$traceability_vars,
        compare = sources[[i]]$traceability_vars,
        list_name = "lstalvdt_source()",
        i = i
      )
    }
    if (!quo_is_null(sources[[i]]$filter)) {
      add_data[[i]] <- sources[[i]]$dataset %>%
        filter(!!(sources[[i]]$filter))
    }
    else {
      add_data[[i]] <- sources[[i]]$dataset
    }
    date_var <- quo_get_expr(sources[[i]]$date_var)
    add_data[[i]] <- filter_extreme(add_data[[i]],
                                    order = vars(!!date_var),
                                    by_vars = subject_keys,
                                    mode = "last",
                                    check_type = "none")
    if (is.Date(add_data[[i]][[as_string(date_var)]])) {
      add_data[[i]] <- transmute(add_data[[i]],
                                 !!!subject_keys,
                                 !!!sources[[i]]$traceability_vars,
                                 LSTALVDT = !!date_var)
    }
    else if (is.instant(add_data[[i]][[as_string(date_var)]])) {
      add_data[[i]] <- transmute(add_data[[i]],
                                 !!!subject_keys,
                                 !!!sources[[i]]$traceability_vars,
                                 LSTALVDT = date(!!date_var))
    }
    else {
      add_data[[i]] <- transmute(add_data[[i]],
                                 !!!subject_keys,
                                 !!!sources[[i]]$traceability_vars,
                                 LSTALVDT = convert_dtc_to_dt(
                                   impute_dtc(!!date_var,
                                              date_imputation = sources[[i]]$date_imputation)
                                 ))
    }
  }

  all_data <- bind_rows(add_data) %>%
    filter(!is.na(LSTALVDT)) %>%
    filter_extreme(
      by_vars = subject_keys,
      order = vars(LSTALVDT),
      mode = "last",
      check_type = "none"
    )

  left_join(dataset, all_data, by = vars2chr(subject_keys))
}

#' Create an `lstalvdt_source` object
#'
#' @param dataset A data.frame containing a source dataset.
#'
#' @param filter An unquoted condition for filtering `dataset`.
#'
#' @param date_var A variable providing a date where the patient was known to be
#'   alive. A date, a datetime, or a character variable containing ISO 8601
#'   dates can be specified. An unquoted symbol is expected.
#'
#' @param date_imputation A string defining the date imputation for `date_var`.
#'   See `date_imputation` parameter of `impute_dtc()` for valid values.
#'
#' @param traceability_vars A named list returned by `vars()` defining the
#'   traceability variables, e.g. `vars(LALVDOM = "AE", LALVSEQ = AESEQ, LALVVAR
#'   = "AESTDTC")`. The values must be a symbol, a character string, or `NA`.
#'
#' @author Stefan Bundfuss
#'
#' @keywords source_specifications
#'
#' @export
#'
#' @return An object of class "lstalvdt_source".
lstalvdt_source <- function(dataset,
                            filter = NULL,
                            date_var,
                            date_imputation = NULL,
                            traceability_vars = NULL) {
  out <- list(
    dataset = dataset,
    filter = enquo(filter),
    date_var = enquo(date_var),
    date_imputation = date_imputation,
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
  if (!is.data.frame(values$dataset)) {
    abort(paste0("`dataset` must be a data frame.\n",
                 "A ", typeof(values$dataset), " was supplied."))
  }
  assert_that(quo_is_null(values$filter) || is.language(quo_get_expr(values$filter)))
  if (!(quo_is_symbol(values$date_var))) {
    abort(paste0("`date_var` must be a symbol.\n",
                 "A ", typeof(quo_get_expr(values$date_var)), " was supplied."))
  }
  date_imputation <- values$date_imputation
  if (!is.null(date_imputation)) {
    assert_that(is_valid_date_entry(date_imputation))
  }
  if (!is.null(values$traceability_vars)) {
    traceability_vars <- values$traceability_vars
    assert_that(is_varval_list(traceability_vars))
  }
  obj
}
