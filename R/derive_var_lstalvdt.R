#' Derive Last Known Alive Date
#'
#' Add the last known alive date (`LSTALVDT`) to the dataset.
#'
#' @param dataset_name Input dataset
#'
#'   The `USUBJID` variable is expected.
#'
#' @param ... Source of known alive dates. A `lstalvdt_source()` object is
#'   expected.
#'
#' @param subject_keys Variables to uniquely identify a subject
#'
#'   A list of quosures where the expressions are symbols as returned by
#'   `vars()` is expected.
#'
#' @details The following steps are performed to create the output dataset:
#'
#'   \enumerate{ \item For each source dataset the observations as specified by
#'   the `filter` element are selected. Then for each patient the last
#'   observation (with respect to `date`) is selected.
#'
#'   \item The `LSTALVDT` variable is set to the variable specified by the
#'   \code{date} element. If the date variable is a datetime variable, only
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
#' library(dplyr, warn.conflicts = FALSE)
#' data("dm")
#' data("ae")
#' data("lb")
#' data("adsl")
#'
#' ae_start <- lstalvdt_source(
#'   dataset_name = ae,
#'   date = AESTDTC,
#'   date_imputation = "first"
#' )
#' ae_end <- lstalvdt_source(
#'   dataset_name = ae,
#'   date = AEENDTC,
#'   date_imputation = "first"
#' )
#' lb_date <- lstalvdt_source(
#'   dataset_name = lb,
#'   date = LBDTC,
#'   filter = nchar(LBDTC) >= 10
#' )
#' adsl_date <- lstalvdt_source(dataset_name = adsl, date = TRTEDT)
#'
#' dm %>%
#'   derive_var_lstalvdt(ae_start, ae_end, lb_date, adsl_date) %>%
#'   select(USUBJID, LSTALVDT)
#'
#' # derive last alive date and traceability variables
#' ae_start <- lstalvdt_source(
#'   dataset_name = ae,
#'   date = AESTDTC,
#'   date_imputation = "first",
#'   traceability_vars = vars(
#'     LALVDOM = "AE",
#'     LALVSEQ = AESEQ,
#'     LALVVAR = "AESTDTC"
#'   )
#' )
#'
#' ae_end <- lstalvdt_source(
#'   dataset_name = ae,
#'   date = AEENDTC,
#'   date_imputation = "first",
#'   traceability_vars = vars(
#'     LALVDOM = "AE",
#'     LALVSEQ = AESEQ,
#'     LALVVAR = "AEENDTC"
#'   )
#' )
#' lb_date <- lstalvdt_source(
#'   dataset_name = lb,
#'   date = LBDTC,
#'   filter = nchar(LBDTC) >= 10,
#'   traceability_vars = vars(
#'     LALVDOM = "LB",
#'     LALVSEQ = LBSEQ,
#'     LALVVAR = "LBDTC"
#'   )
#' )
#'
#' adsl_date <- lstalvdt_source(
#'   dataset_name = adsl,
#'   date = TRTEDT,
#'   traceability_vars = vars(
#'     LALVDOM = "ADSL",
#'     LALVSEQ = NA_integer_,
#'     LALVVAR = "TRTEDTM"
#'   )
#' )
#'
#' dm %>%
#'   derive_var_lstalvdt(ae_start, ae_end, lb_date, adsl_date) %>%
#'   select(USUBJID, LSTALVDT, LALVDOM, LALVSEQ, LALVVAR)
derive_var_lstalvdt <- function(dataset_name,
                                source_datasets,
                                event_conditions,
                                ...,
                                subject_keys = vars(STUDYID, USUBJID)) {
  assert_data_frame(dataset_name)
  assert_vars(subject_keys)

  sources <- list(...)
  assert_list_of(sources, "lstalvdt_source")

  source_names <- names(source_datasets)
  assert_list_element(
    list = event_conditions,
    element = "dataset_name",
    condition = dataset_name %in% source_names,
    source_names = source_names,
    message_text = paste0(
      "The dataset names must be included in the list specified for the ",
      "`source_datasets` parameter.\n",
      "Following names were provided by `source_datasets`:\n",
      enumerate(source_names, quote_fun = squote)
    )
  )

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
    date <- quo_get_expr(sources[[i]]$date)
    add_data[[i]] <- sources[[i]]$dataset_name %>%
      filter_if(sources[[i]]$filter) %>%
      filter_extreme(
        order = vars(!!date),
        by_vars = subject_keys,
        mode = "last",
        check_type = "none"
      )
    if (is.Date(add_data[[i]][[as_string(date)]])) {
      add_data[[i]] <- transmute(
        add_data[[i]],
        !!!subject_keys,
        !!!sources[[i]]$traceability_vars,
        LSTALVDT = !!date
      )
    } else if (is.instant(add_data[[i]][[as_string(date)]])) {
      add_data[[i]] <- transmute(
        add_data[[i]],
        !!!subject_keys,
        !!!sources[[i]]$traceability_vars,
        LSTALVDT = date(!!date)
      )
    } else {
      add_data[[i]] <- transmute(
        add_data[[i]],
        !!!subject_keys,
        !!!sources[[i]]$traceability_vars,
        LSTALVDT = convert_dtc_to_dt(
          !!date,
          date_imputation = sources[[i]]$date_imputation
        )
      )
    }
  }

  all_data <- add_data %>%
    bind_rows() %>%
    filter(!is.na(LSTALVDT)) %>%
    filter_extreme(
      by_vars = subject_keys,
      order = vars(LSTALVDT),
      mode = "last",
      check_type = "none"
    )

  left_join(dataset_name, all_data, by = vars2chr(subject_keys))
}

#' Create an `lstalvdt_source` object
#'
#' @param dataset_name A data.frame containing a source dataset.
#'
#' @param filter An unquoted condition for filtering `dataset`.
#'
#' @param date A variable providing a date where the patient was known to be
#'   alive. A date, a datetime, or a character variable containing ISO 8601
#'   dates can be specified. An unquoted symbol is expected.
#'
#' @param date_imputation A string defining the date imputation for `date`.
#'   See `date_imputation` parameter of `impute_dtc()` for valid values.
#'
#' @param traceability_vars A named list returned by `vars()` defining the
#'   traceability variables, e.g. `vars(LALVDOM = "AE", LALVSEQ = AESEQ, LALVVAR
#'   = "AESTDTC")`. The values must be a symbol, a character string, or `NA`.
#'
#' @param date_var Deprecated, please use `date` instead.
#'
#' @author Stefan Bundfuss
#'
#' @keywords source_specifications
#'
#' @export
#'
#' @return An object of class "lstalvdt_source".
lstalvdt_source <- function(dataset_name,
                            filter = NULL,
                            date,
                            date_imputation = NULL,
                            traceability_vars = NULL,
                            date_var = deprecated()) {

  ### BEGIN DEPRECIATION
  if (!missing(date_var)) {
    deprecate_warn("0.3.0", "lstalvdt_source(date_var = )", "lstalvdt_source(date = )")
    date <- enquo(date_var)
  }
  ### END DEPRECIATION

  if (!is.null(date_imputation)) {
    assert_that(is_valid_date_entry(date_imputation))
  }
  out <- list(
    dataset_name = assert_data_frame(dataset_name),
    filter = assert_filter_cond(enquo(filter), optional = TRUE),
    date = assert_symbol(enquo(date)),
    date_imputation = date_imputation,
    traceability_vars = assert_varval_list(traceability_vars, optional = TRUE)
  )
  class(out) <- c("lstalvdt_source", "list")
  out
}
