#' Derive Datetime of Any Extreme Dates
#'
#' Derives datetime of any extreme dates
#'     e.g. Treatment, Phase or Period start or end dates
#'
#' @param subject_keys Variables to uniquely identify a subject
#'
#'   The variables specified by the `subject_keys` parameter are expected.
#'
#'   Default: vars(STUDYID,USUBJID)
#'
#' @param dataset_date dataset that contains the date values e.g. `ex`
#'
#'   The variables `--ENDTC`, `--SEQ`, and those specified by the `filter_dataset_date`,
#'   `dtc` and `subject_keys` parameters are expected.
#'
#' @param filter_dataset_date Filter condition for the dataset with dates
#'
#'   Only observations of the dataset which fulfill the specified condition
#'   are considered for the treatment start date.
#'
#'   Default: NULL
#'
#'   Permitted Values: logical expression
#'
#' A list of quosures where the expressions are symbols as returned by
#' `vars()` is expected.
#'
#' @param dtc name of variable from which the extreme date is derived e.g. `EXENDTC`
#'
#' @param  new_var_prefix name of the extreme date derived
#'          e.g.Treatment start and end date `TRTST`, `TRT01ED`
#'              Phase start and end date `AP01ST`, `AP01ED`
#'
#' @inheritParams convert_dtc_to_dtm
#'
#' @inheritParams impute_dtc
#'
#' @inheritParams filter_extreme
#'
#' @details For each group (with respect to the variables specified for the
#'   `by_vars` parameter) the first observation (with respect to the order
#'   specified for the `order` parameter) is included in the output dataset.
#'
#' @author Teckla Akinyi
#'
#' @return The input dataset with extreme datetime variables (and date or time imputation flags)
#'  added e.g. Treatment start datetime and associated flags `TRTSTDTM`, `TRTSTDTF`, `TRTSTTMF`
#'              Phase start datetime and associated flags `APO1STDTM`, `APO1STDTF`, `APO1STTMF`
#'
#' @export
#'
#' @keywords adsl timing derivation
#'
#' @examples
#' library(dplyr, warn.conflicts = FALSE)
#' data("ex")
#' data("dm")
#'
#'dm %>%
#'  derive_vars_extreme_date(subject_keys=vars(STUDYID,USUBJID),
#'                           dataset_date = ex,
#'                           filter_dataset_date = EXDOSE>0,
#'                           dtc = EXSTDTC,
#'                           order = vars(EXSTDTC,EXSEQ),
#'                           new_var_prefix=TRTEDTM,
#'                           date_imputation = "first",
#'                           time_imputation = "first",
#'                           flag_imputation = "AUTO",
#'                           min_dates = NULL,
#'                           max_dates = NULL,
#'                           mode="LAST") %>%
#'   select(USUBJID, TRTEDTM)
derive_vars_extreme_dtm <- function(dataset,
                                     subject_keys =vars(STUDYID, USUBJID),
                                     dataset_date,
                                     filter_dataset_date = NULL,
                                     order,
                                     new_var_prefix,
                                     dtc,
                                     date_imputation = NULL,
                                     time_imputation = NULL,
                                     flag_imputation = "AUTO",
                                     min_dates = NULL,
                                     max_dates = NULL,
                                     mode) {

  dtc <- assert_symbol(enquo(dtc))
  assert_data_frame(dataset, subject_keys)
  assert_data_frame(dataset_date, required_vars = quo_c(subject_keys, order))
  assert_character_scalar(mode, values = c("first", "last"), case_sensitive = F)

  filter_dataset_date <- assert_filter_cond(enquo(filter_dataset_date), optional = TRUE)

  flag_imputation <- assert_character_scalar(flag_imputation,
                                             values = c("auto", "both", "date", "time", "none"),
                                             case_sensitive = FALSE)
  assert_character_scalar(new_var_prefix)
  dtm <- paste0(new_var_prefix, "DTM")

  # Issue a warning if --DTM already exists
  warn_if_vars_exist(dataset, dtm)

  # use derive_vars_dtm
  add <- dataset_date %>%
         mutate(imputed_dtc = convert_dtc_to_dtm(!!dtc,
                                            date_imputation = date_imputation,
                                            time_imputation = time_imputation,
                                            min_dates = min_dates,
                                            max_dates = max_dates)
         )

  #filter EX data based on study requirements
  if (quo_is_null(filter_dataset_date)) {
    add
    msg <- sprintf(
      "Input EX dataset has no been filtered. Check study requirement if ex needs to be filtered.")  # nolint
    inform(msg)
  } else {
    add <- add %>% filter(!!filter_dataset_date)
  }

  add <- add %>%
    filter_extreme(order = order,
                   by_vars = subject_keys,
                   mode = mode) %>%
    mutate(!!sym(dtm) := imputed_dtc) %>%
    select(!!!subject_keys, !!sym(dtm), !!dtc, imputed_dtc)

  #add imputation flags
  if (flag_imputation %in% c("both", "date") ||
      flag_imputation == "auto" && !is.null(date_imputation)) {
    # add --DTF if not there already
    dtf <- paste0(new_var_prefix, "DTF")
    dtf_exist <- dtf %in% colnames(add)
    if (!dtf_exist) {
      add <- add %>%
        mutate(!!sym(dtf) := compute_dtf(dtc = !!dtc, dt = imputed_dtc))
    } else {
      msg <- sprintf(
        "The %s variable is already present in the input dataset and will not be re-derived.",
        dtf
      )
      inform(msg)
    }
  }

  if (flag_imputation %in% c("both", "time") ||
      flag_imputation == "auto" && !is.null(time_imputation)) {
    # add --TMF variable
    tmf <- paste0(new_var_prefix, "TMF")
    warn_if_vars_exist(dataset, tmf)
    add <- add %>%
      mutate(!!sym(tmf) := compute_tmf(dtc = !!dtc, dtm = imputed_dtc))
  }
  add <- add %>% select(-!!dtc, -imputed_dtc)

  left_join(dataset, add, by = vars2chr(subject_keys))
}





#' Derive Date of Any Extreme Dates
#'
#' Derives date of any extreme dates
#'     e.g. Treatment, Phase or Period start or end dates
#'
#' @param subject_keys Variables to uniquely identify a subject
#'
#'   The variables specified by the `subject_keys` parameter are expected.
#'
#'   Default: vars(STUDYID,USUBJID)
#'
#' @param dataset_date dataset that contains the date values e.g. `ex`
#'
#'   The variables `--ENDTC`, `--SEQ`, and those specified by the `filter_dataset_date`,
#'   `dtc` and `subject_keys` parameters are expected.
#'
#' @param filter_dataset_date Filter condition for the dataset with dates
#'
#'   Only observations of the dataset which fulfill the specified condition
#'   are considered for the treatment start date.
#'
#'   Default: NULL
#'
#'   Permitted Values: logical expression
#'
#' A list of quosures where the expressions are symbols as returned by
#' `vars()` is expected.
#'
#' @param dtc name of variable from which the extreme date is derived e.g. `EXENDTC`
#'
#' @param  new_var_prefix name of the extreme date derived
#'          e.g.Treatment start and end date `TRTST`, `TRT01ED`
#'              Phase start and end date `AP01ST`, `AP01ED`
#'
#' @inheritParams convert_dtc_to_dt
#'
#' @inheritParams impute_dtc
#'
#' @inheritParams filter_extreme
#'
#' @details For each group (with respect to the variables specified for the
#'   `by_vars` parameter) the first observation (with respect to the order
#'   specified for the `order` parameter) is included in the output dataset.
#'
#' @author Teckla Akinyi
#'
#' @return The input dataset with extreme datetime variables (and date imputation flags) added
#'         e.g. Treatment start datetime and associated flags `TRTSTDT`, `TRTSTDTF`
#'              Phase start datetime and associated flags `APO1STDTM`, `APO1STDTF`
#'
#' @export
#'
#' @keywords adsl timing derivation
#'
#' @examples
#' library(dplyr, warn.conflicts = FALSE)
#' data("ex")
#' data("dm")
#'
#'dm %>%
#'  derive_vars_extreme_date(subject_keys=vars(STUDYID,USUBJID),
#'                           dataset_date = ex,
#'                           filter_dataset_date = EXDOSE>0,
#'                           dtc = EXSTDTC,
#'                           order = vars(EXSTDTC,EXSEQ),
#'                           new_var_prefix=TRTEDTM,
#'                           date_imputation = "first",
#'                           time_imputation = "first",
#'                           flag_imputation = "AUTO",
#'                           min_dates = NULL,
#'                           max_dates = NULL,
#'                           mode="LAST") %>%
#'   select(USUBJID, TRTEDTM)
derive_vars_extreme_dt <- function(dataset,
                                    subject_keys =vars(STUDYID, USUBJID),
                                    dataset_date,
                                    filter_dataset_date = NULL,
                                    order,
                                    new_var_prefix,
                                    dtc,
                                    date_imputation = NULL,
                                    flag_imputation = "AUTO",
                                    min_dates = NULL,
                                    max_dates = NULL,
                                    mode) {

  dtc <- assert_symbol(enquo(dtc))
  assert_data_frame(dataset, subject_keys)
  assert_data_frame(dataset_date, required_vars = quo_c(subject_keys, order))
  assert_character_scalar(mode, values = c("first", "last"), case_sensitive = F)

  filter_dataset_date <- assert_filter_cond(enquo(filter_dataset_date), optional = TRUE)

  flag_imputation <- assert_character_scalar(flag_imputation,
                                             values = c("auto", "both", "date", "time", "none"),
                                             case_sensitive = FALSE)
  assert_character_scalar(new_var_prefix)
  dt <- paste0(new_var_prefix, "DT")

  # Issue a warning if --DTM already exists
  warn_if_vars_exist(dataset, dt)

  # use derive_vars_dtm
  add <- dataset_date %>%
    mutate(imputed_dtc = convert_dtc_to_dt(!!dtc,
                                            date_imputation = date_imputation,
                                            min_dates = min_dates,
                                            max_dates = max_dates)
    )

  #filter EX data based on study requirements
  if (quo_is_null(filter_dataset_date)) {
    add
    msg <- sprintf(
      "Input EX dataset has no been filtered. Check study requirement if ex needs to be filtered.")  # nolint
    inform(msg)
  } else {
    add <- add %>% filter(!!filter_dataset_date)
  }

  add <- add %>%
    filter_extreme(order = order,
                   by_vars = subject_keys,
                   mode = mode) %>%
    mutate(!!sym(dt) := imputed_dtc) %>%
    select(!!!subject_keys, !!sym(dt), !!dtc, imputed_dtc)

  #add imputation flags
  if (flag_imputation %in% c("both", "date") ||
      flag_imputation == "auto" && !is.null(date_imputation)) {
    # add --DTF if not there already
    dtf <- paste0(new_var_prefix, "DTF")
    dtf_exist <- dtf %in% colnames(add)
    if (!dtf_exist) {
      add <- add %>%
        mutate(!!sym(dtf) := compute_dtf(dtc = !!dtc, dt = imputed_dtc))
    } else {
      msg <- sprintf(
        "The %s variable is already present in the input dataset and will not be re-derived.",
        dtf
      )
      inform(msg)
    }
  }

  add <- add %>% select(-!!dtc, -imputed_dtc)

  left_join(dataset, add, by = vars2chr(subject_keys))
}
