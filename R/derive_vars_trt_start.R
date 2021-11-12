#' Derive Datetime of Any Treatment Start Date
#'
#' Derives datetime of first exposure to treatment e.g. (`TRTSDTM`)
#'
#' @param dataset Input dataset
#'
#'   The variables specified by the `subject_keys` parameter are expected.
#'
#' @param dataset_ex `ex` dataset
#'
#'   The variables `EXSTDTC`, `EXSEQ`, and those specified by the `filter_ex`
#'   `dtc` and `subject_keys`parameters are expected.
#'
#' @param filter_ex Filter condition for the ex dataset
#'
#'   Only observations of the ex dataset which fulfill the specified condition
#'   are considered for the treatment start date.
#'
#'   Default: NULL
#'
#'   Permitted Values: logical expression
#'
#' @param subject_keys Variables to uniquely identify a subject
#'
#' A list of quosures where the expressions are symbols as returned by
#' `vars()` is expected.
#'
#' @param new_var name of treatment start date e.g. `TRTSDTM`
#'
#' @param dtc name of variable from which start date is derived e.g. `EXSTDTC`
#'
#' @inheritParams impute_dtc
#'
#' @details For each group (with respect to the variables specified for the
#'   `by_vars` parameter) the first observation (with respect to the order
#'   specified for the `order` parameter) is included in the output dataset.
#'
#' @author Stefan Bundfuss, Teckla Akinyi
#'
#' @return The input dataset with start date variable added e.g. `TRTSDTM`
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
#'  derive_vars_trt_start(dataset_ex = ex,
#'                        new_var=TRTSDTM,
#'                        dtc=EXSTDTC,
#'                        date_imputation = "first",
#'                        time_imputation="first",
#'                        min_dates = NULL,
#'                        max_dates = NULL,
#'                        order = vars(EXSTDTC,EXSEQ),
#'                        subject_keys=vars(STUDYID,USUBJID),
#'                        flag_imputation = "AUTO",
#'                        filter_ex = EXDOSE>0 ,
#'                        ord_filter="first") %>%
#'  select(USUBJID, TRTSDTM)
derive_vars_trt_start <- function(dataset,
                                  dataset_ex,
                                  filter_ex = NULL,
                                  subject_keys,
                                  new_var,
                                  dtc,
                                  date_imputation,
                                  time_imputation,
                                  min_dates = NULL,
                                  max_dates = NULL,
                                  flag_imputation = "AUTO",
                                  order,
                                  ord_filter) {

  new_var <- assert_symbol(enquo(new_var))
  dtc <- assert_symbol(enquo(dtc))
  assert_data_frame(dataset, subject_keys)
  assert_data_frame(dataset_ex, required_vars = quo_c(subject_keys, order))
  assert_character_scalar(ord_filter, values = c("first", "last"), case_sensitive = F)

  filter_ex <- assert_filter_cond(enquo(filter_ex), optional = TRUE)
  filter_ex <- enquo(filter_ex)

  assert_character_scalar(flag_imputation, values = c("auto", "both", "date", "time", "none"),
                          case_sensitive = FALSE)

  #Issue warning that new_var should be DTM
  if (str_detect(rlang::quo_name(new_var), "DTM", negate = T)) {
    msg <- sprintf(
      "This function derives the treatment datetime, new_var must end in DTM.")
    warn(msg)
  }

  add <- dataset_ex %>%
    mutate(imputed_dtc = convert_dtc_to_dtm(!!dtc,
                                            date_imputation = date_imputation,
                                            time_imputation = time_imputation,
                                            min_dates = min_dates,
                                            max_dates = max_dates))

  #filter EX data based on study requirements
  if (quo_is_null(filter_ex)) {
    add
    msg <- sprintf(
      "Input EX dataset has not been filtered. Check study requirement if ex needs to be filtered.")  # nolint
    warn(msg)
  } else {
    add <- add %>% filter(!!filter_ex)
  }

  #pick last or first date in sequence
  add <- add %>%
    filter_extreme(order = order,
                   by_vars = subject_keys,
                   mode = ord_filter) %>%
    mutate(!!new_var := imputed_dtc) %>%
    select(!!!subject_keys, !!new_var, !!dtc, imputed_dtc)

  #add imputation flags
  if (flag_imputation %in% c("both", "date") ||
      flag_imputation == "auto" && !is.null(date_imputation)) {
    # add --DTF if not there already
    dtf <- paste0("TRT", "DTF")
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
    tmf <- paste0("TRT", "TMF")
    warn_if_vars_exist(dataset, tmf)
    add <- add %>%
      mutate(!!sym(tmf) := compute_tmf(dtc = !!dtc, dtm = imputed_dtc))
  }
  add <- add %>% select(-!!dtc, -imputed_dtc)

  left_join(dataset, add, by = vars2chr(subject_keys))
}
