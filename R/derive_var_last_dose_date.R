#' Derive Last Dose Date-Time
#'
#' Add a variable for the dose date or datetime of the last dose to the input dataset.
#'
#' @inheritParams derive_vars_last_dose
#'
#' @param new_var The new date or datetime variable added to `dataset`.
#'
#' @param output_datetime  Display `new_var` as datetime or as date only. Defaults to `TRUE`.
#'
#' @details The last dose date is derived as the maximum dose date where the
#'   `dose_date` is lower to or equal to the `analysis_date` per `by_vars` for
#'   each observation in `dataset`. When `output_datetime` is `TRUE` and time is
#'   missing, then the last dose date time is imputed to `00:00:00`. However, if
#'   date is missing, then no imputation is done.
#'
#'  If dose information is aggregated (i.e. is a dosing frequency other than `"ONCE"`
#'  over a period defined by a start and end date) the function
#'  `create_single_dose_dataset()` can be used to generate single doses from
#'  aggregate dose information and satisfy `single_dose_condition`.
#'
#' @return Input dataset with additional column `new_var`.
#'
#' @author Ben Straub
#'
#' @keywords adam derivation
#'
#' @seealso [derive_vars_last_dose()], [create_single_dose_dataset()]
#'
#' @export
#'
#' @examples
#' library(dplyr, warn.conflicts = FALSE)
#' library(admiral.test)
#' data(admiral_ae)
#' data(ex_single)
#'
#' admiral_ae %>%
#'   head(100) %>%
#'   derive_var_last_dose_date(
#'     head(ex_single, 100),
#'     filter_ex = (EXDOSE > 0 | (EXDOSE == 0 & grepl("PLACEBO", EXTRT))) &
#'       nchar(EXENDTC) >= 10,
#'     dose_date = EXENDTC,
#'     analysis_date = AESTDTC,
#'     single_dose_condition = (EXSTDTC == EXENDTC),
#'     new_var = LDOSEDTM,
#'     traceability_vars = dplyr::vars(LDOSEDOM = "EX", LDOSESEQ = EXSEQ, LDOSEVAR = "EXDOSE")
#'   ) %>%
#'   select(STUDYID, USUBJID, AESEQ, AESTDTC, LDOSEDOM, LDOSESEQ, LDOSEVAR, LDOSEDTM)
derive_var_last_dose_date <- function(dataset,
                                      dataset_ex,
                                      filter_ex = NULL,
                                      by_vars = vars(STUDYID, USUBJID),
                                      dose_id = vars(),
                                      dose_date,
                                      analysis_date,
                                      single_dose_condition = (EXDOSFRQ == "ONCE"),
                                      new_var,
                                      output_datetime = TRUE,
                                      traceability_vars = NULL) {
  filter_ex <- assert_filter_cond(enquo(filter_ex), optional = TRUE)
  by_vars <- assert_vars(by_vars)
  dose_id <- assert_vars(dose_id)
  dose_date <- assert_symbol(enquo(dose_date))
  analysis_date <- assert_symbol(enquo(analysis_date))
  single_dose_condition <- assert_filter_cond(enquo(single_dose_condition))
  new_var <- assert_symbol(enquo(new_var))
  assert_logical_scalar(output_datetime)

  res <- derive_vars_last_dose(
    dataset = dataset,
    dataset_ex = dataset_ex,
    filter_ex = !!filter_ex,
    by_vars = by_vars,
    dose_id = dose_id,
    dose_date = !!dose_date,
    analysis_date = !!analysis_date,
    single_dose_condition = !!single_dose_condition,
    new_vars = vars(!!new_var := !!dose_date),
    traceability_vars = traceability_vars
  )

  # return either date or date-time variable
  if (!output_datetime) {
    res %>% mutate(!!new_var := as.Date(!!new_var))
  } else {
    res %>% mutate(!!new_var := as.POSIXct(as.character(!!new_var), tz = "UTC"))
  }
}
