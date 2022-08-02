#' Derive Last Dose Amount
#'
#' Add a variable for dose amount from the last dose to the input dataset.
#'
#' @inheritParams derive_vars_last_dose
#' @param new_var The new variable added to `dataset`.
#' @param dose_var The EX source dose amount variable. Defaults to `EXDOSE`.
#'
#' @details The last dose amount is derived as the dose amount where the maximum `dose_date` is
#' lower to or equal to the `analysis_date` per `by_vars` for each observation in `dataset`.
#'
#' If dose information is aggregated (i.e. is a dosing frequency other than `"ONCE"`
#' over a period defined by a start and end date) the function
#' `create_single_dose_dataset()` can be used to generate single doses from
#' aggregate dose information and satisfy `single_dose_condition`.
#'
#' @return Input dataset with additional column `new_var`.
#'
#' @author Annie Yang
#'
#' @family der_gen
#' @keywords der_gen
#'
#' @export
#'
#' @seealso [derive_vars_last_dose()], [create_single_dose_dataset()]
#'
#' @examples
#' library(dplyr, warn.conflicts = FALSE)
#' library(admiral.test)
#' data(admiral_ae)
#' data(ex_single)
#'
#' admiral_ae %>%
#'   head(100) %>%
#'   derive_var_last_dose_amt(
#'     head(ex_single, 100),
#'     filter_ex = (EXDOSE > 0 | (EXDOSE == 0 & grepl("PLACEBO", EXTRT))) &
#'       nchar(EXENDTC) >= 10,
#'     dose_date = EXENDTC,
#'     analysis_date = AESTDTC,
#'     single_dose_condition = (EXSTDTC == EXENDTC),
#'     new_var = LDOSE,
#'     dose_var = EXDOSE
#'   ) %>%
#'   select(STUDYID, USUBJID, AESEQ, AESTDTC, LDOSE)
#'
#' # or with traceability variables
#' admiral_ae %>%
#'   head(100) %>%
#'   derive_var_last_dose_amt(
#'     head(ex_single, 100),
#'     filter_ex = (EXDOSE > 0 | (EXDOSE == 0 & grepl("PLACEBO", EXTRT))) &
#'       nchar(EXENDTC) >= 10,
#'     dose_date = EXENDTC,
#'     analysis_date = AESTDTC,
#'     single_dose_condition = (EXSTDTC == EXENDTC),
#'     new_var = LDOSE,
#'     dose_var = EXDOSE,
#'     traceability_vars = dplyr::vars(LDOSEDOM = "EX", LDOSESEQ = EXSEQ, LDOSEVAR = "EXDOSE")
#'   ) %>%
#'   select(STUDYID, USUBJID, AESEQ, AESTDTC, LDOSEDOM, LDOSESEQ, LDOSEVAR, LDOSE)
derive_var_last_dose_amt <- function(dataset,
                                     dataset_ex,
                                     filter_ex = NULL,
                                     by_vars = vars(STUDYID, USUBJID),
                                     dose_id = vars(),
                                     dose_date,
                                     analysis_date,
                                     single_dose_condition = (EXDOSFRQ == "ONCE"),
                                     new_var,
                                     dose_var = EXDOSE,
                                     traceability_vars = NULL) {
  filter_ex <- assert_filter_cond(enquo(filter_ex), optional = TRUE)
  by_vars <- assert_vars(by_vars)
  dose_id <- assert_vars(dose_id)
  dose_date <- assert_symbol(enquo(dose_date))
  analysis_date <- assert_symbol(enquo(analysis_date))
  single_dose_condition <- assert_filter_cond(enquo(single_dose_condition))
  new_var <- assert_symbol(enquo(new_var))
  dose_var <- assert_symbol(enquo(dose_var))

  derive_vars_last_dose(
    dataset = dataset,
    dataset_ex = dataset_ex,
    filter_ex = !!filter_ex,
    by_vars = by_vars,
    dose_id = dose_id,
    dose_date = !!dose_date,
    analysis_date = !!analysis_date,
    single_dose_condition = !!single_dose_condition,
    new_vars = vars(!!new_var := !!dose_var),
    traceability_vars = traceability_vars
  )
}
