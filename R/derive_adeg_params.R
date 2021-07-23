#' Adds a parameter for corrected QT using Bazett's formula
#'
#' @param dataset Input dataset
#'
#'   The variables specified by the `by_vars` parameter are expected.
#'
#' @param new_param Parameter code to add
#'
#'   For the new observations `PARAMCD` is set to the specified value.
#'
#'   Permitted Values: character value
#'
#' @param qt_code QT parameter code
#'
#'   The observations where `PARAMCD` equals the specified value are considered
#'   as the QT interval assessments.
#'
#'   Permitted Values: character value
#'
#' @param rr_code RR parameter code
#'
#'   The observations where `PARAMCD` equals the specified value are considered
#'   as the RR interval assessments.
#'
#'   Permitted Values: character value
#'
#' @param by_vars Grouping variables
#'
#'   Permitted Values: list of variables
#'
#' @inheritParams derive_derived_param
#'
#' @author Stefan Bundfuss
#'
#' @return The input dataset with the new flag variable added
#'
#' @keywords derivation adeg
#'
#' @export
#'
derive_param_qtcb <- function(dataset,
                              filter = NULL,
                              new_param = "QTCBR",
                              qt_code = "QT",
                              rr_code = "RR",
                              by_vars,
                              drop_values_from = vars(ends_with("RESU"))) {
  assert_character_scalar(new_param)
  assert_character_scalar(qt_code)
  assert_character_scalar(rr_code)
  assert_vars(by_vars)
  filter <- assert_filter_cond(enquo(filter), optional = TRUE)
  assert_data_frame(dataset, required_vars = by_vars)

  derive_derived_param(dataset,
                       filter = !!filter,
                       parameters = c(qt_code, rr_code),
                       by_vars = by_vars,
                       analysis_value = !!sym(paste0("AVAL.", qt_code)) / sqrt(!!sym(paste0("AVAL.", rr_code))/1000),
                       set_values_to = vars(PARAMCD = !!new_param),
                       drop_values_from = drop_values_from
                       )
}

#' Adds a parameter for corrected QT using Fridericia's formula
#'
#' @param dataset Input dataset
#'
#'   The variables specified by the `by_vars` parameter are expected.
#'
#' @param new_param Parameter code to add
#'
#'   For the new observations `PARAMCD` is set to the specified value.
#'
#'   Permitted Values: character value
#'
#' @param qt_code QT parameter code
#'
#'   The observations where `PARAMCD` equals the specified value are considered
#'   as the QT interval assessments.
#'
#'   Permitted Values: character value
#'
#' @param rr_code RR parameter code
#'
#'   The observations where `PARAMCD` equals the specified value are considered
#'   as the RR interval assessments.
#'
#'   Permitted Values: character value
#'
#' @param by_vars Grouping variables
#'
#'   Permitted Values: list of variables
#'
#' @inheritParams derive_derived_param
#'
#' @author Stefan Bundfuss
#'
#' @return The input dataset with the new flag variable added
#'
#' @keywords derivation adeg
#'
#' @export
#'
derive_param_qtcf <- function(dataset,
                              filter = NULL,
                              new_param = "QTCFR",
                              qt_code = "QT",
                              rr_code = "RR",
                              by_vars,
                              drop_values_from = vars(ends_with("RESU"))) {
  assert_character_scalar(new_param)
  assert_character_scalar(qt_code)
  assert_character_scalar(rr_code)
  assert_vars(by_vars)
  filter <- assert_filter_cond(enquo(filter), optional = TRUE)
  assert_data_frame(dataset, required_vars = by_vars)

  derive_derived_param(dataset,
                       filter = !!filter,
                       parameters = c(qt_code, rr_code),
                       by_vars = by_vars,
                       analysis_value = !!sym(paste0("AVAL.", qt_code)) / (!!sym(paste0("AVAL.", rr_code))/1000)^(1/3),
                       set_values_to = vars(PARAMCD = !!new_param),
                       drop_values_from = drop_values_from
  )
}
#' Adds a parameter for corrected QT using Sagie's formula
#'
#' @param dataset Input dataset
#'
#'   The variables specified by the `by_vars` parameter are expected.
#'
#' @param new_param Parameter code to add
#'
#'   For the new observations `PARAMCD` is set to the specified value.
#'
#'   Permitted Values: character value
#'
#' @param qt_code QT parameter code
#'
#'   The observations where `PARAMCD` equals the specified value are considered
#'   as the QT interval assessments.
#'
#'   Permitted Values: character value
#'
#' @param rr_code RR parameter code
#'
#'   The observations where `PARAMCD` equals the specified value are considered
#'   as the RR interval assessments.
#'
#'   Permitted Values: character value
#'
#' @param by_vars Grouping variables
#'
#'   Permitted Values: list of variables
#'
#' @inheritParams derive_derived_param
#'
#' @author Stefan Bundfuss
#'
#' @return The input dataset with the new flag variable added
#'
#' @keywords derivation adeg
#'
#' @export
#'
derive_param_qtlc <- function(dataset,
                              filter = NULL,
                              new_param = "QTLCR",
                              qt_code = "QT",
                              rr_code = "RR",
                              by_vars,
                              drop_values_from = vars(ends_with("RESU"))) {
  assert_character_scalar(new_param)
  assert_character_scalar(qt_code)
  assert_character_scalar(rr_code)
  assert_vars(by_vars)
  filter <- assert_filter_cond(enquo(filter), optional = TRUE)
  assert_data_frame(dataset, required_vars = by_vars)

  derive_derived_param(dataset,
                       filter = !!filter,
                       parameters = c(qt_code, rr_code),
                       by_vars = by_vars,
                       analysis_value = 1000*(!!sym(paste0("AVAL.", qt_code)) / 1000 + 0.154*(1 - !!sym(paste0("AVAL.", rr_code))/1000)),
                       set_values_to = vars(PARAMCD = !!new_param),
                       drop_values_from = drop_values_from
  )
}
#' Adds a parameter for derived RR
#'
#' @param dataset Input dataset
#'
#'   The variables specified by the `by_vars` parameter are expected.
#'
#' @param new_param Parameter code to add
#'
#'   For the new observations `PARAMCD` is set to the specified value.
#'
#'   Permitted Values: character value
#'
#' @param hr_code HR parameter code
#'
#'   The observations where `PARAMCD` equals the specified value are considered
#'   as the heart rate assessments.
#'
#'   Permitted Values: character value
#'
#' @param by_vars Grouping variables
#'
#'   Permitted Values: list of variables
#'
#' @inheritParams derive_derived_param
#'
#' @author Stefan Bundfuss
#'
#' @return The input dataset with the new flag variable added
#'
#' @keywords derivation adeg
#'
#' @export
#'
derive_param_rr <- function(dataset,
                            filter = NULL,
                            new_param = "RRR",
                            hr_code = "HR",
                            by_vars,
                            drop_values_from = vars(EGSEQ, starts_with("EGTEST"), contains("RES"))) {
  assert_character_scalar(new_param)
  assert_character_scalar(hr_code)
  assert_vars(by_vars)
  filter <- assert_filter_cond(enquo(filter), optional = TRUE)
  assert_data_frame(dataset, required_vars = by_vars)

  derive_derived_param(dataset,
                       filter = !!filter,
                       parameters = c(hr_code),
                       by_vars = by_vars,
                       analysis_value = 60000/!!sym(paste0("AVAL.", hr_code)),
                       set_values_to = vars(PARAMCD = !!new_param,
                                            AVALU = "msec"),
                       drop_values_from = drop_values_from
  )
}
