#' Adds a Parameter for Corrected QT Using Bazett's Formula
#'
#' Adds a record for corrected QT using Bazett's formula for each by group
#' (e.g., subject and visit) where the source parameters are available.
#'
#' The analysis value of the new parameter is derived as
#' \deqn{\frac{QT}{\sqrt{\frac{RR}{1000}}}}{QT/\sqrt(RR/1000)}
#'
#' @param dataset Input dataset
#'
#'   The variables specified by the `by_vars` and the `unit_var` parameter,
#'   `PARAMCD`, and `AVAL` are expected.
#'
#'   The variable specified by `by_vars` and `PARAMCD` must be a unique key of
#'   the input dataset after restricting it by the filter condition (`filter`
#'   parameter) and to the parameters specified by `qt_code` and `rr_code`.
#'
#' @param qt_code QT parameter code
#'
#'   The observations where `PARAMCD` equals the specified value are considered
#'   as the QT interval assessments. It is expected that QT is measured in msec.
#'
#'   Permitted Values: character value
#'
#' @param rr_code RR parameter code
#'
#'   The observations where `PARAMCD` equals the specified value are considered
#'   as the RR interval assessments. It is expected that RR is measured in msec.
#'
#'   Permitted Values: character value
#'
#' @param by_vars Grouping variables
#'
#'   Permitted Values: list of variables
#'
#' @param unit_var Variable providing the unit of the parameter
#'
#'   The variable is used to check the units of the input parameters and it is
#'   set to `"msec"` for the new parameter.
#'
#'   Permitted Values: A variable of the input dataset
#'
#' @inheritParams derive_derived_param
#'
#' @author Stefan Bundfuss
#'
#' @return The input dataset with the new parameter added
#'
#' @keywords derivation adeg
#'
#' @export
#'
#' @examples
#' adeg <- tibble::tribble(
#'   ~USUBJID,      ~PARAMCD, ~PARAM,        ~AVAL, ~AVALU,      ~VISIT,
#'   "01-701-1015", "HR",     "Heart Rate",  70.14, "beats/min", "BASELINE",
#'   "01-701-1015", "QT",     "QT Duration", 370,   "msec",      "WEEK 2",
#'   "01-701-1015", "HR",     "Heart Rate",  62.66, "beats/min", "WEEK 1",
#'   "01-701-1015", "RR",     "RR Duration", 710,   "msec",      "WEEK 2",
#'   "01-701-1028", "HR",     "Heart Rate",  85.45, "beats/min", "BASELINE",
#'   "01-701-1028", "QT",     "QT Duration", 480,   "msec",      "WEEK 2",
#'   "01-701-1028", "QT",     "QT Duration", 350,   "msec",      "WEEK 3",
#'   "01-701-1028", "HR",     "Heart Rate",  56.54, "beats/min", "WEEK 3",
#'   "01-701-1028", "RR",     "RR Duration", 842,   "msec",      "WEEK 2",
#' )
#' derive_param_qtcb(
#'   adeg,
#'   by_vars = vars(USUBJID, VISIT),
#'   set_values_to = vars(
#'     PARAMCD = "QTCBR",
#'     PARAM = "QTcB - Bazett's Correction Formula Rederived (msec)",
#'     AVALU = "msec"
#'   )
#' )
derive_param_qtcb <- function(dataset,
                              by_vars,
                              set_values_to = vars(PARAMCD = "QTCBR"),
                              qt_code = "QT",
                              rr_code = "RR",
                              unit_var = NULL,
                              filter = NULL) {
  assert_character_scalar(qt_code)
  assert_character_scalar(rr_code)
  assert_vars(by_vars)
  unit_var <- assert_symbol(enquo(unit_var), optional = TRUE)
  filter <- assert_filter_cond(enquo(filter), optional = TRUE)
  assert_data_frame(
    dataset,
    required_vars = quo_c(by_vars, vars(PARAMCD, AVAL), unit_var)
  )
  assert_varval_list(set_values_to, required_elements = "PARAMCD", optional = TRUE)
  assert_param_does_not_exist(dataset, quo_get_expr(set_values_to$PARAMCD))

  if (!quo_is_null(unit_var)) {
    assert_unit(
      dataset,
      param = qt_code,
      unit = "msec",
      unit_var = !!unit_var
    )
    assert_unit(
      dataset,
      param = rr_code,
      unit = "msec",
      unit_var = !!unit_var
    )
    set_unit_var <- vars(!!unit_var := "msec")
  } else {
    set_unit_var <- NULL
  }

  derive_derived_param(
    dataset,
    filter = !!filter,
    parameters = c(qt_code, rr_code),
    by_vars = by_vars,
    analysis_value = compute_qtcb(
      qt = !!sym(paste0("AVAL.", qt_code)),
      rr = !!sym(paste0("AVAL.", rr_code))
    ),
    set_values_to = vars(!!!set_unit_var, !!!set_values_to)
  )
}

#' Compute Corrected QT Using Bazett's Formula
#'
#' Computes corrected QT using Bazett's formula.
#'
#' @param qt QT interval
#'
#'   A numeric vector is expected. It is expected that QT is measured in msec.
#'
#' @param rr RR interval
#'
#'   A numeric vector is expected. It is expected that RR is measured in msec.
#'
#' @author Stefan Bundfuss
#'
#' @return QT interval in msec:
#' \deqn{\frac{QT}{\sqrt{\frac{RR}{1000}}}}{QT/\sqrt(RR/1000)}
#'
#' @keywords computation adeg
#'
#' @export
#'
#' @examples
#' compute_qtcb(qt = 350, rr = 56.54)
compute_qtcb <- function(qt, rr) {
  assert_numeric_vector(qt)
  assert_numeric_vector(rr)
  qt / sqrt(rr / 1000)
}

#' Adds a parameter for corrected QT using Fridericia's formula
#'
#' Adds a record for corrected QT using Fridericia's formula for each by group
#' (e.g., subject and visit) where the source parameters are available.
#'
#' The analysis value of the new parameter is derived as
#' \deqn{\frac{QT}{\sqrt[3]{\frac{RR}{1000}}}}{QT/(RR/1000)^(1/3)}
#'
#' @inheritParams derive_param_qtcb
#'
#' @inheritParams derive_derived_param
#'
#' @author Stefan Bundfuss
#'
#' @return The input dataset with the new parameter added
#'
#' @keywords derivation adeg
#'
#' @export
#'
#' @examples
#' adeg <- tibble::tribble(
#'   ~USUBJID,      ~PARAMCD, ~PARAM,        ~AVAL, ~AVALU,      ~VISIT,
#'   "01-701-1015", "HR",     "Heart Rate",  70.14, "beats/min", "BASELINE",
#'   "01-701-1015", "QT",     "QT Duration", 370,   "msec",      "WEEK 2",
#'   "01-701-1015", "HR",     "Heart Rate",  62.66, "beats/min", "WEEK 1",
#'   "01-701-1015", "RR",     "RR Duration", 710,   "msec",      "WEEK 2",
#'   "01-701-1028", "HR",     "Heart Rate",  85.45, "beats/min", "BASELINE",
#'   "01-701-1028", "QT",     "QT Duration", 480,   "msec",      "WEEK 2",
#'   "01-701-1028", "QT",     "QT Duration", 350,   "msec",      "WEEK 3",
#'   "01-701-1028", "HR",     "Heart Rate",  56.54, "beats/min", "WEEK 3",
#'   "01-701-1028", "RR",     "RR Duration", 842,   "msec",      "WEEK 2",
#' )
#' derive_param_qtcf(
#'   adeg,
#'   by_vars = vars(USUBJID, VISIT),
#'   set_values_to = vars(
#'     PARAMCD = "QTCFR",
#'     PARAM = "QTcF - Fridericia's Correction Formula Rederived (msec)",
#'     AVALU = "msec"
#'   )
#' )
derive_param_qtcf <- function(dataset,
                              by_vars,
                              set_values_to = vars(PARAMCD = "QTCFR"),
                              qt_code = "QT",
                              rr_code = "RR",
                              unit_var = NULL,
                              filter = NULL) {
  assert_character_scalar(qt_code)
  assert_character_scalar(rr_code)
  assert_vars(by_vars)
  unit_var <- assert_symbol(enquo(unit_var), optional = TRUE)
  filter <- assert_filter_cond(enquo(filter), optional = TRUE)
  assert_data_frame(
    dataset,
    required_vars = quo_c(by_vars, vars(PARAMCD, AVAL), unit_var)
  )
  assert_varval_list(set_values_to, required_elements = "PARAMCD", optional = TRUE)
  assert_param_does_not_exist(dataset, quo_get_expr(set_values_to$PARAMCD))

  if (!quo_is_null(unit_var)) {
    assert_unit(
      dataset,
      param = qt_code,
      unit = "msec",
      unit_var = !!unit_var
    )
    assert_unit(
      dataset,
      param = rr_code,
      unit = "msec",
      unit_var = !!unit_var
    )
    set_unit_var <- vars(!!unit_var := "msec")
  } else {
    set_unit_var <- NULL
  }

  derive_derived_param(
    dataset,
    filter = !!filter,
    parameters = c(qt_code, rr_code),
    by_vars = by_vars,
    analysis_value = compute_qtcf(
      qt = !!sym(paste0("AVAL.", qt_code)),
      rr = !!sym(paste0("AVAL.", rr_code))
    ),
    set_values_to = vars(!!!set_unit_var, !!!set_values_to)
  )
}
#' Compute Corrected QT Using Fridericia's Formula
#'
#' Computes corrected QT using Fridericia's formula.
#'
#' @inheritParams compute_qtcb
#'
#' @author Stefan Bundfuss
#'
#' @return QT interval in msec:
#' \deqn{\frac{QT}{\sqrt[3]{\frac{RR}{1000}}}}{QT/(RR/1000)^(1/3)}
#'
#' @keywords computation adeg
#'
#' @export
#'
#' @examples
#' compute_qtcf(qt = 350, rr = 56.54)
compute_qtcf <- function(qt, rr) {
  assert_numeric_vector(qt)
  assert_numeric_vector(rr)
  qt / (rr / 1000)^ (1 / 3)
}

#' Adds a parameter for corrected QT using Sagie's formula
#'
#' Adds a record for corrected QT using Sagie's formula for each by group (e.g.,
#' subject and visit) where the source parameters are available.
#'
#' The analysis value of the new parameter is derived as
#' \deqn{1000\left(\frac{QT}{1000} + 0.154\left(1 - \frac{RR}{1000}\right)\right)}{
#' 1000(QT/1000 + 0.154(1 - RR/1000))}
#'
#' @inheritParams derive_derived_param
#'
#' @inheritParams derive_param_qtcb
#'
#' @author Stefan Bundfuss
#'
#' @return The input dataset with the new parameter added
#'
#' @keywords derivation adeg
#'
#' @export
#'
#' @examples
#' adeg <- tibble::tribble(
#'   ~USUBJID,      ~PARAMCD, ~PARAM,        ~AVAL, ~AVALU,      ~VISIT,
#'   "01-701-1015", "HR",     "Heart Rate",  70.14, "beats/min", "BASELINE",
#'   "01-701-1015", "QT",     "QT Duration", 370,   "msec",      "WEEK 2",
#'   "01-701-1015", "HR",     "Heart Rate",  62.66, "beats/min", "WEEK 1",
#'   "01-701-1015", "RR",     "RR Duration", 710,   "msec",      "WEEK 2",
#'   "01-701-1028", "HR",     "Heart Rate",  85.45, "beats/min", "BASELINE",
#'   "01-701-1028", "QT",     "QT Duration", 480,   "msec",      "WEEK 2",
#'   "01-701-1028", "QT",     "QT Duration", 350,   "msec",      "WEEK 3",
#'   "01-701-1028", "HR",     "Heart Rate",  56.54, "beats/min", "WEEK 3",
#'   "01-701-1028", "RR",     "RR Duration", 842,   "msec",      "WEEK 2",
#' )
#' derive_param_qtlc(
#'   adeg,
#'   by_vars = vars(USUBJID, VISIT),
#'   set_values_to = vars(
#'     PARAMCD = "QTLCR",
#'     PARAM = "QTlc - Sagie's Correction Formula Rederived (msec)",
#'     AVALU = "msec"
#'   )
#' )
derive_param_qtlc <- function(dataset,
                              by_vars,
                              set_values_to = vars(PARAMCD = "QTLCR"),
                              qt_code = "QT",
                              rr_code = "RR",
                              unit_var = NULL,
                              filter = NULL) {
  assert_character_scalar(qt_code)
  assert_character_scalar(rr_code)
  assert_vars(by_vars)
  unit_var <- assert_symbol(enquo(unit_var), optional = TRUE)
  filter <- assert_filter_cond(enquo(filter), optional = TRUE)
  assert_data_frame(
    dataset,
    required_vars = quo_c(by_vars, vars(PARAMCD, AVAL), unit_var)
  )
  assert_varval_list(set_values_to, required_elements = "PARAMCD", optional = TRUE)
  assert_param_does_not_exist(dataset, quo_get_expr(set_values_to$PARAMCD))

  if (!quo_is_null(unit_var)) {
    assert_unit(
      dataset,
      param = qt_code,
      unit = "msec",
      unit_var = !!unit_var
    )
    assert_unit(
      dataset,
      param = rr_code,
      unit = "msec",
      unit_var = !!unit_var
    )
    set_unit_var <- vars(!!unit_var := "msec")
  } else {
    set_unit_var <- NULL
  }

  derive_derived_param(
    dataset,
    filter = !!filter,
    parameters = c(qt_code, rr_code),
    by_vars = by_vars,
    analysis_value = compute_qtlc(
      qt = !!sym(paste0("AVAL.", qt_code)),
      rr = !!sym(paste0("AVAL.", rr_code))
    ),
    set_values_to = vars(!!!set_unit_var, !!!set_values_to)
  )
}

#' Compute Corrected QT Using Sagie's Formula
#'
#' Computes corrected QT using Sagie's formula.
#'
#' @inheritParams compute_qtcb
#'
#' @author Stefan Bundfuss
#'
#' @return QT interval in msec:
#' \deqn{1000\left(\frac{QT}{1000} + 0.154\left(1 - \frac{RR}{1000}\right)\right)}{
#' 1000(QT/1000 + 0.154(1 - RR/1000))}
#'
#' @keywords computation adeg
#'
#' @export
#'
#' @examples
#' compute_qtlc(qt = 350, rr = 56.54)
compute_qtlc <- function(qt, rr) {
  assert_numeric_vector(qt)
  assert_numeric_vector(rr)
  1000 * (qt / 1000 + 0.154 * (1 - rr / 1000))
}

#' Adds a parameter for derived RR
#'
#' Adds a record for derived RR based on heart rate for each by group (e.g.,
#' subject and visit) where the source parameters are available.
#'
#' The analysis value of the new parameter is derived as
#' \deqn{\frac{60000}{HR}}{60000 / HR}
#'
#' @param dataset Input dataset
#'
#'   The variables specified by the `by_vars` parameter, `PARAMCD`, and `AVAL`
#'   are expected.
#'
#'   The variable specified by `by_vars` and `PARAMCD` must be a unique key of
#'   the input dataset after restricting it by the filter condition (`filter`
#'   parameter) and to the parameters specified by `hr_code`.
#'
#' @param hr_code HR parameter code
#'
#'   The observations where `PARAMCD` equals the specified value are considered
#'   as the heart rate assessments.
#'
#'   Permitted Values: character value
#'
#' @inheritParams derive_derived_param
#'
#' @inheritParams derive_param_qtcb
#'
#' @author Stefan Bundfuss
#'
#' @return The input dataset with the new parameter added
#'
#' @keywords derivation adeg
#'
#' @export
#'
#' @examples
#' adeg <- tibble::tribble(
#'   ~USUBJID, ~PARAMCD, ~PARAM, ~AVAL, ~AVALU, ~VISIT,
#'   "01-701-1015", "HR", "Heart Rate", 70.14, "beats/min", "BASELINE",
#'   "01-701-1015", "QT", "QT Duration", 370, "msec", "WEEK 2",
#'   "01-701-1015", "HR", "Heart Rate", 62.66, "beats/min", "WEEK 1",
#'   "01-701-1015", "RR", "RR Duration", 710, "msec", "WEEK 2",
#'   "01-701-1028", "HR", "Heart Rate", 85.45, "beats/min", "BASELINE",
#'   "01-701-1028", "QT", "QT Duration", 480, "msec", "WEEK 2",
#'   "01-701-1028", "QT", "QT Duration", 350, "msec", "WEEK 3",
#'   "01-701-1028", "HR", "Heart Rate", 56.54, "beats/min", "WEEK 3",
#'   "01-701-1028", "RR", "RR Duration", 842, "msec", "WEEK 2",
#' )
#' derive_param_rr(
#'   adeg,
#'   by_vars = vars(USUBJID, VISIT),
#'   set_values_to = vars(
#'     PARAMCD = "RRR",
#'     PARAM = "RR Duration Rederived (msec)",
#'     AVALU = "msec"
#'   )
#' )
derive_param_rr <- function(dataset,
                            by_vars,
                            set_values_to = vars(PARAMCD = "RRR"),
                            hr_code = "HR",
                            unit_var = NULL,
                            filter = NULL) {
  assert_character_scalar(hr_code)
  assert_vars(by_vars)
  unit_var <- assert_symbol(enquo(unit_var), optional = TRUE)
  filter <- assert_filter_cond(enquo(filter), optional = TRUE)
  assert_data_frame(
    dataset,
    required_vars = quo_c(by_vars, vars(PARAMCD, AVAL), unit_var)
  )
  assert_varval_list(set_values_to, required_elements = "PARAMCD", optional = TRUE)
  assert_param_does_not_exist(dataset, quo_get_expr(set_values_to$PARAMCD))

  if (!quo_is_null(unit_var)) {
    assert_unit(
      dataset,
      param = hr_code,
      unit = "beats/min",
      unit_var = !!unit_var
    )
    set_unit_var <- vars(!!unit_var := "msec")
  } else {
    set_unit_var <- NULL
  }

  derive_derived_param(
    dataset,
    filter = !!filter,
    parameters = c(hr_code),
    by_vars = by_vars,
    analysis_value = compute_rr(!!sym(paste0("AVAL.", hr_code))),
    set_values_to = vars(!!!set_unit_var, !!!set_values_to)
  )
}

#' Compute RR Interval From Heart Rate
#'
#' Computes RR interval from heart rate.
#'
#' @param hr Heart rate
#'
#'   A numeric vector is expected. It is expected that heart rate is measured in
#'   beats/min.
#'
#' @author Stefan Bundfuss
#'
#' @return RR interval in msec:
#' \deqn{\frac{60000}{HR}}{60000 / HR}
#'
#' @keywords computation adeg
#'
#' @export
#'
#' @examples
#' compute_rr(hr = 70.14)
compute_rr <- function(hr) {
  assert_numeric_vector(hr)
  60000 / hr
}
