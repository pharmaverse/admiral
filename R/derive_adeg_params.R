#' Adds a Parameter for Corrected QT (an ECG measurement)
#'
#' Adds a record for corrected QT using either Bazett's, Fridericia's or Sagie's
#' formula for each by group (e.g., subject and visit) where the source parameters
#' are available.
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
#'
#' @param by_vars Grouping variables
#'
#'   Only variables specified in `by_vars` will be populated
#'   in the newly created records.
#'
#'   Permitted Values: list of variables
#'
#' @param method Method used to QT correction
#'
#'   Permitted Values: `"Bazett"`, `"Fridericia"`, `"Sagie"`
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
#' @param get_unit_expr An expression providing the unit of the parameter
#'
#'   The result is used to check the units of the input parameters.
#'
#'   Permitted Values: A variable of the input dataset or a function call
#'
#' @inheritParams derive_param_computed
#'
#' @seealso [compute_qtc()]
#'
#' @author Stefan Bundfuss
#'
#' @return The input dataset with the new parameter added. Note, a variable will only
#'    be populated in the new parameter rows if it is specified in `by_vars`.
#'
#' @family der_prm_bds_findings
#'
#' @keywords der_prm_bds_findings
#'
#' @export
#'
#' @examples
#' library(tibble)
#'
#' adeg <- tribble(
#'   ~USUBJID, ~PARAMCD, ~PARAM, ~AVAL, ~AVALU, ~VISIT,
#'   "01-701-1015", "HR", "Heart Rate (beats/min)", 70.14, "beats/min", "BASELINE",
#'   "01-701-1015", "QT", "QT Duration (msec)", 370, "msec", "WEEK 2",
#'   "01-701-1015", "HR", "Heart Rate (beats/min)", 62.66, "beats/min", "WEEK 1",
#'   "01-701-1015", "RR", "RR Duration (msec)", 710, "msec", "WEEK 2",
#'   "01-701-1028", "HR", "Heart Rate (beats/min)", 85.45, "beats/min", "BASELINE",
#'   "01-701-1028", "QT", "QT Duration (msec)", 480, "msec", "WEEK 2",
#'   "01-701-1028", "QT", "QT Duration (msec)", 350, "msec", "WEEK 3",
#'   "01-701-1028", "HR", "Heart Rate (beats/min)", 56.54, "beats/min", "WEEK 3",
#'   "01-701-1028", "RR", "RR Duration (msec)", 842, "msec", "WEEK 2",
#' )
#'
#' derive_param_qtc(
#'   adeg,
#'   by_vars = vars(USUBJID, VISIT),
#'   method = "Bazett",
#'   set_values_to = vars(
#'     PARAMCD = "QTCBR",
#'     PARAM = "QTcB - Bazett's Correction Formula Rederived (msec)",
#'     AVALU = "msec"
#'   ),
#'   get_unit_expr = AVALU
#' )
#'
#' derive_param_qtc(
#'   adeg,
#'   by_vars = vars(USUBJID, VISIT),
#'   method = "Fridericia",
#'   set_values_to = vars(
#'     PARAMCD = "QTCFR",
#'     PARAM = "QTcF - Fridericia's Correction Formula Rederived (msec)",
#'     AVALU = "msec"
#'   ),
#'   get_unit_expr = extract_unit(PARAM)
#' )
#'
#' derive_param_qtc(
#'   adeg,
#'   by_vars = vars(USUBJID, VISIT),
#'   method = "Sagie",
#'   set_values_to = vars(
#'     PARAMCD = "QTLCR",
#'     PARAM = "QTlc - Sagie's Correction Formula Rederived (msec)",
#'     AVALU = "msec"
#'   ),
#'   get_unit_expr = extract_unit(PARAM)
#' )
derive_param_qtc <- function(dataset,
                             by_vars,
                             method,
                             set_values_to = default_qtc_paramcd(method),
                             qt_code = "QT",
                             rr_code = "RR",
                             get_unit_expr,
                             filter = NULL) {
  assert_vars(by_vars)
  assert_data_frame(
    dataset,
    required_vars = vars(!!!by_vars, PARAMCD, AVAL)
  )
  assert_character_scalar(method, values = c("Bazett", "Fridericia", "Sagie"))
  assert_varval_list(set_values_to, required_elements = "PARAMCD")
  assert_param_does_not_exist(dataset, quo_get_expr(set_values_to$PARAMCD))
  assert_character_scalar(qt_code)
  assert_character_scalar(rr_code)
  get_unit_expr <- assert_expr(enquo(get_unit_expr))
  filter <- assert_filter_cond(enquo(filter), optional = TRUE)

  assert_unit(
    dataset,
    param = qt_code,
    required_unit = "msec",
    get_unit_expr = !!get_unit_expr
  )
  assert_unit(
    dataset,
    param = rr_code,
    required_unit = "msec",
    get_unit_expr = !!get_unit_expr
  )

  derive_param_computed(
    dataset,
    filter = !!filter,
    parameters = c(qt_code, rr_code),
    by_vars = by_vars,
    analysis_value = compute_qtc(
      qt = !!sym(paste0("AVAL.", qt_code)),
      rr = !!sym(paste0("AVAL.", rr_code)),
      method = method
    ),
    set_values_to = set_values_to
  )
}

#' Get Default Parameter Code for Corrected QT
#'
#' @param method Method used to QT correction
#'
#'   Permitted Values: `"Bazett"`, `"Fridericia"`, `"Sagie"`
#'
#' @return
#' `"QTCBR"` if `method` is `"Bazett"`, `"QTCFR"` if it's `"Fridericia"` or
#' `"QTLCR"` if it's `"Sagie"`. An error otherwise.
#'
#' @author Thomas Neitmann
#'
#' @export
#'
#' @family der_prm_bds_findings
#' @keywords der_prm_bds_findings
#'
#' @examples
#' default_qtc_paramcd("Sagie")
default_qtc_paramcd <- function(method) {
  assert_character_scalar(method, values = c("Bazett", "Fridericia", "Sagie"))
  paramcd <- c(Bazett = "QTCBR", Fridericia = "QTCFR", Sagie = "QTLCR")
  vars(PARAMCD = !!paramcd[[method]])
}

#' Compute Corrected QT
#'
#' Computes corrected QT using Bazett's, Fridericia's or Sagie's formula.
#'
#' @param qt QT interval
#'
#'   A numeric vector is expected. It is expected that QT is measured in msec.
#'
#' @param rr RR interval
#'
#'   A numeric vector is expected. It is expected that RR is measured in msec.
#'
#' @inheritParams derive_param_qtc
#'
#' @author Stefan Bundfuss
#'
#' @return QT interval in msec
#'
#' @details
#' Depending on the chosen `method` one of the following formulae is used.
#'
#' *Bazett*: \deqn{\frac{QT}{\sqrt{\frac{RR}{1000}}}}{QT/\sqrt(RR/1000)}
#'
#' *Fridericia*: \deqn{\frac{QT}{\sqrt[3]{\frac{RR}{1000}}}}{QT/(RR/1000)^(1/3)}
#'
#' *Sagie*: \deqn{1000\left(\frac{QT}{1000} + 0.154\left(1 - \frac{RR}{1000}\right)\right)}{
#' 1000(QT/1000 + 0.154(1 - RR/1000))}
#'
#' Usually this computation function can not be used with `%>%`.
#'
#' @family com_bds_findings
#'
#' @keywords com_bds_findings
#'
#' @export
#'
#' @examples
#' compute_qtc(qt = 350, rr = 56.54, method = "Bazett")
#'
#' compute_qtc(qt = 350, rr = 56.54, method = "Fridericia")
#'
#' compute_qtc(qt = 350, rr = 56.54, method = "Sagie")
compute_qtc <- function(qt, rr, method) {
  assert_numeric_vector(qt)
  assert_numeric_vector(rr)
  assert_character_scalar(method, values = c("Bazett", "Fridericia", "Sagie"))

  formulae <- alist(
    Bazett = qt / sqrt(rr / 1000),
    Fridericia = qt / (rr / 1000)^(1 / 3), # nolint
    Sagie = 1000 * (qt / 1000 + 0.154 * (1 - rr / 1000))
  )
  eval(formulae[[method]])
}

#' Adds a Parameter for Derived RR (an ECG measurement)
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
#' @inheritParams derive_param_computed
#'
#' @inheritParams derive_param_qtc
#'
#' @author Stefan Bundfuss
#'
#' @return The input dataset with the new parameter added. Note, a variable will only
#'    be populated in the new parameter rows if it is specified in `by_vars`.
#'
#' @family der_prm_bds_findings
#' @keywords der_prm_bds_findings
#'
#' @export
#'
#' @examples
#' library(tibble)
#'
#' adeg <- tribble(
#'   ~USUBJID, ~PARAMCD, ~PARAM, ~AVAL, ~AVALU, ~VISIT,
#'   "01-701-1015", "HR", "Heart Rate", 70.14, "beats/min", "BASELINE",
#'   "01-701-1015", "QT", "QT Duration", 370, "msec", "WEEK 2",
#'   "01-701-1015", "HR", "Heart Rate", 62.66, "beats/min", "WEEK 1",
#'   "01-701-1015", "RR", "RR Duration", 710, "msec", "WEEK 2",
#'   "01-701-1028", "HR", "Heart Rate", 85.45, "beats/min", "BASELINE",
#'   "01-701-1028", "QT", "QT Duration", 480, "msec", "WEEK 2",
#'   "01-701-1028", "QT", "QT Duration", 350, "msec", "WEEK 3",
#'   "01-701-1028", "HR", "Heart Rate", 56.54, "beats/min", "WEEK 3",
#'   "01-701-1028", "RR", "RR Duration", 842, "msec", "WEEK 2"
#' )
#'
#' derive_param_rr(
#'   adeg,
#'   by_vars = vars(USUBJID, VISIT),
#'   set_values_to = vars(
#'     PARAMCD = "RRR",
#'     PARAM = "RR Duration Rederived (msec)",
#'     AVALU = "msec"
#'   ),
#'   get_unit_expr = AVALU
#' )
derive_param_rr <- function(dataset,
                            by_vars,
                            set_values_to = vars(PARAMCD = "RRR"),
                            hr_code = "HR",
                            get_unit_expr,
                            filter = NULL) {
  assert_vars(by_vars)
  assert_data_frame(
    dataset,
    required_vars = vars(!!!by_vars, PARAMCD, AVAL)
  )
  assert_varval_list(set_values_to, required_elements = "PARAMCD", optional = TRUE)
  assert_param_does_not_exist(dataset, quo_get_expr(set_values_to$PARAMCD))
  assert_character_scalar(hr_code)
  get_unit_expr <- assert_expr(enquo(get_unit_expr))
  filter <- assert_filter_cond(enquo(filter), optional = TRUE)

  assert_unit(
    dataset,
    param = hr_code,
    required_unit = "beats/min",
    get_unit_expr = !!get_unit_expr
  )

  derive_param_computed(
    dataset,
    filter = !!filter,
    parameters = c(hr_code),
    by_vars = by_vars,
    analysis_value = compute_rr(!!sym(paste0("AVAL.", hr_code))),
    set_values_to = set_values_to
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
#' @details Usually this computation function can not be used with `%>%`.
#'
#' @return RR interval in msec:
#' \deqn{\frac{60000}{HR}}{60000 / HR}
#'
#' @family com_bds_findings
#'
#' @keywords com_bds_findings
#'
#' @export
#'
#' @examples
#' compute_rr(hr = 70.14)
compute_rr <- function(hr) {
  assert_numeric_vector(hr)
  60000 / hr
}
