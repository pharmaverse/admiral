#' Adds a Parameter for Corrected QT (an ECG measurement)
#'
#' @description Adds a record for corrected QT using either Bazett's, Fridericia's or Sagie's
#' formula for each by group (e.g., subject and visit) where the source parameters
#' are available.
#'
#' **Note:** This is a wrapper function for the more generic `derive_param_computed()`.
#'
#' @param dataset
#'   `r roxygen_param_dataset(expected_vars = c("by_vars", "get_unit_expr"))`
#'   `PARAMCD`, and `AVAL` are expected as well.
#'
#'   The variable specified by `by_vars` and `PARAMCD` must be a unique key of
#'   the input dataset after restricting it by the filter condition (`filter`
#'   argument) and to the parameters specified by `qt_code` and `rr_code`.
#'
#'
#' @param by_vars Grouping variables
#'
#'   Only variables specified in `by_vars` will be populated
#'   in the newly created records.
#'
#'   `r roxygen_param_by_vars()`
#'
#' @permitted [var_list]
#'
#' @param method Method used to QT correction
#'
#'   See [compute_qtc()] for details.
#'
#' @permitted `"Bazett"`, `"Fridericia"`, `"Sagie"`
#'
#' @param qt_code QT parameter code
#'
#'   The observations where `PARAMCD` equals the specified value are considered
#'   as the QT interval assessments. It is expected that QT is measured in ms or
#'   msec.
#'
#' @permitted character value
#'
#' @param rr_code RR parameter code
#'
#'   The observations where `PARAMCD` equals the specified value are considered
#'   as the RR interval assessments. It is expected that RR is measured in ms or
#'   msec.
#'
#' @permitted character value
#'
#' @param get_unit_expr An expression providing the unit of the parameter
#'
#'   The result is used to check the units of the input parameters.
#'
#' @permitted An expression which is evaluable in the input dataset
#'   and results in a character value
#'
#' @inheritParams derive_param_map
#'
#' @inheritParams derive_param_computed
#'
#' @return The input dataset with the new parameter added. Note, a variable will only
#'    be populated in the new parameter rows if it is specified in `by_vars`.
#'
#' @seealso [compute_qtc()]
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
#'   ~USUBJID,      ~PARAMCD, ~PARAM,                   ~AVAL, ~AVALU,      ~VISIT,
#'   "01-701-1015", "HR",     "Heart Rate (beats/min)", 70.14, "beats/min", "BASELINE",
#'   "01-701-1015", "QT",     "QT Duration (ms)",         370, "ms",        "WEEK 2",
#'   "01-701-1015", "HR",     "Heart Rate (beats/min)", 62.66, "beats/min", "WEEK 1",
#'   "01-701-1015", "RR",     "RR Duration (ms)",         710, "ms",        "WEEK 2",
#'   "01-701-1028", "HR",     "Heart Rate (beats/min)", 85.45, "beats/min", "BASELINE",
#'   "01-701-1028", "QT",     "QT Duration (ms)",         480, "ms",        "WEEK 2",
#'   "01-701-1028", "QT",     "QT Duration (ms)",         350, "ms",        "WEEK 3",
#'   "01-701-1028", "HR",     "Heart Rate (beats/min)", 56.54, "beats/min", "WEEK 3",
#'   "01-701-1028", "RR",     "RR Duration (ms)",         842, "ms",        "WEEK 2"
#' )
#'
#' derive_param_qtc(
#'   adeg,
#'   by_vars = exprs(USUBJID, VISIT),
#'   method = "Bazett",
#'   set_values_to = exprs(
#'     PARAMCD = "QTCBR",
#'     PARAM = "QTcB - Bazett's Correction Formula Rederived (ms)",
#'     AVALU = "ms"
#'   ),
#'   get_unit_expr = AVALU
#' )
#'
#' derive_param_qtc(
#'   adeg,
#'   by_vars = exprs(USUBJID, VISIT),
#'   method = "Fridericia",
#'   set_values_to = exprs(
#'     PARAMCD = "QTCFR",
#'     PARAM = "QTcF - Fridericia's Correction Formula Rederived (ms)",
#'     AVALU = "ms"
#'   ),
#'   get_unit_expr = extract_unit(PARAM)
#' )
#'
#' derive_param_qtc(
#'   adeg,
#'   by_vars = exprs(USUBJID, VISIT),
#'   method = "Sagie",
#'   set_values_to = exprs(
#'     PARAMCD = "QTLCR",
#'     PARAM = "QTlc - Sagie's Correction Formula Rederived (ms)",
#'     AVALU = "ms"
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
    required_vars = exprs(!!!by_vars, PARAMCD, AVAL)
  )
  assert_character_scalar(method, values = c("Bazett", "Fridericia", "Sagie"))
  assert_varval_list(set_values_to, required_elements = "PARAMCD")
  assert_param_does_not_exist(dataset, set_values_to$PARAMCD)
  assert_character_scalar(qt_code)
  assert_character_scalar(rr_code)
  get_unit_expr <- assert_expr(enexpr(get_unit_expr))
  filter <- assert_filter_cond(enexpr(filter), optional = TRUE)

  assert_unit(
    dataset,
    param = qt_code,
    required_unit = c("msec", "ms"),
    get_unit_expr = !!get_unit_expr
  )
  assert_unit(
    dataset,
    param = rr_code,
    required_unit = c("msec", "ms"),
    get_unit_expr = !!get_unit_expr
  )

  withCallingHandlers(
    derive_param_computed(
      dataset,
      filter = !!filter,
      parameters = c(qt_code, rr_code),
      by_vars = by_vars,
      set_values_to = exprs(
        AVAL = compute_qtc(
          qt = !!sym(paste0("AVAL.", qt_code)),
          rr = !!sym(paste0("AVAL.", rr_code)),
          method = !!method
        ),
        !!!set_values_to
      )
    ),
    derive_param_computed_all_na = function(cnd) {
      cli_inform(
        c(
          paste(
            "No computed records were added because for all potential computed",
            "records at least one of the contributing values was {.val {NA}}."
          ),
          "If this is not expected, please check the input data."
        ),
        class = class(cnd)
      )
      cnd_muffle(cnd)
    }
  )
}

#' Get Default Parameter Code for Corrected QT
#'
#' @param method Method used to QT correction
#'
#' @permitted `"Bazett"`, `"Fridericia"`, `"Sagie"`
#'
#' @return
#' `"QTCBR"` if `method` is `"Bazett"`, `"QTCFR"` if it's `"Fridericia"` or
#' `"QTLCR"` if it's `"Sagie"`. An error otherwise.
#'
#' @export
#'
#' @family der_prm_bds_findings
#' @keywords der_prm_bds_findings
#'
#' @seealso [derive_param_qtc()]
#'
#' @examples
#' default_qtc_paramcd("Sagie")
default_qtc_paramcd <- function(method) {
  assert_character_scalar(method, values = c("Bazett", "Fridericia", "Sagie"))
  paramcd <- c(Bazett = "QTCBR", Fridericia = "QTCFR", Sagie = "QTLCR")
  exprs(PARAMCD = !!paramcd[[method]])
}

#' Compute Corrected QT
#'
#' Computes corrected QT using Bazett's, Fridericia's or Sagie's formula.
#'
#' @param qt QT interval
#'
#'   A numeric vector is expected. It is expected that QT is measured in ms or
#'   msec.
#'
#' @param rr RR interval
#'
#'   A numeric vector is expected. It is expected that RR is measured in ms or
#'   msec.
#'
#' @param method Method used to QT correction
#'
#' @permitted `"Bazett"`, `"Fridericia"`, `"Sagie"`
#'
#' @return QT interval in ms
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
#' @seealso [derive_param_qtc()]
#'
#' @examples
#' compute_qtc(qt = 350, rr = 857, method = "Bazett")
#'
#' compute_qtc(qt = 350, rr = 857, method = "Fridericia")
#'
#' compute_qtc(qt = 350, rr = 857, method = "Sagie")
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
