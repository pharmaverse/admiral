#' Adds a Parameter for Derived RR (an ECG measurement)
#'
#' @description Adds a record for derived RR based on heart rate for each by group (e.g.,
#' subject and visit) where the source parameters are available.
#'
#' **Note:** This is a wrapper function for the more generic `derive_param_computed()`.
#'
#' The analysis value of the new parameter is derived as
#' \deqn{\frac{60000}{HR}}{60000 / HR}
#'
#' @param dataset
#'   `r roxygen_param_dataset(expected_vars = c("by_vars"))`
#'   `PARAMCD`, and `AVAL` are expected as well.
#'
#'   The variable specified by `by_vars` and `PARAMCD` must be a unique key of
#'   the input dataset after restricting it by the filter condition (`filter`
#'   argument) and to the parameters specified by `hr_code`.
#'
#' @param hr_code HR parameter code
#'
#'   The observations where `PARAMCD` equals the specified value are considered
#'   as the heart rate assessments.
#'
#' @permitted character value
#'
#' @inheritParams derive_param_map
#'
#' @inheritParams derive_param_computed
#'
#' @inheritParams derive_param_qtc
#'
#'
#' @return The input dataset with the new parameter added. Note, a variable will only
#'    be populated in the new parameter rows if it is specified in `by_vars`.
#'
#' @family der_prm_bds_findings
#' @keywords der_prm_bds_findings
#'
#' @export
#'
#' @seealso [compute_rr()]
#'
#' @examples
#' library(tibble)
#'
#' adeg <- tribble(
#'   ~USUBJID, ~PARAMCD, ~PARAM, ~AVAL, ~AVALU, ~VISIT,
#'   "01-701-1015", "HR", "Heart Rate", 70.14, "beats/min", "BASELINE",
#'   "01-701-1015", "QT", "QT Duration", 370, "ms", "WEEK 2",
#'   "01-701-1015", "HR", "Heart Rate", 62.66, "beats/min", "WEEK 1",
#'   "01-701-1015", "RR", "RR Duration", 710, "ms", "WEEK 2",
#'   "01-701-1028", "HR", "Heart Rate", 85.45, "beats/min", "BASELINE",
#'   "01-701-1028", "QT", "QT Duration", 480, "ms", "WEEK 2",
#'   "01-701-1028", "QT", "QT Duration", 350, "ms", "WEEK 3",
#'   "01-701-1028", "HR", "Heart Rate", 56.54, "beats/min", "WEEK 3",
#'   "01-701-1028", "RR", "RR Duration", 842, "ms", "WEEK 2"
#' )
#'
#' derive_param_rr(
#'   adeg,
#'   by_vars = exprs(USUBJID, VISIT),
#'   set_values_to = exprs(
#'     PARAMCD = "RRR",
#'     PARAM = "RR Duration Rederived (ms)",
#'     AVALU = "ms"
#'   ),
#'   get_unit_expr = AVALU
#' )
derive_param_rr <- function(dataset,
                            by_vars,
                            set_values_to = exprs(PARAMCD = "RRR"),
                            hr_code = "HR",
                            get_unit_expr,
                            filter = NULL) {
  assert_vars(by_vars)
  assert_data_frame(
    dataset,
    required_vars = exprs(!!!by_vars, PARAMCD, AVAL)
  )
  assert_varval_list(set_values_to, required_elements = "PARAMCD", optional = TRUE)
  assert_param_does_not_exist(dataset, set_values_to$PARAMCD)
  assert_character_scalar(hr_code)
  get_unit_expr <- assert_expr(enexpr(get_unit_expr))
  filter <- assert_filter_cond(enexpr(filter), optional = TRUE)

  assert_unit(
    dataset,
    param = hr_code,
    required_unit = "beats/min",
    get_unit_expr = !!get_unit_expr
  )

  withCallingHandlers(
    derive_param_computed(
      dataset,
      filter = !!filter,
      parameters = c(hr_code),
      by_vars = by_vars,
      set_values_to = exprs(
        AVAL = compute_rr(!!sym(paste0("AVAL.", hr_code))),
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

#' Compute RR Interval From Heart Rate
#'
#' Computes RR interval from heart rate.
#'
#' @param hr Heart rate
#'
#'   A numeric vector is expected. It is expected that heart rate is measured in
#'   beats/min.
#'
#'
#' @details Usually this computation function can not be used with `%>%`.
#'
#' @return RR interval in ms:
#' \deqn{\frac{60000}{HR}}{60000 / HR}
#'
#' @family com_bds_findings
#'
#' @keywords com_bds_findings
#'
#' @export
#'
#' @seealso [derive_param_rr()]
#'
#' @examples
#' compute_rr(hr = 70.14)
compute_rr <- function(hr) {
  assert_numeric_vector(hr)
  60000 / hr
}
