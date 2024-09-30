#' Adds a Ratio Parameter Computed from the Analysis Value of Other Parameters
#'
#' @description Adds a record for a generic Ratio parameter using two existing parameter
#' (dividend and divisor) each by group (e.g., subject and visit) where the source parameters
#' are available.
#'
#' **Note:** This is a wrapper function for the more generic \code{derive_param_computed()}
#'
#' @param dataset Input dataset
#'
#'   The variables specified by the \code{by_vars} argument are expected to be in the dataset.
#'   \code{PARAMCD}, and \code{AVAL} are expected as well.
#'
#'   The variable specified by \code{by_vars} and \code{PARAMCD} must be a unique key of
#'   the input dataset after restricting it by the filter condition (\code{filter}
#'   parameter) and to the parameters specified by \code{dividend_code} and \code{divisor_code}.
#'
#' @param dividend_code Dividend parameter code
#'
#'   The observations where \code{PARAMCD} equals the specified value are considered
#'   as the dividend
#'
#'   *Permitted Values:* character value
#'
#' @param divisor_code Divisor parameter code
#'
#'   The observations where \code{PARAMCD} equals the specified value are considered
#'   as the divisor
#'
#'   *Permitted Values:* character value
#'
#' @param constant_dividend Is dividend parameter constant?
#'
#'   It is expected that the parameter code (PARAMCD) specified in \code{dividend_code}
#'   which is required to derive the new parameter is measured only once. For example,
#'   if Height-to-Weight Ratio should be derived and height is measured only once
#'   while Weight is measured at each visit. Height could be specified
#'   in the \code{dividend_code} argument and \code{constant_dividend} is to be set to \code{TRUE}.
#'
#'   *Permitted Values:* logical scalar
#'
#' @param constant_divisor Is divisor parameter constant?
#'
#'   It is expected that the parameter code (PARAMCD) specified in \code{dividend_code}
#'   which is required to derive the new parameter is measured only once. For example,
#'   if Waist-to-Height Ratio should be derived and height is measured only once
#'   while Waist Circumference is measured at each visit. Height could be specified
#'   in the \code{divisor_code} argument and \code{constant_divisor} is to be set to \code{TRUE}.
#'
#'   *Permitted Values:* logical scalar
#'
#' @param constant_by_vars By variables for when dividend and/or divisor is constant
#'
#'   When dividend and/or divisor is constant, the parameters (measured only once) are merged
#'   to the other parameters using the specified variables.
#'
#'   If dividend and/or divisor is constant (e.g. only measured once at screening or baseline)
#'   then use \code{constant_by_vars} to select the subject-level variable to merge on
#'   (e.g. \code{USUBJID}). This will produce a generic Ratio parameter at all visits where dividend
#'   and/or divisor is measured. Otherwise it will only be calculated at visits with both dividend
#'   and divisor parameters collected.
#'
#'   *Permitted Values*: list of variables created by \code{exprs()}
#'   e.g. \code{exprs(USUBJID, VISIT)}
#'
#' @inheritParams derive_param_bmi
#'
#' @details
#' The analysis value of the new parameter is derived as
#' \deqn{RATIO = \frac{DIVIDENT}{DIVISOR}}
#'
#'
#' @return The input dataset with the new parameter added. Note, a variable will only
#'    be populated in the new parameter rows if it is specified in \code{by_vars}.
#'
#' @family der_prm_bds_findings
#'
#' @keywords der_prm_bds_findings
#'
#' @export
#'
#' @seealso \code{\link[=compute_ratio]{compute_ratio()}}
#'
#' @examples
#'
#' # Example 1: Derive Waist-to-Hip Ratio where both source parameters are measured multiple times
#'
#' advs <- tribble(
#'   ~USUBJID, ~PARAMCD, ~PARAM, ~AVAL, ~AVALU, ~VISIT,
#'   "01-101-1001", "HIPCIR", "Hip Circumference (cm)", 125, "cm", "SCREENING",
#'   "01-101-1001", "HIPCIR", "Hip Circumference (cm)", 124, "cm", "WEEK 2",
#'   "01-101-1001", "HIPCIR", "Hip Circumference (cm)", 123, "cm", "WEEK 3",
#'   "01-101-1001", "WSTCIR", "Waist Circumference (cm)", 110, "cm", "SCREENING",
#'   "01-101-1001", "WSTCIR", "Waist Circumference (cm)", 108, "cm", "WEEK 2",
#'   "01-101-1001", "WSTCIR", "Waist Circumference (cm)", 107, "cm", "WEEK 3",
#'   "01-101-1002", "HIPCIR", "Hip Circumference (cm)", 135, "cm", "SCREENING",
#'   "01-101-1002", "HIPCIR", "Hip Circumference (cm)", 133, "cm", "WEEK 2",
#'   "01-101-1002", "HIPCIR", "Hip Circumference (cm)", 132, "cm", "WEEK 3",
#'   "01-101-1002", "WSTCIR", "Waist Circumference (cm)", 120, "cm", "SCREENING",
#'   "01-101-1002", "WSTCIR", "Waist Circumference (cm)", 118, "cm", "WEEK 2",
#'   "01-101-1002", "WSTCIR", "Waist Circumference (cm)", 117, "cm", "WEEK 3"
#' )
#'
#' derive_param_ratio(
#'   advs,
#'   by_vars = exprs(USUBJID, VISIT),
#'   dividend_code = "WSTCIR",
#'   divisor_code = "HIPCIR",
#'   set_values_to = exprs(
#'     PARAMCD = "WAISTHIP",
#'     PARAM = "Waist-to-Hip Ratio"
#'   )
#' )
#'
#' # Example 2: Derive Waist-to-Height Ratio where Height is measured only once
#'
#' advs <- tribble(
#'   ~USUBJID, ~PARAMCD, ~PARAM, ~AVAL, ~AVALU, ~VISIT,
#'   "01-101-1001", "HEIGHT", "Height (cm)", 147, "cm", "SCREENING",
#'   "01-101-1001", "WSTCIR", "Waist Circumference (cm)", 110, "cm", "SCREENING",
#'   "01-101-1001", "WSTCIR", "Waist Circumference (cm)", 108, "cm", "WEEK 2",
#'   "01-101-1001", "WSTCIR", "Waist Circumference (cm)", 107, "cm", "WEEK 3",
#'   "01-101-1002", "HEIGHT", "Height (cm)", 163, "cm", "SCREENING",
#'   "01-101-1002", "WSTCIR", "Waist Circumference (cm)", 120, "cm", "SCREENING",
#'   "01-101-1002", "WSTCIR", "Waist Circumference (cm)", 118, "cm", "WEEK 2",
#'   "01-101-1002", "WSTCIR", "Waist Circumference (cm)", 117, "cm", "WEEK 3",
#' )
#'
#' derive_param_ratio(
#'   advs,
#'   by_vars = exprs(USUBJID, VISIT),
#'   dividend_code = "WSTCIR",
#'   divisor_code = "HEIGHT",
#'   set_values_to = exprs(
#'     PARAMCD = "WAISTHGT",
#'     PARAM = "Waist-to-Height Ratio"
#'   ),
#'   constant_divisor = TRUE,
#'   constant_by_vars = exprs(USUBJID)
#' )
derive_param_ratio <- function(dataset,
                               by_vars,
                               dividend_code,
                               divisor_code,
                               set_values_to,
                               constant_dividend = FALSE,
                               constant_divisor = FALSE,
                               filter = NULL,
                               constant_by_vars = NULL) {
  assert_vars(by_vars)
  assert_data_frame(dataset, required_vars = exprs(!!!by_vars, PARAMCD, AVAL))
  assert_character_scalar(dividend_code)
  assert_character_scalar(divisor_code)
  assert_varval_list(set_values_to, required_elements = "PARAMCD")
  assert_param_does_not_exist(dataset, set_values_to$PARAMCD)
  filter <- assert_filter_cond(enexpr(filter), optional = TRUE)
  assert_logical_scalar(constant_dividend)
  assert_logical_scalar(constant_divisor)
  assert_vars(constant_by_vars, optional = TRUE)

  ratio_formula <- expr(
    compute_ratio(
      x = !!sym(paste0("AVAL.", dividend_code)),
      y = !!sym(paste0("AVAL.", divisor_code))
    )
  )

  parameters <- c(dividend_code, divisor_code)
  constant_parameters <- NULL

  if (constant_dividend) {
    constant_parameters <- c(constant_parameters, dividend_code)

    parameters <- parameters %>%
      .[!. == dividend_code]
  }

  if (constant_divisor) {
    constant_parameters <- c(constant_parameters, divisor_code)

    parameters <- parameters %>%
      .[!. == divisor_code] %>%
      ifelse(length(.) == 0, NULL, .)
  }

  derive_param_computed(
    dataset,
    filter = !!filter,
    parameters = parameters,
    by_vars = by_vars,
    set_values_to = exprs(
      AVAL = !!ratio_formula,
      !!!set_values_to
    ),
    constant_parameters = constant_parameters,
    constant_by_vars = constant_by_vars
  )
}

#' Compute ratio
#'
#' Computes ratio as dividend (x) divided by divisor (y)
#'
#' @param x Dividend
#'
#'   *Permitted Values:* numeric vector
#'
#' @param y Divisor
#'
#'   *Permitted Values:* numeric vector
#'
#'
#' @details
#' \code{x} and \code{y} must be of the same length.
#'
#' If the divisor equals 0, the result will be NA.
#'
#' Usually this computation function can not be used with \code{%>%}.
#'
#' @return The ratio computed as x / y.
#'
#' @family com_bds_findings
#'
#' @keywords com_bds_findings
#'
#' @export
#'
#' @seealso \code{\link[=derive_param_ratio]{derive_param_ratio()}}
#'
#' @examples
#'
#' # Example of Waist-to-Height Ratio
#'
#' wstcir <- c(110, 108, 107)
#' hipcir <- c(125, 124, 123)
#'
#' compute_ratio(wstcir, hipcir)
compute_ratio <- function(x, y) {
  assert_numeric_vector(x)
  assert_numeric_vector(y)

  if (length(x) != length(y)) {
    cli_abort(
      c("Argument {.arg x} must be a vector of the same length as {.arg y}",
        i = paste(
          "There are {.val {length(x)}} value{?s} in {.arg x} and",
          "{.val {length(y)}} value{?s} in {.arg y}."
        )
      )
    )
  }

  if_else(y == 0, NA_real_, x / y)
}
