#' Derive Criterion Flag Variables `CRITy`, `CRITyFL`, and `CRITyFN`
#'
#' @description
#'
#' The function derives ADaM compliant criterion flags, e.g., to facilitate
#' subgroup analyses.
#'
#' If a criterion flag can't be derived with this function, the derivation is
#' not ADaM compliant. It helps to ensure that
#' - the condition of the criterion depends only on variables of the same row,
#' - the `CRITyFL` is populated with valid values, i.e, either `"Y"` and `NA` or
#'  `"Y"`, `"N"`, and `NA`,
#' - the `CRITy` variable is populated correctly, i.e.,
#'   - set to a constant value within a parameter if `CRITyFL` is populated with
#'   `"Y"`, `"N"`, and `NA` and
#'   - set to a constant value within a parameter if the criterion condition is
#'   fulfilled and to `NA` otherwise if `CRITyFL` is populated with `"Y"`,  and
#'   `NA`
#'
#' @param dataset Input dataset
#' @param crit_nr The criterion number, i.e., the `y` in `CRITy`
#'
#'   *Permitted Values*: a positive integer
#' @param condition Condition for flagging records
#'
#'   See description of the `values_yn` argument for details on how the
#'   `CRITyFL` variable is populated.
#'
#'   *Permitted Values*: an unquoted expression which evaluates to a logical (in
#'    `dataset`)
#' @param description The description of the criterion
#'
#'   The `CRITy` variable is set to the specified value.
#'
#'   An expression can be specified to set the value depending on the parameter.
#'   Please note that the value must be constant within a parameter.
#'
#'   *Permitted Values*: an unquoted expression which evaluates to a character
#'    (in `dataset`)
#' @param values_yn Should `"Y"` and `"N"` be used for `CRITyFL`?
#'
#'   If set to `TRUE`, the `CRITyFL` variable is set to `"Y"` if the condition
#'   (`condition`) evaluates to `TRUE`, it is set to `"N"` if the condition
#'   evaluate to `FALSE`, and to `NA` if it evaluates to `NA`.
#'
#'   Otherwise, the `CRITyFL` variable is set to `"Y"` if the condition
#'   (`condition`) evaluates to `TRUE`, and to `NA` otherwise.
#'
#'   *Permitted Values*: `TRUE`, `FALSE`
#' @param create_numeric_flag Create a numeric flag?
#'
#'   If set to `TRUE`, the `CRITyFN` variable is created. It is set to `1` if
#'   `CRITyFL == "Y"`, it set to `0` if `CRITyFL == "N"`, and to `NA` otherwise.
#'
#'   *Permitted Values*: `TRUE`, `FALSE`
#' @return The input dataset with the variables `CRITy`, `CRITyFL`, and
#'   optionally `CRITyFN` added.
#'
#' @family der_bds_findings
#' @keywords der_bds_findings
#'
#' @export
#'
#' @examples
#' library(tibble)
#' adbds <- tribble(
#'   ~PARAMCD, ~AVAL,
#'   "AST",    42,
#'   "AST",    52,
#'   "AST",    NA_real_,
#'   "ALT",    33,
#'   "ALT",    51
#' )
#'
#' # Create a criterion flag with values "Y" and NA
#' derive_vars_crit_flag(
#'   adbds,
#'   condition = AVAL > 50,
#'   description = "Absolute value > 50"
#' )
#'
#' # Create criterion flag with values "Y", "N", and NA and parameter dependent
#' # criterion description
#' derive_vars_crit_flag(
#'   adbds,
#'   crit_nr = 2,
#'   condition = AVAL > 50,
#'   description = paste(PARAMCD, "> 50"),
#'   values_yn = TRUE,
#'   create_numeric_flag = TRUE
#' )
derive_vars_crit_flag <- function(dataset,
                                  crit_nr = 1,
                                  condition,
                                  description,
                                  values_yn = FALSE,
                                  create_numeric_flag = FALSE) {
  assert_data_frame(dataset)
  assert_integer_scalar(crit_nr, subset = "positive")
  condition <- assert_filter_cond(enexpr(condition))
  description <- assert_expr(enexpr(description))
  assert_logical_scalar(values_yn)
  assert_logical_scalar(create_numeric_flag)

  new_critvar <- paste0("CRIT", as.character(crit_nr))
  new_critflvar <- paste0("CRIT", as.character(crit_nr), "FL")

  warn_if_vars_exist(dataset, new_critvar)
  warn_if_vars_exist(dataset, new_critflvar)

  if (values_yn) {
    crityfl_no <- "N"
  } else {
    crityfl_no <- NA_character_
  }

  tryCatch(
    dataset <- dataset %>% mutate(!!new_critflvar := if_else(!!condition, "Y", crityfl_no)),
    error = function(cnd) {
      cli_abort(
        c(
          "Evaluating {.arg condition} ({.code {as_label(condition)}}) in {.arg dataset} failed:",
          ` ` = cnd$parent$message
        ),
        call = parent.frame(n = 4)
      )
    }
  )

  tryCatch(
    {
      if (values_yn) {
        dataset <- dataset %>% mutate(!!new_critvar := !!description)
      } else {
        dataset <- dataset %>% mutate(
          !!new_critvar := if_else(!!sym(new_critflvar) == "Y", !!description, NA_character_)
        )
      }
    },
    error = function(cnd) {
      cli_abort(
        c(
          paste(
            "Evaluating {.arg description} ({.code {as_label(description)}}) in",
            "{.arg dataset} failed:"
          ),
          ` ` = cnd$parent$message
        ),
        call = parent.frame(n = 4)
      )
    }
  )

  if (create_numeric_flag) {
    new_critfnvar <- paste0("CRIT", as.character(crit_nr), "FN")
    dataset <- dataset %>% mutate(
      !!new_critfnvar := as.integer(yn_to_numeric(!!sym(new_critflvar)))
    )
  }
  dataset
}
