#' Derive Criterion Flag Variables `CRITy`, `CRITyFL`, and `CRITyFLN`
#'
#' Derives criterion flags based on the provided dataset, condition, and
#' description.
#'
#' @param dataset Input dataset
#' @param crit_nr The criterion number
#' @param condition The condition for deriving the criterion flag.
#' @param description The description of the criterion
#'
#'   The `CRITy` variable is set to the specified value.
#' @param values_yn Should `"Y"` and `"N"` be used for `CRITyFL`?
#'
#'   If set to `TRUE`, the `CRITyFL` variable is set to `"Y"` if the condition
#'   (`condition`) evaluates to `TRUE`, it is set to `"N"` if the condition
#'   evaluate to `FALSE`, and to `NA` if it evaluates to `NA`.
#'
#'   *Permitted Values*: `TRUE`, `FALSE`
#' @param create_numeric_flag Create a numeric flag?
#'
#'   If set to `TRUE`, the `CRITyFLN` variable is created. It is set to `1` if
#'   `CRITyFL == "Y"`, it set to `0` if `CRITyFL == "N"`, and to `NA` otherwise.
#'
#'   *Permitted Values*: `TRUE`, `FALSE`
#' @return The modified dataset with the derived criterion flag variable.
#'
#' @family der_gen
#' @keywords der_gen
#'
#' @export
#'
#' @examples
#' dataset <- data.frame(x = c(1, 2, 3), y = c(4, 5, 6))
#' derive_vars_crit_flag(dataset, crit_nr = 1, condition = x > 2, description = "Flag 1")
derive_vars_crit_flag <- function(dataset,
                                  crit_nr = 1,
                                  condition ,
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
    crityfl_no = "N"
    crityfln_no = 0L
  } else {
    crityfl_no = NA_character_
    crityfln_no = NA_integer_
  }

  tryCatch(
    dataset <- dataset %>% mutate(!!new_critflvar := if_else(!!condition, "Y", crityfl_no)),
    error = function(cnd) {
      cli_abort(c(
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
      cli_abort(c(
        "Evaluating {.arg description} ({.code {as_label(description)}}) in {.arg dataset} failed:",
        ` ` = cnd$parent$message
      ),
      call = parent.frame(n = 4))
    }
  )

  if (create_numeric_flag) {
    new_critflnvar <- paste0("CRIT", as.character(crit_nr), "FLN")
    dataset <- dataset %>% mutate(
       !!new_critflnvar := yn_to_numeric(!!sym(new_critflvar))
    )
  }
  dataset
}
