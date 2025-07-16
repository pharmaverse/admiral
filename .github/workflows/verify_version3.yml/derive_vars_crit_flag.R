#' Derive Criterion Flag Variables `CRITy`, `CRITyFL`, and `CRITyFN`
#'
#' @description
#'
#' The function derives ADaM compliant criterion flags, e.g., to facilitate
#' subgroup analyses.
#'
#' If a criterion flag can't be derived with this function, the derivation is
#' not ADaM compliant. It helps to ensure that:
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
#'
#' @permitted [dataset]
#'
#' @param crit_nr The criterion number, i.e., the `y` in `CRITy`
#'
#' @permitted [pos_int]
#'
#' @param condition Condition for flagging records
#'
#'   See description of the `values_yn` argument for details on how the
#'   `CRITyFL` variable is populated.
#'
#' @permitted [condition]
#'
#' @param description The description of the criterion
#'
#'   The `CRITy` variable is set to the specified value.
#'
#'   An expression can be specified to set the value depending on the parameter.
#'   Please note that the value must be constant within a parameter.
#'
#' @permitted an unquoted expression which evaluates to a character
#'    (in `dataset`)
#'
#' @param values_yn Should `"Y"` and `"N"` be used for `CRITyFL`?
#'
#'   If set to `TRUE`, the `CRITyFL` variable is set to `"Y"` if the condition
#'   (`condition`) evaluates to `TRUE`, it is set to `"N"` if the condition
#'   evaluate to `FALSE`, and to `NA` if it evaluates to `NA`.
#'
#'   Otherwise, the `CRITyFL` variable is set to `"Y"` if the condition
#'   (`condition`) evaluates to `TRUE`, and to `NA` otherwise.
#'
#' @permitted [boolean]
#'
#' @param create_numeric_flag Create a numeric flag?
#'
#'   If set to `TRUE`, the `CRITyFN` variable is created. It is set to `1` if
#'   `CRITyFL == "Y"`, it set to `0` if `CRITyFL == "N"`, and to `NA` otherwise.
#'
#' @permitted [boolean]
#'
#' @return The input dataset with the variables `CRITy`, `CRITyFL`, and
#'   optionally `CRITyFN` added.
#'
#' @family der_bds_findings
#' @keywords der_bds_findings
#'
#' @export
#'
#' @examplesx
#'
#' @caption Data setup
#'
#' @info The following examples use the BDS dataset below as a basis.
#'
#' @code
#' library(tibble, warn.conflicts = FALSE)
#'
#' adbds <- tribble(
#'   ~PARAMCD, ~AVAL,
#'   "AST",    42,
#'   "AST",    52,
#'   "AST",    NA_real_,
#'   "ALT",    33,
#'   "ALT",    51
#' )
#'
#' @caption Creating a simple criterion flag with values `"Y"` and `NA`
#'   (`condition`, `description`)
#'
#' @info The following call is a simple application of `derive_vars_crit_flag()`
#'   to derive a criterion flag/variable pair in a BDS dataset.
#'
#'   - The new variables are named `CRIT1`/`CRIT1FL` because the argument
#'     `crit_nr` has not been passed.
#'   - Since the argument `values_yn` has also not been passed and thus is
#'     set to its default of `FALSE`, `CRIT1FL` is set to `Y` only if
#'     `condition` evaluates to `TRUE`. For example, in both the
#'     first and third records, where `condition` is respectively `FALSE`
#'     and `NA`, we set `CRIT1FL = NA_character_`. The fourth record also
#'     exhibits this behavior. Also, as per CDISC standards, in this case
#'     `CRIT1` is populated only for records where `condition` evaluates
#'     to `TRUE`.
#'
#' @code
#' derive_vars_crit_flag(
#'   adbds,
#'   condition = AVAL > 50,
#'   description = "Absolute value > 50"
#' )
#'
#' @info The `description` argument also accepts expressions which depend
#'   on other variables in the input dataset. This can be useful to
#'   dynamically populate `CRITx`, for instance in the case below where
#'   we improve the `CRIT1` text because the same flag/variable pair is
#'   actually being used for multiple parameters.
#'
#' @code
#' derive_vars_crit_flag(
#'   adbds,
#'   condition = AVAL > 50,
#'   description = paste(PARAMCD, "> 50"),
#' )
#'
#' @caption Creating a criterion flag with values `"Y"`, `"N"` and `NA`
#'   (`values_yn`)
#'
#' @info The next call builds on the previous example by using
#'  `value_yn = TRUE` to distinguish between the cases
#'   where `condition` is `FALSE` and those where it is
#'   not evaluable at all.
#'
#'   - As compared to the previous example, for the first record `condition`
#'     evaluates to `FALSE` and so we set `CRIT1FL = "N"`, whereas for the
#'     third record, `condition` evaluates to `NA` because `AVAL` is
#'     missing and so we set `CRIT1FL` to `NA`.
#'   - Note also that because we are using the values `"Y"`, `"N"` and `NA`
#'     for the flag, as per CDISC standards `CRIT1` is now
#'     populated for all records rather than just for the `"Y"` records.
#'
#' @code
#' derive_vars_crit_flag(
#'   adbds,
#'   condition = AVAL > 50,
#'   description = paste(PARAMCD, "> 50"),
#'   values_yn = TRUE
#' )
#'
#' @info If the user wishes to set the criterion flag to `"N"` whenever
#'   the condition is not fulfilled, `condition` can be updated using
#'   an `if_else` call, where the third argument determines the behavior
#'   when the condition is not evaluable.
#'
#' @code
#' derive_vars_crit_flag(
#'   adbds,
#'   condition = if_else(AVAL > 50, TRUE, FALSE, FALSE),
#'   description = paste(PARAMCD, "> 50"),
#'   values_yn = TRUE
#' )
#'
#' @caption Specifying the criterion variable/flag number and creating
#'   a numeric flag (`crit_nr`, `create_numeric_flag`).
#'
#' @info The user can manually specify the criterion variable/flag number
#'   to use to name `CRITy`/`CRITyFL` by passing the `crit_nr` argument. This
#'   may be necessary if, for instance, other criterion flags already exist
#'   in the input dataset.
#'
#'   The user can also choose to create an additional, equivalent numeric
#'   flag `CRITyFN` by setting `create_numeric_flag` to `TRUE`.
#'
#' @code
#' derive_vars_crit_flag(
#'   adbds,
#'   condition = AVAL > 50,
#'   description = paste(PARAMCD, "> 50"),
#'   values_yn = TRUE,
#'   crit_nr = 2,
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
