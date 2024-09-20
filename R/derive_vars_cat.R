#' Derive pair of variables
#' @param dataset
#' `r roxygen_param_dataset(expected_vars = c("by_vars", "definition"))`
#' @param definition List of expressions created by: `exprs()`.
#' Must be in the format of a tribble.
#' Must contain:
#'  - the column `condition` which evaluates to a logic in `dataset`.
#'  - at least one additional column with the new column name and the category value.
#'  - the column specified in `by_vars` (if `by_vars` is specified)
#'
#' e.g. if `by_vars` is not specified:
#' `exprs(~condition, ~AVALCAT1, ~AVALCA1N
#'        AVAL >= 140, ">=140 cm", 1,
#'        AVAL < 140, "<140 cm", 2)`
#'
#' e.g. if `by_vars` is specified as `exprs(VSTEST)`:
#' `exprs(~VSTEST, ~condition, ~AVALCAT1, ~AVALCA1N
#'        "Height", AVAL >= 140, ">=140 cm", 1,
#'        "Height", AVAL < 140, "<140 cm", 2)`
#'
#' @param by_vars list of expressions with one element. `NULL` by default.
#' Allows for specifying by groups, e.g. `exprs(PARAMCD)`.
#' Variable must be present in both `dataset` and `definition`.
#' The conditions in `definition` are applied only to those records that match `by_vars`.
#' The categorization variables are set to NA for records not matching any of the by groups in `definition`.
#'
#'
#' @details
#' If conditions are overlapping, the row order of `definitions` must be carefully considered.
#' The **first** match will determine the category.
#' i.e. if
#'
#' `AVAL = 155`
#'
#' and the `definition` is:
#'
#' definition <- exprs(
#'   ~VSTEST, ~condition, ~AVALCAT1, ~AVALCA1N,
#'   "Height", AVAL > 170, ">170 cm", 1,
#'   "Height", AVAL <= 170, "<=170 cm", 2,
#'   "Height", AVAL <= 160, "<=160 cm", 3,
#' )
#'
#' then `AVALCAT1` will be `"<=170 cm"`, as this is the first match for `AVAL`.
#' If you specify:
#'
#' definition <- exprs(
#'   ~VSTEST, ~condition, ~AVALCAT1, ~AVALCA1N,
#'   "Height", AVAL <= 160, "<=160 cm", 3,
#'   "Height", AVAL <= 170, "<=170 cm", 2,
#'   "Height", AVAL > 170, ">170 cm", 1,
#' )
#'
#' Then `AVAL <= 160` will lead to `AVALCAT1 == "<=160 cm"`,
#' `AVAL` inbetween `160` and `170` will lead to `AVALCAT1 == "<=170 cm"`,
#' and `AVAL <= 170` will lead to `AVALCAT1 == ">170 cm"`.
#'
#' However, we suggest to be more explicit when defining the `condition`, to avoid overlap.
#' In this case, the middle condition should be:
#' `AVAL <= 170 & AVAL > 160`
#'
#' @return data frame
#' @family der_gen
#' @keywords der_gen
#' @export
#'
#' @examples
#'
#' advs <- tibble::tribble(
#'   ~USUBJID, ~VSTEST, ~AVAL,
#'   "01-701-1015", "Height", 147.32,
#'   "01-701-1015", "Weight", 53.98,
#'   "01-701-1023", "Height", 162.56,
#'   "01-701-1023", "Weight", 78.47,
#'   "01-701-1028", "Height", 177.8,
#'   "01-701-1028", "Weight", 98.88,
#'   "01-701-1033", "Height", 175.26,
#'   "01-701-1033", "Weight", 88.45,
#'   "01-701-1034", "Height", 154.94,
#'   "01-701-1034", "Weight", 63.5,
#'   "01-701-1047", "Height", 148.59,
#'   "01-701-1047", "Weight", 66.23,
#'   "01-701-1097", "Height", 168.91,
#'   "01-701-1097", "Weight", 78.02,
#'   "01-701-1111", "Height", 158.24,
#'   "01-701-1111", "Weight", 60.33,
#'   "01-701-1115", "Height", 181.61,
#'   "01-701-1115", "Weight", 78.7,
#'   "01-701-1118", "Height", 180.34,
#'   "01-701-1118", "Weight", 71.67
#' )
#'
#' definition <- exprs(
#'   ~condition, ~AVALCAT1, ~AVALCA1N, ~NEWCOL,
#'   VSTEST == "Height" & AVAL > 140, ">140 cm", 1, "extra1",
#'   VSTEST == "Height" & AVAL <= 140, "<=140 cm", 2, "extra2"
#' )
#' derive_vars_cat(
#'   dataset = advs %>% dplyr::filter(VSTEST == "Height"),
#'   definition = definition
#' )
#' # using by_vars:
#' definition2 <- exprs(
#'   ~VSTEST, ~condition, ~AVALCAT1, ~AVALCA1N,
#'   "Height", AVAL > 140, ">140 cm", 1,
#'   "Height", AVAL <= 140, "<=140 cm", 2,
#'   "Weight", AVAL > 70, ">70 kg", 1,
#'   "Weight", AVAL <= 70, "<=70 kg", 2
#' )
#'
#' derive_vars_cat(
#'   dataset = advs,
#'   definition = definition2,
#'   by_vars = exprs(VSTEST)
#' )
#'
#' # With three conditions:
#' definition3 <- exprs(
#'   ~VSTEST, ~condition, ~AVALCAT1, ~AVALCA1N,
#'   "Height", AVAL > 170, ">170 cm", 1,
#'   "Height", AVAL <= 170 & AVAL > 160, "<=170 cm", 2,
#'   "Height", AVAL <= 160, "<=160 cm", 3,
#' )
#'
#' derive_vars_cat(
#'   dataset = advs,
#'   definition = definition3,
#'   by_vars = exprs(VSTEST)
#' )

derive_vars_cat <- function(dataset,
                            definition,
                            by_vars = NULL) {
  # assertions
  assert_data_frame(dataset)
  assert_expr_list(definition)
  if(!is.null(by_vars)){
    assert_vars(by_vars)
    assert_data_frame(dataset, required_vars = c(admiraldev::extract_vars(definition) %>% unique(), by_vars))
  } else{
    assert_data_frame(dataset, required_vars = admiraldev::extract_vars(definition) %>% unique())
  }

  # transform definition to tibble
  names(definition) <- NULL
  definition <- tibble::tribble(!!!definition)
  assert_data_frame(definition, required_vars = c(exprs(condition), by_vars))
  if(!is.null(by_vars)){
    # add condition
    definition <- definition %>% mutate(
      condition = extend_condition(as.character(condition), as.character(by_vars), is = !!sym(as.character(by_vars))) %>% parse_exprs()
    ) %>% select(-by_vars[[1]])
  }

  # extract new variable names and conditions
  new_col_names <- names(definition)[!names(definition) == "condition"]
  condition <- definition[["condition"]] # could also be outside of the function


  # (re)apply the function for each new variable name and iteratively derive the categories
  new_dataset <- reduce(new_col_names, function(.data, col_name) {
    # extract conditions
    values <- definition[[col_name]]
    # extract values

    .data %>%
      mutate(!!sym(col_name) := eval(rlang::call2(
        "case_when",
        !!!map2(condition, values, ~ expr(!!.x ~ !!.y))
      )))
  }, .init = dataset)

  return(new_dataset)
}

## helper

definition
extend_condition <- function(cond, var, is) {
  paste(cond, " & ", var, " == '", is, "'", sep = "")
}
