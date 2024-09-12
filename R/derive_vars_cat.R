#' Derive pair of variables
#' @param dataset data frame containing the variables specified in the conditions
#' @param definition of condition and values defined as `exprs()` object.
#' e.g. exprs(~condition, ~AVALCAT1, ~AVALCA1N
#'            AVAL >= 140, ">=140 cm", 1,
#'            AVAL < 140, "<140 cm", 2)
#' @param by_vars list of expressions with one element. `NULL` by default.
#' Allows for specifying a single category, e.g. `exprs(PARAMCD)`.
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
#'
#' definition2 <- exprs(
#'   ~VSTEST, ~condition, ~AVALCAT1, ~AVALCA1N,
#'   "Height", AVAL > 140, ">140 cm", 1,
#'   "Height", AVAL <= 140, "<=140 cm", 2
#' )
#'
#' derive_vars_cat(
#'   dataset = advs,
#'   definition = definition2,
#'   by_vars = exprs(VSTEST)
#' )

derive_vars_cat <- function(dataset,
                            definition,
                            by_vars = NULL) {
  # assertions
  assert_data_frame(dataset)
  assert_expr_list(definition)
  if(!is.null(by_vars)){
    assert_expr_list(by_vars)
    assert_data_frame(dataset, required_vars = c(admiraldev::extract_vars(definition) %>% unique(), extract_vars(by_vars)))

  } else{
    assert_data_frame(dataset, required_vars = admiraldev::extract_vars(definition) %>% unique())
  }


  # transform definition to tibble
  names(definition) <- NULL
  definition <- tibble::tribble(!!!definition)
  if(!is.null(by_vars)){
    assert_data_frame(definition, required_vars = c(exprs(condition), by_vars))
    # add condition
    definition <- definition %>% mutate(
      condition = extend_condition(as.character(condition), as.character(by_vars), is = !!sym(as.character(by_vars))) %>% parse_exprs()
    )
  } else assert_data_frame(definition, required_vars = exprs(condition))


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
