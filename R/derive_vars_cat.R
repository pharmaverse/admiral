#' Derive pair of variables
#' @param dataset data frame containing the variables specified in the conditions
#' @param definition of condition and values defined as `exprs()` object.
#' e.g. exprs(~condition, ~AVALCAT1, ~AVALCA1N
#'            AVAL >= 140, ">=140 cm", 1,
#'            AVAL < 140, "<140 cm", 2)
#'
#' @return data frame
#' @export
#'
#' @examples
#'
#' advs <- pharmaverseadam::advs
#'
#' definition <- exprs(
#'   ~condition,  ~AVALCAT1,  ~AVALCA1N, ~NEWCOL,
#'   VSTESTCD == "Height" & AVAL > 140,  ">140 cm",          1,  "extra1",
#'   VSTESTCD == "Height" & AVAL <= 140, "<=140 cm",         2,  "extra2"
#' )
#'
#' derive_vars_cat(dataset = advs %>% filter(VSTESTCD == "Height"),
#'                 definition = definition) %>%
#'   select(USUBJID, VSTESTCD, AVAL, AVALCA1N, AVALCAT1)

derive_vars_cat <- function(dataset,
                            definition) {
  # assertions
  assert_data_frame(dataset)
  assert_expr_list(definition)
  assert_data_frame(dataset, required_vars = admiraldev::extract_vars(definition) %>% unique())

  # transform definition to tibble
  names(definition) <- NULL
  definition <- tibble::tribble(!!!definition)
  assert_data_frame(definition, required_vars = exprs(condition))


  # extract new variable names and conditions
  new_col_names <- names(definition) [!names(definition) == "condition"]
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
