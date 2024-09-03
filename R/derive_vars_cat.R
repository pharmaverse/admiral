#' Derive pair of variables
#'
#'
#' @param dataset
#' @param definition
#'
#' @return
#' @export
#'
#' @examples
#'
#' advs <- tibble::tribble(
#'   ~USUBJID, ~PARAMCD, ~AVAL,
#'   "01",     "HEIGHT", 150,
#'   "02",     "HEIGHT", 135,
#'   "03",     "HEIGHT", 145
#' )
#'
#' define_tibble <- tribble(
#'   ~condition,  ~AVALCAT1,  ~AVALCA1N, ~NEWCOL,
#'   expr(AVAL > 140),  ">140 cm",          1,  "extra1",
#'   expr(AVAL <= 140), "<=140 cm",         2,  "extra2"
#' )
#'
#'
#' derive_vars_cat(dataset = advs, definition = define_tibble) %>% print()

derive_vars_cat <- function(dataset,
                            definition) {
  # assertions (waiting for feedback first)

  # extract new variable names
  new_col_names <- names(expr_list)[!names(expr_list)=="condition"]

  # (re)apply the function for each new variable name and iteratively derive the categories
  new_dataset <- reduce(new_col_names, function(.data, col_name) {
    # extract conditions
    condition <- definition[["condition"]] # could also be outside of the function
    # extract values
    values <- definition[[col_name]]

    .data %>%
      mutate(!!sym(col_name) := eval(rlang::call2(
        "case_when",
        !!!map2(condition, values, ~ expr(!!.x ~ !!.y))
      )))
  }, .init = dataset)

  return(new_dataset)
}
