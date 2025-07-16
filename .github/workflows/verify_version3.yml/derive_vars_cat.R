#' Derive Categorization Variables Like `AVALCATy` and `AVALCAyN`
#' @param dataset
#' `r roxygen_param_dataset(expected_vars = c("by_vars", "definition"))`
#' @param definition List of expressions created by `exprs()`.
#' Must be in rectangular format and specified using the same syntax as when creating
#' a `tibble` using the `tribble()` function.
#' The `definition` object will be converted to a `tibble` using `tribble()` inside this function.
#'
#' Must contain:
#'  - the column `condition` which will be converted to a logical expression and
#'  will be used on the `dataset` input.
#'  - at least one additional column with the new column name and
#'  the category value(s) used by the logical expression.
#'  - the column specified in `by_vars` (if `by_vars` is specified)
#'
#' e.g. if `by_vars` is not specified:
#'
#' ```{r}
#' #| eval: false
#' exprs(~condition,   ~AVALCAT1, ~AVALCA1N,
#'       AVAL >= 140, ">=140 cm",         1,
#'       AVAL < 140,   "<140 cm",         2)
#' ```
#'
#' e.g. if `by_vars` is specified as `exprs(VSTEST)`:
#'
#' ```{r}
#' #| eval: false
#' exprs(~VSTEST,   ~condition,  ~AVALCAT1, ~AVALCA1N,
#'       "Height", AVAL >= 140, ">=140 cm",         1,
#'       "Height",  AVAL < 140,  "<140 cm",         2)
#' ```
#'
#' @param by_vars list of expressions with one element. `NULL` by default.
#' Allows for specifying by groups, e.g. `exprs(PARAMCD)`.
#' Variable must be present in both `dataset` and `definition`.
#' The conditions in `definition` are applied only to those records that match `by_vars`.
#' The categorization variables are set to `NA` for records
#' not matching any of the by groups in `definition`.
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
#' ```{r}
#' #| eval: false
#' definition <- exprs(
#'   ~VSTEST,   ~condition,  ~AVALCAT1, ~AVALCA1N,
#'   "Height",  AVAL > 170,  ">170 cm",         1,
#'   "Height", AVAL <= 170, "<=170 cm",         2,
#'   "Height", AVAL <= 160, "<=160 cm",         3
#' )
#' ```
#' then `AVALCAT1` will be `"<=170 cm"`, as this is the first match for `AVAL`.
#' If you specify:
#'
#' ```{r}
#' #| eval: false
#' definition <- exprs(
#'   ~VSTEST,   ~condition,  ~AVALCAT1, ~AVALCA1N,
#'   "Height", AVAL <= 160, "<=160 cm",         3,
#'   "Height", AVAL <= 170, "<=170 cm",         2,
#'   "Height",  AVAL > 170,  ">170 cm",         1
#' )
#' ```
#'
#' Then `AVAL <= 160` will lead to `AVALCAT1 == "<=160 cm"`,
#' `AVAL` in-between `160` and `170` will lead to `AVALCAT1 == "<=170 cm"`,
#' and `AVAL <= 170` will lead to `AVALCAT1 == ">170 cm"`.
#'
#' However, we suggest to be more explicit when defining the `condition`, to avoid overlap.
#' In this case, the middle condition should be:
#' `AVAL <= 170 & AVAL > 160`
#'
#' @return The input dataset with the new variables defined in `definition` added
#' @family der_gen
#' @keywords der_gen
#' @export
#'
#' @examples
#' library(dplyr)
#' library(tibble)
#'
#' advs <- tibble::tribble(
#'   ~USUBJID,       ~VSTEST,  ~AVAL,
#'   "01-701-1015", "Height", 147.32,
#'   "01-701-1015", "Weight",  53.98,
#'   "01-701-1023", "Height", 162.56,
#'   "01-701-1023", "Weight",     NA,
#'   "01-701-1028", "Height",     NA,
#'   "01-701-1028", "Weight",     NA,
#'   "01-701-1033", "Height", 175.26,
#'   "01-701-1033", "Weight",  88.45
#' )
#'
#' definition <- exprs(
#'   ~condition,                        ~AVALCAT1, ~AVALCA1N,  ~NEWCOL,
#'   VSTEST == "Height" & AVAL > 160,   ">160 cm",         1, "extra1",
#'   VSTEST == "Height" & AVAL <= 160, "<=160 cm",         2, "extra2"
#' )
#' derive_vars_cat(
#'   dataset = advs,
#'   definition = definition
#' )
#'
#' # Using by_vars:
#' definition2 <- exprs(
#'   ~VSTEST,   ~condition,  ~AVALCAT1, ~AVALCA1N,
#'   "Height",  AVAL > 160,  ">160 cm",         1,
#'   "Height", AVAL <= 160, "<=160 cm",         2,
#'   "Weight",   AVAL > 70,   ">70 kg",         1,
#'   "Weight",  AVAL <= 70,  "<=70 kg",         2
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
#'   ~VSTEST,                ~condition,  ~AVALCAT1, ~AVALCA1N,
#'   "Height",               AVAL > 170,  ">170 cm",         1,
#'   "Height", AVAL <= 170 & AVAL > 160, "<=170 cm",         2,
#'   "Height",              AVAL <= 160, "<=160 cm",         3
#' )
#'
#' derive_vars_cat(
#'   dataset = advs,
#'   definition = definition3,
#'   by_vars = exprs(VSTEST)
#' )
#'
#' # Let's derive both the MCRITyML and the MCRITyMN variables
#' adlb <- tibble::tribble(
#'   ~USUBJID,     ~PARAM, ~AVAL, ~AVALU,  ~ANRHI,
#'   "01-701-1015", "ALT",   150,  "U/L",      40,
#'   "01-701-1023", "ALT",    70,  "U/L",      40,
#'   "01-701-1036", "ALT",   130,  "U/L",      40,
#'   "01-701-1048", "ALT",    30,  "U/L",      40,
#'   "01-701-1015", "AST",    50,  "U/L",      35
#' )
#'
#' definition_mcrit <- exprs(
#'   ~PARAM,                      ~condition,    ~MCRIT1ML, ~MCRIT1MN,
#'   "ALT",                    AVAL <= ANRHI,    "<=ANRHI",         1,
#'   "ALT", ANRHI < AVAL & AVAL <= 3 * ANRHI, ">1-3*ANRHI",         2,
#'   "ALT",                 3 * ANRHI < AVAL,   ">3*ANRHI",         3
#' )
#'
#' adlb %>%
#'   derive_vars_cat(
#'     definition = definition_mcrit,
#'     by_vars = exprs(PARAM)
#'   )
derive_vars_cat <- function(dataset,
                            definition,
                            by_vars = NULL) {
  assert_expr_list(definition)
  assert_vars(by_vars, optional = TRUE)
  if (length(by_vars) > 1) {
    cli_abort("{.arg by_vars} must contain just one variable, e.g. {.code exprs(PARAMCD)}")
  }

  assert_data_frame(dataset,
    required_vars = c(
      admiraldev::extract_vars(definition) %>% unique(),
      by_vars
    )
  )

  # transform definition to tibble
  names(definition) <- NULL
  definition <- tryCatch(
    {
      tibble::tribble(!!!definition)
    },
    error = function(e) {
      # Catch the error and append your own message
      cli_abort(
        c(
          paste(
            "Failed to convert {.arg definition} to {.cls tibble}.",
            "{.arg definition} should be specified similarly to how you would",
            "specify a {.cls tibble} using the {.fun tibble::tribble} function so it",
            "can be converted to {.cls tibble} using {.fun tibble::tribble}."
          ),
          e$message
        )
      )
    }
  )
  assert_data_frame(definition, required_vars = c(exprs(condition), by_vars))
  if (!is.null(by_vars)) {
    # add condition
    definition <- definition %>%
      mutate(
        condition = extend_condition(as.character(condition),
          as.character(by_vars),
          is = !!sym(as.character(by_vars))
        ) %>%
          parse_exprs()
      ) %>%
      select(-by_vars[[1]])
  }

  # extract new variable names and conditions
  new_col_names <- names(definition)[!names(definition) == "condition"]
  condition <- definition[["condition"]]

  # warn if new variables already exist
  if (any(new_col_names %in% names(dataset))) {
    cli_warn(paste(
      "Column(s) in {.arg definition} already exist in {.arg dataset}.",
      "Did you forget to specify {.arg by_vars},",
      "or are you rerunning your code?"
    ))
  }

  # (re)apply the function for each new variable name and iteratively derive the categories
  new_dataset <- reduce(new_col_names, function(.data, col_name) {
    # extract conditions
    values <- definition[[col_name]]

    .data %>%
      mutate(!!sym(col_name) := eval(rlang::call2(
        "case_when",
        !!!map2(condition, values, ~ expr(!!.x ~ !!.y))
      )))
  }, .init = dataset)

  return(new_dataset)
}

#' Extend a condition string by adding a new condition based on a variable and its value
#'
#' This internal helper function extends a condition string by appending a new condition
#' that checks if a variable equals a specific value.
#'
#' @param cond A character string representing an existing condition.
#' @param var A character string representing the name of the variable to check.
#' @param is A character string representing the value the variable should be equal to.
#'
#' @return A character string representing the extended condition.
#' @examples
#' # Extend an existing condition to include a check for 'AGE == "30"'
#' admiral:::extend_condition("SEX == 'M'", "AGE", "30")
#' @keywords internal
extend_condition <- function(cond, var, is) {
  paste(cond, " & ", var, " == '", is, "'", sep = "")
}
