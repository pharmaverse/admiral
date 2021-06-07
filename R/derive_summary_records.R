#' Add new record within by groups using aggregation functions
#'
#' @description
#' It is not uncommon to have an analysis need whereby one needs to derive an
#' analysis value (`AVAL`) from multiple records. The ADaM basic dataset
#' structure variable `DTYPE` is available to indicate when a new derived
#' records has been added to a dataset.
#'
#' @details
#' When all records have same values within `by_vars` then this function will
#' retain those common values in the newly derived records. Otherwise new value
#' will be set to `NA`.
#'
#' @param dataset A data frame.
#' @param by_vars Variables to consider for generation of groupwise summary
#'   records. Providing the names of variables in [vars()] will create a
#'   groupwise summary and generate summary records for the specified groups.
#' @param fns List of formulas specifying variable to use for aggregations.
#'   This can include base functions like `mean`, `min`, `max`, `median`, `sd`,
#'   or `sum` or any other user-defined aggregation function.
#'   For example,
#'
#'   + When a summary function is same for one or more analysis variable, use
#'   `fns = list(vars(AVAL, CHG) ~ mean`).
#'   + If different summary function is required for each analysis variable,
#'   use `fns = list(AVAL ~ mean, CHG ~ sum(., na.rm = TRUE))`.
#'
#'   In general,
#'
#'   + LHS refer to the one or more variable to use for summarizing.
#'   + RHS refer to a **single** summary function.
#'
#'   In the formula representation e.g., `CHG ~ sum(., na.rm = TRUE)`, a `.`
#'   serves as the data to be summarized which refers to the variable `CHG`.
#' @param set_values_to A list of variable name-value pairs. Use this argument
#'   if you need to change the values of any newly derived records. Always new
#'   values in `set_values_to` should be equal to the length of analysis
#'   variable used in the `fns` argument. For example,
#'
#'   + On using single analysis variable in `fns = list(AVAL ~ mean)` and
#'   `set_values_to = vars(AVISITN = 9999, AVISIT = "Endpoint")` would change
#'   the value of `AVISITN` to `9999` and `AVISIT` to `Endpoint` instead of
#'   retaining.
#'   + Multiple analysis variables in `fns = list(vars(AVAL, CHG) ~ mean)` and
#'   `set_values_to = vars(AVISITN = c(9998, 9999))` would change `AVISITN` to
#'   `9998` for `AVAL` and `AVISITN` to `9999` for `CHG`.
#' @param drop_values_from Providing the names of variables in [vars()]
#'   will drop values and set as missing.
#'
#' @family row summary
#'
#' @return A data frame with derived records appended to original dataset.
#'
#' @export
#'
#' @examples
#' # Sample ADEG dataset ---
#' adeg <- tibble::tribble(
#'   ~USUBJID,   ~EGSEQ, ~EGREPNUM, ~PARAM,                 ~VISIT,      ~AVISIT,    ~EGDTC,                ~AVAL, ~TRTA,         ~SAFFL,
#'   "XYZ-1001",      1,         1,  "QTcF Interval (msec)", "SCREENING", "Baseline", "2016-02-24T07:50:16",  385,  "",            "Y",
#'   "XYZ-1001",      2,         2,  "QTcF Interval (msec)", "SCREENING", "Baseline", "2016-02-24T07:52:59",  399,  "",            "Y",
#'   "XYZ-1001",      3,         3,  "QTcF Interval (msec)", "SCREENING", "Baseline", "2016-02-24T07:56:07",  396,  "",            "Y",
#'   "XYZ-1001",      4,         1,  "QTcF Interval (msec)", "VISIT 2",   "Visit 2",  "2016-03-08T09:45:11",  384,  "Placebo",     "Y",
#'   "XYZ-1001",      5,         2,  "QTcF Interval (msec)", "VISIT 2",   "Visit 2",  "2016-03-08T09:48:07",  393,  "Placebo",     "Y",
#'   "XYZ-1001",      6,         3,  "QTcF Interval (msec)", "VISIT 2",   "Visit 2",  "2016-03-08T09:51:04",  388,  "Placebo",     "Y",
#'   "XYZ-1001",      7,         1,  "QTcF Interval (msec)", "VISIT 3",   "Visit 3",  "2016-03-22T10:45:03",  385,  "Placebo",     "Y",
#'   "XYZ-1001",      8,         2,  "QTcF Interval (msec)", "VISIT 3",   "Visit 3",  "2016-03-22T10:48:07",  394,  "Placebo",     "Y",
#'   "XYZ-1001",      9,         3,  "QTcF Interval (msec)", "VISIT 3",   "Visit 3",  "2016-03-22T10:51:05",  402,  "Placebo",     "Y",
#'   "XYZ-1002",      1,         1,  "QTcF Interval (msec)", "SCREENING", "Baseline", "2016-02-22T07:58:05",  399,  "",            "Y",
#'   "XYZ-1002",      2,         2,  "QTcF Interval (msec)", "SCREENING", "Baseline", "2016-02-22T07:58:05",  410,  "",            "Y",
#'   "XYZ-1002",      3,         3,  "QTcF Interval (msec)", "SCREENING", "Baseline", "2016-02-22T08:01:06",  392,  "",            "Y",
#'   "XYZ-1002",      4,         1,  "QTcF Interval (msec)", "VISIT 2",   "Visit 2",  "2016-03-06T09:50:04",  401,  "Active 20mg", "Y",
#'   "XYZ-1002",      5,         2,  "QTcF Interval (msec)", "VISIT 2",   "Visit 2",  "2016-03-06T09:53:51",  407,  "Active 20mg", "Y",
#'   "XYZ-1002",      6,         3,  "QTcF Interval (msec)", "VISIT 2",   "Visit 2",  "2016-03-06T09:56:21",  400,  "Active 20mg", "Y",
#'   "XYZ-1002",      7,         1,  "QTcF Interval (msec)", "VISIT 3",   "Visit 3",  "2016-03-24T10:50:07",  412,  "Active 20mg", "Y",
#'   "XYZ-1002",      8,         2,  "QTcF Interval (msec)", "VISIT 3",   "Visit 3",  "2016-03-24T10:53:08",  414,  "Active 20mg", "Y",
#'   "XYZ-1002",      9,         3,  "QTcF Interval (msec)", "VISIT 3",   "Visit 3",  "2016-03-24T10:56:05",  402,  "Active 20mg", "Y",
#' )
#'
#' # Summarize the average of the triplicate ECG interval values (AVAL)
#' derive_summary_records(adeg,
#'                    by_vars = vars(USUBJID, PARAM, AVISIT),
#'                    fns = list(AVAL ~ mean(., na.rm = TRUE)),
#'                    set_values_to = vars(DTYPE = "AVERAGE"))
#'
#' # Sample ADVS dataset ---
#' advs <- tibble::tribble(
#'   ~STUDYID, ~USUBJID,      ~VSSEQ, ~PARAMCD, ~PARAM,    ~AVAL, ~VSSTRESU, ~VISITNUM, ~VISIT,      ~VSDTC,
#'   "XYZ",    "XYZ-001-001",   1164, "WEIGHT",  "Weight",    99, "kg",              1, "Screening", "2018-03-19",
#'   "XYZ",    "XYZ-001-001",   1165, "WEIGHT",  "Weight",   101, "kg",              2, "Run-In"   , "2018-03-26",
#'   "XYZ",    "XYZ-001-001",   1166, "WEIGHT",  "Weight",   100, "kg",              3, "Baseline" , "2018-04-16",
#'   "XYZ",    "XYZ-001-001",   1167, "WEIGHT",  "Weight",    94, "kg",              4, "Week 24"  , "2018-09-30",
#'   "XYZ",    "XYZ-001-001",   1168, "WEIGHT",  "Weight",    92, "kg",              5, "Week 48"  , "2019-03-17",
#'   "XYZ",    "XYZ-001-001",   1169, "WEIGHT",  "Weight",    95, "kg",              6, "Week 52"  , "2019-04-14",
#' )
#'
#' # set new values to any variable, DTYPE = MAXIMUM refer to `MAX()` record and DTYPE = AVERAGE refer `MEAN()` record
#' # `set_values_to` must be of same length as new records
#' derive_summary_records(advs,
#'   by_vars = vars(STUDYID, USUBJID, PARAMCD),
#'   fns = list(AVAL ~ max, AVAL ~ mean),
#'   set_values_to = vars(DTYPE = c("MAXIMUM", "AVERAGE"))
#' )
#'
#' # drop retained value of VSSTRESU in the derived record
#' derive_summary_records(advs,
#'   by_vars = vars(STUDYID, USUBJID, PARAMCD),
#'   fns = list(AVAL ~ mean),
#'   set_values_to = vars(DTYPE = "MAXIMUM"),
#'   drop_values_from = vars(VSSTRESU)
#' )
#'
derive_summary_records <- function(dataset,
                                   by_vars,
                                   fns,
                                   set_values_to = NULL,
                                   drop_values_from = NULL) {

  # Is quosure?
  assert_that(is_vars(by_vars))
  if (!is.null(drop_values_from)) assert_that(is_vars(drop_values_from))

  by_vars <- vars2chr(by_vars)
  drop_values_from <- vars2chr(drop_values_from)

  # Assert variable names from input dataset
  assert_has_variables(dataset, c(by_vars, drop_values_from))

  # Manipulate functions as direct call for each analysis variable
  # Returns: A list of function call with attributes "variable" and "stats"
  funs <- manip_fun(dataset, fns, .env = caller_env())

  # Get input analysis variable
  summarise_vars <- map_chr(funs, attr, "variable")

  # Get variable values that are constant across the grouped records,
  # that do not change
  keep_vars <- setdiff(names(dataset),
                       union(
                         drop_values_from(dataset, by_vars),
                         drop_values_from))
  # For mutate input
  set_values <- NULL

  # Transpose into list-of-list
  if (!is.null(set_values_to)) {
    assert_that(is_quosures(set_values_to),
                msg = str_glue("`set_values_to` must be a `vars()` object, \\
                             not {friendly_type(typeof(set_values_to))}."))

    set_values <- set_values_to %>%
      map(quo_get_expr) %>%
      map(eval_tidy) %>%
      transpose()

    # Check new values and derived records have same length
    assert_that(
      are_records_same(set_values, funs, x_arg = "set_values_to", y_arg = "fns"))
  }

  # Summaries the analysis value and bind to the original dataset
  summary_data <- bind_rows(
    lapply(
      seq_along(funs),
      function(.x) {
        # Get unique values by grouping variable and do a left join with
        # summarised data
        dataset %>%
          distinct(!!! syms(by_vars), .keep_all = TRUE) %>%
          select(!!! syms(keep_vars)) %>%
          left_join(
            dataset %>%
              group_by(!!! syms(by_vars)) %>%
              summarise(!!! funs[.x]) %>%
              mutate(!!! set_values[[.x]]),
            by = by_vars
          )
      }
    )) %>%
    bind_rows(dataset, .) %>%
    arrange(!!! syms(by_vars))

  summary_data
}

#' Creates an anonymous function without giving it a name
#'
#' @param funs A function expression.
#' @param env The environment in which to evaluate the expression.
#'
#' @family row summary
#'
#' @return An anonymous function as `body`.
#'
#' ```
#' body(function(x = 4) g(x) + h(x))
#' #> g(x) + h(x)
#' ```
#' @noRd
as_inlined_function <- function(funs, env) {
  funs <- new_formula(NULL, f_rhs(funs), env = env)
  # Process unquote operator at inlining time
  funs <- expr_interp(funs)
  # Transform to a purrr-like lambda
  fn <- as_function(funs, env = env)

  body(fn) <- expr({
    # Transform the lambda body into a maskable quosure inheriting
    # from the execution environment
    `_quo` <- quo(!!body(fn))

    # Evaluate the quosure in the mask
    rlang::eval_bare(`_quo`, base::parent.frame())
  })

  structure(
    fn,
    class = "inline_function",
    stats = call_name(funs),
    formula = funs
  )
}

#' Utility to extract a function from an object of class formula
#'
#' @param .funs A two sided formula as a list (e.g. `list(A ~ mean)`).
#'
#'  + LHS refer to the one or more variable to use for summarizing.
#'  + RHS refer to a **single** summary function.
#' @param .env The environment in which to evaluate the expression.
#'
#' @family row summary
#'
#' @return A list of anonymous functions for each formula.
#'
#' @noRd
#'
#' @examples
#' fm <- list(cyl ~ max, hp ~ mean(., na.rm = TRUE))
#' as_fun_list(fm, caller_env())
as_fun_list <- function(.funs, .env) {

  if (is_formula(.funs)) {
    .funs <- list(.funs)
  }

  if (!purrr::every(.funs, is_formula)) {
    abort(
      c("Problem with input `fns`, expecting a two sided formula.",
        i = "LHS = an analysis variable and RHS = a function, or a function name.",
        i = "e.g `fns = list(aval ~ mean)`.")
    )
  }

  funs <- map(.funs, function(.x) {
    rhs <- f_rhs(.x)
    if (is_call(rhs)) {
      if (grepl("^list", deparse(rhs))) {
        abort("The LHS of `fns` must be a string or a function")
      }
      .x <- as_inlined_function(.x, env = .env)
    } else if (is_character(rhs) || is_symbol(rhs)) {
      fn <- get(rhs, .env, mode = "function")
      .x <- structure(fn, stats = as_string(rhs))
    } else {
      abort("Problem with handling RHS of `fns` argument.")
    }
    .x
  })

  funs
}

#' Utility to extract a variable from an object of class formula
#'
#' @inheritParams as_fun_list
#' @param dataset A dataset to conform input variables.
#'
#' @family row summary
#'
#' @return Variable as list
#'
#' @noRd
#'
#' @examples
#' fm <- list(cyl ~ max, hp ~ mean(., na.rm = TRUE))
#' get_fun_vars(fm, mtcars)
get_fun_vars <- function(f, dataset) {
  x <- map_if(f, is_quosure, quo_squash)
  lhs <- map_if(x, is_bare_formula, f_lhs)
  lhs <- map(lhs, function(.x) {
    if (is_symbol(.x)) as_string(.x)
    else eval_bare(.x, caller_env())
  })
  lhs <- map_if(lhs, is_quosures, vars2chr)
  walk(lhs, ~assert_has_variables(dataset, .))
  lhs
}

#' Reconstruct a call object from its components using [rlang::call2()]
#'
#' @param funs Return value of [as_fun_list]
#' @param vars Return value of [get_fun_vars]
#'
#' @family row summary
#'
#' @return An anonymous function with call.
#'
#' @noRd
as_call_list <- function(funs, vars) {
  out <- vector("list", length(funs) * max(lengths(vars)))
  dim(out) <- c(length(funs), max(lengths(vars)))

  for (i in seq_along(funs)) {
    isymvars <- syms(vars[[i]])
    for (j in seq_along(isymvars)) {
      out[[i, j]] <- call2(funs[[i]], isymvars[[j]])
      attr(out[[i, j]], "variable") <- as_string(isymvars[[j]])
      attr(out[[i, j]], "stats") <- attr(funs[[i]], "stats")
    }
  }
  dim(out) <- NULL
  out <- keep(out, ~!is.null(.))
  names(out) <- map_chr(out, attr, "variable")
  out
}

#' Utility function to manipulate `fns` argument of [derive_summary_records]
#'
#' @param dataset A data frame.
#' @param fns A two sided formula as a list (e.g. `list(A ~ mean)`).
#'
#'  + LHS refer to the one or more variable to use for summarizing.
#'  + RHS refer to a **single** summary function.
#' @param .env The environment in which to evaluate the expression.
#'
#' @family row summary
#'
#' @return An anonymous function with call.
#'
#' @noRd
manip_fun <- function(dataset, fns, .env) {
  fn_vars <- get_fun_vars(fns, dataset)
  fn_list <- as_fun_list(fns, .env = .env)
  as_call_list(fn_list, fn_vars)
}

#' Utility to find variable that have constant values across grouped variables
#'
#' `drop_values_from` help [derive_summary_records] to find all variable values that
#' are constant across the original records, that do not change.
#'
#' @param dataset A data frame.
#' @param group_vars Variables to group across records.
#'
#' @family row summary
#'
#' @return Variable vector.
#'
#' @noRd
drop_values_from <- function(dataset, group_vars) {
  non_group_vars <- setdiff(names(dataset), group_vars)

  # Get unique values within each group by variables
  unique_count <- dataset %>%
    group_by(!!! syms(group_vars)) %>%
    summarise_at(vars(!! non_group_vars), n_distinct) %>%
    ungroup() %>%
    select(!!! syms(non_group_vars))

  # If unique values > 1, should be dropped; else retain
  lgl_vars <- unique_count %>%
    map_lgl(~ all(!.x > 1)) %>%
    which() %>%
    names()

  setdiff(non_group_vars, lgl_vars)
}
