#' Add New Records Within By Groups Using Aggregation Functions
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
#'   This can include base functions like `mean()`, `min()`, `max()`, `median()`,
#'    `sd()`, or `sum()` or any other user-defined aggregation function.
#'   For example,
#'
#'   + When a summary function is same for one or more analysis variable, use
#'   `fns = list(vars(AVAL, CHG) ~ mean`).
#'   + If a different summary function is required for each analysis variable,
#'   use `fns = list(AVAL ~ mean, CHG ~ sum(., na.rm = TRUE))`.
#'
#'   In general,
#'
#'   + LHS refer to the one or more variable to use for summarizing.
#'   + RHS refer to a **single** summary function.
#'
#'   In the formula representation e.g., `CHG ~ sum(., na.rm = TRUE)`, a `.`
#'   serves as the data to be summarized which refers to the variable `CHG`.
#' @param filter_rows Filter condition as logical expression to apply during
#'   summary calculation. By default, filtering expressions are computed within
#'   `by_vars` as this will help when an aggregating, lagging, or ranking
#'   function is involved.
#'
#'   For example,
#'
#'   + `filter_rows = (AVAL > mean(AVAL, na.rm = TRUE))` will filter all AVAL
#'   values greater than mean of AVAL with in `by_vars`.
#'   + `filter_rows = (dplyr::n() > 2)` will filter n count of `by_vars` greater
#'   than 2.
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
#' @author Vignesh Thanikachalam
#'
#' @return A data frame with derived records appended to original dataset.
#'
#' @export
#'
#' @examples
#' adeg <- tibble::tribble(
#'   ~USUBJID, ~EGSEQ, ~PARAM,             ~AVISIT,    ~EGDTC,            ~AVAL, ~TRTA,
#'   "XYZ-1001",    1, "QTcF Int. (msec)", "Baseline", "2016-02-24T07:50",  385, "",
#'   "XYZ-1001",    2, "QTcF Int. (msec)", "Baseline", "2016-02-24T07:52",  399, "",
#'   "XYZ-1001",    3, "QTcF Int. (msec)", "Baseline", "2016-02-24T07:56",  396, "",
#'   "XYZ-1001",    4, "QTcF Int. (msec)", "Visit 2",  "2016-03-08T09:45",  384, "Placebo",
#'   "XYZ-1001",    5, "QTcF Int. (msec)", "Visit 2",  "2016-03-08T09:48",  393, "Placebo",
#'   "XYZ-1001",    6, "QTcF Int. (msec)", "Visit 2",  "2016-03-08T09:51",  388, "Placebo",
#'   "XYZ-1001",    7, "QTcF Int. (msec)", "Visit 3",  "2016-03-22T10:45",  385, "Placebo",
#'   "XYZ-1001",    8, "QTcF Int. (msec)", "Visit 3",  "2016-03-22T10:48",  394, "Placebo",
#'   "XYZ-1001",    9, "QTcF Int. (msec)", "Visit 3",  "2016-03-22T10:51",  402, "Placebo",
#'   "XYZ-1002",    1, "QTcF Int. (msec)", "Baseline", "2016-02-22T07:58",  399, "",
#'   "XYZ-1002",    2, "QTcF Int. (msec)", "Baseline", "2016-02-22T07:58",  410, "",
#'   "XYZ-1002",    3, "QTcF Int. (msec)", "Baseline", "2016-02-22T08:01",  392, "",
#'   "XYZ-1002",    4, "QTcF Int. (msec)", "Visit 2",  "2016-03-06T09:50",  401, "Active 20mg",
#'   "XYZ-1002",    5, "QTcF Int. (msec)", "Visit 2",  "2016-03-06T09:53",  407, "Active 20mg",
#'   "XYZ-1002",    6, "QTcF Int. (msec)", "Visit 2",  "2016-03-06T09:56",  400, "Active 20mg",
#'   "XYZ-1002",    7, "QTcF Int. (msec)", "Visit 3",  "2016-03-24T10:50",  412, "Active 20mg",
#'   "XYZ-1002",    8, "QTcF Int. (msec)", "Visit 3",  "2016-03-24T10:53",  414, "Active 20mg",
#'   "XYZ-1002",    9, "QTcF Int. (msec)", "Visit 3",  "2016-03-24T10:56",  402, "Active 20mg",
#')
#'
#' # Summarize the average of the triplicate ECG interval values (AVAL)
#' derive_summary_records(
#'   adeg,
#'   by_vars = vars(USUBJID, PARAM, AVISIT),
#'   fns = list(AVAL ~ mean(., na.rm = TRUE)),
#'   set_values_to = vars(DTYPE = "AVERAGE")
#' )
#'
#' advs <- tibble::tribble(
#'   ~USUBJID,     ~VSSEQ, ~PARAM,  ~AVAL, ~VSSTRESU, ~VISIT,      ~VSDTC,
#'   "XYZ-001-001",  1164, "Weight",   99, "kg",      "Screening", "2018-03-19",
#'   "XYZ-001-001",  1165, "Weight",  101, "kg",      "Run-In"   , "2018-03-26",
#'   "XYZ-001-001",  1166, "Weight",  100, "kg",      "Baseline" , "2018-04-16",
#'   "XYZ-001-001",  1167, "Weight",   94, "kg",      "Week 24"  , "2018-09-30",
#'   "XYZ-001-001",  1168, "Weight",   92, "kg",      "Week 48"  , "2019-03-17",
#'   "XYZ-001-001",  1169, "Weight",   95, "kg",      "Week 52"  , "2019-04-14",
#' )
#'
#' # Set new values to any variable. Here, `DTYPE = MAXIMUM` refers to `max()` records
#' # and `DTYPE = AVERAGE` refers to `mean()` records.
#' # `set_values_to` must be of the same length as `fns`
#' derive_summary_records(
#'   advs,
#'   by_vars = vars(USUBJID, PARAM),
#'   fns = list(AVAL ~ max, AVAL ~ mean),
#'   set_values_to = vars(DTYPE = c("MAXIMUM", "AVERAGE"))
#' )
#'
#' # Drop retained value of VSSTRESU in the derived record
#' derive_summary_records(
#'   advs,
#'   by_vars = vars(USUBJID, PARAM),
#'   fns = list(AVAL ~ mean),
#'   set_values_to = vars(DTYPE = "MAXIMUM"),
#'   drop_values_from = vars(VSSTRESU)
#' )
#'
#' # Sample ADEG dataset with triplicate record for only AVISIT = 'Baseline' ---
#' adeg <- tibble::tribble(
#'   ~USUBJID, ~EGSEQ, ~PARAM,             ~AVISIT,    ~EGDTC,            ~AVAL, ~TRTA,
#'   "XYZ-1001",    1, "QTcF Int. (msec)", "Baseline", "2016-02-24T07:50",  385, "",
#'   "XYZ-1001",    2, "QTcF Int. (msec)", "Baseline", "2016-02-24T07:52",  399, "",
#'   "XYZ-1001",    3, "QTcF Int. (msec)", "Baseline", "2016-02-24T07:56",  396, "",
#'   "XYZ-1001",    4, "QTcF Int. (msec)", "Visit 2",  "2016-03-08T09:48",  393, "Placebo",
#'   "XYZ-1001",    5, "QTcF Int. (msec)", "Visit 2",  "2016-03-08T09:51",  388, "Placebo",
#'   "XYZ-1001",    6, "QTcF Int. (msec)", "Visit 3",  "2016-03-22T10:48",  394, "Placebo",
#'   "XYZ-1001",    7, "QTcF Int. (msec)", "Visit 3",  "2016-03-22T10:51",  402, "Placebo",
#'   "XYZ-1002",    1, "QTcF Int. (msec)", "Baseline", "2016-02-22T07:58",  399, "",
#'   "XYZ-1002",    2, "QTcF Int. (msec)", "Baseline", "2016-02-22T07:58",  410, "",
#'   "XYZ-1002",    3, "QTcF Int. (msec)", "Baseline", "2016-02-22T08:01",  392, "",
#'   "XYZ-1002",    4, "QTcF Int. (msec)", "Visit 2",  "2016-03-06T09:53",  407, "Active 20mg",
#'   "XYZ-1002",    5, "QTcF Int. (msec)", "Visit 2",  "2016-03-06T09:56",  400, "Active 20mg",
#'   "XYZ-1002",    6, "QTcF Int. (msec)", "Visit 3",  "2016-03-24T10:53",  414, "Active 20mg",
#'   "XYZ-1002",    7, "QTcF Int. (msec)", "Visit 3",  "2016-03-24T10:56",  402, "Active 20mg",
#')
#'
#' # Summarize the average of AVAL for AVISIT records greater than 2
#' derive_summary_records(
#'   adeg,
#'   by_vars = vars(USUBJID, PARAM, AVISIT),
#'   fns = list(AVAL ~ mean(., na.rm = TRUE)),
#'   filter_rows = dplyr::n() > 2,
#'   set_values_to = vars(DTYPE = "AVERAGE")
#' )
derive_summary_records <- function(dataset,
                                   by_vars,
                                   fns,
                                   filter_rows = NULL,
                                   set_values_to = NULL,
                                   drop_values_from = NULL) {
  assert_vars(by_vars)
  assert_vars(drop_values_from, optional = TRUE)
  assert_data_frame(dataset, required_vars = quo_c(by_vars, drop_values_from))

  by_vars <- vars2chr(by_vars)
  drop_values_from <- vars2chr(drop_values_from)

  # Manipulate functions as direct call for each analysis variable
  # Returns: A list of function call with attributes "variable" and "stats"
  funs <- manip_fun(dataset, fns, .env = caller_env())

  # Get input analysis variable
  summarise_vars <- map_chr(funs, attr, "variable")

  # Get variable values that are constant across the grouped records,
  # that do not change
  keep_vars <- setdiff(
    names(dataset),
    union(drop_values_from(dataset, by_vars), drop_values_from)
  )
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
      are_records_same(set_values, funs, x_arg = "set_values_to", y_arg = "fns")
    )
  }

  filter_rows <- assert_filter_cond(enquo(filter_rows), optional = TRUE)

  if (!quo_is_null(filter_rows)) {
    subset_ds <- dataset %>%
      group_by(!!! syms(by_vars)) %>%
      filter(!! filter_rows)
  } else {
    subset_ds <- dataset
  }

  # Summaries the analysis value and bind to the original dataset
  summary_data <- bind_rows(
    lapply(
      seq_along(funs),
      function(.x) {
        # Get unique values by grouping variable and do a left join with
        # summarised data
        subset_ds %>%
          distinct(!!! syms(by_vars), .keep_all = TRUE) %>%
          select(!!! syms(keep_vars)) %>%
          left_join(
            subset_ds %>%
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
    `_quo` <- rlang::quo(!!body(fn)) # nolint

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
  funs <- map(.funs, function(.x) {
    rhs <- f_rhs(.x)
    if (is_call(rhs)) {
      if (is_call(rhs, c("list", "vars"))) {
        abort("The LHS of `fns` must be a string or a function")
      }
      .x <- as_inlined_function(.x, env = .env)
    } else if (is_character(rhs) || is_symbol(rhs)) {
      rhs <- as_string(rhs)
      fn <- rlang::as_closure(rhs)
      .x <- structure(fn, stats = as_string(rhs))
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
  if (is_formula(fns)) {
    fns <- list(fns)
  }

  if (!every(fns, is_formula)) {
    abort(
      str_glue(
        "Problem with input `fns`, expecting a two sided formula.
        LHS = an analysis variable and RHS = a function, or a function name.
        e.g `fns = list(aval ~ mean)`."))
  }

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
