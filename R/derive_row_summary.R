#' Add new row within by groups using aggregation functions
#'
#' It is not uncommon to have an analysis need whereby one needs to derive an
#' analysis value (`AVAL`) from multiple rows. The ADaM basic dataset structure
#' variable `DTYPE` is available to indicate when a new derived row has been
#' added to a dataset.
#'
#' @param dataset A data frame.
#' @param by_vars Variables to consider for generation of groupwise summary
#'   rows. Providing the names of variables in [c()] will create a groupwise
#'   summary and generate summary rows for the specified groups.
#' @param fns List of formulas specifying variable to use for aggregations.
#'   This can include base functions like `mean`, `min`, `max`, `median`, `sd`,
#'   or `sum` or any other user-defined aggregation function. For example,
#'   `fns = list(AVAL ~ mean, CHG ~ sum(., na.rm = TRUE)`. In the formula
#'   representation, a `.` serves as the data to be summarized
#'   (e.g., `sum(CHG, na.rm = TRUE)`).
#'
#'  + LHS refer to the one or more variable to use for summarizing.
#'  + RHS refer to a **single** summary function.
#'
#' @param set_values_to A list of variable name-value pairs. Use this argument
#'   if you need to change the values of any newly derived rows. For example,
#'   `values_set = list(AVISITN = 9999, AVISIT= "Endpoint")` would change of
#'   value of AVISITN to 9999 and AVISIT to Endpoint instead of retaining.
#' @param drop_values_from Providing the names of variables in [dplyr::vars()]
#'   will drop values and set as missing.
#'
#' @family row summary
#'
#' @return A data frame with derived records appended to original dataset.
#'
#' @export
#'
#' @examples
#' adeg <- tibble::tribble(
#' ~USUBJID,~EGSEQ,~EGREPNUM,~PARAM,~VISIT,~AVISIT,~EGDTC,~AVAL,~TRTA,~SAFFL,
#' "XYZ-1001",1,1,"QTcF Interval (msec)","SCREENING","Baseline","2016-02-24T07:50:16",385,"","Y",
#' "XYZ-1001",2,2,"QTcF Interval (msec)","SCREENING","Baseline","2016-02-24T07:52:59",399,"","Y",
#' "XYZ-1001",3,3,"QTcF Interval (msec)","SCREENING","Baseline","2016-02-24T07:56:07",396,"","Y",
#' "XYZ-1001",4,1,"QTcF Interval (msec)","VISIT 2","Visit 2","2016-03-08T09:45:11",384,"Placebo","Y",
#' "XYZ-1001",5,2,"QTcF Interval (msec)","VISIT 2","Visit 2","2016-03-08T09:48:07",393,"Placebo","Y",
#' "XYZ-1001",6,3,"QTcF Interval (msec)","VISIT 2","Visit 2","2016-03-08T09:51:04",388,"Placebo","Y",
#' "XYZ-1001",7,1,"QTcF Interval (msec)","VISIT 3","Visit 3","2016-03-22T10:45:03",385,"Placebo","Y",
#' "XYZ-1001",8,2,"QTcF Interval (msec)","VISIT 3","Visit 3","2016-03-22T10:48:07",394,"Placebo","Y",
#' "XYZ-1001",9,3,"QTcF Interval (msec)","VISIT 3","Visit 3","2016-03-22T10:51:05",402,"Placebo","Y",
#' "XYZ-1002",1,1,"QTcF Interval (msec)","SCREENING","Baseline","2016-02-22T07:58:05",399,"","Y",
#' "XYZ-1002",2,2,"QTcF Interval (msec)","SCREENING","Baseline","2016-02-22T07:58:05",410,"","Y",
#' "XYZ-1002",3,3,"QTcF Interval (msec)","SCREENING","Baseline","2016-02-22T08:01:06",392,"","Y",
#' "XYZ-1002",4,1,"QTcF Interval (msec)","VISIT 2","Visit 2","2016-03-06T09:50:04",401,"Active 20mg","Y",
#' "XYZ-1002",5,2,"QTcF Interval (msec)","VISIT 2","Visit 2","2016-03-06T09:53:51",407,"Active 20mg","Y",
#' "XYZ-1002",6,3,"QTcF Interval (msec)","VISIT 2","Visit 2","2016-03-06T09:56:21",400,"Active 20mg","Y",
#' "XYZ-1002",7,1,"QTcF Interval (msec)","VISIT 3","Visit 3","2016-03-24T10:50:07",412,"Active 20mg","Y",
#' "XYZ-1002",8,2,"QTcF Interval (msec)","VISIT 3","Visit 3","2016-03-24T10:53:08",414,"Active 20mg","Y",
#' "XYZ-1002",9,3,"QTcF Interval (msec)","VISIT 3","Visit 3","2016-03-24T10:56:05",402,"Active 20mg","Y",
#' )
#'
#' # summarize the average of the triplicate ECG interval values (AVAL)
#' derive_row_summary(adeg,
#'                    by_vars = vars(USUBJID, PARAM, AVISIT),
#'                    fns = list(AVAL ~ mean(., na.rm = TRUE)),
#'                    set_values_to = list(DTYPE = "AVERAGE"))
derive_row_summary <- function(dataset,
                               by_vars,
                               fns,
                               set_values_to = list(),
                               drop_values_from = NULL) {

  # Assert variable names from input dataset
  by_vars <- assert_summary_vars(dataset,
                                 enexpr(by_vars),
                                 .vars_arg = "by_vars")
  if (!is.null(drop_values_from)) {
    drop_values_from <-
      assert_summary_vars(dataset,
                          enexpr(drop_values_from),
                          .vars_arg = "drop_values_from")
  }

  if (!is.list(set_values_to)) {
    abort(
      sprintf(
        "`set_values_to` must be a list, not %s.",
        friendly_type(typeof(enexpr(set_values_to))) ))
  }

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

  # Transpose into list-of-list
  set_values <- transpose(set_values_to)

  if (length(funs) != length(set_values)) {
    abort(
      c("`set_values_to` must have consistent length to `funs`:",
        x = str_glue("`set_values_to` has length {length(set_values)}"),
        x = str_glue("`funs` has length {length(funs)}")
      )
    )
  }

  funs <- list(funs) %>%
    transpose() %>%
    map2(names(funs), ~ set_names(.x, .y))

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
              summarise(!!! funs[[.x]]) %>%
              mutate(!!! set_values[[.x]]),
            by = by_vars
          )
      }
    )) %>%
    bind_rows(dataset, .) %>%
    arrange(!!! syms(by_vars))

  summary_data
}

#' Assert variables within a dataset
#'
#' `assert_summary_vars()` throws an error when the input variable doesn't
#' conform within the given dataset.
#'
#' @param dataset A data frame.
#' @param vars Input variable name as a character/numeric vector, a symbol or
#'   a [dplyr::vars()] object.
#' @param .env The environment in which to evaluate the expression.
#' @param .vars_arg Name of argument being checked. This is used in error
#'   messages.
#'
#' @family row summary
#'
#' @return Variable as vector
#'
#' @noRd
#'
#' @examples
#' # Variable as symbol
#' assert_summary_vars(mtcars, expr(cyl))
#'
#' # Variable as character
#' assert_summary_vars(mtcars, "cyl")
#' assert_summary_vars(mtcars, c("cyl", "hp"))
#'
#' # Numeric selection
#' assert_summary_vars(mtcars, 1:5)
#'
#' # Variables as `vars()` object
#' assert_summary_vars(mtcars, vars(cyl, hp, wt))
assert_summary_vars <- function(dataset, vars, .env = caller_env(), .vars_arg = "vars") {
  if (is_symbol(vars)) {
    vars <- as_string(vars)
  } else if (is_quosures(vars) || is_character(vars)) {
    vars <- vars
  } else if (is_integerish(vars) || is_call(vars)) {
    vars <- eval_bare(vars, .env)
  } else {
    abort(
      sprintf(
        "`%s` must be a character/numeric vector, a symbol or a `vars()` object, not %s.",
        .vars_arg, friendly_type(typeof(vars))
      )
    )
  }
  out <- select(dataset, !!! vars) %>% names()
  out
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
  vars <- map(lhs, function(x) assert_summary_vars(dataset, x))
  vars
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

#' Utility function to manipulate `fns` argument of [derive_row_summary]
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
#' `drop_values_from` help [derive_row_summary] to find all variable values that
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
    map_lgl(function(.x) {
      res <- ifelse(.x > 1, FALSE, TRUE)
      if (all(res)) TRUE else FALSE
    }) %>%
    which() %>%
    names()

  setdiff(non_group_vars, lgl_vars)
}
