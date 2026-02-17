#' Add New Variable(s) to the Input Dataset Based on Variables from Another
#' Dataset
#'
#' Add new variable(s) to the input dataset based on variables from another
#' dataset. The observations to merge can be selected by a condition
#' (`filter_add` argument) and/or selecting the first or last observation for
#' each by group (`order` and `mode` argument).
#'
#' @param dataset
#'
#' `r roxygen_param_dataset(expected_vars = c("by_vars"))`
#'
#' @permitted [dataset]
#'
#' @param dataset_add Additional dataset
#'
#'   The variables specified by the `by_vars`, the `new_vars`, and the `order`
#'   argument are expected.
#'
#' @permitted [dataset]
#'
#' @param by_vars Grouping variables
#'
#'   The input dataset and the selected observations from the additional dataset
#'   are merged by the specified variables.
#'
#'   `r roxygen_param_by_vars(rename = TRUE)`
#'
#' @permitted [var_list]
#'
#' @param order Sort order
#'
#'   If the argument is set to a non-null value, for each by group the first or
#'   last observation from the additional dataset is selected with respect to the
#'   specified order.
#'
#'   Variables defined by the `new_vars` argument can be used in the sort order.
#'
#'   `r roxygen_order_na_handling()`
#'
#' @permitted [var_list]
#'
#' @param new_vars Variables to add
#'
#'   The specified variables from the additional dataset are added to the output
#'   dataset. Variables can be renamed by naming the element, i.e., `new_vars =
#'   exprs(<new name> = <old name>)`.
#'
#'   For example `new_vars = exprs(var1, var2)` adds variables `var1` and `var2`
#'   from `dataset_add` to the input dataset.
#'
#'   And `new_vars = exprs(var1, new_var2 = old_var2)` takes `var1` and
#'   `old_var2` from `dataset_add` and adds them to the input dataset renaming
#'   `old_var2` to `new_var2`.
#'
#'   Values of the added variables can be modified by specifying an expression.
#'   For example, `new_vars = LASTRSP = exprs(str_to_upper(AVALC))` adds the
#'   variable `LASTRSP` to the dataset and sets it to the upper case value of
#'   `AVALC`.
#'
#'   If the argument is not specified or set to `NULL`, all variables from the
#'   additional dataset (`dataset_add`) are added. In the case when a variable
#'   exists in both datasets, an error is issued to ensure the user either adds
#'   to `by_vars`, removes or renames.
#'
#' @permitted [var_list]
#'
#' @param filter_add Filter for additional dataset (`dataset_add`)
#'
#'   Only observations fulfilling the specified condition are taken into account
#'   for merging. If the argument is not specified, all observations are
#'   considered.
#'
#'   Variables defined by the `new_vars` argument can be used in the filter
#'   condition.
#'
#' @permitted [condition]
#'
#' @param mode Selection mode
#'
#'   Determines if the first or last observation is selected. If the `order`
#'   argument is specified, `mode` must be non-null.
#'
#'   If the `order` argument is not specified, the `mode` argument is ignored.
#'
#' @permitted [mode]
#'
#' @param exist_flag Exist flag
#'
#'   If the argument is specified (e.g., `exist_flag = FLAG`), the specified
#'   variable (e.g., `FLAG`) is added to the input dataset. This variable will
#'   be the value provided in `true_value` for all selected records from `dataset_add`
#'   which are merged into the input dataset, and the value provided in `false_value` otherwise.
#'
#' @permitted [var]
#'
#' @param true_value True value
#'
#'   The value for the specified variable `exist_flag`, applicable to
#'   the first or last observation (depending on the mode) of each by group.
#'
#' @permitted [char_scalar]
#'
#' @param false_value False value
#'
#'   The value for the specified variable `exist_flag`, NOT applicable to
#'   the first or last observation (depending on the mode) of each by group.
#'
#' @permitted [char_scalar]
#'
#' @param missing_values Values for non-matching observations
#'
#'   For observations of the input dataset (`dataset`) which do not have a
#'   matching observation in the additional dataset (`dataset_add`) the values
#'   of the specified variables are set to the specified value. Only variables
#'   specified for `new_vars` can be specified for `missing_values`.
#'
#' @permitted [expr_list_formula]
#'
#' @param check_type Check uniqueness?
#'
#'   If `"warning"`, `"message"`, or `"error"` is specified, the specified message is issued
#'   if the observations of the (restricted) additional dataset are not unique
#'   with respect to the by variables and the order.
#'
#'   If the `order` argument is not specified, the `check_type` argument is ignored:
#'   if the observations of the (restricted) additional dataset are not unique with respect
#'   to the by variables, an error is issued.
#'
#' @permitted [msg_type]
#'
#' @param duplicate_msg Message of unique check
#'
#'   If the uniqueness check fails, the specified message is displayed.
#'
#' @default
#'
#'   ```{r echo=TRUE, eval=FALSE}
#'   paste(
#'     "Dataset {.arg dataset_add} contains duplicate records with respect to",
#'     "{.var {vars2chr(by_vars)}}."
#'   )
#'   ```
#'
#' @permitted [msg]
#'
#' @param relationship Expected merge-relationship between the `by_vars`
#'   variable(s) in `dataset` (input dataset) and the `dataset_add` (additional dataset)
#'    containing the additional `new_vars`.
#'
#'   This argument is passed to the `dplyr::left_join()` function. See
#'   <https://dplyr.tidyverse.org/reference/mutate-joins.html#arguments> for
#'   more details.
#'
#' @permitted [merge_rel]
#'
#' @return The output dataset contains all observations and variables of the
#'   input dataset and additionally the variables specified for `new_vars` from
#'   the additional dataset (`dataset_add`).
#'
#' @details
#'
#'   1. The new variables (`new_vars`) are added to the additional dataset
#'   (`dataset_add`).
#'
#'   1. The records from the additional dataset (`dataset_add`) are restricted
#'   to those matching the `filter_add` condition.
#'
#'   1. If `order` is specified, for each by group the first or last observation
#'   (depending on `mode`) is selected.
#'
#'   1. The variables specified for `new_vars` are merged to the input dataset
#'   using `left_join()`. I.e., the output dataset contains all observations
#'   from the input dataset. For observations without a matching observation in
#'   the additional dataset the new variables are set as specified by
#'   `missing_values` (or to `NA` for variables not in `missing_values`).
#'   Observations in the additional dataset which have no matching observation
#'   in the input dataset are ignored.
#'
#' @family der_gen
#' @keywords der_gen
#'
#' @export
#'
#' @examplesx
#'
#' @caption Note on usage versus `derive_vars_joined()`
#' @info The question between using `derive_vars_merged()` or the more powerful
#'   `derive_vars_joined()` comes down to how you need to select the observations
#'   to be merged.
#'
#' - If the observations from `dataset_add` to merge can be selected
#'   by a condition (`filter_add`) using *only* variables from `dataset_add`, then
#'   always use `derive_vars_merged()` as it requires less resources (time and
#'   memory). A common example of this would be a randomization date in `ADSL`,
#'   where you are simply merging on a date from `DS` according to a certain
#'   `DSDECOD` condition such as `DSDECOD == "RANDOMIZATION"`.
#' - However, if the selection of the observations from `dataset_add` can depend
#'   on variables from *both* datasets, then use `derive_vars_joined()`. An
#'   example of this would be assigning period variables from `ADSL` to an `ADAE`,
#'   where you now need to check each adverse event start date against the period
#'   start and end dates to decide which period value to join.
#' @caption Basic merge of a full dataset
#' @info Merge all demographic variables onto a vital signs dataset.
#'
#' - The variable `DOMAIN` exists in both datasets so note the use of
#'   `select(dm, -DOMAIN)` in the `dataset_add` argument. Without this an error
#'   would be issued to notify the user.
#' @code
#' library(tibble)
#' library(dplyr, warn.conflicts = FALSE)
#' vs <- tribble(
#'   ~DOMAIN,  ~USUBJID, ~VSTESTCD, ~VISIT,      ~VSSTRESN, ~VSDTC,
#'   "VS",     "01",     "HEIGHT",  "SCREENING",     178.0, "2013-08-20",
#'   "VS",     "01",     "WEIGHT",  "SCREENING",      81.9, "2013-08-20",
#'   "VS",     "01",     "WEIGHT",  "BASELINE",       82.1, "2013-08-29",
#'   "VS",     "01",     "WEIGHT",  "WEEK 2",         81.9, "2013-09-15",
#'   "VS",     "01",     "WEIGHT",  "WEEK 4",         82.6, "2013-09-24",
#'   "VS",     "02",     "WEIGHT",  "BASELINE",       58.6, "2014-01-11"
#' ) %>%
#'   mutate(STUDYID = "AB42")
#'
#' dm <- tribble(
#'   ~DOMAIN, ~USUBJID, ~AGE, ~AGEU,
#'   "DM",    "01",       61, "YEARS",
#'   "DM",    "02",       64, "YEARS",
#'   "DM",    "03",       85, "YEARS"
#' ) %>%
#'   mutate(STUDYID = "AB42")
#'
#' derive_vars_merged(
#'   vs,
#'   dataset_add = select(dm, -DOMAIN),
#'   by_vars = exprs(STUDYID, USUBJID)
#' ) %>%
#'   select(USUBJID, VSTESTCD, VISIT, VSSTRESN, AGE, AGEU)
#'
#' @caption Merge only the first/last value (`order` and `mode`)
#' @info Merge the last occurring weight for each subject to the demographics dataset.
#'
#' - To enable sorting by visit date `convert_dtc_to_dtm()` is used to convert
#'   to a datetime, within the `order` argument.
#' - Then the `mode` argument is set to `"last"` to ensure the last sorted value
#'   is taken. Be cautious if `NA` values are possible in the `order` variables -
#'   see [Sort Order](https://pharmaverse.github.io/admiral/articles/generic.html#sort_order).
#' - The `filter_add` argument is used to restrict the vital signs records only
#'   to weight assessments.
#' @code
#' derive_vars_merged(
#'   dm,
#'   dataset_add = vs,
#'   by_vars = exprs(STUDYID, USUBJID),
#'   order = exprs(convert_dtc_to_dtm(VSDTC)),
#'   mode = "last",
#'   new_vars = exprs(LSTWT = VSSTRESN),
#'   filter_add = VSTESTCD == "WEIGHT"
#' ) %>%
#'   select(USUBJID, AGE, AGEU, LSTWT)
#'
#' @caption Handling duplicates (`check_type`)
#' @info The source records are checked regarding duplicates with respect to the
#'   by variables and the order specified.
#'   By default, a warning is issued if any duplicates are found.
#'   Note the results here with a new vital signs dataset containing a
#'   duplicate last weight assessment date.
#' @code [expected_cnds = "duplicate_records"]
#' vs_dup <- tribble(
#'   ~DOMAIN,  ~USUBJID, ~VSTESTCD, ~VISIT,      ~VSSTRESN, ~VSDTC,
#'   "VS",     "01",     "WEIGHT",  "WEEK 2",        81.1, "2013-09-24",
#'   "VS",     "01",     "WEIGHT",  "WEEK 4",        82.6, "2013-09-24"
#' ) %>%
#'   mutate(STUDYID = "AB42")
#'
#' derive_vars_merged(
#'   dm,
#'   dataset_add = vs_dup,
#'   by_vars = exprs(STUDYID, USUBJID),
#'   order = exprs(convert_dtc_to_dtm(VSDTC)),
#'   mode = "last",
#'   new_vars = exprs(LSTWT = VSSTRESN),
#'   filter_add = VSTESTCD == "WEIGHT"
#' ) %>%
#'   select(USUBJID, AGE, AGEU, LSTWT)
#'
#' @info For investigating the issue, the dataset of the duplicate source records
#'   can be obtained by calling `get_duplicates_dataset()`:
#' @code
#' get_duplicates_dataset()
#'
#' @info Common options to solve the issue:
#' - Specifying additional variables for `order` - this is the most common approach,
#'   adding something like a sequence variable.
#' - Restricting the source records by specifying/updating the `filter_add` argument.
#' - Setting `check_type = "none"` to ignore any duplicates, but then in this case
#'   the last occurring record would be chosen according to the sort order of the
#'   input `dataset_add`. This is not often advisable, unless the order has no impact
#'   on the result, as the temporary sort order can be prone to variation across
#'   an ADaM script.
#'
#' @caption Modify values dependent on the merge (`new_vars` and `missing_values`)
#' @info For the last occurring weight for each subject, add a categorization of
#'   which visit it occurred at to the demographics dataset.
#'
#' - In the `new_vars` argument, other functions can be utilized to modify the
#'   merged values. For example, in the below case we want to categorize the
#'   visit as `"BASELINE"` or `"POST-BASELINE"` using `if_else()`.
#' - The `missing_values` argument assigns a specific value for subjects with
#'   no matching observations - see subject `"03"` in the below example.
#' @code
#' derive_vars_merged(
#'   dm,
#'   dataset_add = vs,
#'   by_vars = exprs(STUDYID, USUBJID),
#'   order = exprs(convert_dtc_to_dtm(VSDTC)),
#'   mode = "last",
#'   new_vars = exprs(
#'     LSTWTCAT = if_else(VISIT == "BASELINE", "BASELINE", "POST-BASELINE")
#'   ),
#'   filter_add = VSTESTCD == "WEIGHT",
#'   missing_values = exprs(LSTWTCAT = "MISSING")
#' ) %>%
#'   select(USUBJID, AGE, AGEU, LSTWTCAT)
#'
#' @caption Check existence of records to merge (`exist_flag`, `true_value` and `false_value`)
#' @info Similar to the above example, now we prefer to have a separate flag
#'    variable to show whether a selected record was merged.
#'
#' - The name of the new variable is set with the `exist_flag` argument.
#' - The values of this new variable are assigned via the `true_value` and
#'   `false_value` arguments.
#' @code
#' derive_vars_merged(
#'   dm,
#'   dataset_add = vs,
#'   by_vars = exprs(STUDYID, USUBJID),
#'   order = exprs(convert_dtc_to_dtm(VSDTC)),
#'   mode = "last",
#'   new_vars = exprs(
#'     LSTWTCAT = if_else(VISIT == "BASELINE", "BASELINE", "POST-BASELINE")
#'   ),
#'   filter_add = VSTESTCD == "WEIGHT",
#'   exist_flag = WTCHECK,
#'   true_value = "Y",
#'   false_value = "MISSING"
#' ) %>%
#'   select(USUBJID, AGE, AGEU, LSTWTCAT, WTCHECK)
#'
#' @caption Creating more than one variable from the merge (`new_vars`)
#' @info Derive treatment start datetime and associated imputation flags.
#'
#' - In this example we first impute exposure datetime and associated flag
#'   variables as a separate first step to be used in the `order` argument.
#' - In the `new_vars` arguments, you can see how both datetime and the date and
#'   time imputation flags are all merged in one call.
#' @code
#' ex <- tribble(
#'   ~DOMAIN, ~USUBJID, ~EXSTDTC,
#'   "EX",    "01",     "2013-08-29",
#'   "EX",    "01",     "2013-09-16",
#'   "EX",    "02",     "2014-01-11",
#'   "EX",    "02",     "2014-01-25"
#' ) %>%
#'   mutate(STUDYID = "AB42")
#'
#' ex_ext <- derive_vars_dtm(
#'   ex,
#'   dtc = EXSTDTC,
#'   new_vars_prefix = "EXST",
#'   highest_imputation = "M"
#' )
#'
#' derive_vars_merged(
#'   dm,
#'   dataset_add = ex_ext,
#'   by_vars = exprs(STUDYID, USUBJID),
#'   new_vars = exprs(TRTSDTM = EXSTDTM, TRTSDTF = EXSTDTF, TRTSTMF = EXSTTMF),
#'   order = exprs(EXSTDTM),
#'   mode = "first"
#' ) %>%
#'   select(USUBJID, TRTSDTM, TRTSDTF, TRTSTMF)
#'
#' @caption Further examples
#' @info Further example usages of this function can be found in the
#'   `vignette("generic")`.
derive_vars_merged <- function(dataset,
                               dataset_add,
                               by_vars,
                               order = NULL,
                               new_vars = NULL,
                               filter_add = NULL,
                               mode = NULL,
                               exist_flag = NULL,
                               true_value = "Y",
                               false_value = NA_character_,
                               missing_values = NULL,
                               check_type = "warning",
                               duplicate_msg = NULL,
                               relationship = NULL) {
  filter_add <- assert_filter_cond(enexpr(filter_add), optional = TRUE)
  assert_vars(by_vars)
  by_vars_left <- replace_values_by_names(by_vars)
  by_vars_right <- chr2vars(paste(vars2chr(by_vars)))
  assert_expr_list(order, optional = TRUE)
  assert_expr_list(new_vars, optional = TRUE)
  assert_data_frame(dataset, required_vars = by_vars_left)
  assert_data_frame(
    dataset_add,
    required_vars = expr_c(
      by_vars_right,
      setdiff(extract_vars(order), replace_values_by_names(new_vars)),
      extract_vars(new_vars)
    )
  )

  exist_flag <- assert_symbol(enexpr(exist_flag), optional = TRUE)
  assert_atomic_vector(true_value, optional = TRUE)
  assert_atomic_vector(false_value, optional = TRUE)
  assert_expr_list(missing_values, named = TRUE, optional = TRUE)
  if (!is.null(missing_values) && !is.null(new_vars)) {
    invalid_vars <- setdiff(
      names(missing_values),
      vars2chr(replace_values_by_names(new_vars))
    )
    if (length(invalid_vars) > 0) {
      cli_abort(paste(
        "The variable{?s} {.var {invalid_vars}} w{?as/ere} specified for {.arg missing_values}",
        "but not for {.arg new_vars}."
      ))
    }
  }
  relationship <- assert_character_scalar(
    relationship,
    values = c("one-to-one", "many-to-one"),
    case_sensitive = TRUE,
    optional = TRUE
  )


  add_data <- dataset_add %>%
    mutate(!!!new_vars) %>%
    filter_if(filter_add)

  if (!is.null(order)) {
    add_data <- filter_extreme(
      add_data,
      by_vars = by_vars_right,
      order = order,
      mode = mode,
      check_type = check_type
    )
  } else {
    if (is.null(duplicate_msg)) {
      duplicate_msg <- paste(
        "Dataset {.arg dataset_add} contains duplicate records with respect to",
        "{.var {vars2chr(by_vars)}}."
      )
    }
    signal_duplicate_records(
      add_data,
      by_vars = by_vars_right,
      msg = duplicate_msg
    )
  }
  if (!is.null(new_vars)) {
    add_data <- add_data %>%
      select(!!!by_vars_right, !!!replace_values_by_names(new_vars))
  }

  if (!is.null(missing_values) || !is.null(exist_flag)) {
    match_flag_var <- get_new_tmp_var(add_data, prefix = "tmp_match_flag")
  } else {
    match_flag_var <- NULL
  }
  if (!is.null(match_flag_var)) {
    add_data <- mutate(
      add_data,
      !!match_flag_var := TRUE
    )
  }

  # check if there are any variables in both datasets which are not by vars
  # in this case an error is issued to avoid renaming of varibles by left_join()
  common_vars <-
    setdiff(intersect(names(dataset), names(add_data)), vars2chr(by_vars))
  if (length(common_vars) > 0L) {
    cli_abort(
      c(
        "The variable{?s} {.var {common_vars}} {?is/are} contained in both datasets.",
        i = "Please add them to {.arg by_vars} or remove or rename them in one of the datasets."
      )
    )
  }

  tryCatch(
    dataset <- left_join(
      dataset,
      add_data,
      by = vars2chr(by_vars),
      relationship = relationship
    ),
    "dplyr_error_join_relationship_one_to_one" = function(cnd) {
      cli_abort(
        message = c(
          str_replace(
            str_replace(
              cnd$message, "`x`", "`dataset`"
            ), "`y`", "`dataset_add`"
          ),
          i = str_replace(
            str_replace(
              cnd$body, "`x`", "`dataset`"
            ), "`y`", "`dataset_add`"
          )
        ),
        call = parent.frame(n = 4)
      )
    }
  )

  if (!is.null(match_flag_var)) {
    update_missings <- map2(
      syms(names(missing_values)),
      missing_values,
      ~ expr(if_else(is.na(!!match_flag_var), !!.y, !!.x))
    )
    names(update_missings) <- names(missing_values)
    dataset <- dataset %>%
      mutate(!!!update_missings)
  }

  if (!is.null(exist_flag)) {
    dataset <- dataset %>%
      mutate(!!exist_flag := if_else(is.na(!!match_flag_var), false_value, true_value))
  }

  dataset %>%
    remove_tmp_vars()
}


#' Merge an Existence Flag
#'
#' @description Adds a flag variable to the input dataset which indicates if
#' there exists at least one observation in another dataset fulfilling a certain
#' condition.
#'
#' **Note:** This is a wrapper function for the more generic `derive_vars_merged()`.
#'
#' @param dataset_add Additional dataset
#'
#'   The variables specified by the `by_vars` argument are expected.
#'
#' @permitted [dataset]
#'
#' @param by_vars Grouping variables
#'
#' `r roxygen_param_by_vars()`
#'
#' @permitted [var_list]
#'
#' @param new_var New variable
#'
#'   The specified variable is added to the input dataset.
#'
#' @permitted [var]
#'
#' @param condition Condition
#'
#'   The condition is evaluated at the additional dataset (`dataset_add`). For
#'   all by groups where it evaluates as `TRUE` at least once the new variable
#'   is set to the true value (`true_value`). For all by groups where it
#'   evaluates as `FALSE` or `NA` for all observations the new variable is set
#'   to the false value (`false_value`). The new variable is set to the missing
#'   value (`missing_value`) for by groups not present in the additional
#'   dataset.
#'
#' @permitted [condition]
#'
#' @param true_value True value
#'
#' @permitted [char_scalar]
#'
#' @param false_value False value
#'
#' @permitted [char_scalar]
#'
#' @param missing_value Value used for missing information
#'
#'   The new variable is set to the specified value for all by groups without
#'   observations in the additional dataset.
#'
#' @permitted [char_scalar]
#'
#' @param filter_add Filter for additional data
#'
#'   Only observations fulfilling the specified condition are taken into account
#'   for flagging. If the argument is not specified, all observations are
#'   considered.
#'
#' @permitted [condition]
#'
#' @inheritParams derive_vars_merged
#'
#' @return The output dataset contains all observations and variables of the
#'   input dataset and additionally the variable specified for `new_var` derived
#'   from the additional dataset (`dataset_add`).
#'
#' @details
#'
#'   1. The additional dataset is restricted to the observations matching the
#'   `filter_add` condition.
#'
#'   1. The new variable is added to the input dataset and set to the true value
#'   (`true_value`) if for the by group at least one observation exists in the
#'   (restricted) additional dataset where the condition evaluates to `TRUE`. It
#'   is set to the false value (`false_value`) if for the by group at least one
#'   observation exists and for all observations the condition evaluates to
#'   `FALSE` or `NA`. Otherwise, it is set to the missing value
#'   (`missing_value`).
#'
#' @family der_gen
#' @keywords der_gen
#'
#' @export
#'
#' @examples
#'
#' library(dplyr, warn.conflicts = FALSE)
#'
#' dm <- tribble(
#'   ~STUDYID,  ~DOMAIN,  ~USUBJID, ~AGE,   ~AGEU,
#'   "PILOT01",    "DM", "01-1028",   71, "YEARS",
#'   "PILOT01",    "DM", "04-1127",   84, "YEARS",
#'   "PILOT01",    "DM", "06-1049",   60, "YEARS"
#' )
#'
#' ae <- tribble(
#'   ~STUDYID,  ~DOMAIN,  ~USUBJID,    ~AETERM,     ~AEREL,
#'   "PILOT01",    "AE", "01-1028", "ERYTHEMA", "POSSIBLE",
#'   "PILOT01",    "AE", "01-1028", "PRURITUS", "PROBABLE",
#'   "PILOT01",    "AE", "06-1049",  "SYNCOPE", "POSSIBLE",
#'   "PILOT01",    "AE", "06-1049",  "SYNCOPE", "PROBABLE"
#' )
#'
#'
#' derive_var_merged_exist_flag(
#'   dm,
#'   dataset_add = ae,
#'   by_vars = exprs(STUDYID, USUBJID),
#'   new_var = AERELFL,
#'   condition = AEREL == "PROBABLE"
#' ) %>%
#'   select(STUDYID, USUBJID, AGE, AGEU, AERELFL)
#'
#' vs <- tribble(
#'   ~STUDYID,  ~DOMAIN,  ~USUBJID,      ~VISIT, ~VSTESTCD, ~VSSTRESN, ~VSBLFL,
#'   "PILOT01",    "VS", "01-1028", "SCREENING",  "HEIGHT",     177.8,      NA,
#'   "PILOT01",    "VS", "01-1028", "SCREENING",  "WEIGHT",     98.88,      NA,
#'   "PILOT01",    "VS", "01-1028",  "BASELINE",  "WEIGHT",     99.34,     "Y",
#'   "PILOT01",    "VS", "01-1028",    "WEEK 4",  "WEIGHT",     98.88,      NA,
#'   "PILOT01",    "VS", "04-1127", "SCREENING",  "HEIGHT",     165.1,      NA,
#'   "PILOT01",    "VS", "04-1127", "SCREENING",  "WEIGHT",     42.87,      NA,
#'   "PILOT01",    "VS", "04-1127",  "BASELINE",  "WEIGHT",     41.05,     "Y",
#'   "PILOT01",    "VS", "04-1127",    "WEEK 4",  "WEIGHT",     41.73,      NA,
#'   "PILOT01",    "VS", "06-1049", "SCREENING",  "HEIGHT",    167.64,      NA,
#'   "PILOT01",    "VS", "06-1049", "SCREENING",  "WEIGHT",     57.61,      NA,
#'   "PILOT01",    "VS", "06-1049",  "BASELINE",  "WEIGHT",     57.83,     "Y",
#'   "PILOT01",    "VS", "06-1049",    "WEEK 4",  "WEIGHT",     58.97,      NA
#' )
#' derive_var_merged_exist_flag(
#'   dm,
#'   dataset_add = vs,
#'   by_vars = exprs(STUDYID, USUBJID),
#'   filter_add = VSTESTCD == "WEIGHT" & VSBLFL == "Y",
#'   new_var = WTBLHIFL,
#'   condition = VSSTRESN > 90,
#'   false_value = "N",
#'   missing_value = "M"
#' ) %>%
#'   select(STUDYID, USUBJID, AGE, AGEU, WTBLHIFL)
derive_var_merged_exist_flag <- function(dataset,
                                         dataset_add,
                                         by_vars,
                                         new_var,
                                         condition,
                                         true_value = "Y",
                                         false_value = NA_character_,
                                         missing_value = NA_character_,
                                         filter_add = NULL) {
  condition <- assert_filter_cond(enexpr(condition))
  new_var <- assert_symbol(enexpr(new_var))
  filter_add <- assert_filter_cond(enexpr(filter_add), optional = TRUE)
  add_data <- get_flagged_records(dataset_add,
    new_var = !!new_var,
    condition = !!condition,
    !!filter_add
  )

  derive_vars_merged(
    dataset,
    dataset_add = add_data,
    by_vars = by_vars,
    new_vars = exprs(!!new_var),
    order = exprs(!!new_var),
    check_type = "none",
    mode = "last"
  ) %>%
    mutate(!!new_var := if_else(!!new_var == 1, true_value, false_value, missing_value))
}

#' Merge Lookup Table with Source Dataset
#'
#' Merge user-defined lookup table with the input dataset. Optionally print a
#' list of records from the input dataset that do not have corresponding
#' mapping from the lookup table.
#'
#' @param dataset_add Lookup table
#'
#'   The variables specified by the `by_vars` argument are expected.
#'
#' @permitted [dataset]
#'
#' @param print_not_mapped Print a list of unique `by_vars` values that do not
#' have corresponding records from the lookup table?
#'
#' @permitted [boolean]
#'
#' @inheritParams derive_vars_merged
#'
#' @return The output dataset contains all observations and variables of the
#' input dataset, and add the variables specified in `new_vars` from the lookup
#' table specified in `dataset_add`. Optionally prints a list of unique
#' `by_vars` values that do not have corresponding records
#' from the lookup table (by specifying `print_not_mapped = TRUE`).
#'
#' @keywords der_gen
#' @family der_gen
#'
#' @export
#'
#' @examples
#' library(dplyr, warn.conflicts = FALSE)
#' vs <- tribble(
#'   ~STUDYID,  ~DOMAIN,  ~USUBJID,        ~VISIT, ~VSTESTCD,       ~VSTEST,
#'   "PILOT01",    "VS", "01-1028",   "SCREENING",  "HEIGHT",      "Height",
#'   "PILOT01",    "VS", "01-1028",   "SCREENING",    "TEMP", "Temperature",
#'   "PILOT01",    "VS", "01-1028",    "BASELINE",    "TEMP", "Temperature",
#'   "PILOT01",    "VS", "01-1028",      "WEEK 4",    "TEMP", "Temperature",
#'   "PILOT01",    "VS", "01-1028", "SCREENING 1",  "WEIGHT",      "Weight",
#'   "PILOT01",    "VS", "01-1028",    "BASELINE",  "WEIGHT",      "Weight",
#'   "PILOT01",    "VS", "01-1028",      "WEEK 4",  "WEIGHT",      "Weight",
#'   "PILOT01",    "VS", "04-1325",   "SCREENING",  "HEIGHT",      "Height",
#'   "PILOT01",    "VS", "04-1325",   "SCREENING",    "TEMP", "Temperature",
#'   "PILOT01",    "VS", "04-1325",    "BASELINE",    "TEMP", "Temperature",
#'   "PILOT01",    "VS", "04-1325",      "WEEK 4",    "TEMP", "Temperature",
#'   "PILOT01",    "VS", "04-1325", "SCREENING 1",  "WEIGHT",      "Weight",
#'   "PILOT01",    "VS", "04-1325",    "BASELINE",  "WEIGHT",      "Weight",
#'   "PILOT01",    "VS", "04-1325",      "WEEK 4",  "WEIGHT",      "Weight",
#'   "PILOT01",    "VS", "10-1027",   "SCREENING",  "HEIGHT",      "Height",
#'   "PILOT01",    "VS", "10-1027",   "SCREENING",    "TEMP", "Temperature",
#'   "PILOT01",    "VS", "10-1027",    "BASELINE",    "TEMP", "Temperature",
#'   "PILOT01",    "VS", "10-1027",      "WEEK 4",    "TEMP", "Temperature",
#'   "PILOT01",    "VS", "10-1027", "SCREENING 1",  "WEIGHT",      "Weight",
#'   "PILOT01",    "VS", "10-1027",    "BASELINE",  "WEIGHT",      "Weight",
#'   "PILOT01",    "VS", "10-1027",      "WEEK 4",  "WEIGHT",      "Weight"
#' )
#'
#' param_lookup <- tribble(
#'   ~VSTESTCD,                 ~VSTEST, ~PARAMCD,                       ~PARAM,
#'   "SYSBP", "Systolic Blood Pressure",  "SYSBP", "Syst Blood Pressure (mmHg)",
#'   "WEIGHT",                 "Weight", "WEIGHT",                "Weight (kg)",
#'   "HEIGHT",                 "Height", "HEIGHT",                "Height (cm)",
#'   "TEMP",              "Temperature",   "TEMP",            "Temperature (C)",
#'   "MAP",    "Mean Arterial Pressure",    "MAP",   "Mean Art Pressure (mmHg)",
#'   "BMI",           "Body Mass Index",    "BMI",    "Body Mass Index(kg/m^2)",
#'   "BSA",         "Body Surface Area",    "BSA",     "Body Surface Area(m^2)"
#' )
#'
#' derive_vars_merged_lookup(
#'   dataset = vs,
#'   dataset_add = param_lookup,
#'   by_vars = exprs(VSTESTCD),
#'   new_vars = exprs(PARAMCD, PARAM),
#'   print_not_mapped = TRUE
#' )
derive_vars_merged_lookup <- function(dataset,
                                      dataset_add,
                                      by_vars,
                                      order = NULL,
                                      new_vars = NULL,
                                      mode = NULL,
                                      filter_add = NULL,
                                      check_type = "warning",
                                      duplicate_msg = NULL,
                                      print_not_mapped = TRUE) {
  by_vars_left <- replace_values_by_names(by_vars)
  assert_logical_scalar(print_not_mapped)
  filter_add <- assert_filter_cond(enexpr(filter_add), optional = TRUE)

  tmp_lookup_flag <- get_new_tmp_var(dataset_add, prefix = "tmp_lookup_flag")

  res <- derive_vars_merged(
    dataset,
    dataset_add,
    by_vars = by_vars,
    order = order,
    new_vars = new_vars,
    mode = mode,
    filter_add = !!filter_add,
    exist_flag = !!tmp_lookup_flag,
    check_type = check_type,
    duplicate_msg = duplicate_msg
  )

  if (print_not_mapped) {
    temp_not_mapped <- res %>%
      filter(is.na(!!tmp_lookup_flag)) %>%
      distinct(!!!by_vars_left)

    if (nrow(temp_not_mapped) > 0) {
      # nolint start: undesirable_function_linter
      admiral_environment$nmap <- structure(
        temp_not_mapped,
        class = union("nmap", class(temp_not_mapped)),
        by_vars = vars2chr(by_vars_left)
      )
      # nolint end

      cli_inform(
        c("List of {.var {vars2chr(by_vars_left)}} not mapped:",
          capture.output(temp_not_mapped),
          i = "Run {.run admiral::get_not_mapped()} to access the full list."
        )
      )
    } else if (nrow(temp_not_mapped) == 0) {
      cli_inform(
        "All {.var {vars2chr(by_vars_left)}} are mapped."
      )
    }
  }

  res %>% remove_tmp_vars()
}

#' Get list of records not mapped from the lookup table.
#'
#' @export
#'
#' @return A `data.frame` or `NULL`
#'
#' @keywords utils_help
#' @family utils_help
get_not_mapped <- function() {
  admiral_environment$nmap
}

#' Merge Summary Variables
#'
#' @description Merge a summary variable from a dataset to the input dataset.
#'
#' @param dataset
#'
#' `r roxygen_param_dataset(expected_vars = c("by_vars"))`
#'
#' @permitted [dataset]
#'
#' @param dataset_add Additional dataset
#'
#'   The variables specified by the `by_vars` and the variables used on the left
#'   hand sides of the `new_vars` arguments are expected.
#'
#' @permitted [dataset]
#'
#' @param by_vars Grouping variables
#'
#'   The expressions on the left hand sides of `new_vars` are evaluated by the
#'   specified *variables*. Then the resulting values are merged to the input
#'   dataset (`dataset`) by the specified *variables*.
#'
#'   `r roxygen_param_by_vars()`
#'
#' @permitted [var_list]
#'
#' @param new_vars New variables to add
#'
#'   The specified variables are added to the input dataset.
#'
#'   A named list of expressions is expected:
#'   + LHS refer to a variable.
#'   + RHS refers to the values to set to the variable. This can be a string, a
#'   symbol, a numeric value, an expression or NA. If summary functions are
#'   used, the values are summarized by the variables specified for `by_vars`.
#'
#'   For example:
#'   ```
#'     new_vars = exprs(
#'       DOSESUM = sum(AVAL),
#'       DOSEMEAN = mean(AVAL)
#'     )
#'   ```
#'
#' @permitted [var_list]
#'
#' @param filter_add Filter for additional dataset (`dataset_add`)
#'
#'   Only observations fulfilling the specified condition are taken into account
#'   for summarizing. If the argument is not specified, all observations are
#'   considered.
#'
#' @permitted [condition]
#'
#' @inheritParams derive_vars_merged
#'
#' @details
#'
#'   1. The records from the additional dataset (`dataset_add`) are restricted
#'   to those matching the `filter_add` condition.
#'
#'   1. The new variables (`new_vars`) are created for each by group (`by_vars`)
#'   in the additional dataset (`dataset_add`) by calling `summarize()`. I.e.,
#'   all observations of a by group are summarized to a single observation.
#'
#'   1. The new variables are merged to the input dataset. For observations
#'   without a matching observation in the additional dataset the new variables
#'   are set to `NA`. Observations in the additional dataset which have no
#'   matching observation in the input dataset are ignored.
#'
#' @return The output dataset contains all observations and variables of the
#'   input dataset and additionally the variables specified for `new_vars`.
#'
#' @family der_gen
#' @keywords der_gen
#'
#' @seealso [derive_summary_records()], [get_summary_records()]
#'
#' @export
#'
#' @examples
#' library(tibble)
#'
#' # Add a variable for the mean of AVAL within each visit
#' adbds <- tribble(
#'   ~USUBJID,  ~AVISIT,  ~ASEQ, ~AVAL,
#'   "1",      "WEEK 1",      1,    10,
#'   "1",      "WEEK 1",      2,    NA,
#'   "1",      "WEEK 2",      3,    NA,
#'   "1",      "WEEK 3",      4,    42,
#'   "1",      "WEEK 4",      5,    12,
#'   "1",      "WEEK 4",      6,    12,
#'   "1",      "WEEK 4",      7,    15,
#'   "2",      "WEEK 1",      1,    21,
#'   "2",      "WEEK 4",      2,    22
#' )
#'
#' derive_vars_merged_summary(
#'   adbds,
#'   dataset_add = adbds,
#'   by_vars = exprs(USUBJID, AVISIT),
#'   new_vars = exprs(
#'     MEANVIS = mean(AVAL, na.rm = TRUE),
#'     MAXVIS = max(AVAL, na.rm = TRUE)
#'   )
#' )
#'
#' # Add a variable listing the lesion ids at baseline
#' adsl <- tribble(
#'   ~USUBJID,
#'   "1",
#'   "2",
#'   "3"
#' )
#'
#' adtr <- tribble(
#'   ~USUBJID,     ~AVISIT, ~LESIONID,
#'   "1",       "BASELINE",  "INV-T1",
#'   "1",       "BASELINE",  "INV-T2",
#'   "1",       "BASELINE",  "INV-T3",
#'   "1",       "BASELINE",  "INV-T4",
#'   "1",         "WEEK 1",  "INV-T1",
#'   "1",         "WEEK 1",  "INV-T2",
#'   "1",         "WEEK 1",  "INV-T4",
#'   "2",       "BASELINE",  "INV-T1",
#'   "2",       "BASELINE",  "INV-T2",
#'   "2",       "BASELINE",  "INV-T3",
#'   "2",         "WEEK 1",  "INV-T1",
#'   "2",         "WEEK 1",  "INV-N1"
#' )
#'
#' derive_vars_merged_summary(
#'   adsl,
#'   dataset_add = adtr,
#'   by_vars = exprs(USUBJID),
#'   filter_add = AVISIT == "BASELINE",
#'   new_vars = exprs(LESIONSBL = paste(LESIONID, collapse = ", "))
#' )
#'
derive_vars_merged_summary <- function(dataset,
                                       dataset_add,
                                       by_vars,
                                       new_vars = NULL,
                                       filter_add = NULL,
                                       missing_values = NULL) {
  assert_vars(by_vars)
  by_vars_left <- replace_values_by_names(by_vars)
  by_vars_right <- chr2vars(paste(vars2chr(by_vars)))
  # once new_var is removed new_vars should be mandatory
  assert_expr_list(new_vars, named = TRUE, optional = TRUE)
  filter_add <-
    assert_filter_cond(enexpr(filter_add), optional = TRUE)
  assert_data_frame(
    dataset,
    required_vars = by_vars_left
  )
  assert_data_frame(
    dataset_add,
    required_vars = expr_c(by_vars_right, extract_vars(new_vars))
  )

  # Summarise the analysis value and merge to the original dataset
  # If for one of the new variables no summary function is used, i.e., more than
  # one record is created per by group, the error from signal_duplicates_records()
  # need to be updated and the warning from dplyr needs to be suppressed as it
  # is misleading.
  tryCatch(
    derive_vars_merged(
      dataset,
      dataset_add = derive_summary_records(
        dataset_add = dataset_add,
        by_vars = by_vars_right,
        filter_add = !!filter_add,
        set_values_to = new_vars,
      ) %>%
        select(!!!by_vars_right, names(new_vars)),
      by_vars = by_vars,
      missing_values = missing_values
    ),
    multiple_summary_records = function(cnd) {
      cli_abort(
        c(
          paste(
            "After summarising, the dataset contains multiple records with",
            "respect to {.var {cnd$by_vars}}."
          ),
          paste(
            "Please check the {.arg new_vars} argument if summary functions",
            "like {.fun mean}, {.fun sum}, ... are used on the right hand side."
          ),
          i = "Run {.run admiral::get_duplicates_dataset()} to access the duplicate records"
        ),
        class = cnd$class,
        call = NULL,
        by_vars = cnd$by_vars
      )
    }
  )
}

#' Merge Summary Variables
#'
#' @description
#' `r lifecycle::badge("deprecated")` The `derive_var_merged_summary()`
#' function has been deprecated in favor of `derive_vars_merged_summary()`.
#'
#' @param dataset
#'
#' `r roxygen_param_dataset(expected_vars = c("by_vars"))`
#'
#' @permitted [dataset]
#'
#' @param dataset_add Additional dataset
#'
#'   The variables specified by the `by_vars` and the variables used on the left
#'   hand sides of the `new_vars` arguments are expected.
#'
#' @permitted [dataset]
#'
#' @param by_vars Grouping variables
#'
#'   The expressions on the left hand sides of `new_vars` are evaluated by the
#'   specified *variables*. Then the resulting values are merged to the input
#'   dataset (`dataset`) by the specified *variables*.
#'
#'   `r roxygen_param_by_vars()`
#'
#' @permitted [var_list]
#'
#' @param new_vars New variables to add
#'
#'   The specified variables are added to the input dataset.
#'
#'   A named list of expressions is expected:
#'   + LHS refer to a variable.
#'   + RHS refers to the values to set to the variable. This can be a string, a
#'   symbol, a numeric value, an expression or NA. If summary functions are
#'   used, the values are summarized by the variables specified for `by_vars`.
#'   Any expression on the RHS must result in a single value per by group.
#'
#'   For example:
#'   ```
#'     new_vars = exprs(
#'       DOSESUM = sum(AVAL),
#'       DOSEMEAN = mean(AVAL)
#'     )
#'   ```
#'
#' @permitted [var_list]
#'
#' @param filter_add Filter for additional dataset (`dataset_add`)
#'
#'   Only observations fulfilling the specified condition are taken into account
#'   for summarizing. If the argument is not specified, all observations are
#'   considered.
#'
#' @permitted [condition]
#'
#' @inheritParams derive_vars_merged
#'
#' @details
#'
#'   1. The records from the additional dataset (`dataset_add`) are restricted
#'   to those matching the `filter_add` condition.
#'
#'   1. The new variables (`new_vars`) are created for each by group (`by_vars`)
#'   in the additional dataset (`dataset_add`) by calling `summarize()`. I.e.,
#'   all observations of a by group are summarized to a single observation.
#'
#'   1. The new variables are merged to the input dataset. For observations
#'   without a matching observation in the additional dataset the new variables
#'   are set to `NA`. Observations in the additional dataset which have no
#'   matching observation in the input dataset are ignored.
#'
#' @return The output dataset contains all observations and variables of the
#'   input dataset and additionally the variables specified for `new_vars`.
#'
#' @family deprecated
#' @keywords deprecated
#'
#' @seealso [derive_summary_records()], [get_summary_records()]
#'
#' @export
#'
derive_var_merged_summary <- function(dataset,
                                      dataset_add,
                                      by_vars,
                                      new_vars = NULL,
                                      filter_add = NULL,
                                      missing_values = NULL) {
  deprecate_inform(
    when = "1.4",
    what = "derive_var_merged_summary()",
    with = "derive_vars_merged_summary()",
    details = c(
      x = "Function is brought inline with our programming strategy - warning
      will be issued in January 2027",
      i = "See admiral's deprecation guidance:
      https://pharmaverse.github.io/admiraldev/dev/articles/programming_strategy.html#deprecation"
    )
  )

  # Function code now replaced with call to new function derive_vars_merged_summary(),
  # minimal checks are executed prior to call to ensure all required arguments are present

  assert_data_frame(dataset)
  assert_data_frame(dataset_add)
  assert_vars(by_vars)

  derive_vars_merged_summary(
    dataset = dataset,
    dataset_add = dataset_add,
    by_vars = by_vars,
    new_vars = new_vars,
    filter_add = !!enexpr(filter_add),
    missing_values = missing_values
  )
}
