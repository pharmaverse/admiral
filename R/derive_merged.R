#' Add New Variable(s) to the Input Dataset Based on Variables from Another
#' Dataset
#'
#' Add new variable(s) to the input dataset based on variables from another
#' dataset. The observations to merge can be selected by a condition
#' (`filter_add` argument) and/or selecting the first or last observation for
#' each by group (`order` and `mode` argument).
#'
#' @param dataset
#' `r roxygen_param_dataset(expected_vars = c("by_vars"))`
#'
#' @param dataset_add Additional dataset
#'
#'   The variables specified by the `by_vars`, the `new_vars`, and the `order`
#'   argument are expected.
#'
#' @param by_vars Grouping variables
#'
#'   The input dataset and the selected observations from the additional dataset
#'   are merged by the specified variables.
#'
#'   `r roxygen_param_by_vars(rename = TRUE)`
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
#'   *Permitted Values*: list of expressions created by `exprs()`, e.g.,
#'   `exprs(ADT, desc(AVAL))` or `NULL`
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
#'   additional dataset (`dataset_add`) are added.
#'
#'   *Permitted Values*: list of variables or named expressions created by `exprs()`
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
#'   *Permitted Values*: a condition
#'
#' @param mode Selection mode
#'
#'   Determines if the first or last observation is selected. If the `order`
#'   argument is specified, `mode` must be non-null.
#'
#'   If the `order` argument is not specified, the `mode` argument is ignored.
#'
#'   *Permitted Values*: `"first"`, `"last"`, `NULL`
#'
#' @param match_flag Match flag
#'
#'  `r lifecycle::badge("deprecated")` Please use `exist_flag` instead.
#'
#'   If the argument is specified (e.g., `match_flag = FLAG`), the specified
#'   variable (e.g., `FLAG`) is added to the input dataset. This variable will
#'   be `TRUE` for all selected records from `dataset_add` which are merged into
#'   the input dataset, and `NA` otherwise.
#'
#'   *Permitted Values*: Variable name
#'
#' @param exist_flag Exist flag
#'
#'   If the argument is specified (e.g., `exist_flag = FLAG`), the specified
#'   variable (e.g., `FLAG`) is added to the input dataset. This variable will
#'   be the value provided in `true_value` for all selected records from `dataset_add`
#'   which are merged into the input dataset, and the value provided in `false_value` otherwise.
#'
#'   *Permitted Values*: Variable name
#'
#' @param true_value True value
#'
#'   The value for the specified variable `exist_flag`, applicable to
#'   the first or last observation (depending on the mode) of each by group.
#'
#'   Permitted Values: An atomic scalar
#'
#' @param false_value False value
#'
#'   The value for the specified variable `exist_flag`, NOT applicable to
#'   the first or last observation (depending on the mode) of each by group.
#'
#'   Permitted Values: An atomic scalar
#'
#' @param missing_values Values for non-matching observations
#'
#'   For observations of the input dataset (`dataset`) which do not have a
#'   matching observation in the additional dataset (`dataset_add`) the values
#'   of the specified variables are set to the specified value. Only variables
#'   specified for `new_vars` can be specified for `missing_values`.
#'
#'   *Permitted Values*: named list of expressions, e.g.,
#'   `exprs(BASEC = "MISSING", BASE = -1)`
#'
#' @param check_type Check uniqueness?
#'
#'   If `"warning"` or `"error"` is specified, the specified message is issued
#'   if the observations of the (restricted) additional dataset are not unique
#'   with respect to the by variables and the order.
#'
#'   If the `order` argument is not specified, the `check_type` argument is ignored:
#'    if the observations of the (restricted) additional dataset are not unique with respect
#'    to the by variables, an error is issued.
#'
#'   *Permitted Values*: `"none"`, `"warning"`, `"error"`
#'
#' @param duplicate_msg Message of unique check
#'
#'   If the uniqueness check fails, the specified message is displayed.
#'
#'   *Default*:
#'
#'   ```{r echo=TRUE, eval=FALSE}
#'   paste(
#'     "Dataset {.arg dataset_add} contains duplicate records with respect to",
#'     "{.var {vars2chr(by_vars)}}."
#'   )
#'   ```
#'
#' @param relationship Expected merge-relationship between the `by_vars`
#'   variable(s) in `dataset` (input dataset) and the `dataset_add` (additional dataset)
#'    containing the additional `new_vars`.
#'
#'   This argument is passed to the `dplyr::left_join()` function. See
#'   <https://dplyr.tidyverse.org/reference/mutate-joins.html#arguments> for
#'   more details.
#'
#'   **Permitted Values:** `"one-to-one"`, `"many-to-one"`, `NULL`.
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
#' @examples
#' library(dplyr, warn.conflicts = FALSE)
#' vs <- tribble(
#'   ~STUDYID,  ~DOMAIN,  ~USUBJID, ~VSTESTCD,      ~VISIT, ~VSSTRESN, ~VSSTRESU,       ~VSDTC,
#'   "PILOT01",    "VS", "01-1302",  "HEIGHT", "SCREENING",     177.8,      "cm", "2013-08-20",
#'   "PILOT01",    "VS", "01-1302",  "WEIGHT", "SCREENING",     81.19,      "kg", "2013-08-20",
#'   "PILOT01",    "VS", "01-1302",  "WEIGHT",  "BASELINE",      82.1,      "kg", "2013-08-29",
#'   "PILOT01",    "VS", "01-1302",  "WEIGHT",    "WEEK 2",     81.19,      "kg", "2013-09-15",
#'   "PILOT01",    "VS", "01-1302",  "WEIGHT",    "WEEK 4",     82.56,      "kg", "2013-09-24",
#'   "PILOT01",    "VS", "01-1302",  "WEIGHT",    "WEEK 6",     80.74,      "kg", "2013-10-08",
#'   "PILOT01",    "VS", "01-1302",  "WEIGHT",    "WEEK 8",      82.1,      "kg", "2013-10-22",
#'   "PILOT01",    "VS", "01-1302",  "WEIGHT",   "WEEK 12",      82.1,      "kg", "2013-11-05",
#'   "PILOT01",    "VS", "17-1344",  "HEIGHT", "SCREENING",     163.5,      "cm", "2014-01-01",
#'   "PILOT01",    "VS", "17-1344",  "WEIGHT", "SCREENING",     58.06,      "kg", "2014-01-01",
#'   "PILOT01",    "VS", "17-1344",  "WEIGHT",  "BASELINE",     58.06,      "kg", "2014-01-11",
#'   "PILOT01",    "VS", "17-1344",  "WEIGHT",    "WEEK 2",     58.97,      "kg", "2014-01-24",
#'   "PILOT01",    "VS", "17-1344",  "WEIGHT",    "WEEK 4",     57.97,      "kg", "2014-02-07",
#'   "PILOT01",    "VS", "17-1344",  "WEIGHT",    "WEEK 6",     58.97,      "kg", "2014-02-19",
#'   "PILOT01",    "VS", "17-1344",  "WEIGHT",    "WEEK 8",     57.79,      "kg", "2014-03-14"
#' )
#'
#' dm <- tribble(
#'   ~STUDYID,  ~DOMAIN,  ~USUBJID, ~AGE,   ~AGEU,
#'   "PILOT01",    "DM", "01-1302",   61, "YEARS",
#'   "PILOT01",    "DM", "17-1344",   64, "YEARS"
#' )
#'
#'
#' # Merging all dm variables to vs
#' derive_vars_merged(
#'   vs,
#'   dataset_add = select(dm, -DOMAIN),
#'   by_vars = exprs(STUDYID, USUBJID)
#' ) %>%
#'   select(STUDYID, USUBJID, VSTESTCD, VISIT, VSSTRESN, AGE, AGEU)
#'
#'
#' # Merge last weight to adsl
#' adsl <- tribble(
#'   ~STUDYID,   ~USUBJID, ~AGE,   ~AGEU,
#'   "PILOT01", "01-1302",   61, "YEARS",
#'   "PILOT01", "17-1344",   64, "YEARS"
#' )
#'
#'
#' derive_vars_merged(
#'   adsl,
#'   dataset_add = vs,
#'   by_vars = exprs(STUDYID, USUBJID),
#'   order = exprs(convert_dtc_to_dtm(VSDTC)),
#'   mode = "last",
#'   new_vars = exprs(LASTWGT = VSSTRESN, LASTWGTU = VSSTRESU),
#'   filter_add = VSTESTCD == "WEIGHT",
#'   exist_flag = vsdatafl
#' ) %>%
#'   select(STUDYID, USUBJID, AGE, AGEU, LASTWGT, LASTWGTU, vsdatafl)
#'
#'
#' # Derive treatment start datetime (TRTSDTM)
#' ex <- tribble(
#'   ~STUDYID,  ~DOMAIN,  ~USUBJID, ~EXSTDY, ~EXENDY,     ~EXSTDTC,     ~EXENDTC,
#'   "PILOT01",    "EX", "01-1302",       1,      18, "2013-08-29", "2013-09-15",
#'   "PILOT01",    "EX", "01-1302",      19,      69, "2013-09-16", "2013-11-05",
#'   "PILOT01",    "EX", "17-1344",       1,      14, "2014-01-11", "2014-01-24",
#'   "PILOT01",    "EX", "17-1344",      15,      63, "2014-01-25", "2014-03-14"
#' )
#' ## Impute exposure start date to first date/time
#' ex_ext <- derive_vars_dtm(
#'   ex,
#'   dtc = EXSTDTC,
#'   new_vars_prefix = "EXST",
#'   highest_imputation = "M",
#' )
#' ## Add first exposure datetime and imputation flags to adsl
#' derive_vars_merged(
#'   select(dm, STUDYID, USUBJID),
#'   dataset_add = ex_ext,
#'   by_vars = exprs(STUDYID, USUBJID),
#'   new_vars = exprs(TRTSDTM = EXSTDTM, TRTSDTF = EXSTDTF, TRTSTMF = EXSTTMF),
#'   order = exprs(EXSTDTM),
#'   mode = "first"
#' )
#'
#' # Derive treatment end datetime (TRTEDTM)
#' ## Impute exposure end datetime to last time, no date imputation
#' ex_ext <- derive_vars_dtm(
#'   ex,
#'   dtc = EXENDTC,
#'   new_vars_prefix = "EXEN",
#'   time_imputation = "last",
#' )
#' ## Add last exposure datetime and imputation flag to adsl
#' derive_vars_merged(
#'   select(adsl, STUDYID, USUBJID),
#'   dataset_add = ex_ext,
#'   filter_add = !is.na(EXENDTM),
#'   by_vars = exprs(STUDYID, USUBJID),
#'   new_vars = exprs(TRTEDTM = EXENDTM, TRTETMF = EXENTMF),
#'   order = exprs(EXENDTM),
#'   mode = "last"
#' )
#' # Modify merged values and set value for non matching observations
#' adsl <- tribble(
#'   ~USUBJID, ~SEX, ~COUNTRY,
#'   "ST42-1", "F",  "AUT",
#'   "ST42-2", "M",  "MWI",
#'   "ST42-3", "M",  "NOR",
#'   "ST42-4", "F",  "UGA"
#' )
#'
#' advs <- tribble(
#'   ~USUBJID, ~PARAMCD, ~AVISIT,    ~AVISITN, ~AVAL,
#'   "ST42-1", "WEIGHT", "BASELINE",        0,    66,
#'   "ST42-1", "WEIGHT", "WEEK 2",          1,    68,
#'   "ST42-2", "WEIGHT", "BASELINE",        0,    88,
#'   "ST42-3", "WEIGHT", "WEEK 2",          1,    55,
#'   "ST42-3", "WEIGHT", "WEEK 4",          2,    50
#' )
#'
#' derive_vars_merged(
#'   adsl,
#'   dataset_add = advs,
#'   by_vars = exprs(USUBJID),
#'   new_vars = exprs(
#'     LSTVSCAT = if_else(AVISIT == "BASELINE", "BASELINE", "POST-BASELINE")
#'   ),
#'   order = exprs(AVISITN),
#'   mode = "last",
#'   missing_values = exprs(LSTVSCAT = "MISSING")
#' )
derive_vars_merged <- function(dataset,
                               dataset_add,
                               by_vars,
                               order = NULL,
                               new_vars = NULL,
                               filter_add = NULL,
                               mode = NULL,
                               match_flag,
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
  if (!is_missing(enexpr(match_flag))) {
    deprecate_stop(
      "1.1.0",
      "derive_vars_merged(match_flag =)",
      "derive_vars_merged(exist_flag =)"
    )
    exist_flag <- assert_symbol(enexpr(match_flag), optional = TRUE)
  }
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
    },
    "dplyr_error_join_relationship_many_to_one" = function(cnd) {
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
      mutate(!!exist_flag := ifelse(is.na(!!match_flag_var), false_value, true_value))
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
#' @param by_vars Grouping variables
#'
#' `r roxygen_param_by_vars()`
#'
#' @param new_var New variable
#'
#'   The specified variable is added to the input dataset.
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
#' @param true_value True value
#'
#' @param false_value False value
#'
#' @param missing_value Values used for missing information
#'
#'   The new variable is set to the specified value for all by groups without
#'   observations in the additional dataset.
#'
#'   *Permitted Value*: A character scalar
#'
#' @param filter_add Filter for additional data
#'
#'   Only observations fulfilling the specified condition are taken into account
#'   for flagging. If the argument is not specified, all observations are
#'   considered.
#'
#'   *Permitted Values*: a condition
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
#' The variables specified by the `by_vars` argument are expected.
#'
#' @param print_not_mapped Print a list of unique `by_vars` values that do not
#' have corresponding records from the lookup table?
#'
#' *Default*: `TRUE`
#'
#' *Permitted Values*: `TRUE`, `FALSE`
#'
#' @inheritParams derive_vars_merged
#'
#'
#' @return The output dataset contains all observations and variables of the
#' input dataset, and add the variables specified in `new_vars` from the lookup
#' table specified in `dataset_add`. Optionally prints a list of unique
#' `by_vars` values that do not have corresponding records
#' from the lookup table (by specifying `print_not_mapped = TRUE`).
#'
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
#' `r roxygen_param_dataset(expected_vars = c("by_vars"))`
#'
#' @param dataset_add Additional dataset
#'
#'   The variables specified by the `by_vars` and the variables used on the left
#'   hand sides of the `new_vars` arguments are expected.
#'
#' @param new_var Variable to add
#'
#'  `r lifecycle::badge("deprecated")` Please use `new_vars` instead.
#'
#'   The specified variable is added to the input dataset (`dataset`) and set to
#'   the summarized values.
#'
#' @param by_vars Grouping variables
#'
#'   The expressions on the left hand sides of `new_vars` are evaluated by the
#'   specified *variables*. Then the resulting values are merged to the input
#'   dataset (`dataset`) by the specified *variables*.
#'
#'   `r roxygen_param_by_vars()`
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
#' @param filter_add Filter for additional dataset (`dataset_add`)
#'
#'   Only observations fulfilling the specified condition are taken into account
#'   for summarizing. If the argument is not specified, all observations are
#'   considered.
#'
#'   *Permitted Values*: a condition
#'
#' @param analysis_var Analysis variable
#'
#'  `r lifecycle::badge("deprecated")` Please use `new_vars` instead.
#'
#'   The values of the specified variable are summarized by the function
#'   specified for `summary_fun`.
#'
#' @param summary_fun Summary function
#'
#'  `r lifecycle::badge("deprecated")` Please use `new_vars` instead.
#'
#'   The specified function that takes as input `analysis_var` and performs the
#'   calculation. This can include built-in functions as well as user defined
#'   functions, for example `mean` or `function(x) mean(x, na.rm = TRUE)`.
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
#' derive_var_merged_summary(
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
#' derive_var_merged_summary(
#'   adsl,
#'   dataset_add = adtr,
#'   by_vars = exprs(USUBJID),
#'   filter_add = AVISIT == "BASELINE",
#'   new_vars = exprs(LESIONSBL = paste(LESIONID, collapse = ", "))
#' )
#'
derive_var_merged_summary <- function(dataset,
                                      dataset_add,
                                      by_vars,
                                      new_vars = NULL,
                                      new_var,
                                      filter_add = NULL,
                                      missing_values = NULL,
                                      analysis_var,
                                      summary_fun) {
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

  if (!missing(new_var) || !missing(analysis_var) || !missing(summary_fun)) {
    deprecate_stop(
      "1.1.0",
      I("derive_var_merged_summary(new_var = , anaylsis_var = , summary_fun = )"),
      "derive_var_merged_summary(new_vars = )"
    )
    new_var <- assert_symbol(enexpr(new_var))
    analysis_var <- assert_symbol(enexpr(analysis_var))
    assert_s3_class(summary_fun, "function")
    new_vars <- exprs(!!new_var := {{ summary_fun }}(!!analysis_var), !!!new_vars)
  }

  # Summarise the analysis value and merge to the original dataset
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
  )
}
