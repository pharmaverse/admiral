#' Add New Variable(s) to the Input Dataset Based on Variables from Another
#' Dataset
#'
#' Add new variable(s) to the input dataset based on variables from another
#' dataset. The observations to merge can be selected by a condition
#' (`filter_add` argument) and/or selecting the first or last observation for
#' each by group (`order` and `mode` argument).
#'
#' @param dataset Input dataset
#'
#'   The variables specified by the `by_vars` argument are expected.
#'
#' @param dataset_add Additional dataset
#'
#'   The variables specified by the `by_vars`, the `new_vars`, and the `order`
#'   argument are expected.
#'
#' @param by_vars Grouping variables
#'
#'   The input dataset and the selected observations from the additional dataset
#'   are merged by the specified by variables. The by variables must be a unique
#'   key of the selected observations. Variables from the additional dataset can
#'   be renamed by naming the element, i.e., `by_vars =
#'   exprs(<name in input dataset> = <name in additional dataset>)`, similar to
#'   the dplyr joins.
#'
#'   *Permitted Values*: list of variables created by `exprs()`
#'
#' @param order Sort order
#'
#'   If the argument is set to a non-null value, for each by group the first or
#'   last observation from the additional dataset is selected with respect to the
#'   specified order.
#'
#'   Variables defined by the `new_vars` argument can be used in the sort order.
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
#'   If the argument is specified (e.g., `match_flag = FLAG`), the specified
#'   variable (e.g., `FLAG`) is added to the input dataset. This variable will
#'   be `TRUE` for all selected records from `dataset_add` which are merged into
#'   the input dataset, and `NA` otherwise.
#'
#'   *Permitted Values*: Variable name
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
#'   *Permitted Values*: `"none"`, `"warning"`, `"error"`
#'
#' @param duplicate_msg Message of unique check
#'
#'   If the uniqueness check fails, the specified message is displayed.
#'
#'   *Default*:
#'
#'   ```{r echo=TRUE, eval=FALSE}
#'   paste("Dataset `dataset_add` contains duplicate records with respect to",
#'         enumerate(vars2chr(by_vars)))
#'   ```
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
#'   match_flag = vsdatafl
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
                               match_flag = NULL,
                               missing_values = NULL,
                               check_type = "warning",
                               duplicate_msg = NULL) {
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
  match_flag <- assert_symbol(enexpr(match_flag), optional = TRUE)
  assert_expr_list(missing_values, named = TRUE, optional = TRUE)
  if (!is.null(missing_values)) {
    invalid_vars <- setdiff(
      names(missing_values),
      vars2chr(replace_values_by_names(new_vars))
    )
    if (length(invalid_vars) > 0) {
      abort(paste(
        "The variables",
        enumerate(invalid_vars),
        "were specified for `missing_values` but not for `new_vars`."
      ))
    }
  }

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
        "Dataset `dataset_add` contains duplicate records with respect to",
        enumerate(vars2chr(by_vars))
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

  if (!is.null(missing_values)) {
    match_flag_var <- get_new_tmp_var(add_data, prefix = "tmp_match_flag")
  } else {
    match_flag_var <- match_flag
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
    abort(if_else(
      length(common_vars) == 1L,
      paste0(
        "The variable ",
        common_vars[[1]],
        " is contained in both datasets.\n",
        "Please add it to `by_vars` or remove or rename it in one of the datasets."
      ),
      paste0(
        "The variables ",
        enumerate(common_vars),
        " are contained in both datasets.\n",
        "Please add them to `by_vars` or remove or rename them in one of the datasets."
      )
    ))
  }
  dataset <- left_join(dataset, add_data, by = vars2chr(by_vars))

  if (!is.null(missing_values)) {
    update_missings <- map2(
      syms(names(missing_values)),
      missing_values,
      ~ expr(if_else(is.na(!!match_flag_var), !!.y, !!.x))
    )
    names(update_missings) <- names(missing_values)
    dataset <- dataset %>%
      mutate(!!!update_missings) %>%
      remove_tmp_vars()
  }
  dataset
}

#' Merge a Categorization Variable
#'
#' @description
#' `r lifecycle::badge("deprecated")`
#'
#' This function is *deprecated*, please use `derive_vars_merged()` instead.
#'
#' Merge a categorization variable from a dataset to the input dataset. The
#' observations to merge can be selected by a condition and/or selecting the
#' first or last observation for each by group.
#'
#' @param dataset_add Additional dataset
#'
#'   The variables specified by the `by_vars`, the `source_var`, and the `order`
#'   argument are expected.
#'
#' @param new_var New variable
#'
#'   The specified variable is added to the additional dataset and set to the
#'   categorized values, i.e., `cat_fun(<source variable>)`.
#'
#' @param source_var Source variable
#'
#' @param cat_fun Categorization function
#'
#'   A function must be specified for this argument which expects the values of
#'   the source variable as input and returns the categorized values.
#'
#' @param missing_value Values used for missing information
#'
#'   The new variable is set to the specified value for all by groups without
#'   observations in the additional dataset.
#'
#'   *Default*: `NA_character_`
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
#'   1. The categorization variable is added to the additional dataset.
#'
#'   1. If `order` is specified, for each by group the first or last observation
#'   (depending on `mode`) is selected.
#'
#'   1. The categorization variable is merged to the input dataset.
#'
#'
#' @family deprecated
#' @keywords deprecated
#'
#' @export
#'
#' @examples
#' library(dplyr, warn.conflicts = FALSE)
#'
#' vs <- tribble(
#'   ~STUDYID,  ~DOMAIN,  ~USUBJID,      ~VISIT, ~VSTESTCD, ~VSSTRESN, ~VSSEQ,       ~VSDTC,
#'   "PILOT01",    "VS", "04-1127", "SCREENING",  "HEIGHT",     165.1,     43, "2013-09-16",
#'   "PILOT01",    "VS", "04-1127", "SCREENING",  "WEIGHT",     42.87,    142, "2013-09-16",
#'   "PILOT01",    "VS", "04-1127",  "BASELINE",  "WEIGHT",     41.05,    143, "2013-10-02",
#'   "PILOT01",    "VS", "04-1127",    "WEEK 2",  "WEIGHT",     42.64,    144, "2013-10-16",
#'   "PILOT01",    "VS", "04-1127",    "WEEK 4",  "WEIGHT",     41.73,    145, "2013-10-30",
#'   "PILOT01",    "VS", "04-1127",   "WEEK 26",  "WEIGHT",     43.09,    152, "2014-03-31",
#'   "PILOT01",    "VS", "06-1049", "SCREENING",  "HEIGHT",    167.64,     28, "2013-04-30",
#'   "PILOT01",    "VS", "06-1049", "SCREENING",  "WEIGHT",     57.61,     92, "2013-04-30",
#'   "PILOT01",    "VS", "06-1049",  "BASELINE",  "WEIGHT",     57.83,     93, "2013-05-14",
#'   "PILOT01",    "VS", "06-1049",    "WEEK 2",  "WEIGHT",     58.29,     94, "2013-05-28",
#'   "PILOT01",    "VS", "06-1049",    "WEEK 4",  "WEIGHT",     58.97,     95, "2013-06-11"
#' )
#'
#' dm <- tribble(
#'   ~STUDYID,  ~DOMAIN,  ~USUBJID, ~AGE,   ~AGEU,
#'   "PILOT01",    "DM", "01-1057",   59, "YEARS",
#'   "PILOT01",    "DM", "04-1127",   84, "YEARS",
#'   "PILOT01",    "DM", "06-1049",   60, "YEARS"
#' )
#' wgt_cat <- function(wgt) {
#'   case_when(
#'     wgt < 50 ~ "low",
#'     wgt > 90 ~ "high",
#'     TRUE ~ "normal"
#'   )
#' }
#'
#' derive_var_merged_cat(
#'   dm,
#'   dataset_add = vs,
#'   by_vars = exprs(STUDYID, USUBJID),
#'   order = exprs(VSDTC, VSSEQ),
#'   filter_add = VSTESTCD == "WEIGHT" & substr(VISIT, 1, 9) == "SCREENING",
#'   new_var = WGTBLCAT,
#'   source_var = VSSTRESN,
#'   cat_fun = wgt_cat,
#'   mode = "last"
#' ) %>%
#'   select(STUDYID, USUBJID, AGE, AGEU, WGTBLCAT)
#'
#'
#'
#' # defining a value for missing VS data
#' derive_var_merged_cat(
#'   dm,
#'   dataset_add = vs,
#'   by_vars = exprs(STUDYID, USUBJID),
#'   order = exprs(VSDTC, VSSEQ),
#'   filter_add = VSTESTCD == "WEIGHT" & substr(VISIT, 1, 9) == "SCREENING",
#'   new_var = WGTBLCAT,
#'   source_var = VSSTRESN,
#'   cat_fun = wgt_cat,
#'   mode = "last",
#'   missing_value = "MISSING"
#' ) %>%
#'   select(STUDYID, USUBJID, AGE, AGEU, WGTBLCAT)
derive_var_merged_cat <- function(dataset,
                                  dataset_add,
                                  by_vars,
                                  order = NULL,
                                  new_var,
                                  source_var,
                                  cat_fun,
                                  filter_add = NULL,
                                  mode = NULL,
                                  missing_value = NA_character_) {
  deprecate_warn("0.11.0", "derive_var_merged_cat()", "derive_vars_merged()")
  new_var <- assert_symbol(enexpr(new_var))
  source_var <- assert_symbol(enexpr(source_var))
  filter_add <- assert_filter_cond(enexpr(filter_add), optional = TRUE)
  assert_data_frame(dataset_add, required_vars = expr_c(by_vars, source_var))

  derive_vars_merged(
    dataset,
    dataset_add = dataset_add,
    filter_add = !!filter_add,
    by_vars = by_vars,
    order = order,
    new_vars = exprs(!!new_var := {{ cat_fun }}(!!source_var)),
    mode = mode,
    missing_values = exprs(!!new_var := !!missing_value)
  )
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
#'   *Permitted Values*: list of variables
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
#'   *Default*: `"Y"`
#'
#' @param false_value False value
#'
#'   *Default*: `NA_character_`
#'
#' @param missing_value Values used for missing information
#'
#'   The new variable is set to the specified value for all by groups without
#'   observations in the additional dataset.
#'
#'   *Default*: `NA_character_`
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
  new_var <- assert_symbol(enexpr(new_var))
  condition <- assert_filter_cond(enexpr(condition))
  filter_add <-
    assert_filter_cond(enexpr(filter_add), optional = TRUE)

  add_data <- filter_if(dataset_add, filter_add) %>%
    mutate(!!new_var := if_else(!!condition, 1, 0, 0))

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

#' Merge a Character Variable
#'
#' @description
#' `r lifecycle::badge("deprecated")`
#'
#' This function is *deprecated*, please use `derive_vars_merged()` instead.
#'
#' Merge a character variable from a dataset to the input dataset. The
#' observations to merge can be selected by a condition and/or selecting the
#' first or last observation for each by group.
#'
#' @param dataset_add Additional dataset
#'
#'   The variables specified by the `by_vars`, the `source_var`, and the `order`
#'   argument are expected.
#'
#' @param new_var New variable
#'
#'   The specified variable is added to the additional dataset and set to the
#'   transformed value with respect to the `case` argument.
#'
#' @param source_var Source variable
#'
#' @param case Change case
#'
#'   Changes the case of the values of the new variable.
#'
#'   *Default*: `NULL`
#'
#'   *Permitted Values*: `NULL`, `"lower"`, `"upper"`, `"title"`
#'
#' @param missing_value Values used for missing information
#'
#'   The new variable is set to the specified value for all by groups without
#'   observations in the additional dataset.
#'
#'   *Default*: `NA_character_`
#'
#'   *Permitted Value*: A character scalar
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
#'   1. The (transformed) character variable is added to the additional dataset.
#'
#'   1. If `order` is specified, for each by group the first or last observation
#'   (depending on `mode`) is selected.
#'
#'   1. The character variable is merged to the input dataset.
#'
#'
#' @family deprecated
#' @keywords deprecated
#'
#' @export
derive_var_merged_character <- function(dataset,
                                        dataset_add,
                                        by_vars,
                                        order = NULL,
                                        new_var,
                                        source_var,
                                        case = NULL,
                                        filter_add = NULL,
                                        mode = NULL,
                                        missing_value = NA_character_) {
  deprecate_warn("0.11.0", "derive_var_merged_character()", "derive_vars_merged()")

  new_var <- assert_symbol(enexpr(new_var))
  source_var <- assert_symbol(enexpr(source_var))
  case <-
    assert_character_scalar(
      case,
      values = c("lower", "upper", "title"),
      case_sensitive = FALSE,
      optional = TRUE
    )
  filter_add <- assert_filter_cond(enexpr(filter_add), optional = TRUE)
  assert_data_frame(dataset_add, required_vars = expr_c(by_vars, source_var))
  assert_character_scalar(missing_value)

  if (is.null(case)) {
    trans <- expr(!!source_var)
  } else if (case == "lower") {
    trans <- expr(str_to_lower(!!source_var))
  } else if (case == "upper") {
    trans <- expr(str_to_upper(!!source_var))
  } else if (case == "title") {
    trans <- expr(str_to_title(!!source_var))
  }
  derive_vars_merged(
    dataset,
    dataset_add = dataset_add,
    by_vars = by_vars,
    order = order,
    new_vars = exprs(!!new_var := !!trans),
    filter_add = !!filter_add,
    mode = mode,
    missing_values = exprs(!!new_var := !!missing_value)
  )
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

  res <- derive_vars_merged(
    dataset,
    dataset_add,
    by_vars = by_vars,
    order = order,
    new_vars = new_vars,
    mode = mode,
    filter_add = !!filter_add,
    match_flag = temp_match_flag,
    check_type = check_type,
    duplicate_msg = duplicate_msg
  )

  if (print_not_mapped) {
    temp_not_mapped <- res %>%
      filter(is.na(temp_match_flag)) %>%
      distinct(!!!by_vars_left)

    if (nrow(temp_not_mapped) > 0) {
      admiral_environment$nmap <- structure(
        temp_not_mapped,
        class = union("nmap", class(temp_not_mapped)),
        by_vars = vars2chr(by_vars_left)
      )

      message(
        "List of ", enumerate(vars2chr(by_vars_left)), " not mapped: ", "\n",
        paste0(capture.output(temp_not_mapped), collapse = "\n"),
        "\nRun `get_not_mapped()` to access the full list"
      )
    } else if (nrow(temp_not_mapped) == 0) {
      message(
        "All ", enumerate(vars2chr(by_vars_left)), " are mapped."
      )
    }
  }

  res %>% select(-temp_match_flag)
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

#' Merge a Summary Variable
#'
#' @description Merge a summary variable from a dataset to the input dataset.
#'
#' **Note:** This is a wrapper function for the more generic `derive_vars_merged`.
#'
#' @param dataset Input dataset
#'
#'   The variables specified by the `by_vars` argument are expected.
#'
#' @param dataset_add Additional dataset
#'
#'   The variables specified by the `by_vars` and the `analysis_var` arguments
#'   are expected.
#'
#' @param new_var Variable to add
#'
#'   The specified variable is added to the input dataset (`dataset`) and set to
#'   the summarized values.
#'
#' @param by_vars Grouping variables
#'
#'   The values of `analysis_var` are summarized by the specified variables. The
#'   summarized values are merged to the input dataset (`dataset`) by the
#'   specified by variables.
#'
#'   *Permitted Values*: list of variables created by `exprs()`
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
#'   The values of the specified variable are summarized by the function
#'   specified for `summary_fun`.
#'
#' @param summary_fun Summary function
#'
#'   The specified function that takes as input `analysis_var` and performs the
#'   calculation. This can include built-in functions as well as user defined
#'   functions, for example `mean` or `function(x) mean(x, na.rm = TRUE)`.
#'
#'
#' @details
#'
#'   1. The records from the additional dataset (`dataset_add`) are restricted
#'   to those matching the `filter_add` condition.
#'
#'   1. The values of the analysis variable (`analysis_var`) are summarized by
#'   the summary function (`summary_fun`) for each by group (`by_vars`) in the
#'   additional dataset (`dataset_add`).
#'
#'   1. The summarized values are merged to the input dataset as a new variable
#'   (`new_var`). For observations without a matching observation in the
#'   additional dataset the new variable is set to `NA`. Observations in the
#'   additional dataset which have no matching observation in the input dataset
#'   are ignored.
#'
#' @return The output dataset contains all observations and variables of the
#'   input dataset and additionally the variable specified for `new_var`.
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
#'   new_var = MEANVIS,
#'   analysis_var = AVAL,
#'   summary_fun = function(x) mean(x, na.rm = TRUE)
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
#'   new_var = LESIONSBL,
#'   analysis_var = LESIONID,
#'   summary_fun = function(x) paste(x, collapse = ", ")
#' )
#'
derive_var_merged_summary <- function(dataset,
                                      dataset_add,
                                      by_vars,
                                      new_var,
                                      filter_add = NULL,
                                      analysis_var,
                                      summary_fun) {
  assert_vars(by_vars)
  by_vars_left <- replace_values_by_names(by_vars)
  by_vars_right <- chr2vars(paste(vars2chr(by_vars)))
  new_var <- assert_symbol(enexpr(new_var))
  analysis_var <- assert_symbol(enexpr(analysis_var))
  filter_add <-
    assert_filter_cond(enexpr(filter_add), optional = TRUE)
  assert_s3_class(summary_fun, "function")
  assert_data_frame(
    dataset,
    required_vars = by_vars_left
  )
  assert_data_frame(
    dataset_add,
    required_vars = expr_c(by_vars_right, analysis_var)
  )

  # Summarise the analysis value and merge to the original dataset
  derive_vars_merged(
    dataset,
    dataset_add = get_summary_records(
      dataset_add,
      by_vars = by_vars_right,
      filter = !!filter_add,
      analysis_var = !!analysis_var,
      summary_fun = summary_fun
    ),
    by_vars = by_vars,
    new_vars = exprs(!!new_var := !!analysis_var)
  )
}
