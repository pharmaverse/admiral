#' Add Variables from an Additional Dataset Based on Conditions from Both
#' Datasets
#'
#' The function adds variables from an additional dataset to the input dataset.
#' The selection of the observations from the additional dataset can depend on
#' variables from both datasets. For example, add the lowest value (nadir)
#' before the current observation.
#'
#' @param dataset Input dataset
#'
#'   The variables specified by `by_vars` are expected.
#'
#' @param dataset_add Additional dataset
#'
#'   The variables specified by the `by_vars`, the `new_vars`, the `join_vars`,
#'   and the `order` argument are expected.
#'
#' @param by_vars Grouping variables
#'
#'   The two datasets are joined by the specified variables. Variables from the
#'   additional dataset can be renamed by naming the element, i.e., `by_vars =
#'   exprs(<name in input dataset> = <name in additional dataset>)`.
#'
#'   *Permitted Values*: list of variables created by `exprs()`
#'
#' @param order Sort order
#'
#'   If the argument is set to a non-null value, for each observation of the
#'   input dataset the first or last observation from the joined dataset is
#'   selected with respect to the specified order. The specified variables are
#'   expected in the additional dataset (`dataset_add`). If a variable is
#'   available in both `dataset` and `dataset_add`, the one from `dataset_add`
#'   is used for the sorting.
#'
#'   If an expression is named, e.g., `exprs(EXSTDT =
#'   convert_dtc_to_dt(EXSTDTC), EXSEQ)`, a corresponding variable (`EXSTDT`) is
#'   added to the additional dataset and can be used in the filter conditions
#'   (`filter_add`, `filter_join`) and for `join_vars` and `new_vars`. The
#'   variable is not included in the output dataset.
#'
#'   *Permitted Values*: list of expressions created by `exprs()`, e.g.,
#'    `exprs(ADT, desc(AVAL))` or `NULL`
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
#' @param join_vars Variables to use from additional dataset
#'
#'   Any extra variables required from the additional dataset for `filter_join`
#'   should be specified for this argument. Variables specified for `new_vars`
#'   do not need to be repeated for `join_vars`. If a specified variable exists
#'   in both the input dataset and the additional dataset, the suffix ".join" is
#'   added to the variable from the additional dataset.
#'
#'   If an expression is named, e.g., `exprs(EXTDT =
#'   convert_dtc_to_dt(EXSTDTC))`, a corresponding variable is added to the
#'   additional dataset and can be used in the filter conditions (`filter_add`,
#'   `filter_join`) and for `new_vars`. The variable is not included in the
#'   output dataset.
#'
#'   The variables are not included in the output dataset.
#'
#'   *Permitted Values*: list of variables or named expressions created by `exprs()`
#'
#' @param filter_add Filter for additional dataset (`dataset_add`)
#'
#'   Only observations from `dataset_add` fulfilling the specified condition are
#'   joined to the input dataset. If the argument is not specified, all
#'   observations are joined.
#'
#'   Variables created by `order` or `new_vars` arguments can be used in the
#'   condition.
#'
#'   *Permitted Values*: a condition
#'
#' @param filter_join Filter for the joined dataset
#'
#'   The specified condition is applied to the joined dataset. Therefore
#'   variables from both datasets `dataset` and `dataset_add` can be used.
#'
#'   Variables created by `order` or `new_vars` arguments can be used in the
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
#' @param check_type Check uniqueness?
#'
#'   If `"warning"` or `"error"` is specified, the specified message is issued
#'   if the observations of the (restricted) joined dataset are not unique with
#'   respect to the by variables and the order.
#'
#'   This argument is ignored if `order` is not specified. In this case an error
#'   is issued independent of `check_type` if the restricted joined dataset
#'   contains more than one observation for any of the observations of the input
#'   dataset.
#'
#'   *Permitted Values*: `"none"`, `"warning"`, `"error"`
#'
#'
#' @details
#'
#' 1. The variables specified by `order` are added to the additional dataset
#' (`dataset_add`).
#'
#' 1. The variables specified by `join_vars` are added to the additional dataset
#' (`dataset_add`).
#'
#' 1. The records from the additional dataset (`dataset_add`) are restricted to
#' those matching the `filter_add` condition.
#'
#' 1. The input dataset and the (restricted) additional dataset are left joined
#' by the grouping variables (`by_vars`). If no grouping variables are
#' specified, a full join is performed.
#'
#' 1. The joined dataset is restricted by the `filter_join` condition.
#'
#' 1. If `order` is specified, for each observation of the input dataset the
#' first or last observation (depending on `mode`) is selected.
#'
#' 1. The variables specified for `new_vars` are created (if requested) and
#' merged to the input dataset. I.e., the output dataset contains all
#' observations from the input dataset. For observations without a matching
#' observation in the joined dataset the new variables are set as specified by
#' `missing_values` (or to `NA` for variables not in `missing_values`).
#' Observations in the additional dataset which have no matching observation in
#' the input dataset are ignored.
#'
#' @inheritParams derive_vars_merged
#'
#' @return The output dataset contains all observations and variables of the
#'   input dataset and additionally the variables specified for `new_vars` from
#'   the additional dataset (`dataset_add`).
#'
#' @keywords der_gen
#' @family der_gen
#'
#' @export
#'
#' @examples
#' library(tibble)
#' library(lubridate)
#' library(dplyr, warn.conflicts = FALSE)
#' library(tidyr)
#'
#' # Add AVISIT (based on time windows), AWLO, and AWHI
#' adbds <- tribble(
#'   ~USUBJID, ~ADY,
#'   "1",       -33,
#'   "1",        -2,
#'   "1",         3,
#'   "1",        24,
#'   "2",        NA,
#' )
#'
#' windows <- tribble(
#'   ~AVISIT,    ~AWLO, ~AWHI,
#'   "BASELINE",   -30,     1,
#'   "WEEK 1",       2,     7,
#'   "WEEK 2",       8,    15,
#'   "WEEK 3",      16,    22,
#'   "WEEK 4",      23,    30
#' )
#'
#' derive_vars_joined(
#'   adbds,
#'   dataset_add = windows,
#'   filter_join = AWLO <= ADY & ADY <= AWHI
#' )
#'
#' # derive the nadir after baseline and before the current observation
#' adbds <- tribble(
#'   ~USUBJID, ~ADY, ~AVAL,
#'   "1",        -7,    10,
#'   "1",         1,    12,
#'   "1",         8,    11,
#'   "1",        15,     9,
#'   "1",        20,    14,
#'   "1",        24,    12,
#'   "2",        13,     8
#' )
#'
#' derive_vars_joined(
#'   adbds,
#'   dataset_add = adbds,
#'   by_vars = exprs(USUBJID),
#'   order = exprs(AVAL),
#'   new_vars = exprs(NADIR = AVAL),
#'   join_vars = exprs(ADY),
#'   filter_add = ADY > 0,
#'   filter_join = ADY.join < ADY,
#'   mode = "first",
#'   check_type = "none"
#' )
#'
#' # add highest hemoglobin value within two weeks before AE,
#' # take earliest if more than one
#' adae <- tribble(
#'   ~USUBJID, ~ASTDY,
#'   "1",           3,
#'   "1",          22,
#'   "2",           2
#' )
#'
#' adlb <- tribble(
#'   ~USUBJID, ~PARAMCD, ~ADY, ~AVAL,
#'   "1",      "HGB",       1,   8.5,
#'   "1",      "HGB",       3,   7.9,
#'   "1",      "HGB",       5,   8.9,
#'   "1",      "HGB",       8,   8.0,
#'   "1",      "HGB",       9,   8.0,
#'   "1",      "HGB",      16,   7.4,
#'   "1",      "HGB",      24,   8.1,
#'   "1",      "ALB",       1,    42,
#' )
#'
#' derive_vars_joined(
#'   adae,
#'   dataset_add = adlb,
#'   by_vars = exprs(USUBJID),
#'   order = exprs(AVAL, desc(ADY)),
#'   new_vars = exprs(HGB_MAX = AVAL, HGB_DY = ADY),
#'   filter_add = PARAMCD == "HGB",
#'   filter_join = ASTDY - 14 <= ADY & ADY <= ASTDY,
#'   mode = "last"
#' )
#'
#' # Add APERIOD, APERIODC based on ADSL
#' adsl <- tribble(
#'   ~USUBJID, ~AP01SDT,     ~AP01EDT,     ~AP02SDT,     ~AP02EDT,
#'   "1",      "2021-01-04", "2021-02-06", "2021-02-07", "2021-03-07",
#'   "2",      "2021-02-02", "2021-03-02", "2021-03-03", "2021-04-01"
#' ) %>%
#'   mutate(across(ends_with("DT"), ymd)) %>%
#'   mutate(STUDYID = "xyz")
#'
#' period_ref <- create_period_dataset(
#'   adsl,
#'   new_vars = exprs(APERSDT = APxxSDT, APEREDT = APxxEDT)
#' )
#'
#' period_ref
#'
#' adae <- tribble(
#'   ~USUBJID, ~ASTDT,
#'   "1",      "2021-01-01",
#'   "1",      "2021-01-05",
#'   "1",      "2021-02-05",
#'   "1",      "2021-03-05",
#'   "1",      "2021-04-05",
#'   "2",      "2021-02-15",
#' ) %>%
#'   mutate(
#'     ASTDT = ymd(ASTDT),
#'     STUDYID = "xyz"
#'   )
#'
#' derive_vars_joined(
#'   adae,
#'   dataset_add = period_ref,
#'   by_vars = exprs(STUDYID, USUBJID),
#'   join_vars = exprs(APERSDT, APEREDT),
#'   filter_join = APERSDT <= ASTDT & ASTDT <= APEREDT
#' )
#'
#' # Add day since last dose (LDRELD)
#' adae <- tribble(
#'   ~USUBJID, ~ASTDT,       ~AESEQ,
#'   "1",      "2020-02-02",      1,
#'   "1",      "2020-02-04",      2
#' ) %>%
#'   mutate(ASTDT = ymd(ASTDT))
#'
#' ex <- tribble(
#'   ~USUBJID, ~EXSDTC,
#'   "1",      "2020-01-10",
#'   "1",      "2020-01",
#'   "1",      "2020-01-20",
#'   "1",      "2020-02-03"
#' )
#'
#' ## Please note that EXSDT is created via the order argument and then used
#' ## for new_vars, filter_add, and filter_join
#' derive_vars_joined(
#'   adae,
#'   dataset_add = ex,
#'   by_vars = exprs(USUBJID),
#'   order = exprs(EXSDT = convert_dtc_to_dt(EXSDTC)),
#'   new_vars = exprs(LDRELD = compute_duration(
#'     start_date = EXSDT, end_date = ASTDT
#'   )),
#'   filter_add = !is.na(EXSDT),
#'   filter_join = EXSDT <= ASTDT,
#'   mode = "last"
#' )
derive_vars_joined <- function(dataset,
                               dataset_add,
                               by_vars = NULL,
                               order = NULL,
                               new_vars = NULL,
                               join_vars = NULL,
                               filter_add = NULL,
                               filter_join = NULL,
                               mode = NULL,
                               missing_values = NULL,
                               check_type = "warning") {
  assert_vars(by_vars, optional = TRUE)
  by_vars_left <- replace_values_by_names(by_vars)
  assert_expr_list(order, optional = TRUE)
  assert_expr_list(new_vars, optional = TRUE)
  assert_expr_list(join_vars, optional = TRUE)
  assert_data_frame(dataset, required_vars = by_vars_left)
  assert_data_frame(
    dataset_add,
    required_vars = expr_c(
      by_vars,
      extract_vars(order),
      setdiff(extract_vars(join_vars), replace_values_by_names(order))
    )
  )

  filter_add <- assert_filter_cond(enexpr(filter_add), optional = TRUE)
  filter_join <- assert_filter_cond(enexpr(filter_join), optional = TRUE)

  if (is.null(new_vars)) {
    new_vars <- chr2vars(colnames(dataset_add))
  }

  # number observations of the input dataset to get a unique key
  # (by_vars and tmp_obs_nr)
  tmp_obs_nr <- get_new_tmp_var(dataset, prefix = "tmp_obs_nr")
  data <- dataset %>%
    derive_var_obs_number(
      new_var = !!tmp_obs_nr,
      by_vars = by_vars_left,
      check_type = "none"
    )

  # prepare right side of the join,
  # by_vars are renamed here, new_vars will be renamed at the end
  data_right <- dataset_add %>%
    mutate(!!!order, !!!join_vars) %>%
    filter_if(filter_add) %>%
    select(
      !!!by_vars,
      !!!chr2vars(names(order)),
      !!!replace_values_by_names(join_vars),
      !!!intersect(unname(extract_vars(new_vars)), chr2vars(colnames(dataset_add)))
    )

  # join dataset (if no by variable, a full join is performed)
  data_joined <- left_join(
    data,
    data_right,
    by = vars2chr(by_vars_left),
    suffix = c("", ".join")
  )

  # select observations for the new variables
  data_return <- filter_if(data_joined, filter_join)

  common_vars <-
    chr2vars(setdiff(intersect(colnames(data), colnames(data_right)), vars2chr(by_vars)))
  if (!is.null(order)) {
    data_return <- filter_extreme(
      data_return,
      by_vars = expr_c(by_vars_left, tmp_obs_nr),
      order = add_suffix_to_vars(
        replace_values_by_names(order),
        vars = common_vars,
        suffix = ".join"
      ),
      mode = mode,
      check_type = check_type
    )
  }

  # merge new variables to the input dataset and rename them
  data %>%
    derive_vars_merged(
      dataset_add = data_return,
      by_vars = exprs(!!!by_vars_left, !!tmp_obs_nr),
      new_vars = add_suffix_to_vars(new_vars, vars = common_vars, suffix = ".join"),
      missing_values = missing_values,
      duplicate_msg = paste(
        paste(
          "After applying `filter_join` the joined dataset contains more",
          "than one observation per observation of the input dataset."
        ),
        paste(
          "Please adjust `filter_add` and/or `filter_join` or specify `order`",
          "and `mode` to select one of the observations."
        ),
        sep = "\n"
      )
    ) %>%
    remove_tmp_vars()
}
