#' Derives a Flag Based on an Existing Flag
#'
#' Derive a flag which depends on other observations of the dataset. For
#' example, flagging events which need to be confirmed by a second event.
#'
#' An example usage might be flagging if a patient received two required
#' medications within a certain timeframe of each other.
#'
#' In the oncology setting, for example, the function could be used to flag if a
#' response value can be confirmed by an other assessment. This is commonly
#' used in endpoints such as best overall response.
#'
#' @param dataset
#'   `r roxygen_param_dataset(expected_vars = c("by_vars", "join_vars"))`
#'
#' @param dataset_add Additional dataset
#'
#'   The variables specified for `by_vars`, `join_vars`, and `order` are
#'   expected.
#'
#' @param by_vars Grouping variables
#'
#'   The specified variables are used for joining the input
#'   dataset (`dataset`) with the additional dataset (`dataset_add`).
#'
#'   `r roxygen_param_by_vars()`
#'
#' @param order Order
#'
#'   The observations are ordered by the specified order.
#'
#'   `r roxygen_order_na_handling()`
#'
#' @param new_var New variable
#'
#'   The specified variable is added to the input dataset.
#'
#' @param tmp_obs_nr_var Temporary observation number
#'
#'   The specified variable is added to the input dataset (`dataset`) and the
#'   additional dataset (`dataset_add`). It is set to the observation number
#'   with respect to `order`. For each by group (`by_vars`) the observation
#'   number starts with `1`. The variable can be used in the conditions
#'   (`filter_join`, `first_cond_upper`, `first_cond_lower`). It is not included
#'   in the output dataset. It can also be used to flag consecutive observations
#'   or the last observation (see last example below).
#'
#' @param join_vars Variables to keep from joined dataset
#'
#'   The variables needed from the other observations should be specified
#'   for this parameter. The specified variables are added to the joined dataset
#'   with suffix ".join". For example to flag all observations with `AVALC ==
#'   "Y"` and `AVALC == "Y"` for at least one subsequent visit `join_vars =
#'   exprs(AVALC, AVISITN)` and `filter_join = AVALC == "Y" & AVALC.join == "Y" &
#'   AVISITN < AVISITN.join` could be specified.
#'
#'   The `*.join` variables are not included in the output dataset.
#'
#' @param first_cond Condition for selecting range of data
#'
#'   `r lifecycle::badge("deprecated")`
#'
#'   This argument is *deprecated*, please use `first_cond_upper` instead.
#'
#'   If this argument is specified, the other observations are restricted up to
#'   the first observation where the specified condition is fulfilled. If the
#'   condition is not fulfilled for any of the other observations, no
#'   observations are considered, i.e., the observation is not flagged.
#'
#'   This parameter should be specified if `filter_join` contains summary
#'   functions which should not apply to all observations but only up to the
#'   confirmation assessment. For an example see the third example below.
#'
#' @param first_cond_lower Condition for selecting range of data (before)
#'
#'   If this argument is specified, the other observations are restricted from
#'   the first observation before the current observation where the specified
#'   condition is fulfilled up to the current observation. If the condition is
#'   not fulfilled for any of the other observations, no observations are
#'   considered, i.e., the observation is not flagged.
#'
#'   This parameter should be specified if `filter_join` contains summary
#'   functions which should not apply to all observations but only from a
#'   certain observation before the current observation up to the current
#'   observation. For an example see the last example below.
#'
#' @param first_cond_upper Condition for selecting range of data (after)
#'
#'   If this argument is specified, the other observations are restricted up to
#'   the first observation where the specified condition is fulfilled. If the
#'   condition is not fulfilled for any of the other observations, no
#'   observations are considered, i.e., the observation is not flagged.
#'
#'   This parameter should be specified if `filter_join` contains summary
#'   functions which should not apply to all observations but only up to the
#'   confirmation assessment. For an example see the third example below.
#'
#' @param filter_join Condition for selecting observations
#'
#'   The filter is applied to the joined dataset for flagging the confirmed
#'   observations. The condition can include summary functions like `all()` or
#'   `any()`. The joined dataset is grouped by the original observations. I.e.,
#'   the summary function are applied to all observations up to the confirmation
#'   observation. For example, `filter_join = AVALC == "CR" & all(AVALC.join
#'   %in% c("CR", "NE")) & count_vals(var = AVALC.join, val = "NE") <= 1`
#'   selects observations with response "CR" and for all observations up to the
#'   confirmation observation the response is "CR" or "NE" and there is at most
#'   one "NE".
#'
#' @param filter Condition for selecting observations
#'
#'   `r lifecycle::badge("deprecated")`
#'
#'   This argument is *deprecated*, please use `filter_join` instead.
#'
#'   The filter is applied to the joined dataset for flagging the confirmed
#'   observations. The condition can include summary functions. The joined
#'   dataset is grouped by the original observations. I.e., the summary function
#'   are applied to all observations up to the confirmation observation. For
#'   example, `filter = AVALC == "CR" & all(AVALC.join %in% c("CR", "NE")) &
#'   count_vals(var = AVALC.join, val = "NE") <= 1` selects observations with
#'   response "CR" and for all observations up to the confirmation observation
#'   the response is "CR" or "NE" and there is at most one "NE".
#'
#' @param check_type Check uniqueness?
#'
#'   If `"warning"` or `"error"` is specified, the specified message is issued
#'   if the observations of the input dataset are not unique with respect to the
#'   by variables and the order.
#'
#'   *Permitted Values:* `"none"`, `"warning"`, `"error"`
#'
#' @param true_value Value of `new_var` for flagged observations
#'
#' @param false_value Value of `new_var` for observations not flagged
#'
#' @inheritParams get_joined_data
#'
#' @details
#'   The following steps are performed to produce the output dataset.
#'
#'   ## Step 1
#'
#'   - The variables specified by `order` are added to the additional dataset
#'   (`dataset_add`).
#'
#'   - The variables specified by `join_vars` are added to the additional dataset
#'   (`dataset_add`).
#'
#'   - The records from the additional dataset (`dataset_add`) are restricted to
#'   those matching the `filter_add` condition.
#'
#'   The input dataset (`dataset`) is joined with the restricted additional
#'   dataset by the variables specified for `by_vars`. From the additional
#'   dataset only the variables specified for `join_vars` are kept. The suffix
#'   ".join" is added to those variables which also exist in the input dataset.
#'
#'   For example, for `by_vars = USUBJID`, `join_vars = exprs(AVISITN, AVALC)`
#'   and input dataset and additional dataset
#'
#'   ```{r eval=FALSE}
#'   # A tibble: 2 x 4
#'   USUBJID AVISITN AVALC  AVAL
#'   <chr>     <dbl> <chr> <dbl>
#'   1             1 Y         1
#'   1             2 N         0
#'   ```
#'
#'   the joined dataset is
#'
#'   ```{r eval=FALSE}
#'   A tibble: 4 x 6
#'   USUBJID AVISITN AVALC  AVAL AVISITN.join AVALC.join
#'   <chr>     <dbl> <chr> <dbl>        <dbl> <chr>
#'   1             1 Y         1            1 Y
#'   1             1 Y         1            2 N
#'   1             2 N         0            1 Y
#'   1             2 N         0            2 N
#'   ```
#'
#'   ## Step 2
#'
#'   The joined dataset is restricted to observations with respect to
#'   `join_type` and `order`.
#'
#'   The dataset from the example in the previous step with `join_type =
#'   "after"` and `order = exprs(AVISITN)` is restricted to
#'
#'   ```{r eval=FALSE}
#'   A tibble: 4 x 6
#'   USUBJID AVISITN AVALC  AVAL AVISITN.join AVALC.join
#'   <chr>     <dbl> <chr> <dbl>        <dbl> <chr>
#'   1             1 Y         1            2 N
#'   ```
#'
#'   ## Step 3
#'
#'   If `first_cond_lower` is specified, for each observation of the input
#'   dataset the joined dataset is restricted to observations from the first
#'   observation where `first_cond_lower` is fulfilled (the observation
#'   fulfilling the condition is included) up to the observation of the input
#'   dataset. If for an observation of the input dataset the condition is not
#'   fulfilled, the observation is removed.
#'
#'   If `first_cond_upper` is specified, for each observation of the input
#'   dataset the joined dataset is restricted to observations up to the first
#'   observation where `first_cond_upper` is fulfilled (the observation
#'   fulfilling the condition is included). If for an observation of the input
#'   dataset the condition is not fulfilled, the observation is removed.
#'
#'   For an example see the last example in the "Examples" section.
#'
#'   ## Step 4
#'
#'   The joined dataset is grouped by the observations from the input dataset
#'   and restricted to the observations fulfilling the condition specified by
#'   `filter_join`.
#'
#'   ## Step 5
#'
#'   The first observation of each group is selected
#'
#'   ## Step 6
#'
#'   The variable specified by `new_var` is added to the input dataset. It is
#'   set to `true_value` for all observations which were selected in the
#'   previous step. For the other observations it is set to `false_value`.
#'
#' @return The input dataset with the variable specified by `new_var` added.
#'
#'
#' @keywords der_gen
#' @family der_gen
#'
#' @seealso [filter_joined()], [derive_vars_joined()]
#'
#' @export
#'
#' @examples
#' library(tibble)
#'
#' # flag observations with a duration longer than 30 and
#' # at, after, or up to 7 days before a COVID AE (ACOVFL == "Y")
#' adae <- tribble(
#'   ~USUBJID, ~ADY, ~ACOVFL, ~ADURN,
#'   "1",        10, "N",          1,
#'   "1",        21, "N",         50,
#'   "1",        23, "Y",         14,
#'   "1",        32, "N",         31,
#'   "1",        42, "N",         20,
#'   "2",        11, "Y",         13,
#'   "2",        23, "N",          2,
#'   "3",        13, "Y",         12,
#'   "4",        14, "N",         32,
#'   "4",        21, "N",         41
#' )
#'
#' derive_var_joined_exist_flag(
#'   adae,
#'   dataset_add = adae,
#'   new_var = ALCOVFL,
#'   by_vars = exprs(USUBJID),
#'   join_vars = exprs(ACOVFL, ADY),
#'   join_type = "all",
#'   order = exprs(ADY),
#'   filter_join = ADURN > 30 & ACOVFL.join == "Y" & ADY >= ADY.join - 7
#' )
#'
#' # flag observations with AVALC == "Y" and AVALC == "Y" at one subsequent visit
#' data <- tribble(
#'   ~USUBJID, ~AVISITN, ~AVALC,
#'   "1",      1,        "Y",
#'   "1",      2,        "N",
#'   "1",      3,        "Y",
#'   "1",      4,        "N",
#'   "2",      1,        "Y",
#'   "2",      2,        "N",
#'   "3",      1,        "Y",
#'   "4",      1,        "N",
#'   "4",      2,        "N",
#' )
#'
#' derive_var_joined_exist_flag(
#'   data,
#'   dataset_add = data,
#'   by_vars = exprs(USUBJID),
#'   new_var = CONFFL,
#'   join_vars = exprs(AVALC, AVISITN),
#'   join_type = "after",
#'   order = exprs(AVISITN),
#'   filter_join = AVALC == "Y" & AVALC.join == "Y" & AVISITN < AVISITN.join
#' )
#'
#' # select observations with AVALC == "CR", AVALC == "CR" at a subsequent visit,
#' # only "CR" or "NE" in between, and at most one "NE" in between
#' data <- tribble(
#'   ~USUBJID, ~AVISITN, ~AVALC,
#'   "1",      1,        "PR",
#'   "1",      2,        "CR",
#'   "1",      3,        "NE",
#'   "1",      4,        "CR",
#'   "1",      5,        "NE",
#'   "2",      1,        "CR",
#'   "2",      2,        "PR",
#'   "2",      3,        "CR",
#'   "3",      1,        "CR",
#'   "4",      1,        "CR",
#'   "4",      2,        "NE",
#'   "4",      3,        "NE",
#'   "4",      4,        "CR",
#'   "4",      5,        "PR"
#' )
#'
#' derive_var_joined_exist_flag(
#'   data,
#'   dataset_add = data,
#'   by_vars = exprs(USUBJID),
#'   join_vars = exprs(AVALC),
#'   join_type = "after",
#'   order = exprs(AVISITN),
#'   new_var = CONFFL,
#'   first_cond_upper = AVALC.join == "CR",
#'   filter_join = AVALC == "CR" & all(AVALC.join %in% c("CR", "NE")) &
#'     count_vals(var = AVALC.join, val = "NE") <= 1
#' )
#'
#' # flag observations with AVALC == "PR", AVALC == "CR" or AVALC == "PR"
#' # at a subsequent visit at least 20 days later, only "CR", "PR", or "NE"
#' # in between, at most one "NE" in between, and "CR" is not followed by "PR"
#' data <- tribble(
#'   ~USUBJID, ~ADY, ~AVALC,
#'   "1",         6, "PR",
#'   "1",        12, "CR",
#'   "1",        24, "NE",
#'   "1",        32, "CR",
#'   "1",        48, "PR",
#'   "2",         3, "PR",
#'   "2",        21, "CR",
#'   "2",        33, "PR",
#'   "3",        11, "PR",
#'   "4",         7, "PR",
#'   "4",        12, "NE",
#'   "4",        24, "NE",
#'   "4",        32, "PR",
#'   "4",        55, "PR"
#' )
#'
#' derive_var_joined_exist_flag(
#'   data,
#'   dataset_add = data,
#'   by_vars = exprs(USUBJID),
#'   join_vars = exprs(AVALC, ADY),
#'   join_type = "after",
#'   order = exprs(ADY),
#'   new_var = CONFFL,
#'   first_cond_upper = AVALC.join %in% c("CR", "PR") & ADY.join - ADY >= 20,
#'   filter_join = AVALC == "PR" &
#'     all(AVALC.join %in% c("CR", "PR", "NE")) &
#'     count_vals(var = AVALC.join, val = "NE") <= 1 &
#'     (
#'       min_cond(var = ADY.join, cond = AVALC.join == "CR") >
#'         max_cond(var = ADY.join, cond = AVALC.join == "PR") |
#'         count_vals(var = AVALC.join, val = "CR") == 0
#'     )
#' )
#'
#' # flag observations with CRIT1FL == "Y" at two consecutive visits or at the last visit
#' data <- tribble(
#'   ~USUBJID, ~AVISITN, ~CRIT1FL,
#'   "1",      1,        "Y",
#'   "1",      2,        "N",
#'   "1",      3,        "Y",
#'   "1",      5,        "N",
#'   "2",      1,        "Y",
#'   "2",      3,        "Y",
#'   "2",      5,        "N",
#'   "3",      1,        "Y",
#'   "4",      1,        "Y",
#'   "4",      2,        "N",
#' )
#'
#' derive_var_joined_exist_flag(
#'   data,
#'   dataset_add = data,
#'   by_vars = exprs(USUBJID),
#'   new_var = CONFFL,
#'   tmp_obs_nr_var = tmp_obs_nr,
#'   join_vars = exprs(CRIT1FL),
#'   join_type = "all",
#'   order = exprs(AVISITN),
#'   filter_join = CRIT1FL == "Y" & CRIT1FL.join == "Y" &
#'     (tmp_obs_nr + 1 == tmp_obs_nr.join | tmp_obs_nr == max(tmp_obs_nr.join))
#' )
#'
#' # first_cond_lower and first_cond_upper argument
#' myd <- tribble(
#'   ~subj, ~day, ~val,
#'   "1",      1, "++",
#'   "1",      2, "-",
#'   "1",      3, "0",
#'   "1",      4, "+",
#'   "1",      5, "++",
#'   "1",      6, "-",
#'   "2",      1, "-",
#'   "2",      2, "++",
#'   "2",      3, "+",
#'   "2",      4, "0",
#'   "2",      5, "-",
#'   "2",      6, "++"
#' )
#'
#' # flag "0" where all results from the first "++" before the "0" up to the "0"
#' # (excluding the "0") are "+" or "++"
#' derive_var_joined_exist_flag(
#'   myd,
#'   dataset_add = myd,
#'   by_vars = exprs(subj),
#'   order = exprs(day),
#'   new_var = flag,
#'   join_vars = exprs(val),
#'   join_type = "before",
#'   first_cond_lower = val.join == "++",
#'   filter_join = val == "0" & all(val.join %in% c("+", "++"))
#' )
#'
#' # flag "0" where all results from the "0" (excluding the "0") up to the first
#' # "++" after the "0" are "+" or "++"
#' derive_var_joined_exist_flag(
#'   myd,
#'   dataset_add = myd,
#'   by_vars = exprs(subj),
#'   order = exprs(day),
#'   new_var = flag,
#'   join_vars = exprs(val),
#'   join_type = "after",
#'   first_cond_upper = val.join == "++",
#'   filter_join = val == "0" & all(val.join %in% c("+", "++"))
#' )
derive_var_joined_exist_flag <- function(dataset,
                                         dataset_add,
                                         by_vars,
                                         order,
                                         new_var,
                                         tmp_obs_nr_var = NULL,
                                         join_vars,
                                         join_type,
                                         first_cond = NULL,
                                         first_cond_lower = NULL,
                                         first_cond_upper = NULL,
                                         filter = NULL,
                                         filter_add = NULL,
                                         filter_join,
                                         true_value = "Y",
                                         false_value = NA_character_,
                                         check_type = "warning") {
  new_var <- assert_symbol(enexpr(new_var))
  tmp_obs_nr_var <- assert_symbol(enexpr(tmp_obs_nr_var), optional = TRUE)
  first_cond_lower <- assert_filter_cond(enexpr(first_cond_lower), optional = TRUE)
  first_cond_upper <- assert_filter_cond(enexpr(first_cond_upper), optional = TRUE)
  if (!missing(first_cond)) {
    deprecate_stop(
      "1.1.0",
      "derive_var_joined_exist_flag(first_cond=)",
      "derive_var_joined_exist_flag(first_cond_upper=)"
    )
    first_cond_upper <- assert_filter_cond(enexpr(first_cond), optional = TRUE)
  }
  filter_add <- assert_filter_cond(enexpr(filter_add), optional = TRUE)
  filter_join <- assert_filter_cond(enexpr(filter_join))
  if (!missing(filter)) {
    deprecate_stop(
      "1.1.0",
      "derive_var_joined_exist_flag(filter=)",
      "derive_var_joined_exist_flag(filter_join=)"
    )
    filter_join <- assert_filter_cond(enexpr(filter))
  }
  assert_data_frame(dataset)

  tmp_obs_nr <- get_new_tmp_var(dataset, prefix = "tmp_obs_nr_")

  data <- derive_var_obs_number(
    dataset,
    new_var = !!tmp_obs_nr
  )

  data_filtered <- filter_joined(
    data,
    dataset_add = dataset_add,
    by_vars = by_vars,
    order = order,
    tmp_obs_nr_var = !!tmp_obs_nr_var,
    join_vars = join_vars,
    join_type = join_type,
    first_cond_lower = !!first_cond_lower,
    first_cond_upper = !!first_cond_upper,
    filter_join = !!filter_join,
    check_type = check_type
  )

  derive_var_merged_exist_flag(
    data,
    dataset_add = data_filtered,
    by_vars = exprs(!!tmp_obs_nr),
    new_var = !!new_var,
    condition = TRUE,
    true_value = true_value,
    false_value = false_value,
    missing_value = false_value
  ) %>%
    remove_tmp_vars()
}
