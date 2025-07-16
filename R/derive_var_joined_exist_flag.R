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
#'
#'   `r roxygen_param_dataset(expected_vars = c("by_vars", "join_vars"))`
#'
#' @permitted [dataset]
#'
#' @param dataset_add Additional dataset
#'
#'   The variables specified for `by_vars`, `join_vars`, and `order` are
#'   expected.
#'
#' @permitted [dataset]
#'
#' @param by_vars Grouping variables
#'
#'   The specified variables are used for joining the input
#'   dataset (`dataset`) with the additional dataset (`dataset_add`).
#'
#'   `r roxygen_param_by_vars()`
#'
#' @permitted [var_list]
#'
#' @param order Order
#'
#'   The observations are ordered by the specified order if `join_type =
#'   "after"`, `join_type = "before"`, `first_cond_lower`, `first_cond_upper`,
#'   or `tmp_obs_nr_var` are specified.
#'
#'   `r roxygen_order_na_handling()`
#'
#' @permitted [var_list]
#'
#' @param new_var New variable
#'
#'   The specified variable is added to the input dataset.
#'
#' @permitted [var]
#'
#' @param tmp_obs_nr_var Temporary observation number
#'
#'   The specified variable is added to the input dataset (`dataset`) and the
#'   additional dataset (`dataset_add`). It is set to the observation number
#'   with respect to `order`. For each by group (`by_vars`) the observation
#'   number starts with `1`. If there is more than one record for specific
#'   values for `by_vars` and `order`, all records get the same observation
#'   number. By default, a warning (see `check_type`) is issued in this case.
#'   The variable can be used in the conditions (`filter_join`,
#'   `first_cond_upper`, `first_cond_lower`). It is not included in the output
#'   dataset. It can also be used to flag consecutive observations or the last
#'   observation (see last example below).
#'
#' @permitted [var]
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
#' @permitted [var_list]
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
#' @permitted [condition]
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
#' @permitted [condition]
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
#' @permitted [condition]
#'
#' @param check_type Check uniqueness?
#'
#'   If `"message"`, `"warning"`, or `"error"` is specified, the specified
#'   message is issued if the observations of the input dataset are not unique
#'   with respect to the by variables and the order.
#'
#' @permitted [msg_type]
#'
#' @param true_value Value of `new_var` for flagged observations
#'
#' @permitted [char_scalar]
#'
#' @param false_value Value of `new_var` for observations not flagged
#'
#' @permitted [char_scalar]
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
#'   For examples see the "Examples" section.
#'
#'   ## Step 4
#'
#'   The joined dataset is grouped by the observations from the input dataset
#'   and restricted to the observations fulfilling the condition specified by
#'   `filter_join`.
#'
#'   ## Step 5
#'
#'   The first observation of each group is selected.
#'
#'   ## Step 6
#'
#'   The variable specified by `new_var` is added to the input dataset. It is
#'   set to `true_value` for all observations which were selected in the
#'   previous step. For the other observations it is set to `false_value`.
#'
#' `r roxygen_save_memory()`
#'
#' @return The input dataset with the variable specified by `new_var` added.
#'
#' @keywords der_gen
#' @family der_gen
#'
#' @seealso [filter_joined()], [derive_vars_joined()]
#'
#' @export
#'
#' @examplesx
#' @caption Flag records considering other records (`filter_join`, `join_vars`)
#' @info In this example, records with a duration longer than 30 and where a
#'   COVID AE (`ACOVFL == "Y"`) occurred before or up to seven days after the
#'   record should be flagged. The condition for flagging the records is
#'   specified by the `filter_join` argument. Variables from the other records
#'   are referenced by variable names with the suffix `.join`. These variables
#'   have to be specified for the `join_vars` argument. As records before _and_
#'   after the current record should be considered, `join_type = "all"` is
#'   specified.
#' @code
#' library(tibble)
#'
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
#'   filter_join = ADURN > 30 & ACOVFL.join == "Y" & ADY.join <= ADY + 7
#' )
#'
#' @caption Considering only records after the current one (`join_type =
#'   "after"`, `true_value`, `false_value`)
#' @info In this example, records with `AVALC == "Y"` and `AVALC == "Y"` at a
#'   subsequent visit should be flagged. `join_type = "after"` is specified to
#'   consider only records after the current one. Please note that the `order`
#'   argument must be specified, as otherwise it is not possible to determine
#'   which records are after the current record.
#'
#'   Please note that a numeric flag is created here by specifying the
#'   `true_value` and the `false_value` argument.
#' @code
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
#'   new_var = CONFFLN,
#'   join_vars = exprs(AVALC, AVISITN),
#'   join_type = "after",
#'   order = exprs(AVISITN),
#'   filter_join = AVALC == "Y" & AVALC.join == "Y",
#'   true_value = 1,
#'   false_value = 0
#' )
#'
#' @caption Considering a range of records only (`first_cond_lower`, `first_cond_upper`)
#' @info Consider the following data.
#' @code
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
#' @info To flag `"0"` where all results from the first `"++"` before the `"0"`
#' up to the `"0"` (excluding the `"0"`) are `"+"` or `"++"` the
#' `first_cond_lower` argument and `join_type = "before"` are specified.
#' @code
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
#' @info To flag `"0"` where all results from the `"0"` (excluding the `"0"`) up
#' to the first `"++"` after the `"0"` are `"+"` or `"++"` the
#' `first_cond_upper` argument and `join_type = "after"` are specified.
#' @code
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
#'
#' @caption Considering only records up to a condition (`first_cond_upper`)
#' @info In this example from deriving confirmed response in oncology, the
#'   records with
#' - `AVALC == "CR"`,
#' - `AVALC == "CR"` at a subsequent visit,
#' - only `"CR"` or `"NE"` in between, and
#' - at most one `"NE"` in between
#'
#' should be flagged. The other records to be considered are restricted to those
#' up to the first occurrence of `"CR"` by specifying the `first_cond_upper`
#' argument. The `count_vals()` function is used to count the `"NE"`s for the
#' last condition.
#' @code
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
#' @caption Considering order of values (`min_cond()`, `max_cond()`)
#' @info In this example from deriving confirmed response in oncology, records
#'   with
#' - `AVALC == "PR"`,
#' - `AVALC == "CR"` or `AVALC == "PR"` at a subsequent visit at least 20 days later,
#' - only `"CR"`, `"PR"`, or `"NE"` in between,
#' - at most one `"NE"` in between, and
#' - `"CR"` is
#'   not followed by `"PR"`
#'
#' should be flagged. The last condition is realized by using `min_cond()` and
#' `max_cond()`, ensuring that the first occurrence of `"CR"` is after the last
#' occurrence of `"PR"`. The second call to `count_vals()` in the condition is
#' required to cover the case of no `"CR"`s (the `min_cond()` call returns `NA`
#' then).
#' @code
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
#' @caption Considering the order of records (`tmp_obs_nr_var`)
#' @info In this example, the records with `CRIT1FL == "Y"` at two consecutive
#'   visits or at the last visit should be flagged. A temporary order variable
#'   is created by specifying the `tmp_obs_nr_var` argument. Then it is used in
#'   `filter_join`. The temporary variable doesn't need to be specified for
#'   `join_vars`.
#' @code
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
#' @caption Flag each dose which is lower than the previous dose
#'   (`tmp_obs_nr_var`)
#' @code
#' ex <- tribble(
#'   ~USUBJID, ~EXSTDTM,           ~EXDOSE,
#'   "1",      "2024-01-01T08:00",       2,
#'   "1",      "2024-01-02T08:00",       4,
#'   "2",      "2024-01-01T08:30",       1,
#'   "2",      "2024-01-02T08:30",       4,
#'   "2",      "2024-01-03T08:30",       3,
#'   "2",      "2024-01-04T08:30",       2,
#'   "2",      "2024-01-05T08:30",       2
#' )
#'
#' derive_var_joined_exist_flag(
#'   ex,
#'   dataset_add = ex,
#'   by_vars = exprs(USUBJID),
#'   order = exprs(EXSTDTM),
#'   new_var = DOSREDFL,
#'   tmp_obs_nr_var = tmp_dose_nr,
#'   join_vars = exprs(EXDOSE),
#'   join_type = "before",
#'   filter_join = (
#'     tmp_dose_nr == tmp_dose_nr.join + 1 # Look only at adjacent doses
#'     & EXDOSE > 0 & EXDOSE.join > 0 # Both doses are valid
#'     & EXDOSE < EXDOSE.join # Dose is lower than previous
#'   )
#' )
#'
#' @caption Derive definitive deterioration flag
#' @info In this example a definitive deterioration flag should be derived as
#'   any deterioration (`CHGCAT1 = "Worsened"`) by parameter that is not
#'   followed by a non-deterioration. Please note that `join_type = "after"`
#'   can't by used here, as otherwise the last record wouldn't be flagged.
#' @code
#' adqs <- tribble(
#'   ~USUBJID, ~PARAMCD, ~ADY, ~CHGCAT1,
#'   "1",      "QS1",      10, "Improved",
#'   "1",      "QS1",      21, "Improved",
#'   "1",      "QS1",      23, "Improved",
#'   "1",      "QS2",      32, "Worsened",
#'   "1",      "QS2",      42, "Improved",
#'   "2",      "QS1",      11, "Worsened",
#'   "2",      "QS1",      24, "Worsened"
#' )
#'
#' derive_var_joined_exist_flag(
#'   adqs,
#'   dataset_add = adqs,
#'   new_var = DDETERFL,
#'   by_vars = exprs(USUBJID, PARAMCD),
#'   join_vars = exprs(CHGCAT1, ADY),
#'   join_type = "all",
#'   filter_join = all(CHGCAT1.join == "Worsened" | ADY > ADY.join)
#' )
#'
#' @caption Handling duplicates (`check_type`)
#' @info If the `order` argument is used, it is checked if the records are
#'   unique with respect to `by_vars` and `order`. Consider for example the
#'   derivation of `CONFFL` which flags records with `AVALC == "Y"` which are
#'   confirmed at a subsequent visit.
#' @code [expected_cnds = "duplicate_records"]
#' data <- tribble(
#'   ~USUBJID, ~AVISITN, ~ADY, ~AVALC,
#'   "1",      1,           1, "Y",
#'   "1",      2,           8, "N",
#'   "1",      3,          15, "Y",
#'   "1",      4,          22, "N",
#'   "2",      1,           1, "Y",
#'   "2",      2,           8, "Y",
#'   "2",      2,          10, "Y"
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
#'   filter_join = AVALC == "Y" & AVALC.join == "Y"
#' )
#' @info The records for `USUBJID == "2"` are not unique with respect to
#'   `USUBJID` and `AVISITN`. Thus a warning is issued. The duplicates can be
#'   accessed by calling `get_duplicates_dataset()`:
#' @code
#' get_duplicates_dataset()
#' @info In this example, confirmation is required at a subsequent _visit_.
#'   Please note that the first record for subject `"2"` at visit `2` is not
#'   flagged. Thus the warning can be suppressed by specifying `check_type =
#'   "none"`.
#' @code
#' derive_var_joined_exist_flag(
#'   data,
#'   dataset_add = data,
#'   by_vars = exprs(USUBJID),
#'   new_var = CONFFL,
#'   join_vars = exprs(AVALC, AVISITN),
#'   join_type = "after",
#'   order = exprs(AVISITN),
#'   filter_join = AVALC == "Y" & AVALC.join == "Y",
#'   check_type = "none"
#' )
derive_var_joined_exist_flag <- function(dataset,
                                         dataset_add,
                                         by_vars,
                                         order = NULL,
                                         new_var,
                                         tmp_obs_nr_var = NULL,
                                         join_vars,
                                         join_type,
                                         first_cond_lower = NULL,
                                         first_cond_upper = NULL,
                                         filter_add = NULL,
                                         filter_join,
                                         true_value = "Y",
                                         false_value = NA_character_,
                                         check_type = "warning") {
  new_var <- assert_symbol(enexpr(new_var))
  tmp_obs_nr_var <- assert_symbol(enexpr(tmp_obs_nr_var), optional = TRUE)
  first_cond_lower <- assert_filter_cond(enexpr(first_cond_lower), optional = TRUE)
  first_cond_upper <- assert_filter_cond(enexpr(first_cond_upper), optional = TRUE)
  filter_add <- assert_filter_cond(enexpr(filter_add), optional = TRUE)
  filter_join <- assert_filter_cond(enexpr(filter_join))
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
