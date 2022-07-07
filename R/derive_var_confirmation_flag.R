#' Derive Confirmation Flag
#'
#' Derive a flag which depends on subsequent observations of the dataset. For
#' example, flagging events which need to be confirmed by a second event.
#'
#' @param dataset Input dataset
#'
#'   The variables specified by the `by_vars` and `join_vars` parameter are
#'   expected.
#'
#' @param by_vars By variables
#'
#'   The specified variables are used as by variables for joining the input
#'   dataset with itself.
#'
#' @param join_vars Variables to keep from joined dataset
#'
#'   The variables needed from the subsequent observations should be specified
#'   for this parameter. The specified variables are added to the joined dataset
#'   with suffix ".join". For example to flag all observations with `AVALC ==
#'   "Y"` and `AVALC == "Y"` for at least one subsequent visit `join_vars =
#'   vars(AVALC, AVISITN)` and `filter = AVALC == "Y" & AVALC.join == "Y" &
#'   AVISITN < AVISITN.join` could be specified.
#'
#'   The `*.join` variables are not included in the output dataset.
#'
#' @param first_cond Condition for selecting range of data
#'
#'   If this argument is specified, the subsequent observations are restricted
#'   up to the first observation where the specified condition is fulfilled. If
#'   the condition is not fulfilled for any of the subsequent observations, no
#'   observations are considered, i.e., the observation is not flagged.
#'
#' @param order Order
#'
#'   The observations are ordered by the specified order.
#'
#' @param filter Condition for selecting observations
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
#'   *Default:* `"none"`
#'
#'   *Permitted Values:* `"none"`, `"warning"`, `"error"`
#'
#' @param true_value Value of `new_var` for flagged observations
#'
#'   *Default*: `"Y"`
#'
#' @param false_value Value of `new_var` for observations not flagged
#'
#'   *Default*: `NA_character_`
#'
#' @details
#'   The following steps are performed to produce the output dataset.
#'
#'   ## Step 1
#'
#'   The input dataset is joined with itself by the variables specified for
#'   `by_vars`. From the right hand side of the join only the variables
#'   specified for `join_vars` are kept. The suffix ".join" is added to these
#'   variables.
#'
#'   For example, for `by_vars = USUBJID`, `join_vars = vars(AVISITN, AVALC)` and input dataset
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
#'   The joined dataset is restricted to observations where the joined variables
#'   are at or after the other variables with respect to `order`.
#'
#'   The dataset from the example in the previous step with `order =
#'   vars(AVISITN)` is restricted to
#'
#'   ```{r eval=FALSE}
#'   A tibble: 4 x 6
#'   USUBJID AVISITN AVALC  AVAL AVISITN.join AVALC.join
#'   <chr>     <dbl> <chr> <dbl>        <dbl> <chr>
#'   1             1 Y         1            1 Y
#'   1             1 Y         1            2 N
#'   1             2 N         0            2 N
#'   ```
#'
#'   ## Step 3
#'
#'   If `first_cond` is specified, for each observation of the input dataset the
#'   joined dataset is restricted to observations up to the first observation
#'   where `first_cond` is fulfilled (the observation fulfilling the condition
#'   is included). If for an observation of the input dataset the condition is
#'   not fulfilled, the observation is removed.
#'
#'   ## Step 4
#'
#'   The joined dataset is grouped by the observations from the input dataset
#'   and restricted to the observations fulfilling the condition specified by
#'   `filter`.
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
#' @author Stefan Bundfuss
#'
#' @keywords adam derivation
#'
#' @seealso [filter_confirmation()]
#'
#' @export
#'
#' @examples
#' library(tibble)
#' library(admiral)
#'
#' # flag observations with a duration longer than 30 and
#' # at, after, or up to 7 days before a COVID AE (ACOVFL == "Y")
#' adae <- tribble(
#'   ~USUBJID, ~ADY, ~ACOVFL, ~ADURN,
#'   "1",        10, "N",          1,
#'   "1",        21, "N",         50,
#'   "1",        23, "Y",         14,
#'   "1",        42, "N",         20,
#'   "2",        11, "Y",         13,
#'   "2",        23, "N",          2,
#'   "3",        13, "Y",         12,
#'   "4",        14, "N",         32,
#'   "4",        21, "N",         41
#' )
#'
#' derive_var_confirmation_flag(
#'   adae,
#'   new_var = ALCOVFL,
#'   by_vars = vars(USUBJID),
#'   join_vars = vars(ACOVFL, ADY),
#'   order = vars(ADY),
#'   filter = ADURN > 30 & ACOVFL.join == "Y" & ADY >= ADY.join - 7
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
#' derive_var_confirmation_flag(
#'   data,
#'   by_vars = vars(USUBJID),
#'   new_var = CONFFL,
#'   join_vars = vars(AVALC, AVISITN),
#'   order = vars(AVISITN),
#'   filter = AVALC == "Y" & AVALC.join == "Y" & AVISITN < AVISITN.join
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
#' derive_var_confirmation_flag(
#'   data,
#'   by_vars = vars(USUBJID),
#'   join_vars = vars(AVALC),
#'   order = vars(AVISITN),
#'   new_var = CONFFL,
#'   first_cond = AVALC.join == "CR",
#'   filter = AVALC == "CR" & all(AVALC.join %in% c("CR", "NE")) &
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
#' derive_var_confirmation_flag(
#'   data,
#'   by_vars = vars(USUBJID),
#'   join_vars = vars(AVALC, ADY),
#'   order = vars(ADY),
#'   new_var = CONFFL,
#'   first_cond = AVALC.join %in% c("CR", "PR") & ADY.join - ADY >= 20,
#'   filter = AVALC == "PR" &
#'     all(AVALC.join %in% c("CR", "PR", "NE")) &
#'     count_vals(var = AVALC.join, val = "NE") <= 1 &
#'     (
#'       min_cond(var = ADY.join, cond = AVALC.join == "CR") >
#'         max_cond(var = ADY.join, cond = AVALC.join == "PR") |
#'         count_vals(var = AVALC.join, val = "CR") == 0
#'     )
#' )
derive_var_confirmation_flag <- function(dataset,
                                         by_vars,
                                         order,
                                         new_var,
                                         join_vars,
                                         first_cond = NULL,
                                         filter,
                                         true_value = "Y",
                                         false_value = NA_character_,
                                         check_type = "warning") {
  new_var <- assert_symbol(enquo(new_var))
  first_cond <- assert_filter_cond(enquo(first_cond), optional = TRUE)
  filter <- assert_filter_cond(enquo(filter))

  data <- derive_var_obs_number(
    dataset,
    new_var = tmp_obs_nr_var_conf_flag
  )

  data_filtered <- admiralonco::filter_confirmation(
    data,
    by_vars = by_vars,
    order = order,
    join_vars = join_vars,
    first_cond = !!first_cond,
    filter = !!filter,
    check_type = check_type
  )

  derive_var_merged_exist_flag(
    data,
    dataset_add = data_filtered,
    by_vars = vars(tmp_obs_nr_var_conf_flag),
    new_var = !!new_var,
    condition = TRUE,
    true_value = true_value,
    false_value = false_value,
    missing_value = false_value
  ) %>%
    select(-tmp_obs_nr_var_conf_flag)
}
