#' Add Variables from an Additional Dataset Based on Conditions from Both
#' Datasets
#'
#' The function adds variables from an additional dataset to the input dataset.
#' The selection of the observations from the additional dataset can depend on
#' variables from both datasets. For example, add the lowest value (nadir)
#' before the current observation.
#'
#' @param dataset
#'
#' `r roxygen_param_dataset(expected_vars = c("by_vars"))`
#'
#' @permitted [dataset]
#'
#' @param dataset_add Additional dataset
#'
#'   The variables specified by the `by_vars`, the `new_vars`, the `join_vars`,
#'   and the `order` argument are expected.
#'
#' @permitted [dataset]
#'
#' @param by_vars Grouping variables
#'
#'   The two datasets are joined by the specified variables.
#'
#'   `r roxygen_param_by_vars(rename = TRUE)`
#'
#' @permitted [var_list]
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
#' @param tmp_obs_nr_var Temporary observation number
#'
#'   The specified variable is added to the input dataset (`dataset`) and the
#'   additional dataset (`dataset_add`). It is set to the observation number
#'   with respect to `order`. For each by group (`by_vars`) the observation
#'   number starts with `1`. If there is more than one record for specific
#'   values for `by_vars` and `order`, all records get the same observation
#'   number. By default, a warning (see `check_type`) is issued in this case.
#'   The variable can be used in the conditions (`filter_join`,
#'   `first_cond_upper`, `first_cond_lower`). It can also be used to select
#'   consecutive observations or the last observation.
#'
#'   The variable is not included in the output dataset. To include it specify
#'   it for `new_vars`.
#'
#' @permitted [var]
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
#' @permitted [var_list]
#'
#' @param first_cond_lower Condition for selecting range of data (before)
#'
#'   If this argument is specified, the other observations are restricted from
#'   the last observation before the current observation where the specified
#'   condition is fulfilled up to the current observation. If the condition is
#'   not fulfilled for any of the other observations, no observations are
#'   considered.
#'
#'   This argument should be specified if `filter_join` contains summary
#'   functions which should not apply to all observations but only from a
#'   certain observation before the current observation up to the current
#'   observation. For an example, see the "Examples" section below.
#'
#' @permitted [condition]
#'
#' @param first_cond_upper Condition for selecting range of data (after)
#'
#'   If this argument is specified, the other observations are restricted up to
#'   the first observation where the specified condition is fulfilled. If the
#'   condition is not fulfilled for any of the other observations, no
#'   observations are considered.
#'
#'   This argument should be specified if `filter_join` contains summary
#'   functions which should not apply to all observations but only up to the
#'   confirmation assessment. For an example, see the "Examples" section below.
#'
#' @permitted [condition]
#'
#' @param filter_join Filter for the joined dataset
#'
#'   The specified condition is applied to the joined dataset. Therefore
#'   variables from both datasets `dataset` and `dataset_add` can be used.
#'
#'   Variables created by `order` or `new_vars` arguments can be used in the
#'   condition.
#'
#'   The condition can include summary functions like `all()` or `any()`. The
#'   joined dataset is grouped by the original observations.
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
#' @param check_type Check uniqueness?
#'
#'   If `"message"`, `"warning"` or `"error"` is specified, the specified
#'   message is issued if the observations of the (restricted) joined dataset
#'   are not unique with respect to the by variables and the order.
#'
#'   This argument is ignored if `order` is not specified. In this case an error
#'   is issued independent of `check_type` if the restricted joined dataset
#'   contains more than one observation for any of the observations of the input
#'   dataset.
#'
#' @permitted [msg_type]
#'
#' @inheritParams get_joined_data
#' @inheritParams derive_vars_merged
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
#' 1. If `first_cond_lower` is specified, for each observation of the input
#'     dataset the joined dataset is restricted to observations from the first
#'     observation where `first_cond_lower` is fulfilled (the observation fulfilling
#'     the condition is included) up to the observation of the input dataset. If for
#'     an observation of the input dataset the condition is not fulfilled, the
#'     observation is removed.
#'
#'     If `first_cond_upper` is specified, for each observation of the input
#'     dataset the joined dataset is restricted to observations up to the first
#'     observation where `first_cond_upper` is fulfilled (the observation
#'     fulfilling the condition is included). If for an observation of the input
#'     dataset the condition is not fulfilled, the observation is removed.
#'
#'     For an example, see the "Examples" section below.
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
#' `r roxygen_save_memory()`
#'
#' @return The output dataset contains all observations and variables of the
#'   input dataset and additionally the variables specified for `new_vars` from
#'   the additional dataset (`dataset_add`).
#'
#' @seealso [derive_var_joined_exist_flag()], [filter_joined()]
#'
#' @keywords der_gen
#' @family der_gen
#'
#' @export
#'
#' @examplesx
#'
#' @caption Note on usage versus `derive_vars_merged()`
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
#' @caption Basic join based on a generic time window (`filter_join`)
#' @info Derive a visit based on where the study day falls according to a
#'   scheduled set of time windows.
#'
#' - The `filter_join` argument here can check conditions using variables from
#'   both the `dataset` and `dataset_add`, so the study day is compared to the
#'   start and end of the time window.
#' - As no grouping variables are assigned using the `by_vars` argument, a full
#'   join is performed keeping all variables from `dataset_add`.
#' @code
#' library(tibble)
#' library(lubridate)
#' library(dplyr, warn.conflicts = FALSE)
#' library(tidyr, warn.conflicts = FALSE)
#'
#' adbds <- tribble(
#'   ~USUBJID, ~ADY, ~AVAL,
#'   "1",       -33,    11,
#'   "1",        -7,    10,
#'   "1",         1,    12,
#'   "1",         8,    12,
#'   "1",        15,     9,
#'   "1",        20,    14,
#'   "1",        24,    12,
#'   "2",        -1,    13,
#'   "2",        13,     8
#' ) %>%
#'   mutate(STUDYID = "AB42")
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
#'   join_type = "all",
#'   filter_join = AWLO <= ADY & ADY <= AWHI
#' ) %>%
#'   select(USUBJID, ADY, AWLO, AWHI, AVISIT)
#'
#' @caption Join only the lowest/highest value occurring within a condition (`filter_join`,
#'   `order` and `mode`)
#' @info Derive the nadir value for each observation (i.e. the lowest value
#'   occurring before) by subject.
#'
#' - Note how `dataset` and `dataset_add` are the same here, so we are joining
#'   a dataset with itself. This enables us to compare records within the dataset
#'   to each other.
#' - Now we use `by_vars` as we only want to perform the join by subject.
#' - To find the lowest value we use the `order` and `mode` arguments.
#' - We subsequently need to check `ADY` to only check assessments occurring
#'   before. As this is not included in `by_vars` or `order`, we have to ensure
#'   it also gets joined by adding to `join_vars`. Then in `filter_join` note
#'   how `ADY.join < ADY` is used as the same variable exists in both datasets,
#'   so the version from `dataset_add` has `.join` added.
#' - According to the `AVAL` sort order used there could be duplicates (e.g. see
#'   subject `"1"` records at day 1 and 8), but given we only need to join `AVAL`
#'   itself here it doesn't actually matter to us which exact record is taken.
#'   So, in this example, we silence the uniqueness check by using
#'   `check_type = "none"`.
#' @code
#' derive_vars_joined(
#'   adbds,
#'   dataset_add = adbds,
#'   by_vars = exprs(STUDYID, USUBJID),
#'   order = exprs(AVAL),
#'   new_vars = exprs(NADIR = AVAL),
#'   join_vars = exprs(ADY),
#'   join_type = "all",
#'   filter_join = ADY.join < ADY,
#'   mode = "first",
#'   check_type = "none"
#' ) %>%
#'   select(USUBJID, ADY, AVAL, NADIR)
#'
#' @caption Filtering which records are joined from the additional dataset (`filter_add`)
#' @info Imagine we wanted to achieve the same as above, but we now want to derive
#'   this allowing only post-baseline values to be possible for the nadir.
#'
#' - The `filter_add` argument can be used here as we only need to restrict the
#'   source data from `dataset_add`.
#' @code
#' derive_vars_joined(
#'   adbds,
#'   dataset_add = adbds,
#'   by_vars = exprs(STUDYID, USUBJID),
#'   order = exprs(AVAL),
#'   new_vars = exprs(NADIR = AVAL),
#'   join_vars = exprs(ADY),
#'   join_type = "all",
#'   filter_add = ADY > 0,
#'   filter_join = ADY.join < ADY,
#'   mode = "first",
#'   check_type = "none"
#' ) %>%
#'   select(USUBJID, ADY, AVAL, NADIR)
#'
#' @caption Combining all of the above examples
#' @info Using all of the arguments demonstrated above, here is a more complex
#'   example to add to `ADAE` the highest hemoglobin value occurring within two weeks
#'   before each adverse event. Also join the day it occurred, taking the earliest
#'   occurrence if more than one assessment with the same value.
#'
#' - Note how we used `mode = "last"` to get the highest lab value, but then as we
#'   wanted the earliest occurrence if more than one it means we need to add
#'   `desc(ADY)` to `order`. i.e. the last day when in descending order is the first.
#' @code
#' adae <- tribble(
#'   ~USUBJID, ~ASTDY,
#'   "1",           3,
#'   "1",          22,
#'   "2",           2
#' ) %>%
#'   mutate(STUDYID = "AB42")
#'
#' adlb <- tribble(
#'   ~USUBJID, ~PARAMCD, ~ADY, ~AVAL,
#'   "1",      "HGB",       1,   8.5,
#'   "1",      "HGB",       3,   7.9,
#'   "1",      "HGB",       5,   8.9,
#'   "1",      "HGB",       8,   8.0,
#'   "1",      "HGB",       9,   8.0,
#'   "1",      "HGB",      16,   7.4,
#'   "1",      "ALB",       1,    42,
#' ) %>%
#'   mutate(STUDYID = "AB42")
#'
#' derive_vars_joined(
#'   adae,
#'   dataset_add = adlb,
#'   by_vars = exprs(STUDYID, USUBJID),
#'   order = exprs(AVAL, desc(ADY)),
#'   new_vars = exprs(HGB_MAX = AVAL, HGB_DY = ADY),
#'   join_type = "all",
#'   filter_add = PARAMCD == "HGB",
#'   filter_join = ASTDY - 14 <= ADY & ADY <= ASTDY,
#'   mode = "last"
#' ) %>%
#'   select(USUBJID, ASTDY, HGB_MAX, HGB_DY)
#'
#' @caption Compute values in `new_vars` and `order`
#' @info Add to `ADAE` the number of days since the last dose of treatment, plus
#'   1 day. If the dose occurs on the same day as the AE then include it as the
#'   last dose.
#'
#' - In the `new_vars` argument, other functions can be utilized to modify the
#'   joined values using variables from both `dataset` and `dataset_add`.
#'   For example, in the below case we want to calculate the number of days
#'   between the AE and the last dose using `compute_duration()`. This function
#'   includes the plus 1 day as default.
#' - Also note how in this example `EXSDT` is created via the `order` argument
#'   and then used for `new_vars`, `filter_add` and `filter_join`.
#' - The reason to use `join_type = "all"` here instead of `"before"` is that we
#'   want to include any dose occurring on the same day as the AE, hence the
#'   `filter_join = EXSDT <= ASTDT`. Whereas using `join_type = "before"`
#'   would have resulted in the condition `EXSDT < ASTDT`. See the next example
#'   instead for `join_type = "before"`.
#' @code
#' adae <- tribble(
#'   ~USUBJID, ~ASTDT,
#'   "1",      "2020-02-02",
#'   "1",      "2020-02-04",
#'   "2",      "2021-01-08"
#' ) %>%
#'   mutate(
#'     ASTDT = ymd(ASTDT),
#'     STUDYID = "AB42"
#'   )
#'
#' ex <- tribble(
#'   ~USUBJID, ~EXSDTC,
#'   "1",      "2020-01-10",
#'   "1",      "2020-01",
#'   "1",      "2020-01-20",
#'   "1",      "2020-02-03",
#'   "2",      "2021-01-05"
#' ) %>%
#'   mutate(STUDYID = "AB42")
#'
#' derive_vars_joined(
#'   adae,
#'   dataset_add = ex,
#'   by_vars = exprs(STUDYID, USUBJID),
#'   order = exprs(EXSDT = convert_dtc_to_dt(EXSDTC)),
#'   join_type = "all",
#'   new_vars = exprs(LDRELD = compute_duration(
#'     start_date = EXSDT, end_date = ASTDT
#'   )),
#'   filter_add = !is.na(EXSDT),
#'   filter_join = EXSDT <= ASTDT,
#'   mode = "last"
#' ) %>%
#'   select(USUBJID, ASTDT, LDRELD)
#'
#' @caption Join records occurring before a condition (`join_type = "before"`)
#' @info In an arbitrary dataset where subjects have values of `"0"`, `"-"`, `"+"`
#'   or `"++"`, for any value of `"0"` derive the last occurring `"++"` day that
#'   occurs before the `"0"`.
#'
#' - The `AVAL.join == "++"` in `filter_join`, along with `order` and `mode`
#'   taking the last day, identifies the target records to join from
#'   `dataset_add` for each observation of `dataset`.
#' - Then `join_type = "before"` is now used instead of `join_type = "all"`.
#'   This is because we only want to join the records occurring before the
#'   current observation in `dataset`. Including `AVAL == "0"` in `filter_join`
#'   ensures here that we only populate the new variable for records with
#'   `AVAL == "0"` in our `dataset`.
#' @code
#' myd <- tribble(
#'   ~USUBJID, ~ADY, ~AVAL,
#'   "1",         1, "++",
#'   "1",         2, "-",
#'   "1",         3, "0",
#'   "1",         4, "+",
#'   "1",         5, "++",
#'   "1",         6, "-",
#'   "2",         1, "-",
#'   "2",         2, "++",
#'   "2",         3, "+",
#'   "2",         4, "0",
#'   "2",         5, "-",
#'   "2",         6, "++",
#'   "2",         7, "0"
#' ) %>%
#'   mutate(STUDYID = "AB42")
#'
#' derive_vars_joined(
#'   myd,
#'   dataset_add = myd,
#'   by_vars = exprs(STUDYID, USUBJID),
#'   order = exprs(ADY),
#'   mode = "last",
#'   new_vars = exprs(PREVPLDY = ADY),
#'   join_vars = exprs(AVAL),
#'   join_type = "before",
#'   filter_join = AVAL == "0" & AVAL.join == "++"
#' ) %>%
#'   select(USUBJID, ADY, AVAL, PREVPLDY)
#'
#' @caption Join records occurring before a condition and checking all values in
#'   between (`first_cond_lower`, `join_type` and `filter_join`)
#' @info In the same example as above, now additionally check that in between the
#'   `"++"` and the `"0"` all results must be either `"+"` or `"++"`.
#'
#' - Firstly, `first_cond_lower = AVAL.join == "++"` is used so that for each
#'   observation of `dataset` the joined records from `dataset_add` are restricted
#'   to only include from the last occurring `"++"` before. This is necessary
#'   because of the use of a summary function in `filter_join` only on a subset
#'   of the joined observations as explained below.
#' - The `filter_join` condition used here now includes `all(AVAL.join %in% c("+", "++"))`
#'   to further restrict the joined records from `dataset_add` to only where all
#'   the values are either `"+"` or `"++"`.
#' - The `order` and `mode` arguments ensure only the day of the `"++"` value
#'   is joined. For example, for subject `"2"` it selects the day 2 record
#'   instead of day 3, by using `"first"`.
#' @code
#' derive_vars_joined(
#'   myd,
#'   dataset_add = myd,
#'   by_vars = exprs(STUDYID, USUBJID),
#'   order = exprs(ADY),
#'   mode = "first",
#'   new_vars = exprs(PREVPLDY = ADY),
#'   join_vars = exprs(AVAL),
#'   join_type = "before",
#'   first_cond_lower = AVAL.join == "++",
#'   filter_join = AVAL == "0" & all(AVAL.join %in% c("+", "++"))
#' ) %>%
#'   select(USUBJID, ADY, AVAL, PREVPLDY)
#'
#' @caption Join records occurring after a condition checking all values in between
#'   (`first_cond_upper`, `join_type` and `filter_join`)
#' @info Similar to the above, now derive the first `"++"` day after any `"0"`
#'   where all results in between are either `"+"` or `"++"`.
#'
#' - Note how the main difference here is the use of `join_type = "after"`,
#'   `mode = "last"` and the `first_cond_upper` argument, instead of
#'   `first_cond_lower`.
#' @code
#' derive_vars_joined(
#'   myd,
#'   dataset_add = myd,
#'   by_vars = exprs(STUDYID, USUBJID),
#'   order = exprs(ADY),
#'   mode = "last",
#'   new_vars = exprs(NEXTPLDY = ADY),
#'   join_vars = exprs(AVAL),
#'   join_type = "after",
#'   first_cond_upper = AVAL.join == "++",
#'   filter_join = AVAL == "0" & all(AVAL.join %in% c("+", "++"))
#' ) %>%
#'   select(USUBJID, ADY, AVAL, NEXTPLDY)
#'
#' @caption Join a value from the next occurring record (`join_type = "after"`)
#' @info Add the value from the next occurring record as a new variable.
#'
#' - The `join_type = "after"` here essentially acts as a lag to join variables from
#'   the next occurring record, and `mode = "first"` selects the first of these.
#' @code
#' derive_vars_joined(
#'   myd,
#'   dataset_add = myd,
#'   by_vars = exprs(STUDYID, USUBJID),
#'   order = exprs(ADY),
#'   mode = "first",
#'   new_vars = exprs(NEXTVAL = AVAL),
#'   join_vars = exprs(AVAL),
#'   join_type = "after"
#' ) %>%
#'   select(USUBJID, ADY, AVAL, NEXTVAL)
#'
#' @caption Join records after a condition occurring in consecutive visits
#'   (`tmp_obs_nr_var`, `join_type` and `filter_join`)
#' @info Find the last occurring value on any of the next 3 unique visit days.
#'
#' - The `tmp_obs_nr_var` argument can be useful as shown here to help pick out
#'   records happening before or after with respect to `order`, as you can see
#'   in the `filter_join`.
#' @code
#' derive_vars_joined(
#'   myd,
#'   dataset_add = myd,
#'   by_vars = exprs(STUDYID, USUBJID),
#'   order = exprs(ADY),
#'   mode = "last",
#'   new_vars = exprs(NEXTVAL = AVAL),
#'   tmp_obs_nr_var = tmp_obs_nr,
#'   join_vars = exprs(AVAL),
#'   join_type = "after",
#'   filter_join = tmp_obs_nr + 3 >= tmp_obs_nr.join
#' ) %>%
#'   select(USUBJID, ADY, AVAL, NEXTVAL)
#'
#' @caption Derive period variables (`APERIOD`, `APERSDT`, `APEREDT`)
#' @info Create a period reference dataset from `ADSL` and join this with `ADAE`
#'   to identify within which period each AE occurred.
#' @code
#' adsl <- tribble(
#'   ~USUBJID, ~AP01SDT,     ~AP01EDT,     ~AP02SDT,     ~AP02EDT,
#'   "1",      "2021-01-04", "2021-02-06", "2021-02-07", "2021-03-07",
#'   "2",      "2021-02-02", "2021-03-02", "2021-03-03", "2021-04-01"
#' ) %>%
#'   mutate(across(ends_with("DT"), ymd)) %>%
#'   mutate(STUDYID = "AB42")
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
#'     STUDYID = "AB42"
#'   )
#'
#' derive_vars_joined(
#'   adae,
#'   dataset_add = period_ref,
#'   by_vars = exprs(STUDYID, USUBJID),
#'   join_vars = exprs(APERSDT, APEREDT),
#'   join_type = "all",
#'   filter_join = APERSDT <= ASTDT & ASTDT <= APEREDT
#' ) %>%
#'   select(USUBJID, ASTDT, APERSDT, APEREDT, APERIOD)
#'
#' @caption Further examples
#' @info Further example usages of this function can be found in the
#'   [Generic Derivations vignette](../articles/generic.html).
#'
#'   Equivalent examples for using the `exist_flag`, `true_value`, `false_value`,
#'   `missing_values` and `check_type` arguments can be found in `derive_vars_merged()`.
derive_vars_joined <- function(dataset,
                               dataset_add,
                               by_vars = NULL,
                               order = NULL,
                               new_vars = NULL,
                               tmp_obs_nr_var = NULL,
                               join_vars = NULL,
                               join_type,
                               filter_add = NULL,
                               first_cond_lower = NULL,
                               first_cond_upper = NULL,
                               filter_join = NULL,
                               mode = NULL,
                               exist_flag = NULL,
                               true_value = "Y",
                               false_value = NA_character_,
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

  tmp_obs_nr_var <- assert_symbol(enexpr(tmp_obs_nr_var), optional = TRUE)
  filter_add <- assert_filter_cond(enexpr(filter_add), optional = TRUE)
  first_cond_lower <- assert_filter_cond(enexpr(first_cond_lower), optional = TRUE)
  first_cond_upper <- assert_filter_cond(enexpr(first_cond_upper), optional = TRUE)
  filter_join <- assert_filter_cond(enexpr(filter_join), optional = TRUE)
  exist_flag <- assert_symbol(enexpr(exist_flag), optional = TRUE)

  if (is.null(new_vars)) {
    new_vars <- setdiff(chr2vars(colnames(dataset_add)), by_vars)
  }
  preexisting_vars <- chr2vars(colnames(dataset))
  preexisting_vars_no_by_vars <- preexisting_vars[which(!(preexisting_vars %in% by_vars))]
  duplicates <- intersect(replace_values_by_names(new_vars), preexisting_vars_no_by_vars)
  if (length(duplicates) > 0) {
    cli_abort(
      paste(
        "The variable{?s} {.var {duplicates}} in {.arg dataset_add} ha{?s/ve} naming",
        "conflicts with {.arg dataset}, please make the appropriate modifications",
        "to {.arg new_vars}."
      )
    )
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

  data_joined <- get_joined_data(
    data,
    dataset_add = dataset_add,
    by_vars = by_vars,
    join_vars = expr_c(
      join_vars,
      intersect(unname(extract_vars(new_vars)), chr2vars(colnames(dataset_add)))
    ),
    join_type = join_type,
    first_cond_lower = !!first_cond_lower,
    first_cond_upper = !!first_cond_upper,
    order = order,
    tmp_obs_nr_var = !!tmp_obs_nr_var,
    filter_add = !!filter_add,
    filter_join = !!filter_join,
    check_type = check_type
  )

  common_vars <-
    chr2vars(setdiff(intersect(colnames(data), colnames(dataset_add)), vars2chr(by_vars)))
  if (!is.null(order)) {
    data_joined <- filter_extreme(
      data_joined,
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
      dataset_add = data_joined,
      by_vars = exprs(!!!by_vars_left, !!tmp_obs_nr),
      new_vars = add_suffix_to_vars(new_vars, vars = common_vars, suffix = ".join"),
      missing_values = missing_values,
      check_type = check_type,
      exist_flag = !!exist_flag,
      true_value = true_value,
      false_value = false_value,
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

#' Join Data for "joined" functions
#'
#' The helper function joins the data for the "joined" functions. All `.join`
#' variables are included in the output dataset.
#'
#' @param dataset
#'
#' `r roxygen_param_dataset(expected_vars = c("by_vars"))`
#'
#' @permitted [dataset]
#'
#' @param dataset_add Additional dataset
#'
#'   The variables specified by the `by_vars`, the `new_vars`, the `join_vars`,
#'   and the `order` argument are expected.
#'
#' @permitted [dataset]
#'
#' @param by_vars Grouping variables
#'
#'   The two datasets are joined by the specified variables.
#'
#'   `r roxygen_param_by_vars(rename = TRUE)`
#'
#' @permitted [var_list]
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
#'   `r roxygen_order_na_handling()`
#'
#' @permitted [var_list]
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
#' @permitted [var_list]
#'
#' @param join_type Observations to keep after joining
#'
#'   The argument determines which of the joined observations are kept with
#'   respect to the original observation. For example, if `join_type = "after"`
#'   is specified all observations after the original observations are kept.
#'
#'   For example for confirmed response or BOR in the oncology setting or
#'   confirmed deterioration in questionnaires the confirmatory assessment must
#'   be after the assessment. Thus `join_type = "after"` could be used.
#'
#'   Whereas, sometimes you might allow for confirmatory observations to occur
#'   prior to the observation. For example, to identify AEs occurring on or
#'   after seven days before a COVID AE. Thus `join_type = "all"` could be used.
#'
#' @permitted [join_type]
#'
#' @param tmp_obs_nr_var Temporary observation number
#'
#'   The specified variable is added to the input dataset (`dataset`) and the
#'   additional dataset (`dataset_add`). It is set to the observation number
#'   with respect to `order`. For each by group (`by_vars`) the observation
#'   number starts with `1`. The variable can be used in the conditions
#'   (`filter_join`, `first_cond_upper`, `first_cond_lower`). It can also be
#'   used to select consecutive observations or the last observation.
#'
#' @permitted [var]
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
#'   The condition can include summary functions like `all()` or `any()`. The
#'   additional dataset is grouped by the by variables (`by_vars`).
#'
#' @permitted [condition]
#'
#' @param first_cond_lower Condition for selecting range of data (before)
#'
#'   If this argument is specified, the other observations are restricted from
#'   the first observation before the current observation where the specified
#'   condition is fulfilled up to the current observation. If the condition is
#'   not fulfilled for any of the other observations, no observations are
#'   considered, i.e., the observation is not flagged.
#'
#'   This argument should be specified if `filter_join` contains summary
#'   functions which should not apply to all observations but only from a
#'   certain observation before the current observation up to the current
#'   observation.
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
#'   This argument should be specified if `filter_join` contains summary
#'   functions which should not apply to all observations but only up to the
#'   confirmation assessment.
#'
#' @permitted [condition]
#'
#' @param filter_join Filter for the joined dataset
#'
#'   The specified condition is applied to the joined dataset. Therefore
#'   variables from both datasets `dataset` and `dataset_add` can be used.
#'
#'   Variables created by `order` or `new_vars` arguments can be used in the
#'   condition.
#'
#'   The condition can include summary functions like `all()` or `any()`. The
#'   joined dataset is grouped by the original observations.
#'
#' @permitted [condition]
#'
#' @param check_type Check uniqueness?
#'
#'   If `"message"`, `"warning"` or `"error"` is specified, the specified
#'   message is issued if the observations of the (restricted) joined dataset
#'   are not unique with respect to the by variables and the order.
#'
#'   This argument is ignored if `order` is not specified. In this case an error
#'   is issued independent of `check_type` if the restricted joined dataset
#'   contains more than one observation for any of the observations of the input
#'   dataset.
#'
#' @permitted [msg_type]
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
#' @keywords internal
get_joined_data <- function(dataset,
                            dataset_add,
                            by_vars = NULL,
                            join_vars = NULL,
                            join_type,
                            first_cond_lower = NULL,
                            first_cond_upper = NULL,
                            order = NULL,
                            tmp_obs_nr_var = NULL,
                            filter_add = NULL,
                            filter_join = NULL,
                            check_type = "warning") {
  # Check input arguments
  assert_vars(by_vars, optional = TRUE)
  by_vars_left <- replace_values_by_names(by_vars)
  assert_expr_list(join_vars, optional = TRUE)
  join_type <-
    assert_character_scalar(
      join_type,
      values = c("before", "after", "all"),
      case_sensitive = FALSE
    )
  first_cond_lower <- assert_filter_cond(enexpr(first_cond_lower), optional = TRUE)
  first_cond_upper <- assert_filter_cond(enexpr(first_cond_upper), optional = TRUE)
  assert_expr_list(order, optional = TRUE)
  tmp_obs_nr_var <- assert_symbol(enexpr(tmp_obs_nr_var), optional = TRUE)
  filter_add <- assert_filter_cond(enexpr(filter_add), optional = TRUE)
  filter_join <- assert_filter_cond(enexpr(filter_join), optional = TRUE)
  check_type <-
    assert_character_scalar(
      check_type,
      values = c("none", "warning", "error"),
      case_sensitive = FALSE
    )
  any_first_cond <- !is.null(first_cond_lower) || !is.null(first_cond_upper)
  if (join_type != "all" || any_first_cond) {
    dataset_order_vars <- extract_vars(order)
  } else {
    dataset_order_vars <- NULL
  }

  assert_data_frame(
    dataset,
    required_vars = expr_c(by_vars_left, dataset_order_vars)
  )

  assert_data_frame(
    dataset_add,
    required_vars = expr_c(
      by_vars,
      extract_vars(order),
      setdiff(extract_vars(join_vars), replace_values_by_names(order))
    )
  )

  # number observations of the input dataset to get a unique key
  # (by_vars and tmp_obs_nr_left), it is used later to apply filter_join
  tmp_obs_nr_left <- get_new_tmp_var(dataset, prefix = "tmp_obs_nr_left")
  data <- dataset %>%
    derive_var_obs_number(
      new_var = !!tmp_obs_nr_left,
      by_vars = by_vars_left,
      check_type = "none"
    )

  # derive variables defined by order and join_vars and restrict the additional
  # dataset
  data_add <- dataset_add %>%
    group_by(!!!by_vars) %>%
    mutate(!!!order, !!!join_vars) %>%
    filter_if(filter_add) %>%
    ungroup()

  # number groups with respect to by_vars and order in the input dataset and the
  # additional dataset for relation of records, e.g., join_type = before|after,
  # first_cond_lower, first_cond_upper
  tmp_obs_nr_var_join <- NULL
  if (join_type != "all" || any_first_cond || !is.null(tmp_obs_nr_var)) {
    if (is.null(tmp_obs_nr_var)) {
      tmp_obs_nr_var <- get_new_tmp_var(dataset, prefix = "tmp_obs_nr")
      tmp_obs_nr_var_join <- paste0(as_name(tmp_obs_nr_var), ".join")
    }

    if (check_type != "none") {
      # check if order results in unique values in dataset and dataset_add
      signal_duplicate_records(
        dataset = data,
        by_vars = c(by_vars, order),
        msg = paste(
          "Dataset {.arg dataset} contains duplicate records with respect to",
          "{.var {replace_values_by_names(by_vars)}}"
        ),
        cnd_type = check_type
      )
      signal_duplicate_records(
        dataset = data_add,
        by_vars = c(by_vars, order),
        msg = paste(
          "Dataset {.arg dataset_add} contains duplicate records with respect to",
          "{.var {replace_values_by_names(by_vars)}}"
        ),
        cnd_type = check_type
      )
    }

    # derive variables defined by order, then these can be used later, e.g., in
    # filter_join
    data <- data %>%
      mutate(!!!order)

    # if order contains unnamed expressions like floor(event_nr), the previous
    # mutate() creates a variable named `floor(event_nr)`. We need to use this
    # variable name instead of the expression in the following calls.
    order_vars <- map(
      replace_values_by_names(order),
      function(x) sym(as_label(x))
    )

    groups <- bind_rows(
      select(data, !!!by_vars, !!!order_vars),
      select(data_add, !!!by_vars, !!!order_vars)
    ) %>%
      distinct() %>%
      derive_var_obs_number(
        new_var = !!tmp_obs_nr_var,
        by_vars = by_vars,
        order = order_vars,
        check_type = "none"
      )

    data <- data %>%
      derive_vars_merged(
        dataset_add = groups,
        by_vars = exprs(!!!by_vars, !!!order_vars),
        check_type = "none"
      )

    data_add <- data_add %>%
      derive_vars_merged(
        dataset_add = groups,
        by_vars = exprs(!!!by_vars, !!!order_vars),
        check_type = "none"
      )
  }

  # join the input dataset with itself such that to each observation of the
  # input dataset all following observations are joined
  data_add_to_join <- select(
    data_add,
    !!!by_vars,
    !!!replace_values_by_names(extract_vars(order)),
    !!!replace_values_by_names(join_vars),
    !!tmp_obs_nr_var
  )

  if (get_admiral_option("save_memory")) {
    # split input dataset into smaller pieces and process them separately
    # This reduces the memory consumption.
    if (is.null(by_vars_left)) {
      # create batches of about 1MB input data
      obs_per_batch <- floor(1000000 / as.numeric(object.size(data) / nrow(data)))
      tmp_batch_nr <- get_new_tmp_var(data, prefix = "tmp_batch_nr")
      data_list <- data %>%
        mutate(!!tmp_batch_nr := ceiling(row_number() / obs_per_batch)) %>%
        group_by(!!tmp_batch_nr) %>%
        group_split(.keep = FALSE)
      data_add_list <- list(data_add_to_join)
    } else {
      data_nest <- nest(data, data = everything(), .by = vars2chr(unname(by_vars_left)))
      data_add_nest <- nest(data_add, data_add = everything(), .by = vars2chr(unname(by_vars_left)))
      data_all_nest <- inner_join(data_nest, data_add_nest, by = vars2chr(by_vars_left))
      data_list <- data_all_nest$data
      data_add_list <- data_all_nest$data_add
    }

    joined_data <- map2(
      data_list,
      data_add_list,
      function(x, y) {
        get_joined_sub_data(
          x,
          y,
          by_vars = by_vars_left,
          tmp_obs_nr_var = tmp_obs_nr_var,
          tmp_obs_nr_left = tmp_obs_nr_left,
          join_type = join_type,
          first_cond_upper = first_cond_upper,
          first_cond_lower = first_cond_lower,
          filter_join = filter_join
        )
      }
    )
  } else {
    joined_data <- get_joined_sub_data(
      data,
      dataset_add = data_add_to_join,
      by_vars = by_vars_left,
      tmp_obs_nr_var = tmp_obs_nr_var,
      tmp_obs_nr_left = tmp_obs_nr_left,
      join_type = join_type,
      first_cond_upper = first_cond_upper,
      first_cond_lower = first_cond_lower,
      filter_join = filter_join
    )
  }

  bind_rows(joined_data) %>%
    remove_tmp_vars() %>%
    select(-!!tmp_obs_nr_var_join)
}

#' Join Data for "joined" functions
#'
#' The helper function joins the data for the "joined" functions. All `.join`
#' variables are included in the output dataset. It is called by
#' `get_joined_data()` to process each by group separately. This reduces the
#' memory consumption.
#'
#' @param tmp_obs_nr_left Temporary observation number for `dataset`
#'
#' The specified variable has to be in the input dataset (`dataset`) and has to
#' be a unique key.
#'
#' @inheritParams get_joined_data
#'
#' @details
#'
#' 1. The input dataset (`dataset`) and the additional dataset (`dataset_add`)
#' are left joined by the grouping variables (`by_vars`). If no grouping
#' variables are specified, a full join is performed.
#'
#' 1. The joined dataset is restricted as specified by arguments `join_type`,
#' `first_cond_upper`, and `first_cond_lower`. See argument descriptions for
#' details.
#'
#' 1. The joined dataset is restricted by the `filter_join` condition.
#'
#' @keywords internal
get_joined_sub_data <- function(dataset,
                                dataset_add,
                                by_vars,
                                tmp_obs_nr_var,
                                tmp_obs_nr_left,
                                join_type,
                                first_cond_upper,
                                first_cond_lower,
                                filter_join) {
  data_joined <-
    left_join(
      dataset,
      dataset_add,
      by = vars2chr(by_vars),
      suffix = c("", ".join"),
      relationship = "many-to-many"
    )

  if (join_type != "all") {
    operator <- c(before = "<", after = ">")

    data_joined <- filter(
      data_joined,
      !!parse_expr(paste0(
        as_name(tmp_obs_nr_var), ".join",
        operator[join_type],
        as_name(tmp_obs_nr_var)
      ))
    )
  }

  if (!is.null(first_cond_upper)) {
    # select all observations up to the first confirmation observation
    data_joined <- filter_relative(
      data_joined,
      by_vars = expr_c(by_vars, tmp_obs_nr_left),
      condition = !!first_cond_upper,
      order = exprs(!!parse_expr(paste0(as_name(tmp_obs_nr_var), ".join"))),
      mode = "first",
      selection = "before",
      inclusive = TRUE,
      keep_no_ref_groups = FALSE
    )
  }

  if (!is.null(first_cond_lower)) {
    # select all observations up to the first confirmation observation
    data_joined <- filter_relative(
      data_joined,
      by_vars = expr_c(by_vars, tmp_obs_nr_left),
      condition = !!first_cond_lower,
      order = exprs(!!parse_expr(paste0("desc(", as_name(tmp_obs_nr_var), ".join)"))),
      mode = "first",
      selection = "before",
      inclusive = TRUE,
      keep_no_ref_groups = FALSE
    )
  }
  # apply confirmation condition, which may include summary functions
  data_joined %>%
    group_by(!!!by_vars, !!tmp_obs_nr_left) %>%
    filter_if(filter_join) %>%
    ungroup()
}
