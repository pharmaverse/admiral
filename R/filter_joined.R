#' Filter Observations Taking Other Observations into Account
#'
#' @description
#'
#' The function filters observation using a condition taking other observations
#' into account. For example, it could select all observations with `AVALC ==
#' "Y"` and `AVALC == "Y"` for at least one subsequent observation. The input
#' dataset is joined with itself to enable conditions taking variables from both
#' the current observation and the other observations into account. The suffix
#' ".join" is added to the variables from the subsequent observations.
#'
#' An example usage might be checking if a patient received two required
#' medications within a certain timeframe of each other.
#'
#' In the oncology setting, for example, we use such processing to check if a
#' response value can be confirmed by a subsequent assessment. This is commonly
#' used in endpoints such as best overall response.
#'
#' @param dataset Input dataset
#'
#'   The variables specified for `by_vars`, `join_vars`, and `order` are
#'   expected.
#'
#' @param by_vars By variables
#'
#'   The specified variables are used as by variables for joining the input
#'   dataset with itself.
#'
#' @param join_vars Variables to keep from joined dataset
#'
#'   The variables needed from the other observations should be specified for
#'   this parameter. The specified variables are added to the joined dataset
#'   with suffix ".join". For example to select all observations with `AVALC ==
#'   "Y"` and `AVALC == "Y"` for at least one subsequent visit `join_vars =
#'   exprs(AVALC, AVISITN)` and `filter = AVALC == "Y" & AVALC.join == "Y" &
#'   AVISITN < AVISITN.join` could be specified.
#'
#'   The `*.join` variables are not included in the output dataset.
#'
#' @param join_type Observations to keep after joining
#'
#'   The argument determines which of the joined observations are kept with
#'   respect to the original observation. For example, if `join_type =
#'   "after"` is specified all observations after the original observations are
#'   kept.
#'
#'   *Permitted Values:* `"before"`, `"after"`, `"all"`
#'
#' @param first_cond Condition for selecting range of data
#'
#'   If this argument is specified, the other observations are restricted up to
#'   the first observation where the specified condition is fulfilled. If the
#'   condition is not fulfilled for any of the subsequent observations, all
#'   observations are removed.
#'
#' @param order Order
#'
#'   The observations are ordered by the specified order.
#'
#'   *Permitted Values:* list of expressions created by `exprs()`, e.g.,
#'   `exprs(ADT, desc(AVAL))`
#'
#' @param tmp_obs_nr_var Temporary observation number
#'
#'   The specified variable is added to the input dataset and set to the
#'   observation number with respect to `order`. For each by group (`by_vars`)
#'   the observation number starts with `1`. The variable can be used in the
#'   conditions (`filter`, `first_cond`). It is not included in the output
#'   dataset. It can be used to select consecutive observations or the last
#'   observation (see last example below).
#'
#' @param filter Condition for selecting observations
#'
#'   The filter is applied to the joined dataset for selecting the confirmed
#'   observations. The condition can include summary functions. The joined
#'   dataset is grouped by the original observations. I.e., the summary function
#'   are applied to all observations up to the confirmation observation. For
#'   example in the oncology setting when using this function for confirmed best
#'   overall response,  `filter = AVALC == "CR" & all(AVALC.join %in% c("CR",
#'   "NE")) & count_vals(var = AVALC.join, val = "NE") <= 1` selects
#'   observations with response "CR" and for all observations up to the
#'   confirmation observation the response is "CR" or "NE" and there is at most
#'   one "NE".
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
#' @details
#'
#'   The following steps are performed to produce the output dataset.
#'
#'   ## Step 1
#'
#'   The input dataset is joined with itself by the variables specified for
#'   `by_vars`. From the right hand side of the join only the variables
#'   specified for `join_vars` are kept. The suffix ".join" is added to these
#'   variables.
#'
#'   For example, for `by_vars = USUBJID`, `join_vars = exprs(AVISITN, AVALC)` and input dataset
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
#'   The first observation of each group is selected and the `*.join` variables
#'   are dropped.
#'
#' @returns A subset of the observations of the input dataset. All variables of
#'   the input dataset are included in the output dataset.
#'
#'
#' @keywords utils_fil
#' @family utils_fil
#'
#' @seealso [count_vals()], [min_cond()], [max_cond()]
#'
#' @export
#'
#' @examples
#'
#' library(tibble)
#' library(admiral)
#'
#' # filter observations with a duration longer than 30 and
#' # on or after 7 days before a COVID AE (ACOVFL == "Y")
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
#' filter_joined(
#'   adae,
#'   by_vars = exprs(USUBJID),
#'   join_vars = exprs(ACOVFL, ADY),
#'   join_type = "all",
#'   order = exprs(ADY),
#'   filter = ADURN > 30 & ACOVFL.join == "Y" & ADY >= ADY.join - 7
#' )
#'
#' # filter observations with AVALC == "Y" and AVALC == "Y" at a subsequent visit
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
#' filter_joined(
#'   data,
#'   by_vars = exprs(USUBJID),
#'   join_vars = exprs(AVALC, AVISITN),
#'   join_type = "after",
#'   order = exprs(AVISITN),
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
#' filter_joined(
#'   data,
#'   by_vars = exprs(USUBJID),
#'   join_vars = exprs(AVALC),
#'   join_type = "after",
#'   order = exprs(AVISITN),
#'   first_cond = AVALC.join == "CR",
#'   filter = AVALC == "CR" & all(AVALC.join %in% c("CR", "NE")) &
#'     count_vals(var = AVALC.join, val = "NE") <= 1
#' )
#'
#' # select observations with AVALC == "PR", AVALC == "CR" or AVALC == "PR"
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
#' filter_joined(
#'   data,
#'   by_vars = exprs(USUBJID),
#'   join_vars = exprs(AVALC, ADY),
#'   join_type = "after",
#'   order = exprs(ADY),
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
#'
#' # select observations with CRIT1FL == "Y" at two consecutive visits or at the last visit
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
#' filter_joined(
#'   data,
#'   by_vars = exprs(USUBJID),
#'   tmp_obs_nr_var = tmp_obs_nr,
#'   join_vars = exprs(CRIT1FL),
#'   join_type = "all",
#'   order = exprs(AVISITN),
#'   filter = CRIT1FL == "Y" & CRIT1FL.join == "Y" &
#'     (tmp_obs_nr + 1 == tmp_obs_nr.join | tmp_obs_nr == max(tmp_obs_nr.join))
#' )
#'
filter_joined <- function(dataset,
                          by_vars,
                          join_vars,
                          join_type,
                          first_cond = NULL,
                          order,
                          tmp_obs_nr_var = NULL,
                          filter,
                          check_type = "warning") {
  # Check input parameters
  assert_vars(by_vars)
  assert_vars(join_vars)
  join_type <-
    assert_character_scalar(
      join_type,
      values = c("before", "after", "all"),
      case_sensitive = FALSE
    )
  first_cond <- assert_filter_cond(enexpr(first_cond), optional = TRUE)
  assert_expr_list(order)
  tmp_obs_nr_var <- assert_symbol(enexpr(tmp_obs_nr_var), optional = TRUE)
  filter <- assert_filter_cond(enexpr(filter))
  check_type <-
    assert_character_scalar(
      check_type,
      values = c("none", "warning", "error"),
      case_sensitive = FALSE
    )
  assert_data_frame(
    dataset,
    required_vars = expr_c(by_vars, join_vars, extract_vars(order))
  )

  # number observations of the input dataset to get a unique key
  # (by_vars and tmp_obs_nr_var)
  if (is.null(tmp_obs_nr_var)) {
    tmp_obs_nr_var <- get_new_tmp_var(dataset, prefix = "tmp_obs_nr_")
  }
  data <- dataset %>%
    derive_var_obs_number(
      new_var = !!tmp_obs_nr_var,
      by_vars = by_vars,
      order = order,
      check_type = check_type
    )
  # join the input dataset with itself such that to each observation of the
  # input dataset all following observations are joined
  data_joined <-
    left_join(
      data,
      select(data, !!!by_vars, !!!join_vars, !!tmp_obs_nr_var),
      by = vars2chr(by_vars),
      suffix = c("", ".join")
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

  if (!is.null(first_cond)) {
    # select all observations up to the first confirmation observation
    data_joined <- filter_relative(
      data_joined,
      by_vars = expr_c(by_vars, tmp_obs_nr_var),
      condition = !!first_cond,
      order = exprs(!!parse_expr(paste0(as_name(tmp_obs_nr_var), ".join"))),
      mode = "first",
      selection = "before",
      inclusive = TRUE,
      keep_no_ref_groups = FALSE
    )
  }

  # apply confirmation condition, which may include summary functions
  data_joined %>%
    group_by(!!!by_vars, !!tmp_obs_nr_var) %>%
    filter(!!filter) %>%
    # select one observation of each group, as the joined variables are removed
    # it doesn't matter which one, so we take just the first one
    slice(1L) %>%
    ungroup() %>%
    select(colnames(dataset))
}

#' Count Number of Observations Where a Variable Equals a Value
#'
#' Count number of observations where a variable equals a value.
#'
#' @param var A vector
#'
#' @param val A value
#'
#'
#' @keywords utils_fil
#' @family utils_fil
#'
#' @export
#'
#' @examples
#'
#' library(tibble)
#' library(dplyr, warn.conflicts = FALSE)
#' library(admiral)
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
#' # add variable providing the number of NEs for each subject
#' group_by(data, USUBJID) %>%
#'   mutate(nr_nes = count_vals(var = AVALC, val = "NE"))
count_vals <- function(var, val) {
  length(var[var == val])
}

#' Minimum Value on a Subset
#'
#' The function derives the minimum value of a vector/column on a subset of
#' entries/observations.
#'
#' @param var A vector
#'
#' @param cond A condition
#'
#'
#' @keywords utils_fil
#' @family utils_fil
#'
#' @export
#'
#' @examples
#'
#' library(tibble)
#' library(dplyr, warn.conflicts = FALSE)
#' library(admiral)
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
#' )
#'
#' # In oncology setting, when needing to check the first time a patient had
#' # a Complete Response (CR) to compare to see if any Partial Response (PR)
#' # occurred after this add variable indicating if PR occurred after CR
#' group_by(data, USUBJID) %>% mutate(
#'   first_cr_vis = min_cond(var = AVISITN, cond = AVALC == "CR"),
#'   last_pr_vis = max_cond(var = AVISITN, cond = AVALC == "PR"),
#'   pr_after_cr = last_pr_vis > first_cr_vis
#' )
min_cond <- function(var, cond) {
  assert_filter_cond(enexpr(cond))
  if (length(var[cond]) == 0) {
    NA
  } else {
    min(var[cond])
  }
}

#' Maximum Value on a Subset
#'
#' The function derives the maximum value of a vector/column on a subset of
#' entries/observations.
#'
#' @param var A vector
#'
#' @param cond A condition
#'
#'
#' @keywords utils_fil
#' @family utils_fil
#'
#' @export
#'
#' @examples
#'
#' library(tibble)
#' library(dplyr, warn.conflicts = FALSE)
#' library(admiral)
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
#' )
#'
#' # In oncology setting, when needing to check the first time a patient had
#' # a Complete Response (CR) to compare to see if any Partial Response (PR)
#' # occurred after this add variable indicating if PR occurred after CR
#' group_by(data, USUBJID) %>% mutate(
#'   first_cr_vis = min_cond(var = AVISITN, cond = AVALC == "CR"),
#'   last_pr_vis = max_cond(var = AVISITN, cond = AVALC == "PR"),
#'   pr_after_cr = last_pr_vis > first_cr_vis
#' )
max_cond <- function(var, cond) {
  assert_filter_cond(enexpr(cond))
  if (length(var[cond]) == 0) {
    NA
  } else {
    max(var[cond])
  }
}
