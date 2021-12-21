#' Merge Variables from a Dataset to the Input Dataset
#'
#' Merge variables from a dataset to the input dataset. The observations to
#' merge can be selected by a condition and/or selecting the first or last
#' observation for each by group.
#'
#' @param dataset Input dataset
#'
#'   The variables specified by the `by_vars` parameter are expected.
#'
#' @param dataset_add Additional dataset
#'
#'   The variables specified by the `by_vars`, the `new_vars`, and the `order`
#'   parameter are expected.
#'
#' @param order Sort order
#'
#'   The first or last observation from the addtional dataset is selected with
#'   respect to the specified order.
#'
#'   *Permitted Values*: list of variables or `desc(<variable>)` function calls
#'
#' @param new_vars Variables to add
#'
#'   The specified variables are added to the output dataset. Variables can be
#'   renamed by naming the element, i.e., `new_vars = vars(<new name> = <old
#'   name>)`.
#'
#'   If the parameter is not specified or set to `NULL`, all variables from the
#'   additional dataset are added.
#'
#'   *Permitted Values*: list of variables created by `vars()`
#'
#' @param mode Selection mode
#'
#'   Determines if the first or last observation is selected.
#'
#'   *Permitted Values*: `"first"`, `"last"`
#'
#' @param by_vars Grouping variables
#'
#'   The input dataset and the selected observations from the additional dataset
#'   are merged by the specified by variables. The by variables must be a unique
#'   key of the selected observations.
#'
#'   *Permitted Values*: list of variables
#'
#' @param filter_add Filter for additional data
#'
#'   Only observations fulfilling the specified condition are taken into account
#'   for merging. If the parameter is not specified, all observations are
#'   considered.
#'
#'   *Permitted Values*: a condition
#'
#' @param check_type Check uniqueness?
#'
#'   If `"warning"` or `"error"` is specified, the specified message is issued
#'   if the observations of the (restricted) additional dataset are not unique
#'   with respect to the by variables and the order.
#'
#'   *Permitted Values*: `"none"`, `"warning"`, `"error"`
#'
#' @details
#'
#'   1. The additional dataset is restricted to the observations matching the
#'   `filter_add` condition.
#'
#'   1. If `order` is specified, for each by group the first or last observation
#'   (depending on `mode`) is selected.
#'
#'   1. The variables specified for `new_vars` are merged to the input dataset.
#'
#' @author Stefan Bundfuss
#'
#' @keywords derivation adam
#'
#' @export
#'
#' @examples
#' library(admiral.test)
#' library(dplyr)
#' data("vs")
#' data("dm")
#'
#' # merging all dm variables to vs
#' derive_vars_merged(
#'   vs,
#'   dataset_add = select(dm, -DOMAIN),
#'   by_vars = vars(STUDYID, USUBJID)) %>%
#' select(STUDYID, USUBJID, VSTESTCD, VISIT, VSTPT, VSSTRESN, AGE, AGEU)
#'
#' # merge last weight to adsl
#' data("adsl")
#' derive_vars_merged(
#'   adsl,
#'   dataset_add = vs,
#'   by_vars = vars(STUDYID, USUBJID),
#'   order = vars(VSDTC),
#'   mode = "last",
#'   new_vars = vars(LASTWGT = VSSTRESN, LASTWGTU = VSSTRESU),
#'   filter_add = VSTESTCD == "WEIGHT"
#' ) %>%
#' select(STUDYID, USUBJID, AGE, AGEU, LASTWGT, LASTWGTU)
derive_vars_merged <- function(dataset,
                               dataset_add,
                               by_vars,
                               order = NULL,
                               new_vars = NULL,
                               mode = NULL,
                               filter_add = NULL,
                               check_type = "warning") {
  filter_add <- assert_filter_cond(enquo(filter_add), optional = TRUE)
  assert_vars(by_vars)
  assert_order_vars(order, optional = TRUE)
  assert_vars(new_vars, optional = TRUE)
  assert_data_frame(dataset, required_vars = by_vars)
  assert_data_frame(dataset_add, required_vars = quo_c(by_vars, extract_vars(order)))

  add_data <- filter_if(dataset_add, filter_add)
  if (!is.null(order)) {
    add_data <- filter_extreme(add_data,
                               by_vars = by_vars,
                               order = order,
                               mode = mode,
                               check_type = check_type)
  } else {
    signal_duplicate_records(
      add_data,
      by_vars = by_vars,
      msg = paste(
        "Dataset `dataset_add` contains duplicate records with respect to",
        enumerate(vars2chr(by_vars))
      )
    )
  }
  if (!is.null(new_vars)) {
    add_data <- select(add_data, !!!by_vars, !!!new_vars)
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
        "The variables",
        enumerate(common_vars),
        "are contained in both datasets.\n",
        "Please add them to `by_vars` or remove or rename them in one of the datasets."
      )
    ))
  }
  left_join(dataset, add_data, by = vars2chr(by_vars))
}

#' Merge a (Imputed) Date Variable
#'
#' Merge a imputed date variable and date imputation  flag from a dataset to the
#' input dataset. The observations to merge can be selected by a condition
#' and/or selecting the first or last observation for each by group.
#'
#' @param dataset_add Additional dataset
#'
#'   The variables specified by the `by_vars`, the `dtc`, and the `order`
#'   parameter are expected.
#'
#' @inheritParams derive_vars_merged
#' @inheritParams derive_vars_dt
#'
#' @details
#'
#'   1. The additional dataset is restricted to the observations matching the
#'   `filter_add` condition.
#'
#'   1. The date variable and if requested, the date imputation flag is added to
#'   the additional dataset.
#'
#'   1. If `order` is specified, for each by group the first or last observation
#'   (depending on `mode`) is selected.
#'
#'   1. The date and flag variables are merged to the input dataset.
#'
#' @author Stefan Bundfuss
#'
#' @keywords derivation adam timing
#'
#' @export
#'
#' @examples
#' library(admiral.test)
#' data("dm")
#' data("ex")
#'
#' # derive treatment start date (TRTSDT)
#' derive_vars_merged_dt(
#'   select(dm, STUDYID, USUBJID),
#'   dataset_add = ex,
#'   by_vars = vars(STUDYID, USUBJID),
#'   new_vars_prefix = "TRTS",
#'   dtc = EXSTDTC,
#'   date_imputation = "first",
#'   order = vars(TRTSDT),
#'   mode = "first"
#' )
#'
#' # derive treatment end date (TRTEDT) (without imputation)
#' derive_vars_merged_dt(
#'   select(dm, STUDYID, USUBJID),
#'   dataset_add = ex,
#'   by_vars = vars(STUDYID, USUBJID),
#'   new_vars_prefix = "TRTE",
#'   dtc = EXENDTC,
#'   order = vars(TRTEDT),
#'   mode = "last"
#' )
derive_vars_merged_dt <- function(dataset,
                                  dataset_add,
                                  by_vars,
                                  order = NULL,
                                  new_vars_prefix,
                                  filter_add = NULL,
                                  mode = NULL,
                                  dtc,
                                  date_imputation = NULL,
                                  flag_imputation = TRUE,
                                  min_dates = NULL,
                                  max_dates = NULL) {
  dtc <- assert_symbol(enquo(dtc))
  filter_add <- assert_filter_cond(enquo(filter_add), optional = TRUE)
  assert_data_frame(dataset_add, required_vars = quo_c(by_vars, dtc))

  old_vars <- names(dataset_add)
  add_data <- filter_if(dataset_add, filter_add) %>%
    derive_vars_dt(
      new_vars_prefix = new_vars_prefix,
      dtc = !!dtc,
      date_imputation = date_imputation,
      flag_imputation = flag_imputation,
      min_dates = min_dates,
      max_dates = max_dates
    )
  new_vars <- quos(!!!syms(setdiff(names(add_data), old_vars)))
  derive_vars_merged(dataset,
                     dataset_add = add_data,
                     by_vars = by_vars,
                     order = order,
                     new_vars = new_vars,
                     mode = mode)
}

#' Merge a (Imputed) Datetime Variable
#'
#' Merge a imputed datetime variable, date imputation  flag, and time imputation
#' flag from a dataset to the input dataset. The observations to merge can be
#' selected by a condition and/or selecting the first or last observation for
#' each by group.
#'
#' @param dataset_add Additional dataset
#'
#'   The variables specified by the `by_vars`, the `dtc`, and the `order`
#'   parameter are expected.
#'
#' @inheritParams derive_vars_merged
#' @inheritParams derive_vars_dtm
#'
#' @details
#'
#'   1. The additional dataset is restricted to the observations matching the
#'   `filter_add` condition.
#'
#'   1. The datetime variable and if requested, the date imputation flag and
#'   time imputation flag is added to the additional dataset.
#'
#'   1. If `order` is specified, for each by group the first or last observation
#'   (depending on `mode`) is selected.
#'
#'   1. The date and flag variables are merged to the input dataset.
#'
#' @author Stefan Bundfuss
#'
#' @keywords derivation adam timing
#'
#' @export
#'
#' @examples
#' library(admiral.test)
#' data("dm")
#' data("ex")
#'
#' # derive treatment start datetime (TRTSDTM)
#' derive_vars_merged_dtm(
#'   select(dm, STUDYID, USUBJID),
#'   dataset_add = ex,
#'   by_vars = vars(STUDYID, USUBJID),
#'   new_vars_prefix = "TRTS",
#'   dtc = EXSTDTC,
#'   date_imputation = "first",
#'   time_imputation = "first",
#'   order = vars(TRTSDTM),
#'   mode = "first"
#' )
#'
#' # derive treatment end datetime (TRTEDTM) (without date imputation)
#' derive_vars_merged_dtm(
#'   select(dm, STUDYID, USUBJID),
#'   dataset_add = ex,
#'   by_vars = vars(STUDYID, USUBJID),
#'   new_vars_prefix = "TRTE",
#'   dtc = EXENDTC,
#'   time_imputation = "last",
#'   order = vars(TRTEDTM),
#'   mode = "last"
#' )
derive_vars_merged_dtm <- function(dataset,
                                   dataset_add,
                                   by_vars,
                                   order = NULL,
                                   new_vars_prefix,
                                   filter_add = NULL,
                                   mode = NULL,
                                   dtc,
                                   date_imputation = NULL,
                                   time_imputation = "00:00:00",
                                   flag_imputation = "auto",
                                   min_dates = NULL,
                                   max_dates = NULL) {
  dtc <- assert_symbol(enquo(dtc))

  filter_add <- assert_filter_cond(enquo(filter_add), optional = TRUE)
  assert_data_frame(dataset_add, required_vars = quo_c(by_vars, dtc))

  old_vars <- names(dataset_add)
  add_data <- filter_if(dataset_add, filter = filter_add) %>%
    derive_vars_dtm(
      new_vars_prefix = new_vars_prefix,
      dtc = !!dtc,
      date_imputation = date_imputation,
      time_imputation = time_imputation,
      flag_imputation = flag_imputation,
      min_dates = min_dates,
      max_dates = max_dates
    )
  new_vars <- quos(!!!syms(setdiff(names(add_data), old_vars)))
  derive_vars_merged(dataset,
                     dataset_add = add_data,
                     by_vars = by_vars,
                     order = order,
                     new_vars = new_vars,
                     mode = mode)
}

#' Merge a Categorization Variable
#'
#' Merge a categorization variable from a dataset to the input dataset. The
#' observations to merge can be selected by a condition and/or selecting the
#' first or last observation for each by group.
#'
#' @param dataset_add Additional dataset
#'
#'   The variables specified by the `by_vars`, the `source_var`, and the `order`
#'   parameter are expected.
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
#'   A function must be specified for this parameter which expects the values of
#'   the source variable as input and returns the categorized values.
#'
#' @inheritParams derive_vars_merged
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
#' @author Stefan Bundfuss
#'
#' @keywords derivation adam
#'
#' @export
#'
#' @examples
#' library(admiral.test)
#' data("dm")
#' data("vs")
#' wgt_cat <- function(wgt) {
#'   case_when(wgt < 50 ~ "low",
#'             wgt > 90 ~ "high",
#'             TRUE ~ "normal")
#' }
#'
#' derive_var_merged_cat(
#'   dm,
#'   dataset_add = vs,
#'   by_vars = vars(STUDYID, USUBJID),
#'   order = vars(VSDTC, VSSEQ),
#'   filter_add = VSTESTCD == "WEIGHT" & substr(VISIT, 1, 9) == "SCREENING" ,
#'   new_var = WGTBLCAT,
#'   source_var = VSSTRESN,
#'   cat_fun = wgt_cat,
#'   mode = "last"
#' )
derive_var_merged_cat <- function(dataset,
                                  dataset_add,
                                  by_vars,
                                  order = NULL,
                                  new_var,
                                  source_var,
                                  cat_fun,
                                  filter_add = NULL,
                                  mode = NULL) {
  new_var <- assert_symbol(enquo(new_var))
  source_var <- assert_symbol(enquo(source_var))
  filter_add <- assert_filter_cond(enquo(filter_add), optional = TRUE)
  assert_data_frame(dataset_add, quo_c(by_vars, source_var))

  add_data <- filter_if(dataset_add, filter_add) %>%
    mutate(!!new_var := cat_fun(!!source_var))
  derive_vars_merged(dataset,
                     dataset_add = add_data,
                     by_vars = by_vars,
                     order = order,
                     new_vars = vars(!!new_var),
                     mode = mode)
}

#' Merge a Existence Flag
#'
#' @export
#'
#' @examples
#'
#' library(admiral.test)
#' data("dm")
#' data("ae")
#' derive_var_merged_exist_flag(
#'   dm,
#'   dataset_add = ae,
#'   by_vars = vars(STUDYID, USUBJID),
#'   new_var = AESERFL,
#'   condition = AESER == "Y"
#' )
derive_var_merged_exist_flag <- function(
  dataset,
  dataset_add,
  by_vars,
  new_var,
  condition,
  true_value = "Y",
  false_value = NA_character_,
  filter_add
) {
  new_var <- assert_symbol(enquo(new_var))
  condition <- assert_filter_cond(enquo(condition), optional = TRUE)

  add_data <- mutate(dataset_add, cond_flag := if_else(!!condition, 1, 0, 0))

  derive_vars_merged(
    dataset,
    dataset_add = add_data,
    by_vars = by_vars,
    new_vars = vars(cond_flag),
    order = vars(cond_flag),
    check_type = "none",
    mode = "last"
  ) %>%
    mutate(!!new_var := if_else(cond_flag == 1, true_value, false_value, false_value)) %>%
    select(-cond_flag)
}

#' Merge a Character Variable
#'
#' @export
#'
#' @examples
#' library(admiral.test)
#' data("dm")
#' data("ds")
#' derive_var_merged_character(
#'   dm,
#'   dataset_add = ds,
#'   by_vars = vars(STUDYID, USUBJID),
#'   new_var = DISPSTAT,
#'   filter_add = DSCAT == "DISPOSITION EVENT",
#'   source_var = DSDECOD,
#'   case = "title"
#' )
derive_var_merged_character <- function(dataset,
                                        dataset_add,
                                        by_vars,
                                        order = NULL,
                                        new_var,
                                        source_var,
                                        case = NULL,
                                        filter_add = NULL,
                                        mode = NULL) {

  new_var <- assert_symbol(enquo(new_var))
  source_var <- assert_symbol(enquo(source_var))
  case <-
    assert_character_scalar(
      case,
      values = c("lower", "upper", "title"),
      case_sensitive = FALSE,
      optional = TRUE
    )
  filter_add <- assert_filter_cond(enquo(filter_add), optional = TRUE)
  assert_data_frame(dataset_add, required_vars = quo_c(by_vars, source_var))

  if (is.null(case)) {
    trans <- expr(!!source_var)
  }
  else if (case == "lower") {
    trans <- expr(str_to_lower(!!source_var))
  }
  else if (case == "upper") {
    trans <- expr(str_to_upper(!!source_var))
  }
  else if (case == "title") {
    trans <- expr(str_to_title(!!source_var))
  }
  add_data <- mutate(dataset_add, !!new_var := !!trans)
  derive_vars_merged(dataset,
                     dataset_add = add_data,
                     by_vars = by_vars,
                     order = order,
                     new_vars = vars(!!new_var),
                     filter_add = !!filter_add,
                     mode = mode)
}
