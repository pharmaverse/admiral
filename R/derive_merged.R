#' Merge Variables
#'
#' @export
#'
#' @examples
#' library(admiral.test)
#' data("vs")
#' data("dm")
#'
#' derive_vars_merged(vs,
#'                    dataset_add = select(dm, -DOMAIN),
#'                    by_vars = vars(STUDYID, USUBJID))
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

#' Merge a (imputed) Date Variable
#'
#' @export
#'
#' @examples
#' library(admiral.test)
#' data("dm")
#' data("ex")
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
#' derive_vars_merged_dt(
#'   select(dm, STUDYID, USUBJID),
#'   dataset_add = ex,
#'   by_vars = vars(STUDYID, USUBJID),
#'   new_vars_prefix = "TRTE",
#'   dtc = EXENDTC,
#'   date_imputation = "last",
#'   order = vars(TRTEDT),
#'   mode = "last"
#' )
derive_vars_merged_dt <- function(dataset,
                                  dataset_add,
                                  by_vars,
                                  order,
                                  new_vars_prefix,
                                  filter_add,
                                  mode,
                                  dtc,
                                  date_imputation = NULL,
                                  flag_imputation = TRUE,
                                  min_dates = NULL,
                                  max_dates = NULL) {
  dtc <- assert_symbol(enquo(dtc))
  assert_data_frame(dataset_add, required_vars = quo_c(by_vars, dtc))

  old_vars <- names(dataset_add)
  add_data <- derive_vars_dt(
    dataset_add,
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

#' Merge a (imputed) Datetime Variable
#'
#' @export
#'
#' @examples
#' library(admiral.test)
#' data("dm")
#' data("ex")
#' derive_vars_merged_dtm(
#'   select(dm, STUDYID, USUBJID),
#'   dataset_add = ex,
#'   by_vars = vars(STUDYID, USUBJID),
#'   new_vars_prefix = "TRTS",
#'   date = EXSTDTC,
#'   date_imputation = "first",
#'   time_imputation = "first",
#'   order = vars(TRTSDTM),
#'   mode = "first"
#' )
#'
#' derive_vars_merged_dtm(
#'   select(dm, STUDYID, USUBJID),
#'   dataset_add = ex,
#'   by_vars = vars(STUDYID, USUBJID),
#'   new_vars_prefix = "TRTE",
#'   date = EXENDTC,
#'   date_imputation = "last",
#'   time_imputation = "last",
#'   order = vars(TRTEDTM),
#'   mode = "last"
#' )

derive_vars_merged_dtm <- function(dataset,
                                  dataset_add,
                                  by_vars,
                                  order,
                                  new_vars_prefix,
                                  filter_add,
                                  mode,
                                  dtc,
                                  date_imputation = NULL,
                                  time_imputation = "00:00:00",
                                  flag_imputation = "auto",
                                  min_dates = NULL,
                                  max_dates = NULL) {
  dtc <- assert_symbol(enquo(dtc))
  assert_data_frame(dataset_add, required_vars = quo_c(by_vars, dtc))

  old_vars <- names(dataset_add)
  add_data <- derive_vars_dtm(
    dataset_add,
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

#' Merge a Categorized Variable
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

  add_data <- mutate(dataset_add,
                     !!new_var := cat_fun(!!source_var))
  derive_vars_merged(dataset,
                     dataset_add = add_data,
                     by_vars = by_vars,
                     order = order,
                     new_vars = vars(!!new_var),
                     filter_add = !!filter_add,
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
