#' Derive a disposition reason at a specific timepoint
#'
#' Derive a disposition reason from the the relevant records in the disposition domain.
#'
#' @param dataset Input dataset.
#'
#' @param dataset_ds Dataset containing the disposition information (e.g.: `ds`).
#'
#' The variable specified in the `status_var` parameter must be in `dataset_ds`.
#'
#' @param new_var Name of the disposition reason variable.
#'
#' a variable name is expected (e.g. `DCSREAS`).
#'
#' @param reason_var The variable used to derive the disposition reason
#'
#'   A character vector is expected (e.g. `DSDECOD`).
#'
#' @param format_new_vars The format used to derive the reason(s)
#'
#' Default: format_reason_default defined as:
#' format_reason_default<-function(x, y=NULL){
#' case_when (
#'   x == "COMPLETED" ~ x,
#'   x %!in% c("COMPLETED") & ! is.na(x)~"DISCONTINUED",
#'   TRUE ~ "ONGOING"
#' )
#' }
#' where x is the status_var.
#'
#' @param filter_ds Filter condition for the disposition data.
#'
#'   Filter used to select the relevant disposition data.
#'
#'   Permitted Values: logical expression.
#'
#' @return the input dataset with the disposition reason (`new_var`) added.
#'
#' @keywords adsl
#'
#' @author Samia Kabi
#'
#' @export
#'
#' @examples
#' library(dplyr)
#' data("dm")
#' data("ds")
derive_disposition_vars <- function(dataset,
                                    dataset_ds,
                                    filter_ds,
                                    disposition_vars
                                    #disposition_date,
                                    #disposition_status,
                                    #disposition_reason,
                                    #disposition_reason_spe
) {
  print("is list")
  print(is.list(disposition_vars))
  print("is list 1st part")
  print(is.list(disposition_vars[[1]]))

  if (!is.list(disposition_vars[[1]])) {
    list<-list(disposition_vars)
  }
  else {
    list<-disposition_vars
  }

  print("check here")
  print(list[[1]][[1]])
  print("enfd of check")
  # Checks
  assert_that(is.data.frame(dataset_ds))

  # Subset Disposition data as needed
  ds_subset <- dataset_ds %>%
    filter(!!enquo(filter_ds))

  # Expect 1 record per subject - issue a warning otherwise
  has_unique_records(
    dataset = ds_subset,
    by_vars = "USUBJID",
    message_type = "warning",
    message = "the filter used for DS results in several records per patient - please check"
  )

  # process each source
  # add_data <- vector("list", length(disposition_vars))

  for (ii in seq_along(list)) {
    print(ii)

    print(list[ii])
    print("class source var")
    print(class(list[[ii]]$source_var))
    print("source var")
    print(list[[ii]]$source_var)
    print(quo_get_expr(list[[ii]]$source_var))
    print(class(list[[ii]])[1])
    # print(as_name(disposition_vars[[ii]]$source_var))
    warn_if_vars_exist(dataset, as_name(list[[ii]]$new_var))
    assert_has_variables(dataset_ds, as_name(list[[ii]]$source_var))

    if (class(list[[ii]])[1] == "disposition_date_source") {
      prefix <- sub("\\DT.*", "", as_name(list[[ii]]$new_var))
      newvar <- paste0(prefix, "DT")
      ds_tmp <- ds_subset %>%
        mutate(datedtc___ = !!list[[ii]]$source_var) %>%
        derive_vars_dt(
          new_vars_prefix = prefix,
          dtc = datedtc___,
          date_imputation = list[[ii]]$date_imputation,
          flag_imputation = FALSE
        ) %>%
        select(STUDYID, USUBJID, !!list[[ii]]$new_var := !!sym(newvar))

      print(head(ds_tmp))

      # add the new dispo date to the input dataset
      dataset <- dataset %>%
        left_join(ds_tmp, by = c("STUDYID", "USUBJID"))
    }
    else {
      print("CLASSSSSS")
      print(class(list[[ii]]))
      if (class(list[[ii]])[1] == "disposition_reason_spe_source") {
        assert_has_variables(dataset_ds, as_name(list[[ii]]$source_var_main))
        ds_tmp <- ds_subset %>%
          select(STUDYID, USUBJID, !!list[[ii]]$source_var, !!list[[ii]]$source_var_main)
        print("i should see ")
        print(head(ds_tmp))

        check<-dataset %>%
          left_join(ds_tmp, by = c("STUDYID", "USUBJID"))
        print(head(check))
        print(list[[ii]]$new_var)
        print(list[[ii]]$format)

        dataset <- dataset %>%
          left_join(ds_tmp, by = c("STUDYID", "USUBJID")) %>%
          mutate(!!list[[ii]]$new_var := !!list[[ii]]$format) %>%
          select(-!!list[[ii]]$source_var, -!!list[[ii]]$source_var_main)
        print(head(dataset))
      }
      else {
        print(list[[ii]]$format)
        ds_tmp <- ds_subset %>%
          select(STUDYID, USUBJID, !!list[[ii]]$source_var)
        print(head(ds_tmp))

        dataset <- dataset %>%
          left_join(ds_tmp, by = c("STUDYID", "USUBJID")) %>%
          mutate(!!list[[ii]]$new_var := !!list[[ii]]$format) %>%
          select(-!!list[[ii]]$source_var)
        print(head(dataset))
      }
    }
  }
  dataset
}

disposition_date_source <- function(source_var=DSSTDTC, new_var, date_imputation = NULL) {
  out <- list(
    source_var = enquo(source_var),
    new_var = enquo(new_var),
    date_imputation = date_imputation
  )
  class(out) <- c("disposition_date_source", "list")
  validate_disposition_date_source(out)
}

validate_disposition_date_source <- function(x) {
  assert_that(inherits(x, "disposition_date_source"))
  values <- unclass(x)
  assert_that(is_expr(values$source_var))
  assert_that(is_expr(values$new_var))
  if (!is.null(values$date_imputation)) {
    assert_that(is_valid_date_entry(values$date_imputation))
  }
  x
}

disposition_source <- function(source_var=DSDECOD, new_var, format) {
  out <- list(
    source_var = enquo(source_var),
    new_var = enquo(new_var),
    format = enquo(format)
  )
  class(out) <- c("disposition_source", "list")
  validate_disposition_source(out)
}

validate_disposition_source <- function(x) {
  assert_that(inherits(x, "disposition_source"))
  values <- unclass(x)
  # assert_that(is_expr(values$source_var))
  assert_that(is_expr(values$new_var))
  assert_that(is_expr(values$format))
  x
}


disposition_reason_spe_source <- function(source_var=DSTERM, source_var_main=DSDECOD, new_var, format) {
  out <- list(
    source_var = enquo(source_var),
    source_var_main = enquo(source_var_main),
    new_var = enquo(new_var),
    format = enquo(format)
  )
  class(out) <- c("disposition_reason_spe_source", "list")
  validate_disposition_reason_spe_source(out)
}

validate_disposition_reason_spe_source <- function(x) {
  assert_that(inherits(x, "disposition_reason_spe_source"))
  values <- unclass(x)
  assert_that(is_expr(values$source_var))
  assert_that(is_expr(values$source_var_main))
  assert_that(is_expr(values$new_var))
  assert_that(is_expr(values$format))
  x
}



