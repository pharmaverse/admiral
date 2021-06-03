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
derive_disposition_reason <- function(dataset,
                                      dataset_ds,
                                      #new_var_reason,
                                      #status_var_reason,
                                      new_var,
                                      add_detail_var = NULL,
                                      #status_var_reason_spe = NULL,
                                      format_new_vars = format_reason_default,
                                      filter_ds) {
  # Checks
  assert_that(is.data.frame(dataset_ds))
  #print("1")
  #print(new_var)
  #print("1.1")
  #print(new_var$new_var)
  #print("2")
  #print(add_detail_var)
  #print("2.1")
  #print(add_detail_var$new_var)
  #print(add_detail_var$source_var)


  warn_if_vars_exist(dataset, as_name(new_var$new_var))

  if (!is.null(add_detail_var)) {
    warn_if_vars_exist(dataset, as_name(add_detail_var$new_var))
    source_var <- c(as_name(new_var$source_var),as_name(add_detail_var$source_var))

  }
  else {
    source_var <- as_name(new_var$source_var)
  }
  print("sopurce")
  print(source_var)
  assert_has_variables(dataset_ds, source_var)

  # Process the disposition data
  ds_subset <- dataset_ds %>%
    filter(!!enquo(filter_ds)) %>%
    select(STUDYID, USUBJID, source_var)
  print(head(ds_subset))

  # Expect 1 record per subject in the subsetted DS - issue a warning otherwise
  has_unique_records(
    dataset = ds_subset,
    by_vars = "USUBJID",
    message_type = "warning",
    message = "The filter used for DS results in several records per patient - please check"
  )
  # Add the status variable, derive the reason...
  dataset <- dataset %>%
    left_join(ds_subset, by = c("STUDYID", "USUBJID")) %>%
    mutate(!!new_var$new_var := format_new_vars(!!new_var$source_var))
  print(head(dataset))


  # ...and add the details (sepcify) if required
  print("check here")
  print(is.null(add_detail_var))
  if (!is.null(add_detail_var)) {
    dataset<-  dataset %>%
        mutate(!!add_detail_var$new_var := format_new_vars(!!new_var$source_var, !!add_detail_var$source_var))
  }

  print(head(dataset))
  dataset %>%
    select(-source_var)
}
#' Default format for the disposition reason
#'
#' Define a function to map the disposition reason
#'
#' @param x the disposition variable used for the mapping of the main  reason (e.g. `DSDECOD`).
#' @param y the disposition variable used for the mapping of the details  reason (e.g. `DSTERM`).
format_reason_default <- function(x, y = NULL) {
  if (is.null(y)) {
    case_when(
      x != "COMPLETED" & !is.na(x) ~ x,
      TRUE ~ NA_character_
    )
  }
  else {
    case_when(
      x != "COMPLETED" & !is.na(x) ~ y,
      TRUE ~ NA_character_
    )
  }
}

reason_source <- function(source_var, new_var) {
  out <- list(
    source_var=enquo(source_var),
    new_var = enquo(new_var)
  )
  class(out) <- c("reason_source", "list")
  validate_reason_source(out)
}

validate_reason_source <- function(x) {
  #print("heya")
  #print(class(x))
  #print(class(reason_source))
  assert_that(inherits(x, "reason_source"))
  values <- unclass(x)
  assert_that(is_expr(values$source_var))
  #print(values$source_var)
  #print(is_expr(values$source_var))
  assert_that(is_expr(values$new_var))
  x
}





