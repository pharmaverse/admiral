#' Derive a date representing a disposition status.
#'
#' Derive dates from the the relevant reecords in the disposition domain.
#'
#' @param dataset Input dataset.
#'
#' @param dataset_ds Datasets containing the disposition information
#' (usually: ds).
#'
#' The variable specified in the status_var parameter must be in dataset_ds.
#'
#'
#' @param new_var Name of the disposition date variable.
#'
#' a variable name is expected.
#'
#' @param status_var The variable used to derive the disposition status.
#'
#'   A character vector is expected:
#'   If status_var equals 'COMPLETED' then the disposition status is set to COMPLETED,
#'   If status_var is not equal to 'COMPLETED' then the disposition status is set to DISCONTINUED,
#'   Otherwise if there is no record available for status_var in dataset_ds, set to ONGOING.
#'
#' @author Samia Kabi
#'
#' @export
#'
#' @examples
#' derive_disposition_eoxxstt(
#'   dataset = dm,
#'   dataset_ds = ds,
#'   new_var = EOSSTT,
#'   status_var = DSDECOD,
#'   filter = expr(DSCAT == "DISPOSITION EVENT" & DSSCAT == "STUDY COMPLETION/EARLY DISCONTINUATION")
#' )
derive_disposition_eoxxstt <- function(dataset,
                                       dataset_ds,
                                       new_var,
                                       status_var,
                                       filter = NULL) {
  # check if dataset_ds exists
  assert_dataset_exist(deparse(substitute(dataset_ds)))
  # Check status_var is present in dataset_ds
  assert_has_variables(dataset_ds, deparse(substitute(status_var)))

  # Warn if the variable to derive already exists in the input dataset
  warn_if_vars_exist(dataset, deparse(substitute(new_var)))

  # if DS needs to be filtered, filter
  if (!is.null(filter)) {
    ds_subset <- dataset_ds %>%
      filter(!!!filter)
  }
  else {
    ds_subset <- dataset_ds
  }

  # only 1 record per subject is expected - issue a warning otherwise
  assert_has_unique_records(
    dataset = ds_subset,
    by_vars = "USUBJID",
    message_type = "warning",
    message = "the filter used for DS results in several records per patient - please check"
  )

  # set status_var in status (resolves in mutate)
  ds_subset <- ds_subset %>%
    mutate(status___ = !!enquo(status_var)) %>%
    select(STUDYID, USUBJID, status___)

  # Add the status var and Derive the new dispo status in the input dataset
  dataset <- dataset %>%
    left_join(ds_subset, by = c("STUDYID", "USUBJID")) %>%
    mutate(
      !!enquo(new_var) := case_when(
        status___ == "COMPLETED" ~ "COMPLETED",
        status___ != "COMPLETED" ~ "DISCONTINUED",
        TRUE ~ "ONGOING"
      )
    ) %>%
    select(-ends_with("___"))
}
