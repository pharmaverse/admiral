#' Returns records that fit into existing by groups in a filtered source dataset
#'
#' Returns all records in the input dataset that belong to by groups that are present
#' in a source dataset, after the source dataset is optionally filtered. For example,
#' this could be used to return ADSL records for subjects that experienced a certain
#' adverse event during the course of the study (as per records in ADAE).
#'
#' @param dataset `r roxygen_param_dataset(expected_vars = c("by_vars"))`
#'
#' @param dataset_add Source dataset
#'
#'   The source dataset, which determines the by groups returned in the input dataset,
#'   based on the groups that exist in this dataset after being subset by `filter_add`.
#'
#'   The variables specified in the `by_vars` and `filter_add` parameters are expected
#'   in this dataset.
#'
#' @param by_vars Grouping variables
#'
#' `r roxygen_param_by_vars()`
#'
#' @param filter_add Filter for the source dataset
#'
#'   The filter condition which will be used to subset the source dataset.
#'   Alternatively, if no filter condition is supplied, no subsetting of the source
#'   dataset will be performed.
#'
#'
#' @details Returns the records in `dataset` which match an existing by group in `dataset_add`,
#'   after being filtered according to `filter_add`. If there are no by groups that exist
#'   in both datasets, an empty dataset will be returned.
#'
#' @return The records in the input dataset which are contained within an existing by group in
#'   the filtered source dataset.
#'
#' @keywords utils_fil
#'
#' @family utils_fil
#'
#' @export
#'
#' @examples
#' # Get demographic information about subjects who have suffered from moderate or
#' # severe fatigue
#'
#' library(tibble)
#'
#' adsl <- tribble(
#'   ~USUBJID,      ~AGE, ~SEX,
#'   "01-701-1015", 63,   "F",
#'   "01-701-1034", 77,   "F",
#'   "01-701-1115", 84,   "M",
#'   "01-701-1146", 75,   "F",
#'   "01-701-1444", 63,   "M"
#' )
#'
#' adae <- tribble(
#'   ~USUBJID,      ~AEDECOD,                    ~AESEV,     ~AESTDTC,
#'   "01-701-1015", "DIARRHOEA",                 "MODERATE", "2014-01-09",
#'   "01-701-1034", "FATIGUE",                   "SEVERE",   "2014-11-02",
#'   "01-701-1034", "APPLICATION SITE PRURITUS", "MODERATE", "2014-08-27",
#'   "01-701-1115", "FATIGUE",                   "MILD",     "2013-01-14",
#'   "01-701-1146", "FATIGUE",                   "MODERATE", "2013-06-03"
#' )
#'
#' filter_exist(
#'   dataset = adsl,
#'   dataset_add = adae,
#'   by_vars = exprs(USUBJID),
#'   filter_add = AEDECOD == "FATIGUE" & AESEV %in% c("MODERATE", "SEVERE")
#' )
#'
filter_exist <- function(dataset,
                         dataset_add,
                         by_vars,
                         filter_add = NULL) {
  assert_vars(by_vars)
  assert_data_frame(
    dataset,
    required_vars = by_vars
  )
  filter_add <- assert_filter_cond(enexpr(filter_add), optional = TRUE)
  assert_data_frame(
    dataset_add,
    required_vars = by_vars
  )

  dataset %>%
    semi_join(
      dataset_add %>%
        filter_if(filter_add),
      by = vars2chr(by_vars)
    )
}

#' Returns records that don't fit into existing by groups in a filtered source dataset
#'
#' Returns all records in the input dataset that belong to by groups that are not
#' present in a source dataset, after the source dataset is optionally filtered. For
#' example, this could be used to return ADSL records for subjects that didn't take certain
#' concomitant medications during the course of the study (as per records in ADCM).
#'
#' @inheritParams filter_exist
#'
#' @param dataset_add Source dataset
#'
#'   The source dataset, which determines the by groups returned in the input dataset,
#'   based on the groups that don't exist in this dataset after being subset by `filter_add`.
#'
#'   The variables specified in the `by_vars` and `filter_add` parameters are expected
#'   in this dataset.
#'
#' @details Returns the records in `dataset` which don't match any existing by groups in
#'   `dataset_add`, after being filtered according to `filter_add`. If all by
#'    groups that exist in `dataset` don't exist in `dataset_add`, an empty dataset will
#'    be returned.
#'
#' @return The records in the input dataset which are not contained within any existing by
#'    group in the filtered source dataset.
#'
#' @keywords utils_fil
#'
#' @family utils_fil
#'
#' @export
#'
#' @examples
#' # Get demographic information about subjects who didn't take vitamin supplements
#' # during the study
#'
#' library(tibble)
#'
#' adsl <- tribble(
#'   ~USUBJID,      ~AGE, ~SEX,
#'   "01-701-1015", 63,   "F",
#'   "01-701-1023", 64,   "M",
#'   "01-701-1034", 77,   "F",
#'   "01-701-1118", 52,   "M"
#' )
#'
#' adcm <- tribble(
#'   ~USUBJID,      ~CMTRT,         ~CMSTDTC,
#'   "01-701-1015", "ASPIRIN",      "2013-05-14",
#'   "01-701-1023", "MYLANTA",      "2014-01-04",
#'   "01-701-1023", "CALCIUM",      "2014-02-25",
#'   "01-701-1034", "VITAMIN C",    "2013-12-12",
#'   "01-701-1034", "CALCIUM",      "2013-03-27",
#'   "01-701-1118", "MULTIVITAMIN", "2013-02-21"
#' )
#'
#' filter_not_exist(
#'   dataset = adsl,
#'   dataset_add = adcm,
#'   by_vars = exprs(USUBJID),
#'   filter_add = str_detect(CMTRT, "VITAMIN")
#' )
#'
filter_not_exist <- function(dataset,
                             dataset_add,
                             by_vars,
                             filter_add = NULL) {
  assert_vars(by_vars)
  assert_data_frame(
    dataset,
    required_vars = by_vars
  )
  filter_add <- assert_filter_cond(enexpr(filter_add), optional = TRUE)
  assert_data_frame(
    dataset_add,
    required_vars = by_vars
  )

  dataset %>%
    anti_join(
      dataset_add %>%
        filter_if(filter_add),
      by = vars2chr(by_vars)
    )
}
