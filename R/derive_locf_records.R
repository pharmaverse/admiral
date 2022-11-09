#' Derive LOCF(Last observation carried forward) records
#'
#' Adds LOCF records as new observations for each 'by group' when the data set
#' does not contain observations for missed visits/time points.
#'
#' @param dataset Input dataset
#'
#'   The columns specified by the `by_vars` and the `order`
#'   parameter are expected.
#'
#' @param dataset_expected_obs Expected observations dataset
#'
#'   Data frame with all the combinations of `PARAMCD`, `PARAM`, `AVISIT`,
#'   `AVISITN`, ... which are expected in the dataset is expected.
#'
#' @param by_vars Grouping variables
#'
#'   For each group defined by `by_vars` an observation is added to the output
#'   dataset. Only variables specified in `by_vars` will be populated
#'   in the newly created records.
#'
#' @param order List of variables for sorting a dataset
#'
#'   The dataset is sorted by `order` before carrying the last
#'   observation (eg. `AVAL`) within each `by_vars` forward.
#'
#' @author G Gayatri
#'
#' @details For each group (with respect to the variables specified for the
#'   `by_vars` parameter) an observation is added to the output dataset from
#'   `dataset_expected_obs`, if the input dataset does not contain observations
#'   for missed visits/time points or if AVAL is missing.
#'
#'   For the new observations, `AVAL` is set to the non-missing `AVAL` of the
#'   previous observation in the input dataset (when sorted by `order`) and
#'   `DTYPE` is set to `LOCF`.
#'
#' @return The input dataset with the new `LOCF` observations added for each
#' `by_vars`. Note, a variable will only be populated in the new parameter rows
#' if it is specified in `by_vars`.
#'
#' @keywords der_bds_findings
#' @family der_bds_findings
#'
#' @export
#'
#' @examples
#'
#' advs <- tibble::tribble(
#'   ~STUDYID,  ~USUBJID,      ~PARAMCD, ~PARAM,                          ~AVAL, ~AVISITN, ~AVISIT,
#'   "CDISC01", "01-701-1015", "DIABP", "Diastolic Blood Pressure (mmHg)", 51,    0,       "BASELINE",
#'   "CDISC01", "01-701-1015", "DIABP", "Diastolic Blood Pressure (mmHg)", 50,    2,       "WEEK 2",
#'   "CDISC01", "01-701-1015", "DIABP", "Diastolic Blood Pressure (mmHg)", 51,    4,       "WEEK 4",
#'   "CDISC01", "01-701-1015", "DIABP", "Diastolic Blood Pressure (mmHg)", 50,    6,       "WEEK 6",
#'   "CDISC01", "01-701-1015", "SYSBP", "Systolic Blood Pressure (mmHg)", 121,    0,       "BASELINE",
#'   "CDISC01", "01-701-1015", "SYSBP", "Systolic Blood Pressure (mmHg)", 121,    2,       "WEEK 2",
#'   "CDISC01", "01-701-1015", "SYSBP", "Systolic Blood Pressure (mmHg)", 121,    4,       "WEEK 4",
#'   "CDISC01", "01-701-1015", "SYSBP", "Systolic Blood Pressure (mmHg)", 121,    6,       "WEEK 6",
#'   "CDISC01", "01-701-1028", "DIABP", "Diastolic Blood Pressure (mmHg)", 79,    0,       "BASELINE",
#'   "CDISC01", "01-701-1028", "DIABP", "Diastolic Blood Pressure (mmHg)", 80,    2,       "WEEK 2",
#'   "CDISC01", "01-701-1028", "DIABP", "Diastolic Blood Pressure (mmHg)", NA,    4,       "WEEK 4",
#'   "CDISC01", "01-701-1028", "DIABP", "Diastolic Blood Pressure (mmHg)", NA,    6,       "WEEK 6",
#'   "CDISC01", "01-701-1028", "SYSBP", "Systolic Blood Pressure (mmHg)", 130,    0,       "BASELINE",
#'   "CDISC01", "01-701-1028", "SYSBP", "Systolic Blood Pressure (mmHg)", 132,    2,       "WEEK 2"
#' )
#'
#'
#' # A dataset with all the combinations of PARAMCD, PARAM, AVISIT, AVISITN, ... which are expected.
#' advs_expected_obsv <- tibble::tribble(
#'  ~PARAMCD, ~PARAM, ~AVISITN, ~AVISIT,
#'   "DIABP", "Diastolic Blood Pressure (mmHg)", 0, "BASELINE",
#'   "DIABP", "Diastolic Blood Pressure (mmHg)", 2, "WEEK 2",
#'   "DIABP", "Diastolic Blood Pressure (mmHg)", 4, "WEEK 4",
#'   "DIABP", "Diastolic Blood Pressure (mmHg)", 6, "WEEK 6",
#'   "SYSBP", "Systolic Blood Pressure (mmHg)",  0, "BASELINE",
#'   "SYSBP", "Systolic Blood Pressure (mmHg)",  2, "WEEK 2",
#'   "SYSBP", "Systolic Blood Pressure (mmHg)",  4, "WEEK 4",
#'   "SYSBP", "Systolic Blood Pressure (mmHg)",  6, "WEEK 6"
#' )
#'
#' # Add a new record for each 'by_vars' when the input dataset does not contain
#' # observations for missed visits/timepoints or when AVAL is missing.
#' # Set DTYPE = 'LOCF' for the newly added records.
#' # Carry forward the AVAL of the previous observation after sorting the dataset
#' # by 'order' parameter.
#' Final_advs <- derive_locf_records(data=advs,
#'                                   dataset_expected_obs=advs_expected_obsv,
#'                                   by_vars=vars(STUDYID, USUBJID, PARAM, PARAMCD),
#'                                   order=vars(STUDYID, USUBJID, PARAM, PARAMCD, AVISITN, AVISIT))



derive_locf_records <- function(dataset,
                                dataset_expected_obs,
                                by_vars,
                                order) {
  library(tidyr)

  #### Input Checking ####

  # Check by_vars and order variables in input dataset
  assert_data_frame(dataset, required_vars = by_vars)
  assert_data_frame(dataset, required_vars = order)

  # Check if input parameters is a valid list of variables
  assert_vars(by_vars, optional = TRUE)
  assert_order_vars(order)



  #### Prepping 'dataset_expected_obs' ####

  # Add temp_dtype='LOCF' for the 'expected observation dataset'
  dataset_expected_obs <- dataset_expected_obs %>%
    mutate(temp_dtype="LOCF")


  # Get the STUDYID and USUBJID from input dataset to set expected
  # observations for each subject
  subjids <- dataset %>%
    select(STUDYID, USUBJID) %>%
    distinct()

  exp_obsv <- subjids %>%
    crossing(dataset_expected_obs)



  ##### Add LOCF records ####

  # Split input dataset into the missing and non-missing AVAL records
  aval_missing <- dataset %>%
    filter(is.na(AVAL))

  aval_not_missing <- dataset %>%
    drop_na(AVAL)


  # Get the variable names from the expected observation dataset
  exp_obs_vars <- exp_obsv %>%
    select(-temp_dtype) %>%
    colnames()


  # Get unique combination of visits/timepoints per parameter per subject
  # from the original input dataset (with non-missing AVAL)
  advs_unique_original <- aval_not_missing %>%
    select(all_of(exp_obs_vars)) %>%
    drop_na() %>%
    distinct()


  # Get all the expected observations that are to be added to the input
  # dataset (with non-missing AVAL)
  advs_exp_obsv3 <- exp_obsv %>%
    anti_join(advs_unique_original, by=c(exp_obs_vars))


  # Merge the expected observations with the input dataset (with non-missing AVAL)
  # Arrange the dataset by 'order' and group it by 'by_vars'
  # Use fill() to fill the AVAL from the previous observation for the newly added records
  aval_not_missing_locf <- aval_not_missing %>%
    full_join(advs_exp_obsv3,  by=c(exp_obs_vars)) %>%
    mutate(DTYPE = case_when(temp_dtype=="LOCF" ~ "LOCF")) %>%
    select(-temp_dtype) %>%
    arrange(!!!order) %>%
    group_by(!!!by_vars) %>%
    fill("AVAL") %>%
    ungroup()


  # Output dataset - merge the AVAL missing with non-missing+newly added LOCF records
   output <- bind_rows(aval_not_missing_locf, aval_missing) %>%
    arrange(!!!order)



}






