#' Map 1:1 EX variables
#'
#' Maps the 1:1 variables from SDTM.EX to ADEX
#'
#' @param dataset Input dataset
#'
#' @param source_vars
#'   A list of the variables from SDTM.EX to be mapped is expected.
#'
#' @param by_vars
#'
#'   A list of variables by which observations are unique in SDTM.EX is expected.
#'
#'   Default:
#'
#' @param new_var
#'
#' @details The age is derived as the integer part of the duration from start to
#'   end date in the specified unit.
#'
#' @author Teckla G Akinyi
#'
#' @return The input dataset with ``AAGE`` and ``AAGEU`` added
#'
#' @export
#'
#'
#' @examples
#'

map_ex <- function(dataset,
                   source_var,
                   by_vars,
                   new_var){
  # Checks
  warn_if_vars_exist(dataset, deparse(substitute(new_var)))
  assert_that(is.data.frame(dataset),
              is.character(by_vars)
  )
  assert_has_variables(dataset, by_vars)
  assert_has_variables(dataset, source_var)

  # Expect 1 record per subject per treatment per time point
  #- issue a warning otherwise
  signal_duplicate_records(
    dataset,
    by_vars = vars(STUDYID, USUBJID, EXTRT,EXSTDTC),
    msg = "The dataset SDTM.EX results in multiple records
           per patient per treatment per timepoint",
    cnd_type="error"
  )

  #MAPPINGS
  #start times
  convert_dtc_to_dt(EXSTDTC) #out is ASTDT
  convert_dtc_to_dtm(EXSTDTC) #out is ASTDTM
  #end times
  convert_dtc_to_dt(EXENDTC) #out is AENDT
  convert_dtc_to_dtm(EXENDTC) #out is AENDTM

  #assign mappings
  dataset <- mutate(dataset,
                    ASTDY = EXSTDY,
                    AENDY = EXENDY,
                    AVAL = EXDOSE,
                    ATPTREF = EXTPTREF,
                    ATPTN = EXTPTNUM,
                    ATPT = EXTPT,
                    AVISIT = VISIT,
                    AVISITN = VISITNUM)
}
