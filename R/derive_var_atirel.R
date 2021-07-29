#' Derive the variable ATIREL
#'
#' Derives the variable ATIREL to CONCOMITANT, PRIOR, PRIOR_CONCOMITANT or NULL based on the relationship
#' of cm Analysis start/end date/times to treatment start date/time
#'
#' @param dataset Input dataset
#'
#'   The columns specified by the source_vars are expected
#'
#' @param source_vars Variables used to assign ATIREL
#'
#'   The variables TRTSDTM,ASTDTM,AENDTM,ASTTMF are expected
#'
#' @param new_var Name of variable to create
#'
#' @details ATIREL is set to:
#'            - null, if Datetime of First Exposure to Treatment is missing,
#'            - "CONCOMITANT", if the Analysis Start Date/Time is greater than or equal to Datetime of
#'               First Exposure to Treatment,
#'            - "PRIOR", if the Analysis End Date/Time is not missing and less than
#'               the Datetime of First Exposure to Treatment,
#'            - "CONCOMITANT" if the date part of Analysis Start Date/Time is equal to the date part of Datetime of
#'               First Exposure to Treatment and the Analysis Start Time Imputation Flag is 'H' or 'M'.
#'            - otherwise it is set to "PRIOR_CONCOMITANT"
#'
#' @author Teckla Akinyi
#'
#' @return The input data set with the ATIREL variable is returned
#'
#' @keywords ADaM Relationship Var ATIREL
#'
#' @export
#'
#' @seealso
#'
#' @examples
#' library(lubridate)
#' STUDYID <- "000-100"
#' USUBJID <- c("8712","8713","8714","8715","8716")
#' TRTSDTM  <- as_datetime(c("2012-02-25 23:00:00",
#'                           "",
#'                           "2017-02-25 23:00:00",
#'                           "2017-02-25 16:00:00",
#'                           "2017-02-25 16:00:00"))
#'
#' ASTDTM <- as_datetime(c("2012-02-28 19:00:00",
#'                        "2012-02-28 19:00:00",
#'                        "2013-02-25 19:00:00",
#'                        "2017-02-25 14:00:00",
#'                        "2017-02-25 14:00:00"))
#'
#' AENDTM  <- as_datetime(c("2012-02-25 23:00:00",
#'                          "",
#'                          "2014-02-25 19:00:00",
#'                          "2017-03-25 23:00:00",
#'                          "2017-04-29 14:00:00"))
#' ASTTMF <-  c("","","","M","")
#'
#' adcm <-data.frame(STUDYID,USUBJID,TRTSDTM,ASTDTM,AENDTM,ASTTMF)
#'
#' derive_var_atirel(dataset = adcm,
#'                   source_vars = vars(STUDYID,USUBJID,TRTSDTM,ASTDTM,AENDTM,ASTTMF),
#'                   new_var = ATIREL)
#'
#'

#function to derive ATIREL
derive_var_atirel <- function(dataset,
                              source_vars,
                              new_var){
  #checks
  dataout <- deparse(substitute(dataset))
  assert_data_frame(dataset, required_vars = source_vars)
  new_var <- assert_symbol(enquo(new_var))
  warn_if_vars_exist(dataset, quo_text(new_var))

  #merge to adsl and logic to create ATIREL
  dataset <-  dataset %>% mutate(!!new_var := case_when(is.na(TRTSDTM) ~ NA_character_,
                                                        ASTDTM >= TRTSDTM ~ "CONCOMITANT",
                                                        !is.na(AENDTM) & AENDTM<TRTSDTM  ~ "PRIOR",
                                                        date(ASTDTM) == date(TRTSDTM) & toupper(ASTTMF) %in% c("H","M") ~ "CONCOMITANT",
                                                        TRUE ~ "PRIOR_CONCOMITANT")
  )
  assign(x=dataout, value=dataset, envir=.GlobalEnv)
}
