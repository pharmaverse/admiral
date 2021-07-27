#' Derive the variable ATIREL
#'
#' Derives the variable ATIREL to CONCOMITANT, PRIOR, PRIOR_CONCOMITANT or NULL based on the relationship
#' of cm Analysis start/end date/times to treatment start date/time
#'
#' @param dataset Input dataset
#'
#'   The columns specified by the by_vars are expected
#'
#' @param dataset Input dataset2
#'   ADSL with treatment start date is expected
#'
#' @param by_var The start date
#'
#'   STUDYID and USUBJID to merge adsl to adcm are expected.
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
#' adsl <-data.frame(STUDYID,USUBJID,TRTSDTM)
#' adcm <-data.frame(STUDYID,USUBJID,ASTDTM,AENDTM,ASTTMF)
#'
#' derive_var_atirel(dataset=adcm,
#'                   dataset_adsl=adsl,
#'                   by_vars = vars(STUDYID,USUBJID)
#'                   new_var = ATIREL)
#'
#'

#function to derive ATIREL
derive_var_atirel <- function(dataset,
                              dataset_adsl,
                              by_vars,
                              new_var){
  # #checks
  assert_data_frame(dataset,vars(STUDYID, USUBJID, ASTDTM,AENDTM,ASTTMF))
  assert_data_frame(dataset_adsl, vars(STUDYID,USUBJID,TRTSDTM))
  by_vars <- assert_vars(by_vars)
  new_var <- assert_symbol(enquo(new_var))
  warn_if_vars_exist(dataset, quo_text(new_var))
  #source_vars <- assert_symbol(enquo(source_vars))

  #merge to adsl and logic to create ATIREL
  dataset_adsl <- dataset_adsl %>% select(STUDYID,USUBJID,TRTSDTM)
  dataset <-  dataset %>% dplyr::left_join(dataset_adsl, by=vars2chr(by_vars))%>%
    mutate(!!new_var :=ifelse(is.na(TRTSDTM),NA_character_,
                          ifelse( ASTDTM >= TRTSDTM,"CONCOMITANT",
                                  ifelse(!is.na(AENDTM) & AENDTM<TRTSDTM ,"PRIOR",
                                         ifelse(date(ASTDTM) == date(TRTSDTM) & toupper(ASTTMF) %in% c("H","M"), "CONCOMITANT",
                                                "PRIOR_CONCOMITANT"))))
          )
}











