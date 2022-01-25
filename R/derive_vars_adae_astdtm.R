#' @title Derive ASTDTM and Flags for ADAE
#' 
#' @description
#' This function derives ADAE Start Date/Time Variable (ADAE.ASTDTM) and imputation flags 
#' (ADAE.ASTDTF and ADAE.ASTTMF). Admiral support funciton derive_vars_dtm is
#' used for initial date/time imputation. \cr
#' The derivation is based on the default specifications for ADAE.ASTDTM from the GDSR
#' (Version 2021-Oct.v1: Shenzou12). \cr
#' The following variables have to be available in the input data: 
#' AESTDTC, TRTSDTM and AENDTM
#'
#' @param dataset Input AE dataset with derived TRTSDTM and AENDTM and source AESTDTC
#'
#' @param drop_empty_flag Flag to indicate whether the date/time imputation flag(s) 
#' should be dropped in case all values are missing. Possible Options are T(RUE) and 
#' F(ALSE). The default is F
#' 
#' @author Sandra Ott
#'
#' @details
#' The presence of 'ASTDTM', 'ASTDTF' and 'ASTTMF' is checked and a warning
#' is issued when any of them is found. The existing variables will be overwritten.
#' 
#' @return
#' The input AE dataset with the datetime 'ASTDTM' and the date/time imputation 
#' flags 'ASTDTF' and 'ASTTMF' added.
#'
#' @keywords adam timing derivation adae
#' 
#' @examples
#' \dontrun{
#' adae <- tibble::tribble(
#'   ~AESTDTC, ~TRTSDTM,  ~AENDTM,
#'   "2020-07-02T22:10:09", ymd_hms("2020-12-06T12:12:12"), ymd_hms("2020-07-08T12:12:12"),
#'   "2020-07-02T22:10", ymd_hms("2020-12-06T12:12:12"), ymd_hms("2020-07-08T12:12:12"),
#'   "2020-12", ymd_hms("2020-12-06T12:12:12"), ymd_hms("2020-12-06T12:12:12"),
#'   "2020", ymd_hms("2020-12-06T12:12:12"), ymd_hms("2020-12-06T12:12:12")
#' )
#' 
#' # Run function with default parameter setting
#' derive_vars_adae_astdtm(adae) 
#' 
#' # Run function and drop empty imputation flags
#' derive_vars_adae_astdtm(adae, drop_empty_flag = T) 
#' }
#' 
#' @export
#' 
derive_vars_adae_astdtm <- 
  function(dataset, drop_empty_flag = F) {
    
  assert_data_frame(dataset, required_vars = vars(AESTDTC, TRTSDTM, AENDTM))
  assert_logical_scalar(drop_empty_flag)
    
  # Warn if ASTDTM, ASTDTF or ASTTMF variables already exist
  check <- warn_if_vars_exist(dataset, c("ASTDTM", "ASTDTF", "ASTTMF"))
  if (!is.null(check)) {
    warning("Existing variables will be overwritten")
    # Drop variables as they will be re-derived in the next step
    dataset <- select(dataset, -any_of(c("ASTDTM", "ASTDTF", "ASTTMF")))
  }
  

  # Derive ASTDTM with standard imputation rules
  # This will also derive the flags
  dataset01 <- dataset %>% derive_vars_dtm(
    dtc = AESTDTC,
    new_vars_prefix = "AST",
    date_imputation = "first",
    time_imputation = "last",
    flag_imputation = "both"
  ) 
  
  # update date/time variables class to POSIXct/ POSIXt
  class(dataset01$ASTDTM) <- c("POSIXct", "POSIXt")
  class(dataset01$TRTSDTM) <- c("POSIXct", "POSIXt")
  class(dataset01$AENDTM) <- c("POSIXct", "POSIXt")
  
  # If missing date part or missing year of Start Date/Time of Adverse Event 
  # [ADAE.AESTDTC], set to earliest of non-missing Datetime of First Exposure to 
  # Treatment [ADSL.TRTSDTM] or the (imputed) Analysis End Date/Time [ADAE.AENDTM]. 
  # Set to null, if Datetime of First Exposure to Treatment [ADSL.TRTSDTM] and 
  # Analysis End Date/Time [ADAE.AENDTM] are is missing.
  dataset02 <- 
    dataset01 %>% 
    mutate(ASTDTM = if_else(AESTDTC == "", pmin(TRTSDTM, AENDTM , na.rm = T), ASTDTM),
           ASTDTF = if_else(AESTDTC == "" & !is.na(ASTDTM), "Y", ASTDTF),
           ASTTMF = if_else(AESTDTC == "" & !is.na(ASTDTM), "H", ASTTMF))
  
  # For missing month update ASTDTM to Datetime of First Exposure to Treatment 
  # [ADSL.TRTSDTM] (with time replaced as 23:59:59) if Datetime of First Exposure to 
  # Treatment [ADSL.TRTSDTM] is not missing and has matching year (and month) part to 
  # Start Date/Time of Adverse Event [ADAE.AESTDTC].
  # Else set to Analysis End Date/Time [ADAE.AENDTM] if Datetime 
  # of First Exposure to Treatment [ADSL.TRTSDTM] is missing or has a year part 
  # greater than the year part of Start Date/Time of Adverse Event [ADAE.AESTDTC],
  # and Analysis End Date/Time [ADAE.AENDTM] has matching year part to Start
  # Date/Time of Adverse Event [ADAE.AESTDTC].
  # If start date imputation takes date past Analysis End Date/Time [ADAE.AENDTM], 
  # then set to Analysis End Date/Time [ADAE.AENDTM].
  
  # Update Time Part of TRTSDTM
  dataset02$TRTSDTM_ORIG <- dataset02$TRTSDTM
  hour(dataset02$TRTSDTM) <-  23
  minute(dataset02$TRTSDTM) <-  59
  second(dataset02$TRTSDTM) <-  59
  
  dataset03 <- 
    dataset02 %>% 
    mutate(
      ASTDTM = case_when(
        (nchar(AESTDTC) == 4 | word(AESTDTC, 2, sep = "-") == "")  & 
          as.numeric(word(AESTDTC, 1, sep = "-")) == year(TRTSDTM) & 
          TRTSDTM <= AENDTM ~ TRTSDTM, 
        (nchar(AESTDTC) == 4 | word(AESTDTC, 2, sep = "-") == "")  &  
          ( as.numeric(word(AESTDTC, 1, sep = "-")) != year(TRTSDTM) | is.na(TRTSDTM) |
            as.numeric(word(AESTDTC, 1, sep = "-")) == year(TRTSDTM)  & TRTSDTM > AENDTM  ) &
          as.numeric(word(AESTDTC, 1, sep = "-")) == year(AENDTM) ~ AENDTM,
       TRUE  ~ ASTDTM
      ))
  
  # For missing day update ASTDTM to the Datetime of First Exposure to Treatment
  # [ADSL.TRTSDTM] (with time replaced as 23:59:59) if Datetime of First 
  # Exposure to Treatment [ADSL.TRTSDTM] is in same month and year as start 
  # date of adverse event [ADAE.AESTDTC] and Analysis End Date/Time [ADAE.AENDTM] 
  # is not before Datetime of First Exposure to Treatment [ADSL.TRTSDTM]. 
  
  dataset04 <- 
    dataset03 %>% 
    mutate(
      ASTDTM = case_when(
        nchar(word(AESTDTC, 1, sep = "T")) == 7 & 
          as.numeric(word(AESTDTC, 1, sep = "-")) == year(TRTSDTM) & 
          as.numeric(word(AESTDTC, 2, sep = "-")) == month(TRTSDTM) &
          AENDTM >= TRTSDTM ~ TRTSDTM,
        TRUE  ~ ASTDTM
      ))

  # If after imputation Analysis Start Date/Time [ADAE.ASTDTM] is on same day as 
  # the Analysis End Date/Time [ADAE.AENDTM] but the start time goes beyond the 
  # end time, then set Analysis Start Date/Time [ADAE.ASTDTM] to Analysis 
  # End Date/Time [ADAE.AENDTM].
  dataset <- 
    dataset04 %>% 
    mutate(
      ASTDTM = case_when(
        date(AENDTM) == date(ASTDTM) &
        AENDTM < ASTDTM ~ AENDTM,
        TRUE  ~ ASTDTM
      ))
  
  # Put TRTSDTM back to original
  dataset$TRTSDTM <-  dataset$TRTSDTM_ORIG
  dataset <- select(dataset, -TRTSDTM_ORIG)
  
  # Drop empty flags if requested 
  if (drop_empty_flag == T) {
    if(suppressWarnings(max(nchar(dataset$ASTDTF), na.rm = T)) == -Inf) {
      dataset <- select(dataset, -ASTDTF)
    }
    if(suppressWarnings(max(nchar(dataset$ASTTMF), na.rm = T)) == -Inf) {
      dataset <- select(dataset, -ASTTMF)
    }
  }
  
  # Return updated dataset
  dataset
  
}
