#' Impute partial date/time portion of a --DTC variable
#'
#' Imputation partial date/time portion of a --DTC variable
#' based on user input
#'
#' @param dtc The --DTC date to impute
#'
#'   A character date is expected in a format like yyy-mm-yy
#'
#' @param impute Whether to impute partial date or no
#'
#'   A logical value is expected.
#'
#'   Default: TRUE
#'
#' @param day The value to impute the day when day is missing
#'
#'   Permitted Values: number from 1 to 31, text 'FIRST', 'MID', 'LAST'
#'   'FIRST' must be used to impute to the first day of the month
#'   'MID' must be used to impute to the 15th of the month
#'   'LAST' must be used to impute to the last day of the month
#'
#'   Default: 1
#'
#' @param month the value to impute the month when month is missing
#'
#'   Permitted Values: number from 1 to 12
#'
#'   Default: 1
#'
#' @param year the value to impute the year when year is missing
#'
#'   Permitted Values: any 4-digit number or 'NONE' if a missing date shoudl not be imputed
#'
#'   Default: 'NONE'
#'
#' @author Samia Kabi
#'
#' @return  a character
#'  when impute = TRUE: with a complete/time  : e.g. 2021-03-13T00:00:00
#'  when impute = FALSE: with a complete date/time if the date is complete,
#'  or a "" if the date is partial
#'
#'
#' @family {general functions}
#'
#' @export
#'
#' @examples
#'
#' date<-"2019-07-18T15:25"
#' compute_imputed_dtc(dtc="2019-07-18T15:25")#returns "2019-07-18T15:25:00"
#' compute_imputed_dtc(dtc="2019-07-18")#returns "2019-07-18T00:00:00"
#' compute_imputed_dtc(dtc="2019-07")#returns "2019-07-01T00:00:00"
#' compute_imputed_dtc(dtc="2019", day="LAST")#returns "2019-01-31T00:00:00"
#' compute_imputed_dtc(dtc="2019", day="LAST", month=12)#returns "2019-12-31T00:00:00"
#' compute_imputed_dtc(dtc="2019-07", day="LAST", hour=23, min=59, sec=59)#returns "2019-07-31T23:59:59"
#' ompute_imputed_dtc(dtc="2019-07", day="FIRST")#returns "2019-07-01T00:00:00"
#' compute_imputed_dtc(dtc="2019-07", day="MID")#returns "2019-07-15T00:00:00"
#' compute_imputed_dtc(dtc="2019")#returns "2019-01-01T00:00:00"
#' compute_imputed_dtc(dtc="2019---07")#returns "2019-01-01T00:00:00"
#'
compute_imputed_dtc<-function(dtc,
                              impute=TRUE,
                              day=01, month=01,year='NONE',
                              hour=0, min=0, sec=0)
  {

  #check inputs are valid
  assert_that(is_valid_sec_min(sec))
  assert_that(is_valid_sec_min(min))
  assert_that(is_valid_hour(hour))
  assert_that(is_valid_day(day))
  assert_that(is_valid_month(month))
  assert_that(is_valid_year(year))
  assert_that(is.logical(impute))

  if(!day %in% c("FIRST", "MID", "LAST")) {
    d<-sprintf("%02d", as.integer(day))}
  else {
    if (day %in% c("FIRST",  "LAST")){
      d<-"01"}
    else {
      if (day =="MID"){
        d<-"15"
      }
    }
  }

  mo<-sprintf("%02d", as.integer(month))
  if (year !='NONE'){y<-sprintf("%04d", as.integer(year))}

  h<-sprintf("%02d", as.integer(hour))
  mi<-sprintf("%02d", as.integer(min))
  s<-sprintf("%02d", as.integer(sec))

  if (impute==TRUE){
    #add missing text in date/time
    tmp__ <- case_when(
         nchar(dtc)==19 ~ dtc,
         nchar(dtc)==16 ~ paste0(dtc,":", s),
         nchar(dtc)==10 ~ paste0(dtc,"T",h,":", mi,":", s),
         nchar(dtc)==9 ~  paste0(substr(dtc,1,4), "-",mo,"-",d,"T",h,":", mi,":", s),#dates like 2021---14 - use only year part
         nchar(dtc)==7 ~  paste0(dtc,"-",d,"T",h,":", mi,":", s),
         nchar(dtc)==4  ~  paste0(dtc, "-",mo,"-",d,"T",h,":", mi,":", s),
         if (year !='NONE'){
           nchar(dtc)<4  ~  paste0(y, "-",mo,"-",d,"T",h,":", mi,":", s)
         }
         else {TRUE ~ ""}
       )
    if (day=="LAST"){
        tmp__<-if_else(nchar(tmp__)>0,
                    paste0(as.character(ceiling_date(as.Date(tmp__, format = "%Y-%m-%dT%H:%M:%S"), "month") - days(1)),"T",h,":", mi,":", s)
                    ,"")
    }
  }
  #No imputation: keep the datetime or teh date with a time set to 0
  else{
    tmp__ <- case_when(
      nchar(dtc)==19 ~ dtc,
      nchar(dtc)==16 ~ paste0(dtc,":", s),
      nchar(dtc)==10 ~ paste0(dtc,"T00:00:00"),
      TRUE ~ "" #as.character(dtc)
    )
  }
  return(tmp__)
}


dtc_dt<-function(dtc){
  case_when(
    nchar(dtc)>=10 ~ymd(substr(dtc,1,10)),
    TRUE ~ymd(NA)
  )

}

dtc_dtm<-function(dtc){

  #note T00:00:00 is not printed in dataframe
  case_when(
    nchar(dtc)==19 ~ymd_hms(dtc),
    TRUE ~ymd_hms(NA)
  )
}


compute_dtf <- function(dtc,dt){

  case_when(
    !is.na(dt) & nchar(dtc)< 4 ~"Y",
    !is.na(dt) & nchar(dtc)==4 ~"M",
    !is.na(dt) & nchar(dtc)==7 ~"D",
    !is.na(dt) & nchar(dtc)==9 ~"M",#dates like "2019---07"
    (!is.na(dt) & nchar(dtc)>=10) | is.na(dt) ~" "

  )
}

compute_tmf <- function(dtc,dt){

  case_when(
    (!is.na(dt) & nchar(dtc)==19)| is.na(dt) ~" ",
    !is.na(dt) & nchar(dtc)==16~"S",
    !is.na(dt) & nchar(dtc)==13~"M",
    (!is.na(dt) & nchar(dtc)==10) | (nchar(dtc)> 0 & nchar(dtc) <10) ~"H"
  )
}

#' derive_vars_dt
#'
#' Add --DT/--DTF variable(s) based on --DTC
#' --DT is imputed basd on user input
#'
#' @param dataset Input dataset
#'
#'   The --DTC variable must be present.
#'
#' @param new_vars_prefix
#'
#' Prefix used for the output variable: e.g new_vars_prefix="AST" will output ASTDT and ASTDTF
#'
#' @param dtc The --DTC date to impute
#'
#'   A character date is expected in a format like yyy-mm-yy
#'
#' @param impute Whether to impute partial date or no
#'
#'   A logical value is expected.
#'
#'   Default: TRUE
#'
#' @param day The value to impute the day when day is missing
#'
#'   Permitted Values: number from 1 to 31, text 'FIRST', 'MID', 'LAST'
#'   'FIRST' must be used to impute to the first day of the month
#'   'MID' must be used to impute to the 15th of the month
#'   'LAST' must be used to impute to the last day of the month
#'
#'   Default: 1
#'
#' @param month the value to impute the month when month is missing
#'
#'   Permitted Values: number from 1 to 12
#'
#'   Default: 1
#'
#' @param year the value to impute the year when year is missing
#'
#'   Permitted Values: any 4-digit number or 'NONE' if a missing date shoudl not be imputed
#'
#'   Default: 'NONE'
#'
#' @param DTF Whether to add the date iomputation flag
#'
#'   A logical value is expected.
#'
#'   Default: TRUE
#'
#' @author Samia Kabi
#'
#' @return The input dataset with --DT (and --DTF if DTF=TRUE)
#'
#'
#' @family {general functions}
#'
#' @export
#'
#' @examples
#' mhdt <- tibble::tribble(
#' ~MHSTDTC,
#' "2019-07-18T15:25:40",
#' "2019-07-18T15:25",
#' "2019-07-18",
#' "2019-02",
#' "2019",
#' "2019---07",
#' "")
##Derive ASTDT and ASDTF based on MHSTDTC - day/month are imputed to 01 by default, missing dates are not imputed
#mhdt1<-derive_vars_dt(mhdt,new_vars_prefix = "AST", dtc=MHSTDTC)
#
##Derive AENDT and AENTF based on MHSTDTC - day/month are imputed are imputed to the last, missing dates are not imputed
#mhdt2<-derive_vars_dt(mhdt,new_vars_prefix = "AST", dtc=MHSTDTC, day="LAST", month=12)
#


derive_vars_dt<- function(
  dataset,
  new_vars_prefix,
  dtc,
  impute=TRUE,
  day=01,#use "FIRST", "LAST", "MID"  to have the first/last/15th day of teh month
  month=01, year='NONE',#2020,
  DTF=TRUE
){
  #Check DTC is present in input dataset
  assert_has_variables(dataset, deparse(substitute(dtc)))

  tmp<-dataset
  warn_if_vars_exist(dataset, paste0(deparse(substitute(new_vars_prefix)),"DT") )
  tmp<-tmp %>%
    mutate(dt=dtc_dt(dtc=compute_imputed_dtc(dtc=!!enquo(dtc),
                                             impute=impute,
                                             year=year, month=month, day=day,
                                             hour=0, min=0, sec=0)
                     )
    )

  if(DTF==TRUE){
    warn_if_vars_exist(dataset, paste0(deparse(substitute(new_vars_prefix)),"DTF") )
    tmp<-tmp %>%
        mutate(dtf=compute_dtf(dtc=!!enquo(dtc), dt=dt))#use dt since as.name(paste0(new_vars_prefix,"DT")) does not resolve...
    #rename DTF with relevant prefix
    names(tmp)[names(tmp) == "dtf"] <- paste0(new_vars_prefix,"DTF")
  }
  #rename DT with relevant prefix
  names(tmp)[names(tmp) == "dt"] <- paste0(new_vars_prefix,"DT")

  return(tmp)

}


#' derive_vars_dtm
#'
#' Add --DTM/--TMF variable(s) based on --DTC
#' --DTM is imputed basd on user input
#'
#' @param dataset Input dataset
#'
#'   The --DTC variable must be present.
#'
#' @param new_vars_prefix
#'
#' Prefix used for the output variable: e.g new_vars_prefix="AST" will output ASTDT and ASTDTF
#'
#' @param dtc The --DTC date to impute
#'
#'   A character date is expected in a format like yyy-mm-yy
#'
#' @param impute Whether to impute partial date or no
#'
#'   A logical value is expected.
#'
#'   Default: TRUE
#'
#' @param day The value to impute the day when day is missing
#'
#'   Permitted Values: number from 1 to 31, text 'FIRST', 'MID', 'LAST'
#'   'FIRST' must be used to impute to the first day of the month
#'   'MID' must be used to impute to the 15th of the month
#'   'LAST' must be used to impute to the last day of the month
#'
#'   Default: 1
#'
#' @param month the value to impute the month when month is missing
#'
#'   Permitted Values: number from 1 to 12
#'
#'   Default: 1
#'
#' @param year the value to impute the year when year is missing
#'
#'   Permitted Values: any 4-digit number or 'NONE' if a missing date shoudl not be imputed
#'
#'   Default: 'NONE'
#'
#' @param hour the value to impute the hour when time is missing
#'
#'   Permitted Values: number from 0 to 23
#'
#'   Default: 0
#'
#' @param min the value to impute the minutes when time is missing
#'
#'   Permitted Values: number from 0 to 59
#'
#'   Default: 0
#'
#' @param sec the value to impute the seconds when time is missing
#'
#'   Permitted Values: number from 0 to 59
#'
#'   Default: 0
#'
#' @param TMF Whether to add the time imputation flag
#'
#'   A logical value is expected.
#'
#'   Default: TRUE
#'
#' @author Samia Kabi
#'
#' @return The input dataset with --DTM (and --TMF if DTF=TRUE)
#'
#'
#' @family {general functions}
#'
#' @export
#'
#' @examples
#' mhdt <- tibble::tribble(
#' ~MHSTDTC,
#' "2019-07-18T15:25:40",
#' "2019-07-18T15:25",
#' "2019-07-18",
#' "2019-02",
#' "2019",
#' "2019---07",
#' "")
#' #Derive ASTDTM and ASTTMF based on MHSTDTC - day/month are imputed to 01 by default, time is imputed to 00:00
#' by defaul, missing dates are not imputed
#' mhdt1<-derive_vars_dtm(mhdt,new_vars_prefix = "A", dtc=MHSTDTC)
#' #Derive AENDTM and AENTMF based on MHSTDTC - day/month are imputed are imputed to the last, missing dates are not imputed
#' mhdt2<-derive_vars_dtm(mhdt,new_vars_prefix = "AST", dtc=MHSTDTC, day="LAST", month=12, hour=23, min=59, sec=59)
#'
#'

derive_vars_dtm<- function(
  dataset,
  new_vars_prefix,
  dtc,
  impute=TRUE,
  day=01,#use "FIRST", "LAST", "MID"  to have the first/last/15th day of teh month
  month=01, year=2020,
  hour=0, min=0, sec=0 ,
  TMF=TRUE
){
  #Check DTC is present in input dataset
  assert_has_variables(dataset, deparse(substitute(dtc)))

  tmp<-dataset
  warn_if_vars_exist(dataset, paste0(deparse(substitute(new_vars_prefix)),"DTM") )
  tmp<-tmp %>%
    mutate(dtm=dtc_dtm(dtc=compute_imputed_dtc(dtc=!!enquo(dtc),
                                               impute=impute,
                                               year=year, month=month, day=day,
                                               hour=hour, min=min, sec=sec)
                       )
    )
   if(TMF==TRUE){
    warn_if_vars_exist(dataset, paste0(deparse(substitute(new_vars_prefix)),"TMF") )
    tmp<-tmp %>%
      mutate(tmf=compute_tmf(dtc=!!enquo(dtc), dt=dtm))

    names(tmp)[names(tmp) == "tmf"] <- paste0(new_vars_prefix,"TMF")
  }
  names(tmp)[names(tmp) == "dtm"] <- paste0(new_vars_prefix,"DTM")

  return(tmp)

}


