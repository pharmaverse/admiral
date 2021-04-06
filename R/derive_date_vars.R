#' Impute partial date/time portion of a --DTC variable
#'
#' Imputation partial date/time portion of a --DTC variable
#' based on user input
#'
#' @param dtc The --DTC date to impute
#'
#'   A character date is expected in a format like yyyy-mm-dd or yyyy-mm-ddThh:mm:ss
#'   if the year part is not recorded (missing date), no imputation is done
#'
#' @param date_imputation The value to impute the day/month when a datepart is missing
#'
#'   If NULL: no date imputation is perfomed and partial dates are returned as missing
#'
#'   Otherwise, a character value is expected, either as a
#'   - format with day and month specified as 'dd-mm': e.g. '15-06' for the 15th of June
#'   - or as a keyword: 'FIRST', 'MID', 'LAST' to impute to the first/mid/last day/month
#'
#'   Default is NULL
#'
#' @param time_imputation The value to impute the time when a timepart is missing
#'
#'   A character value is expected, either as a
#'   - format with hour, min and sec specified as 'hh:mm:ss': e.g. '00:00:00' for the start of the day
#'   - or as a keyword: 'FIRST','LAST' to impute to the start/end of a day
#'
#'   Default is '00:00:00'
#'
#' @author Samia Kabi
#'
#' @return  a character vector
#'
#' @family {general functions}
#'
#' @export
#'
#' @examples
#'
#' mhdt <- tibble::tribble(
#' ~MHSTDTC,
#' "2019-07-18T15:25:40",
#' "2019-07-18T15:25",
#' "2019-07-18T15",
#' "2019-07-18",
#' "2019-02",
#' "2019",
#' "2019",
#' "2019---07",
#' "")
#' #No date imputation (date_imputation defaulted to NULL)
#' #Missing time part imputed with 00:00:00 portion by default
#' idtc<-impute_dtc(dtc=mhdt$MHSTDTC)
#'
#' #No date imputation (date_imputation defaulted to NULL)
#' #Missing time part imputed with 23:59:59 portion
#' idtc<-impute_dtc(dtc=mhdt$MHSTDTC,
#'                 time_imputation = "23:59:59")
#' #same as above
#' idtc<-impute_dtc(dtc=mhdt$MHSTDTC,
#'                  time_imputation = "LAST")
#'
#' #impute to first day/month if date is partial
#' #Missing time part imputed with 00:00:00 portion by default
#' idtc<-impute_dtc(dtc=mhdt$MHSTDTC,
#'                  date_imputation="01-01")
#' #same as above
#' idtc<-impute_dtc(dtc=mhdt$MHSTDTC,
#'                  date_imputation="FIRST")
#'
#' #impute to last day/month if date is partial
#' #Missing time part imputed with 23:59:59 portion
#' idtc<-impute_dtc(dtc=mhdt$MHSTDTC,
#'                  date_imputation="LAST",
#'                  time_imputation = "LAST")
#'
#' #impute to Mid day/month if date is partial
#' #Missing time part imputed with 00:00:00 portion by default
#' idtc<-impute_dtc(dtc=mhdt$MHSTDTC,
#'                   date_imputation="MID")
#'

impute_dtc<-function(dtc,
                     date_imputation=NULL,
                     time_imputation="00:00:00"
                     )
  {

  # check format of DTC is as expected
  assert_is_valid_dtc(dtc)

  #date imputation
  if (!is.null(date_imputation)){
    #check input for date_imputation
    assert_that(is_valid_date_entry(date_imputation))
    #Specific setup for FISRT/MID/LAST
    if(date_imputation =="FIRST") {
      d<-"01"
      mo<-"01"
    }
    else if(date_imputation== "MID") {
      d<-"15"
      mo<-"06"
    }
    else if(date_imputation =="LAST") {
      d<-"01"
      mo<-"12"
    }
    #otherwise, use time_imputation input
    else {
      mo__=as.integer(sub(".*-","",date_imputation))
      day__=as.integer(sub("-.*","",date_imputation))
      #check input for day and moth are valid
      assert_that(is_valid_day(day__))
      assert_that(is_valid_month(mo__))

      d<-sprintf("%02d", day__)
      mo<-sprintf("%02d", mo__)
    }
    imputed_date <- case_when(
      nchar(dtc)>=10 ~ substr(dtc,1,10),
      nchar(dtc)==9 ~  paste0(substr(dtc,1,4), "-",mo,"-",d),#dates like 2021---14 - use only year part
      nchar(dtc)==7 ~  paste0(dtc,"-",d),
      nchar(dtc)==4  ~  paste0(dtc, "-",mo,"-",d),
      TRUE ~ ""
    )
    if (date_imputation== "LAST"){
      imputed_date<-case_when(
        nchar(imputed_date)>0 & nchar(dtc)<10 ~ as.character(ceiling_date(as.Date(imputed_date, format = "%Y-%m-%d"), "month") - days(1)),
        TRUE~imputed_date
      )
    }
  }
  #no imputation
  else{
    imputed_date <- case_when(
      nchar(dtc)>=10 ~ substr(dtc,1,10),
      TRUE ~ ""
    )
  }

  #impute time
  assert_that(is_valid_time_entry(time_imputation))
  if(time_imputation =="FIRST") {
    imputed_time<-"00:00:00"
    sec<-":00"
    min<-":00"
    h<-"00"
  }
  else {
    if(time_imputation =="LAST") {
      imputed_time<-"23:59:59"
      sec<-":59"
      min<-":59"
      h<-"23"
    }
    else {
      imputed_time<-time_imputation
      sec<-paste0(":",paste0(substr(dtc,18,19), substr(time_imputation,7,8)))
      min<-paste0(":",paste0(substr(dtc,15,16), substr(time_imputation,4,5)))
      h<-paste0(substr(dtc,12,13),substr(time_imputation,1,2))
    }
  }

  imputed_time <- case_when(
    nchar(dtc)>=19 ~ substr(dtc,12,19),
    nchar(dtc)==16 ~ paste0(substr(dtc,12,16),sec),
    nchar(dtc)==13 ~ paste0(substr(dtc,12,13),min,sec),
    nchar(dtc)==10 ~ paste0(h,min,sec),
    TRUE ~ imputed_time
  )

  imputed_dtc<-if_else(imputed_date!="",paste0(imputed_date,"T", imputed_time),"")

  return(imputed_dtc)

}

#' Convert a --DTC to a --DT
#'
#' Convert a complete character --DTC variable into a date object
#'
#' @param dtc The --DTC date to convert
#'
#'   A character date is expected in a format like yyyy-mm-dd or yyyy-mm-ddThh:mm:ss
#'   a partial date will return a NA date and a warning will be issued: 'All formats failed to parse. No formats found.'
#'   Note: you can use impute_dtc function to build a complete date
#'
#' @author Samia Kabi
#'
#' @return  a date object
#'
#' @family {general functions}
#'
#' @export
#'
#' @examples
#'
#' convert_dtc_dt("2019-07-18")
#' convert_dtc_dt("2019-07")
#'
convert_dtc_to_dt<-function(dtc){
  case_when(
    nchar(dtc)>=10 ~ymd(substr(dtc,1,10)),
    TRUE ~ymd(NA)
  )

}

#' Convert a --DTC to a --DTM
#'
#' Convert a complete character --DTC variable into a datetime object
#'
#' @param dtc The --DTC date to convert
#'
#'   A character date is expected in a format like yyyy-mm-ddThh:mm:ss
#'   a partial datetime will return a NA date and a warning will be issued: 'All formats failed to parse. No formats found.'
#'   Note: you can use impute_dtc function to build a complete datetime
#'
#' @author Samia Kabi
#'
#' @return  a datetime  object
#'
#' @family {general functions}
#'
#' @export
#'
#' @examples
#'
#' convert_dtc_dtm("2019-07-18T15:25:00")
#' convert_dtc_dtm("2019-07-18T00:00:00")#note Time = 00:00:00 is not printed
#' convert_dtc_dtm("2019-07-18")
#'
convert_dtc_dtm<-function(dtc){

  #note T00:00:00 is not printed in dataframe
  case_when(
    nchar(dtc)==19 ~ymd_hms(dtc),
    TRUE ~ymd_hms(NA)
  )
}

#' Derive --DTF
#'
#' Derive --DTF based on --DTC and --DT
#'
#' @param dtc The --DTC used as input for --DT
#'
#'   A character date is expected in a format like yyyy-mm-ddThh:mm:ss (partial or complete)
#'
#' @param dt The --DT resulting from the imputation of --DTC
#'
#'   A date object is expected
#'
#' @author Samia Kabi
#'
#' @return  the date imputation flag --DTF (character value of 'D', 'M' , 'Y' or missing )
#'
#' @family {general functions}
#'
#' @export
#'
#' @examples
#'
#' compute_dtf(dtc= "2019-07",dt=as.Date("2019-07-18"))
#' compute_dtf(dtc= "2019",dt=as.Date("2019-07-18"))
#'
compute_dtf <- function(dtc,dt){

  case_when(
    !is.na(dt) & nchar(dtc)< 4 ~"Y",
    !is.na(dt) & nchar(dtc)==4 ~"M",
    !is.na(dt) & nchar(dtc)==7 ~"D",
    !is.na(dt) & nchar(dtc)==9 ~"M",#dates like "2019---07"
    (!is.na(dt) & nchar(dtc)>=10) | is.na(dt) ~" "
  )
}

#' Derive --TMF
#'
#' Derive --TMF based on --DTC and --DT
#'
#' @param dtc The --DTC used as input for --DT
#'
#'   A character date is expected in a format like yyyy-mm-ddThh:mm:ss (partial or complete)
#'
#' @param dtm The --DTM resulting from the imputation of --DTC
#'
#'   A datetime object is expected
#'
#' @author Samia Kabi
#'
#' @return  the time imputation flag --DTF (character value of 'H', 'M' , 'S' or missing)
#'
#' @family {general functions}
#'
#' @export
#'
#' @examples
#'
#' compute_tmf(dtc= "2019-07-18T15:25",dtm=as.POSIXct("2019-07-18T15:25:00"))
#' compute_tmf(dtc= "2019-07-18T15",dtm=as.POSIXct("2019-07-18T15:25:00"))
#' compute_tmf(dtc= "2019-07-18",dtm=as.POSIXct("2019-07-18"))
#'
compute_tmf <- function(dtc,dtm){

  case_when(
    (!is.na(dtm) & nchar(dtc)>=19)|is.na(dtm) ~" ",
    !is.na(dtm) & nchar(dtc)==16~"S",
    !is.na(dtm) & nchar(dtc)==13~"M",
    (!is.na(dtm) & (nchar(dtc)==10|is_valid_dtc(dtc)==FALSE)) |
      (nchar(dtc)> 0 & nchar(dtc) <10) ~"H"
  )
}


#' Derive --DT (and --DTF)
#'
#' Derive --DT based on --DTC. --DT is imputed based on user input
#' Derive --DTF if needed based on --DTC and --DT
#'
#' @param dataset Input dataset
#'
#'   The --DTC variable must be present.
#'
#' @param new_vars_prefix Prefix used for the output variable(s)
#'
#' a character is expected: e.g new_vars_prefix="AST"
#'
#' @param dtc The --DTC date used to derive/impute --DT
#'
#'   A character date is expected in a format like yyyy-mm-dd or yyyy-mm-ddThh:mm:ss
#'   if the year part is not recorded (missing date), no imputation is done
#'
#' @param date_imputation The value to impute the day/month when a datepart is missing
#'
#'   If NULL: no date imputation is perfomed and partial dates are returned as missing
#'
#'   Otherwise, a character value is expected, either as a
#'   - format with day and month specified as 'dd-mm': e.g. '15-06' for the 15th of June
#'   - or as a keyword: 'FIRST', 'MID', 'LAST' to impute to the first/mid/last day/month
#'
#'   Default is NULL
#'
#' @param flag_imputation Whether the --DTF variable must also be derived
#'
#' A logical value
#'
#' Default: TRUE
#'
#'@examples
#'
#'mhdt <- tibble::tribble(~MHSTDTC,
#'"2019-07-18T15:25:40",
#'"2019-07-18T15:25",
#'"2019-07-18",
#'"2019-02",
#'"2019",
#'"2019---07",
#'"")
#'
#'#Create ASTDT and ASTDTF
#'#no imputation for partial date
#'mhdt1<-derive_vars_dt(mhdt,
#'                      new_vars_prefix = "AST",
#'                      dtc=MHSTDTC)
#'
#'#Create ASTDT and ASTDTF
#'#Impute partial dates to first day/month
#'mhdt2<-derive_vars_dt(mhdt,
#'                      new_vars_prefix = "AST",
#'                      dtc=MHSTDTC,
#'                      date_imputation = "FIRST")
#'#Impute partial dates to 4th of june
#'mhdt2a<-derive_vars_dt(mhdt,
#'                       new_vars_prefix = "AST",
#'                       dtc=MHSTDTC,
#'                       date_imputation = "04-06")
#'#Create AENDT and AENDTF
#'#Impute partial dates to last day/month
#'mhdt3<-derive_vars_dt(mhdt,
#'                      new_vars_prefix = "AEN",
#'                      dtc=MHSTDTC,
#'                      date_imputation="LAST"
#'                      )
#'#Create BIRTHDT
#'#Impute partial dates to 15of june. No DTF
#'mhdt4<-derive_vars_dt(mhdt,
#'                      new_vars_prefix = "BIRTH",
#'                      dtc=MHSTDTC,
#'                      date_imputation="MID",
#'                      flag_imputation=FALSE
#'                      )

derive_vars_dt<- function(
  dataset,
  new_vars_prefix,
  dtc,
  date_imputation=NULL,# "02-01" or "LAST"
  flag_imputation=TRUE
){
  #Check DTC is present in input dataset
  assert_has_variables(dataset, deparse(substitute(dtc)))

  #output varname
  dt <- paste0(new_vars_prefix, "DT")
  warn_if_vars_exist(dataset,dt)

  #derive --DT var
  dataset<-dataset %>%
    mutate(
      idtc__=impute_dtc(dtc=!!enquo(dtc),
                            date_imputation=date_imputation),
      !!sym(dt):=convert_dtc_dt(dtc=idtc__)
    )%>%
    select(-ends_with(("__")))


  #derive DTF
  if(flag_imputation){
    dtf <- paste0(new_vars_prefix, "DTF")
    warn_if_vars_exist(dataset, dtf)
    dataset<-dataset  %>%
        mutate(!!sym(dtf) := compute_dtf(dtc=!!enquo(dtc), dt=!!sym(dt)))
  }
  return(dataset)
}


#' Derive --DTM (and --DTF/--TMF)
#'
#' Derive --DTM based on --DTC. --DTM is imputed based on user input
#' Derive --DTF/--TMF if needed based on --DTC and --DT
#'
#' @param dataset Input dataset
#'
#'   The --DTC variable must be present.
#'
#' @param new_vars_prefix Prefix used for the output variable(s)
#'
#' a character is expected: e.g new_vars_prefix="AST"
#'
#' @param dtc The --DTC date used to derive/impute --DT
#'
#'   A character date is expected in a format like yyyy-mm-dd or yyyy-mm-ddThh:mm:ss
#'   if the year part is not recorded (missing date), no imputation is done
#'
#' @param date_imputation The value to impute the day/month when a datepart is missing
#'
#'   If NULL: no date imputation is perfomed and partial dates are returned as missing
#'
#'   Otherwise, a character value is expected, either as a
#'   - format with day and month specified as 'dd-mm': e.g. '15-06' for the 15th of June
#'   - or as a keyword: 'FIRST', 'MID', 'LAST' to impute to the first/mid/last day/month
#'
#'   Default is NULL
#'
#' @param time_imputation The value to impute the time when a timepart is missing
#'
#'   A character value is expected, either as a
#'   - format with hour, min and sec specified as 'hh:mm:ss': e.g. '00:00:00' for the start of the day
#'   - or as a keyword: 'FIRST','LAST' to impute to the start/end of a day
#'
#'   Default is '00:00:00'
#'

#' @param flag_imputation Whether the --DTF and or --TMF variable must also be derived
#'
#' A logical value
#'
#' Default: TRUE
#'
#' @details:
#'
#' the presence of --DTF is checked and the variable is not derived if it already exists in the input datasets
#' however, if --TMF already exists in the dataset, a warning is issued and --TMF will be overwritten
#'
#'@examples
#'
#'mhdt <- tibble::tribble(~MHSTDTC,
#'"2019-07-18T15:25:40",
#'"2019-07-18T15:25",
#'"2019-07-18",
#'"2019-02",
#'"2019",
#'"2019---07",
#'"")
#'
derive_vars_dtm<- function(
  dataset,
  new_vars_prefix,
  dtc,
  date_imputation=NULL,# "02-01" or "LAST"
  time_imputation="00:00:00",#or 'FIRST' 'LAST'
  flag_imputation=TRUE
){

  #Check DTC is present in input dataset
  assert_has_variables(dataset, deparse(substitute(dtc)))

  dtm <- paste0(new_vars_prefix, "DTM")

  #Issue a warning if --DTM already exists
  warn_if_vars_exist(dataset, dtm)
  dataset<-dataset %>%
    mutate(
      idtc__=impute_dtc(dtc=!!enquo(dtc),
                        date_imputation=date_imputation,
                        time_imputation=time_imputation
                                  ),
      !!sym(dtm):=convert_dtc_dtm(dtc=idtc__)
      )%>%
    select(-ends_with(("__")))

  if(flag_imputation){
    dtf<-paste0(new_vars_prefix, "DTF")
    tmf <- paste0(new_vars_prefix, "TMF")

    #add --DTF if not there already
    dtf_exist<- dtf %in% colnames(dataset)
    if (!dtf_exist){
      dataset<-dataset %>%
        mutate(!!sym(dtf):=compute_dtf(dtc=!!enquo(dtc), dt=!!sym(dtm)))

    }
    else{
      print(paste0(dtf, " was already present in dataset ", dtf , " is not derived."  ))
    }
    #add --TMF
    warn_if_vars_exist(dataset, tmf )
    dataset<-dataset %>%
       mutate(!!sym(tmf):=compute_tmf(dtc=!!enquo(dtc), dtm=!!sym(dtm)))
   }

  return(dataset)

}



