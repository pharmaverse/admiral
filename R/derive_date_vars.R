compute_imputed_dtc<-function(dtc,
                              impute=TRUE,
                              day=01, month=01,year=2020,
                              hour=0, min=0, sec=0)
  {

  #check input for imputed day, month and year
  assert_that(is_valid_sec_min(sec))
  assert_that(is_valid_sec_min(min))
  assert_that(is_valid_hour(hour))
  assert_that(is_valid_day(day))
  assert_that(is_valid_month(month))
  assert_that(is_valid_year(year))
  assert_that(is.logical(impute))

  if(day != "LAST") {
    d<-sprintf("%02d", as.integer(day))}
  else {d<-"01"}
  mo<-sprintf("%02d", as.integer(month))
  h<-sprintf("%02d", as.integer(hour))
  mi<-sprintf("%02d", as.integer(min))
  s<-sprintf("%02d", as.integer(sec))

  if (impute==TRUE){
    #add missing text in date/time
    tmp__ <- case_when(
         nchar(dtc)>10 ~ dtc,
         nchar(dtc)==10 ~ paste0(dtc,"T",h,":", mi,":", s),
         nchar(dtc)==9 ~  paste0(substr(dtc,1,4), "-",mo,"-",d,"T",h,":", mi,":", s),
         nchar(dtc)==7 ~  paste0(dtc,"-",d,"T",h,":", mi,":", s),
         nchar(dtc)==4  ~  paste0(dtc, "-",mo,"-",d,"T",h,":", mi,":", s),
         TRUE ~ as.character(dtc)
       )
    if (day=="LAST"){
      #last day of the month
      tmp__<-paste0(as.character(ceiling_date(as.Date(tmp__), "month") - days(1))
                    ,"T",h,":", mi,":", s)
    }
  }
  #No imputation: keep the datetime or teh date with a time set to 0
  else{
    tmp__ <- case_when(
      nchar(dtc)>10 ~ dtc,
      nchar(dtc)==10 ~ paste0(dtc,"T00:00:00"),
      TRUE ~ "" #as.character(dtc)
    )
  }
  return(tmp__)
}

dtc_dt<-function(dtc){

  ymd(substr(dtc,1,10))
}

dtc_dtm<-function(dtc){

  #note T00:00:00 is not printed in dataframe
  ymd_hms(dtc)
}


compute_dtf <- function(dtc,dt){

  flag<-case_when(
    (!is.na(dt) & nchar(dtc)>=10)| is.na(dt) ~" ",
    !is.na(dt) & nchar(dtc)==7~"D",
    !is.na(dt) & nchar(dtc)==4 ~"M"
  )
  return(flag)
}

compute_tmf <- function(dtc,dt){

  flag<-case_when(
    (!is.na(dt) & nchar(dtc)==19)| is.na(dt) ~" ",
    !is.na(dt) & nchar(dtc)==15~"S",
    !is.na(dt) & nchar(dtc)==13~"M",
    !is.na(dt) & nchar(dtc)==10~"H"
  )
  return(flag)
}

derive_date_vars<- function(
                      dataset,
                      new_vars_prefix,
                      dtc,
                      impute=TRUE,
                      day=01,#use "LAST"  to have the last day of teh month
                      month=01, year=2020,
                      hour=0, min=0, sec=0 ,
                      DT=TRUE, DTM=TRUE, DTF=TRUE, TMF=TRUE
                    ){
  #Check DTC is present in input dataset
  dtcvar<-str_replace_all(deparse(expr(!!enquo(dtc))),'~','')
  assert_has_variables(dataset, dtcvar)

  tmp<-dataset
  if(DT==TRUE){
    warn_if_vars_exist(dataset, paste0(deparse(substitute(new_vars_prefix)),"DT") )
    tmp<-tmp %>%
      mutate(dt=dtc_dt(compute_imputed_dtc(!!enquo(dtc), impute=impute,year=year, month=month, day=day, hour=hour, min=min, sec=sec))
      )
    names(tmp)[names(tmp) == "dt"] <- paste0(new_vars_prefix,"DT")
    }
  if(DTM==TRUE){
    warn_if_vars_exist(dataset, paste0(deparse(substitute(new_vars_prefix)),"DTM") )
    tmp<-tmp %>%
      mutate(dtm = dtc_dtm(compute_imputed_dtc(!!enquo(dtc), impute=impute,year=year, month=month, day=day, hour=hour, min=min, sec=sec))
    )
    names(tmp)[names(tmp) == "dtm"] <- paste0(new_vars_prefix,"DTM")
  }
  if(DTF==TRUE){
   assert_has_a_date_variable(tmp, c( paste0(new_vars_prefix,"DT"), paste0(new_vars_prefix,"DTM")))
   warn_if_vars_exist(dataset, paste0(deparse(substitute(new_vars_prefix)),"DTF") )
   #Use DT var if availabe
   if (DT==TRUE) {
     tmp<-tmp %>%
       mutate(dtf=compute_dtf(!!enquo(dtc), as.name(paste0(new_vars_prefix,"DT"))))
    }
    else {
     if (DTM==TRUE){
       tmp<-tmp %>%
         mutate(dtf=compute_dtf(!!enquo(dtc), as.name(paste0(new_vars_prefix,"DTM"))))
     }
   }
    names(tmp)[names(tmp) == "dtf"] <- paste0(new_vars_prefix,"DTF")
  }
  if(TMF==TRUE){
    assert_has_a_date_variable(tmp,  paste0(new_vars_prefix,"DTM"))
    warn_if_vars_exist(dataset, paste0(deparse(substitute(new_vars_prefix)),"TMF") )
    tmp<-tmp %>%
      mutate(tmf=compute_tmf(!!enquo(dtc), as.name(paste0(new_vars_prefix,"DTM"))))

    names(tmp)[names(tmp) == "tmf"] <- paste0(new_vars_prefix,"TMF")
  }

  return(tmp)

}

