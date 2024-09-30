## ----setup, include = FALSE----------------------------------------------------------------------------------------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

library(admiraldev)


## ----message=FALSE, warning=FALSE----------------------------------------------------------------------------------------------------------------------------------
library(admiral)
library(dplyr, warn.conflicts = FALSE)
library(pharmaversesdtm)
library(lubridate)
library(stringr)

dm <- pharmaversesdtm::dm
ds <- pharmaversesdtm::ds
ex <- pharmaversesdtm::ex
ae <- pharmaversesdtm::ae
lb <- pharmaversesdtm::lb

dm <- convert_blanks_to_na(dm)
ds <- convert_blanks_to_na(ds)
ex <- convert_blanks_to_na(ex)
ae <- convert_blanks_to_na(ae)
lb <- convert_blanks_to_na(lb)


## ----eval=TRUE-----------------------------------------------------------------------------------------------------------------------------------------------------
adsl <- dm %>%
  select(-DOMAIN)


## ----eval=TRUE, echo=FALSE-----------------------------------------------------------------------------------------------------------------------------------------
z =dataset_vignette(
  adsl,
  display_vars = exprs(USUBJID, RFSTDTC, COUNTRY, AGE, SEX, RACE, ETHNIC, ARM, ACTARM)
)


## ----eval=TRUE-----------------------------------------------------------------------------------------------------------------------------------------------------
adsl <- dm %>%
  mutate(TRT01P = ARM, TRT01A = ACTARM)


## ----eval=TRUE-----------------------------------------------------------------------------------------------------------------------------------------------------
# impute start and end time of exposure to first and last respectively,
# do not impute date
ex_ext <- ex %>%
  derive_vars_dtm(
    dtc = EXSTDTC,
    new_vars_prefix = "EXST"
  ) %>%
  derive_vars_dtm(
    dtc = EXENDTC,
    new_vars_prefix = "EXEN",
    time_imputation = "last"
  )

adsl <- adsl %>%
  derive_vars_merged(
    dataset_add = ex_ext,
    filter_add = (EXDOSE > 0 |
      (EXDOSE == 0 &
        str_detect(EXTRT, "PLACEBO"))) & !is.na(EXSTDTM),
    new_vars = exprs(TRTSDTM = EXSTDTM, TRTSTMF = EXSTTMF),
    order = exprs(EXSTDTM, EXSEQ),
    mode = "first",
    by_vars = exprs(STUDYID, USUBJID)
  ) %>%
  derive_vars_merged(
    dataset_add = ex_ext,
    filter_add = (EXDOSE > 0 |
      (EXDOSE == 0 &
        str_detect(EXTRT, "PLACEBO"))) & !is.na(EXENDTM),
    new_vars = exprs(TRTEDTM = EXENDTM, TRTETMF = EXENTMF),
    order = exprs(EXENDTM, EXSEQ),
    mode = "last",
    by_vars = exprs(STUDYID, USUBJID)
  )


## ----eval=TRUE-----------------------------------------------------------------------------------------------------------------------------------------------------
adsl <- adsl %>%
  derive_vars_dtm_to_dt(source_vars = exprs(TRTSDTM, TRTEDTM))


## ----eval=TRUE-----------------------------------------------------------------------------------------------------------------------------------------------------
adsl <- adsl %>%
  derive_var_trtdurd()


## ----eval=TRUE, echo=FALSE-----------------------------------------------------------------------------------------------------------------------------------------
dataset_vignette(
  adsl,
  display_vars = exprs(USUBJID, RFSTDTC, TRTSDTM, TRTSDT, TRTEDTM, TRTEDT, TRTDURD)
)


## ----eval=TRUE-----------------------------------------------------------------------------------------------------------------------------------------------------
# convert character date to numeric date without imputation
ds_ext <- derive_vars_dt(
  ds,
  dtc = DSSTDTC,
  new_vars_prefix = "DSST"
)

adsl <- adsl %>%
  derive_vars_merged(
    dataset_add = ds_ext,
    by_vars = exprs(STUDYID, USUBJID),
    new_vars = exprs(EOSDT = DSSTDT),
    filter_add = DSCAT == "DISPOSITION EVENT" & DSDECOD != "SCREEN FAILURE"
  )


## ----eval=TRUE, echo=FALSE-----------------------------------------------------------------------------------------------------------------------------------------
dataset_vignette(
  ds_ext,
  display_vars = exprs(USUBJID, DSCAT, DSDECOD, DSTERM, DSSTDT, DSSTDTC),
  filter = DSDECOD != "SCREEN FAILURE"
)


## ----eval=TRUE, echo=FALSE-----------------------------------------------------------------------------------------------------------------------------------------
dataset_vignette(adsl, display_vars = exprs(USUBJID, EOSDT))


## ----eval=TRUE-----------------------------------------------------------------------------------------------------------------------------------------------------
format_eosstt <- function(x) {
  case_when(
    x %in% c("COMPLETED") ~ "COMPLETED",
    x %in% c("SCREEN FAILURE") ~ NA_character_,
    TRUE ~ "DISCONTINUED"
  )
}


## ----eval=TRUE-----------------------------------------------------------------------------------------------------------------------------------------------------
adsl <- adsl %>%
  derive_vars_merged(
    dataset_add = ds,
    by_vars = exprs(STUDYID, USUBJID),
    filter_add = DSCAT == "DISPOSITION EVENT",
    new_vars = exprs(EOSSTT = format_eosstt(DSDECOD)),
    missing_values = exprs(EOSSTT = "ONGOING")
  )


## ----eval=TRUE, echo=FALSE-----------------------------------------------------------------------------------------------------------------------------------------
dataset_vignette(adsl, display_vars = exprs(USUBJID, EOSDT, EOSSTT))


## ------------------------------------------------------------------------------------------------------------------------------------------------------------------
adsl <- adsl %>%
  derive_vars_merged(
    dataset_add = ds,
    by_vars = exprs(USUBJID),
    new_vars = exprs(DCSREAS = DSDECOD, DCSREASP = DSTERM),
    filter_add = DSCAT == "DISPOSITION EVENT" &
      !(DSDECOD %in% c("SCREEN FAILURE", "COMPLETED", NA))
  )


## ----eval=TRUE, echo=FALSE-----------------------------------------------------------------------------------------------------------------------------------------
dataset_vignette(adsl, display_vars = exprs(USUBJID, EOSDT, EOSSTT, DCSREAS, DCSREASP))


## ----eval=TRUE, echo=FALSE-----------------------------------------------------------------------------------------------------------------------------------------
adsl <- adsl %>%
  select(-DCSREAS, -DCSREASP)


## ------------------------------------------------------------------------------------------------------------------------------------------------------------------
adsl <- adsl %>%
  derive_vars_merged(
    dataset_add = ds,
    by_vars = exprs(USUBJID),
    new_vars = exprs(DCSREAS = DSDECOD),
    filter_add = DSCAT == "DISPOSITION EVENT" &
      DSDECOD %notin% c("SCREEN FAILURE", "COMPLETED", NA)
  ) %>%
  derive_vars_merged(
    dataset_add = ds,
    by_vars = exprs(USUBJID),
    new_vars = exprs(DCSREASP = DSTERM),
    filter_add = DSCAT == "DISPOSITION EVENT" & DSDECOD %in% "OTHER"
  )


## ----eval=TRUE, echo=FALSE-----------------------------------------------------------------------------------------------------------------------------------------
dataset_vignette(adsl, display_vars = exprs(USUBJID, EOSDT, EOSSTT, DCSREAS, DCSREASP))


## ----eval=TRUE-----------------------------------------------------------------------------------------------------------------------------------------------------
adsl <- adsl %>%
  derive_vars_merged(
    dataset_add = ds_ext,
    filter_add = DSDECOD == "RANDOMIZED",
    by_vars = exprs(STUDYID, USUBJID),
    new_vars = exprs(RANDDT = DSSTDT)
  )


## ----eval=TRUE, echo=FALSE-----------------------------------------------------------------------------------------------------------------------------------------
dataset_vignette(adsl, display_vars = exprs(USUBJID, RANDDT))


## ----eval=TRUE-----------------------------------------------------------------------------------------------------------------------------------------------------
adsl <- adsl %>%
  derive_vars_dt(
    new_vars_prefix = "DTH",
    dtc = DTHDTC
  )


## ----eval=TRUE, echo=FALSE-----------------------------------------------------------------------------------------------------------------------------------------
dataset_vignette(adsl %>% filter(!is.na(DTHDT) | row_number() %% 50 == 0), display_vars = exprs(USUBJID, TRTEDT, DTHDTC, DTHDT, DTHFL))


## ----eval=FALSE----------------------------------------------------------------------------------------------------------------------------------------------------
## adsl <- adsl %>%
##   derive_vars_dt(
##     new_vars_prefix = "DTH",
##     dtc = DTHDTC,
##     date_imputation = "first"
##   )


## ----eval=TRUE-----------------------------------------------------------------------------------------------------------------------------------------------------
adsl <- adsl %>%
  derive_vars_extreme_event(
    by_vars = exprs(STUDYID, USUBJID),
    events = list(
      event(
        dataset_name = "ae",
        condition = AEOUT == "FATAL",
        set_values_to = exprs(DTHCAUS = AEDECOD),
      ),
      event(
        dataset_name = "ds",
        condition = DSDECOD == "DEATH" & grepl("DEATH DUE TO", DSTERM),
        set_values_to = exprs(DTHCAUS = DSTERM),
      )
    ),
    source_datasets = list(ae = ae, ds = ds),
    tmp_event_nr_var = event_nr,
    order = exprs(event_nr),
    mode = "first",
    new_vars = exprs(DTHCAUS)
  )


## ----eval=TRUE, echo=FALSE-----------------------------------------------------------------------------------------------------------------------------------------
dataset_vignette(
  adsl,
  display_vars = exprs(USUBJID, EOSDT, DTHDTC, DTHDT, DTHCAUS),
  filter = DTHFL == "Y"
)


## ----eval=TRUE-----------------------------------------------------------------------------------------------------------------------------------------------------
adsl <- adsl %>%
  select(-DTHCAUS) %>% # remove it before deriving it again
  derive_vars_extreme_event(
    by_vars = exprs(STUDYID, USUBJID),
    events = list(
      event(
        dataset_name = "ae",
        condition = AEOUT == "FATAL",
        set_values_to = exprs(DTHCAUS = AEDECOD, DTHDOM = "AE", DTHSEQ = AESEQ),
      ),
      event(
        dataset_name = "ds",
        condition = DSDECOD == "DEATH" & grepl("DEATH DUE TO", DSTERM),
        set_values_to = exprs(DTHCAUS = DSTERM, DTHDOM = "DS", DTHSEQ = DSSEQ),
      )
    ),
    source_datasets = list(ae = ae, ds = ds),
    tmp_event_nr_var = event_nr,
    order = exprs(event_nr),
    mode = "first",
    new_vars = exprs(DTHCAUS, DTHDOM, DTHSEQ)
  )


## ----eval=TRUE, echo=FALSE-----------------------------------------------------------------------------------------------------------------------------------------
dataset_vignette(
  adsl,
  display_vars = exprs(USUBJID, TRTEDT, DTHDTC, DTHDT, DTHCAUS, DTHDOM, DTHSEQ),
  filter = DTHFL == "Y"
)


## ----eval=TRUE-----------------------------------------------------------------------------------------------------------------------------------------------------
adsl <- adsl %>%
  mutate(DTHCGR1 = case_when(
    is.na(DTHDOM) ~ NA_character_,
    DTHDOM == "AE" ~ "ADVERSE EVENT",
    str_detect(DTHCAUS, "(PROGRESSIVE DISEASE|DISEASE RELAPSE)") ~ "PROGRESSIVE DISEASE",
    TRUE ~ "OTHER"
  ))


## ----eval=TRUE-----------------------------------------------------------------------------------------------------------------------------------------------------
adsl <- adsl %>%
  derive_vars_duration(
    new_var = DTHADY,
    start_date = TRTSDT,
    end_date = DTHDT
  )


## ----eval=TRUE-----------------------------------------------------------------------------------------------------------------------------------------------------
adsl <- adsl %>%
  derive_vars_duration(
    new_var = LDDTHELD,
    start_date = TRTEDT,
    end_date = DTHDT,
    add_one = FALSE
  )


## ----eval=TRUE, echo=FALSE-----------------------------------------------------------------------------------------------------------------------------------------
dataset_vignette(
  adsl,
  display_vars = exprs(USUBJID, TRTEDT, DTHDTC, DTHDT, DTHCAUS, DTHADY, LDDTHELD),
  filter = DTHFL == "Y"
)


## ----eval=TRUE-----------------------------------------------------------------------------------------------------------------------------------------------------
adsl <- adsl %>%
  derive_vars_extreme_event(
    by_vars = exprs(STUDYID, USUBJID),
    events = list(
      event(
        dataset_name = "ae",
        order = exprs(AESTDTC, AESEQ),
        condition = !is.na(AESTDTC),
        set_values_to = exprs(
          LSTALVDT = convert_dtc_to_dt(AESTDTC, highest_imputation = "M"),
          seq = AESEQ
        ),
      ),
      event(
        dataset_name = "ae",
        order = exprs(AEENDTC, AESEQ),
        condition = !is.na(AEENDTC),
        set_values_to = exprs(
          LSTALVDT = convert_dtc_to_dt(AEENDTC, highest_imputation = "M"),
          seq = AESEQ
        ),
      ),
      event(
        dataset_name = "lb",
        order = exprs(LBDTC, LBSEQ),
        condition = !is.na(LBDTC),
        set_values_to = exprs(
          LSTALVDT = convert_dtc_to_dt(LBDTC, highest_imputation = "M"),
          seq = LBSEQ
        ),
      ),
      event(
        dataset_name = "adsl",
        condition = !is.na(TRTEDT),
        set_values_to = exprs(LSTALVDT = TRTEDT, seq = 0),
      )
    ),
    source_datasets = list(ae = ae, lb = lb, adsl = adsl),
    tmp_event_nr_var = event_nr,
    order = exprs(LSTALVDT, seq, event_nr),
    mode = "last",
    new_vars = exprs(LSTALVDT)
  )


## ----eval=TRUE, echo=FALSE-----------------------------------------------------------------------------------------------------------------------------------------
dataset_vignette(
  adsl,
  display_vars = exprs(USUBJID, TRTEDT, DTHDTC, LSTALVDT),
  filter = !is.na(TRTSDT)
)


## ----eval=TRUE-----------------------------------------------------------------------------------------------------------------------------------------------------
adsl <- adsl %>%
  select(-LSTALVDT) %>% # created in the previous call
  derive_vars_extreme_event(
    by_vars = exprs(STUDYID, USUBJID),
    events = list(
      event(
        dataset_name = "ae",
        order = exprs(AESTDTC, AESEQ),
        condition = !is.na(AESTDTC),
        set_values_to = exprs(
          LSTALVDT = convert_dtc_to_dt(AESTDTC, highest_imputation = "M"),
          LALVSEQ = AESEQ,
          LALVDOM = "AE",
          LALVVAR = "AESTDTC"
        ),
      ),
      event(
        dataset_name = "ae",
        order = exprs(AEENDTC, AESEQ),
        condition = !is.na(AEENDTC),
        set_values_to = exprs(
          LSTALVDT = convert_dtc_to_dt(AEENDTC, highest_imputation = "M"),
          LALVSEQ = AESEQ,
          LALVDOM = "AE",
          LALVVAR = "AEENDTC"
        ),
      ),
      event(
        dataset_name = "lb",
        order = exprs(LBDTC, LBSEQ),
        condition = !is.na(LBDTC),
        set_values_to = exprs(
          LSTALVDT = convert_dtc_to_dt(LBDTC, highest_imputation = "M"),
          LALVSEQ = LBSEQ,
          LALVDOM = "LB",
          LALVVAR = "LBDTC"
        ),
      ),
      event(
        dataset_name = "adsl",
        condition = !is.na(TRTEDT),
        set_values_to = exprs(LSTALVDT = TRTEDT, LALVSEQ = NA_integer_, LALVDOM = "ADSL", LALVVAR = "TRTEDTM"),
      )
    ),
    source_datasets = list(ae = ae, lb = lb, adsl = adsl),
    tmp_event_nr_var = event_nr,
    order = exprs(LSTALVDT, LALVSEQ, event_nr),
    mode = "last",
    new_vars = exprs(LSTALVDT, LALVSEQ, LALVDOM, LALVVAR)
  )


## ----eval=TRUE, echo=FALSE-----------------------------------------------------------------------------------------------------------------------------------------
dataset_vignette(
  adsl,
  display_vars = exprs(USUBJID, TRTEDT, DTHDTC, LSTALVDT, LALVDOM, LALVSEQ, LALVVAR),
  filter = !is.na(TRTSDT)
)


## ----eval=TRUE-----------------------------------------------------------------------------------------------------------------------------------------------------
format_agegr1 <- function(var_input) {
  case_when(
    var_input < 18 ~ "<18",
    between(var_input, 18, 64) ~ "18-64",
    var_input > 64 ~ ">64",
    TRUE ~ "Missing"
  )
}

format_region1 <- function(var_input) {
  case_when(
    var_input %in% c("CAN", "USA") ~ "North America",
    !is.na(var_input) ~ "Rest of the World",
    TRUE ~ "Missing"
  )
}


## ----eval=TRUE-----------------------------------------------------------------------------------------------------------------------------------------------------
adsl <- adsl %>%
  mutate(
    AGEGR1 = format_agegr1(AGE),
    REGION1 = format_region1(COUNTRY)
  )


## ----eval=TRUE, echo=FALSE-----------------------------------------------------------------------------------------------------------------------------------------
dataset_vignette(
  adsl,
  display_vars = exprs(USUBJID, AGE, SEX, COUNTRY, AGEGR1, REGION1)
)


## ----eval=TRUE-----------------------------------------------------------------------------------------------------------------------------------------------------
adsl <- adsl %>%
  derive_var_merged_exist_flag(
    dataset_add = ex,
    by_vars = exprs(STUDYID, USUBJID),
    new_var = SAFFL,
    condition = (EXDOSE > 0 | (EXDOSE == 0 & str_detect(EXTRT, "PLACEBO")))
  )


## ----eval=TRUE, echo=FALSE-----------------------------------------------------------------------------------------------------------------------------------------
dataset_vignette(
  adsl,
  display_vars = exprs(USUBJID, TRTSDT, ARM, ACTARM, SAFFL)
)

