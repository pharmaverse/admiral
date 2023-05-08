## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

library(admiraldev)

## ---- warning=FALSE, message=FALSE--------------------------------------------
library(admiral)
library(dplyr, warn.conflicts = FALSE)

## ---- eval = TRUE-------------------------------------------------------------
data(admiral_adlb)
adlb <- admiral_adlb %>%
  filter(PARAMCD %in% c("AST", "ALT", "BILI") & is.na(DTYPE))

## ---- echo = FALSE------------------------------------------------------------
head(adlb) %>%
  dataset_vignette()

## -----------------------------------------------------------------------------
adlb_annotated <- adlb %>%
  mutate(
    CRIT1 = case_when(
      PARAMCD == "AST" ~ "AST >=3xULN",
      PARAMCD == "ALT" ~ "ALT >=3xULN",
      PARAMCD == "BILI" ~ "BILI >=2xULN"
    ),
    CRIT1FL = if_else(
      (AVAL / ANRHI >= 3 & PARAMCD %in% c("AST", "ALT")) |
        (AVAL / ANRHI >= 2 & PARAMCD == "BILI"),
      "Y",
      NA_character_
    )
  ) %>%
  select(STUDYID, USUBJID, TRT01A, PARAMCD, LBSEQ, ADT, AVISIT, ADY, AVAL, ANRHI, CRIT1, CRIT1FL)

## ---- echo = FALSE------------------------------------------------------------
dataset_vignette(adlb_annotated)

## ---- warning = FALSE---------------------------------------------------------
altast_records <- adlb_annotated %>%
  filter(PARAMCD %in% c("AST", "ALT"))

bili_records <- adlb_annotated %>%
  filter(PARAMCD %in% c("BILI"))

hylaw_records <- derive_vars_joined(
  dataset = altast_records,
  dataset_add = bili_records,
  by_vars = exprs(STUDYID, USUBJID),
  order = exprs(ADY),
  filter_join = 0 <= ADT.join - ADT & ADT.join - ADT <= 14 & CRIT1FL == "Y" & CRIT1FL.join == "Y",
  new_vars = exprs(BILI_DT = ADT, BILI_CRITFL = CRIT1FL),
  mode = "first"
)

## ---- echo = FALSE------------------------------------------------------------
hylaw_records %>%
  arrange(desc(BILI_CRITFL), desc(CRIT1FL)) %>%
  dataset_vignette()

## -----------------------------------------------------------------------------
hylaw_records_pts_visits <- hylaw_records %>%
  select(STUDYID, USUBJID, TRT01A) %>% # add AVISIT, ADT for by visit
  distinct()

hylaw_records_fls <- hylaw_records %>%
  select(STUDYID, USUBJID, TRT01A, CRIT1FL, BILI_CRITFL) %>% # add AVISIT, ADT for by visit
  distinct()

hylaw_params <- derive_param_exist_flag(
  dataset_adsl = hylaw_records_pts_visits,
  dataset_add = hylaw_records_fls,
  condition = CRIT1FL == "Y" & BILI_CRITFL == "Y",
  false_value = "N",
  missing_value = "N",
  subject_keys = exprs(STUDYID, USUBJID, TRT01A), # add AVISIT, ADT for by visit
  set_values_to = exprs(
    PARAMCD = "HYSLAW",
    PARAM = "ALT/AST >= 3xULN and BILI >= 2xULN"
  )
)

## ---- echo = FALSE------------------------------------------------------------
hylaw_params %>%
  arrange(desc(AVAL)) %>%
  dataset_vignette()

## -----------------------------------------------------------------------------
adlbhy <- adlb_annotated %>%
  bind_rows(hylaw_params)

## ---- echo = FALSE------------------------------------------------------------
dataset_vignette(adlbhy)

