## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

library(admiraldev)

## ---- warning=FALSE, message=FALSE--------------------------------------------
library(dplyr)
library(tidyr)
library(tibble)
library(admiral)

## -----------------------------------------------------------------------------
data("example_qs")
qs <- example_qs

## ----echo=FALSE---------------------------------------------------------------
dataset_vignette(qs)

## -----------------------------------------------------------------------------
adsl <- tribble(
  ~STUDYID, ~USUBJID, ~SITEID, ~ITTFL,    ~TRTP, ~TRTSDT,                      ~DTHCAUS,
  "STUDYX",  "P0001",     13L,    "Y", "DRUG A", lubridate::ymd("2012-11-16"), NA_character_,
  "STUDYX",  "P0002",     11L,    "Y", "DRUG B", lubridate::ymd("2012-11-16"), NA_character_
)

## ----echo=FALSE---------------------------------------------------------------
dataset_vignette(adsl)

## ----eval=TRUE----------------------------------------------------------------
adqs <- qs %>%
  # Add ADSL variables
  derive_vars_merged(
    dataset_add = adsl,
    new_vars = exprs(TRTSDT, DTHCAUS),
    by_vars = exprs(STUDYID, USUBJID)
  ) %>%
  # Add analysis parameter variables
  mutate(
    PARAMCD = QSTESTCD,
    PARAM = QSTEST,
    PARCAT1 = QSCAT,
    AVALC = QSORRES,
    AVAL = QSSTRESN
  ) %>%
  # Add timing variables
  derive_vars_dt(new_vars_prefix = "A", dtc = QSDTC) %>%
  derive_vars_dy(reference_date = TRTSDT, source_vars = exprs(ADT)) %>%
  mutate(
    AVISIT = if_else(ADT <= TRTSDT, "BASELINE", VISIT),
    AVISITN = if_else(ADT <= TRTSDT, 0, VISITNUM)
  )

## ----echo=FALSE---------------------------------------------------------------
dataset_vignette(
  arrange(adqs, USUBJID, PARCAT1, ADY, PARAMCD),
  display_vars = exprs(USUBJID, PARAMCD, PARAM, PARCAT1, AVALC, AVAL, ADY, AVISIT)
)

## ----eval=TRUE----------------------------------------------------------------
adgad7 <- adqs %>%
  # Select records to keep in the GAD-7 ADaM
  filter(PARCAT1 == "GAD-7 V2") %>%
  derive_summary_records(
    by_vars = exprs(STUDYID, USUBJID, AVISIT, ADT, ADY, TRTSDT),
    analysis_var = AVAL,
    summary_fun = function(x) sum(x, na.rm = TRUE),
    # Select records contributing to total score
    filter = str_detect(PARAMCD, "GAD020[1-7]"),
    set_values_to = exprs(
      PARAMCD = "GAD02TS",
      PARAM = "GAD02-Total Score - Analysis"
    )
  )

## ----echo=FALSE---------------------------------------------------------------
dataset_vignette(
  arrange(adgad7, USUBJID, ADY, PARAMCD),
  display_vars = exprs(USUBJID, PARAMCD, PARAM, AVAL, ADY, AVISIT)
)

## ----eval=TRUE----------------------------------------------------------------
adgdssf <- adqs %>%
  # Select records to keep in the GDS-SF ADaM
  filter(PARCAT1 == "GDS SHORT FORM") %>%
  derive_summary_records(
    by_vars = exprs(STUDYID, USUBJID, AVISIT, ADT, ADY, TRTSDT),
    analysis_var = AVAL,
    summary_fun = function(x) {
      compute_scale(
        x,
        source_range = c(0, 1),
        target_range = c(0, 15),
        min_n = 10
      ) %>%
        ceiling()
    },
    # Select records contributing to total score
    filter = str_detect(PARAMCD, "GDS02[01][0-9]"),
    set_values_to = exprs(
      PARAMCD = "GDS02TS",
      PARAM = "GDS02- Total Score - Analysis"
    )
  )

## ----echo=FALSE---------------------------------------------------------------
dataset_vignette(
  arrange(adgdssf, USUBJID, ADY, PARAMCD),
  display_vars = exprs(USUBJID, PARAMCD, PARAM, AVAL, ADY, AVISIT)
)

## ----eval=TRUE----------------------------------------------------------------
adgdssf <- adgdssf %>%
  # Flag baseline records (last before treatement start)
  restrict_derivation(
    derivation = derive_var_extreme_flag,
    args = params(
      by_vars = exprs(STUDYID, USUBJID, PARAMCD),
      order = exprs(ADT),
      new_var = ABLFL,
      mode = "last"
    ),
    filter = !is.na(AVAL) & ADT <= TRTSDT
  ) %>%
  # Derive baseline and change from baseline variables
  derive_var_base(
    by_vars = exprs(STUDYID, USUBJID, PARAMCD),
    source_var = AVAL,
    new_var = BASE
  ) %>%
  derive_var_chg() %>%
  derive_var_pchg() %>%
  # Derive sequence number
  derive_var_obs_number(
    by_vars = exprs(STUDYID, USUBJID),
    order = exprs(PARAMCD, ADT),
    check_type = "error"
  )

## ----echo=FALSE---------------------------------------------------------------
dataset_vignette(
  adgdssf,
  display_vars = exprs(USUBJID, PARAMCD, PARAM, AVISIT, AVAL, BASE, CHG, PCHG)
)

## -----------------------------------------------------------------------------
adgdssf <- adgdssf %>%
  mutate(
    AVALCAT1 = if_else(
      PARAMCD == "GDS02TS",
      case_when(
        AVAL <= 5 ~ "Normal",
        AVAL <= 10 ~ "Possible Depression",
        AVAL > 10 ~ "Likely Depression"
      ),
      NA_character_
    ),
    AVALCAT1N = if_else(
      PARAMCD == "GDS02TS",
      case_when(
        AVAL <= 5 ~ 0L,
        AVAL <= 10 ~ 1L,
        AVAL > 10 ~ 2L
      ),
      NA_integer_
    )
  ) %>%
  derive_var_base(
    by_vars = exprs(STUDYID, USUBJID, PARAMCD),
    source_var = AVALCAT1,
    new_var = BASECAT1
  ) %>%
  derive_var_base(
    by_vars = exprs(STUDYID, USUBJID, PARAMCD),
    source_var = AVALCAT1N,
    new_var = BASECA1N
  ) %>%
  mutate(
    CHGCAT1 = case_when(
      AVALCAT1N > BASECA1N ~ "WORSENED",
      AVALCAT1N == BASECA1N ~ "NO CHANGE",
      AVALCAT1N < BASECA1N ~ "IMPROVED",
    )
  )

## ----echo=FALSE---------------------------------------------------------------
dataset_vignette(
  arrange(adgdssf, USUBJID, desc(PARAMCD), ADY),
  display_vars = exprs(USUBJID, PARAMCD, PARAM, AVISIT, AVAL, AVALCAT1, CHGCAT1)
)

## -----------------------------------------------------------------------------
# Define event
deterioration_event <- event_source(
  dataset_name = "adqs",
  filter = PARAMCD == "GDS02TS" & CHGCAT1 == "WORSENED",
  date = ADT,
  set_values_to = exprs(
    EVNTDESC = "DEPRESSION WORSENED",
    SRCDOM = "ADQS",
    SRCVAR = "ADT",
    SRCSEQ = ASEQ
  )
)

# Define censoring at last assessment
last_valid_assessment <- censor_source(
  dataset_name = "adqs",
  filter = PARAMCD == "GDS02TS" & !is.na(CHGCAT1),
  date = ADT,
  set_values_to = exprs(
    EVNTDESC = "LAST ASSESSMENT",
    SRCDOM = "ADQS",
    SRCVAR = "ADT",
    SRCSEQ = ASEQ
  )
)

# Define censoring at treatment start (for subjects without assessment)
start <- censor_source(
  dataset_name = "adsl",
  date = TRTSDT,
  set_values_to = exprs(
    EVNTDESC = "TREATMENT START",
    SRCDOM = "ADSL",
    SRCVAR = "TRTSDT"
  )
)

adgdstte <- derive_param_tte(
  dataset_adsl = adsl,
  source_datasets = list(adsl = adsl, adqs = adgdssf),
  start_date = TRTSDT,
  event_conditions = list(deterioration_event),
  censor_conditions = list(last_valid_assessment, start),
  set_values_to = exprs(
    PARAMCD = "TTDEPR",
    PARAM = "Time to depression"
  )
) %>%
  derive_vars_duration(
    new_var = AVAL,
    start_date = STARTDT,
    end_date = ADT
  )

## ----echo=FALSE---------------------------------------------------------------
dataset_vignette(
  adgdstte,
  display_vars = exprs(USUBJID, PARAMCD, PARAM, AVAL, CNSR, EVNTDESC, SRCDOM, SRCVAR)
)

## -----------------------------------------------------------------------------
adgdssf <- adgdssf %>%
  derive_var_joined_exist_flag(
    by_vars = exprs(USUBJID, PARAMCD),
    order = exprs(ADT),
    new_var = CDETFL,
    join_vars = exprs(CHGCAT1, ADY),
    join_type = "after",
    filter = CHGCAT1 == "WORSENED" &
      CHGCAT1.join == "WORSENED" &
      ADY.join >= ADY + 7
  )

## ----echo=FALSE---------------------------------------------------------------
dataset_vignette(
  arrange(adgdssf, USUBJID, desc(PARAMCD), ADY),
  display_vars = exprs(USUBJID, PARAMCD, PARAM, ADY, CHGCAT1, CDETFL)
)

## -----------------------------------------------------------------------------
# Flagging deterioration at two consecutive assessments
adgdssf <- adgdssf %>%
  derive_var_joined_exist_flag(
    by_vars = exprs(USUBJID, PARAMCD),
    order = exprs(ADT),
    new_var = CONDETFL,
    join_vars = exprs(CHGCAT1),
    join_type = "after",
    tmp_obs_nr_var = tmp_obs_nr,
    filter = CHGCAT1 == "WORSENED" &
      CHGCAT1.join == "WORSENED" &
      tmp_obs_nr.join == tmp_obs_nr + 1
  ) %>%
  # Flagging deterioration confirmed by
  # - a second deterioration at least 7 days later or
  # - deterioration at the last assessment and death due to progression
  derive_var_joined_exist_flag(
    by_vars = exprs(USUBJID, PARAMCD),
    order = exprs(ADT),
    new_var = CDTDTHFL,
    join_vars = exprs(CHGCAT1, ADY),
    join_type = "after",
    tmp_obs_nr_var = tmp_obs_nr,
    filter = CHGCAT1 == "WORSENED" & (
      CHGCAT1.join == "WORSENED" & ADY.join >= ADY + 7 |
        tmp_obs_nr == max(tmp_obs_nr.join) & DTHCAUS == "PROGRESSIVE DISEASE")
  )

## ----echo=FALSE---------------------------------------------------------------
dataset_vignette(
  arrange(adgdssf, USUBJID, desc(PARAMCD), ADY),
  display_vars = exprs(USUBJID, PARAMCD, PARAM, ADY, CHGCAT1, CONDETFL, CDTDTHFL)
)

## -----------------------------------------------------------------------------
adgdssf <- adgdssf %>%
  derive_var_joined_exist_flag(
    by_vars = exprs(USUBJID, PARAMCD),
    order = exprs(ADT),
    new_var = DEFDETFL,
    join_vars = exprs(CHGCAT1),
    join_type = "after",
    filter = CHGCAT1 == "WORSENED" & all(CHGCAT1.join == "WORSENED")
  )

## ----echo=FALSE---------------------------------------------------------------
dataset_vignette(
  arrange(adgdssf, USUBJID, desc(PARAMCD), ADY),
  display_vars = exprs(USUBJID, PARAMCD, PARAM, ADY, CHGCAT1, DEFDETFL)
)

## -----------------------------------------------------------------------------
adsp <- adqs %>%
  filter(PARCAT1 == "SLEEPING PROBLEMS") %>%
  derive_extreme_event(
    by_vars = exprs(USUBJID, AVISIT),
    order = exprs(ADY, QSSEQ),
    mode = "first",
    events = list(
      event(
        condition = PARAMCD == "SP0101" & AVALC == "YES",
        set_values_to = exprs(
          AVALC = "No sleep",
          AVAL = 1
        )
      ),
      event(
        condition = PARAMCD == "SP0102" & AVALC == "YES",
        set_values_to = exprs(
          AVALC = "Waking up more than three times",
          AVAL = 2
        )
      ),
      event(
        condition = PARAMCD == "SP0103" & AVALC == "YES",
        set_values_to = exprs(
          AVALC = "More than 30 mins to fall asleep",
          AVAL = 3
        )
      ),
      event(
        condition = all(AVALC == "NO"),
        set_values_to = exprs(
          AVALC = "No sleeping problems",
          AVAL = 4
        )
      ),
      event(
        condition = TRUE,
        set_values_to = exprs(
          AVALC = "Missing",
          AVAL = 99
        )
      )
    ),
    set_values_to = exprs(
      PARAMCD = "SP01WSP",
      PARAM = "Worst Sleeping Problems"
    )
  )

## ----echo=FALSE---------------------------------------------------------------
dataset_vignette(
  arrange(adsp, USUBJID, ADY, PARAMCD),
  display_vars = exprs(USUBJID, PARAMCD, PARAM, AVISIT, AVALC)
)

## -----------------------------------------------------------------------------
adgdssf <- adgdssf %>%
  derive_summary_records(
    filter = str_detect(PARAMCD, "GDS02[01][0-9]"),
    by_vars = exprs(USUBJID, AVISIT),
    analysis_var = AVAL,
    summary_fun = function(x) sum(!is.na(x)) / 15 >= 0.9,
    set_values_to = exprs(
      PARAMCD = "COMPL90P",
      PARAM = "Completed at least 90% of questions?",
    )
  ) %>%
  mutate(
    AVALC = if_else(
      PARAMCD == "COMPL90P",
      if_else(
        AVAL == 1,
        "YES",
        "NO"
      ),
      AVALC
    )
  )

## ----echo=FALSE---------------------------------------------------------------
dataset_vignette(
  arrange(adgdssf, USUBJID, PARAMCD, ADY),
  display_vars = exprs(USUBJID, PARAMCD, PARAM, AVISIT, AVALC)
)

## -----------------------------------------------------------------------------
# Create dataset with expected visits and parameters (GDS0201 - GDS0215)
parm_visit_ref <- crossing(
  tribble(
    ~AVISIT,    ~AVISITN,
    "BASELINE",        0,
    "VISIT 2",         2,
    "VISIT 3",         3,
    "VISIT 4",         4,
    "VISIT 5",         5
  ),
  tibble(PARAMCD = sprintf("GDS02%02d", seq(1, 15)))
)

adgdssf <- adgdssf %>%
  derive_expected_records(
    dataset_expected_obs = parm_visit_ref,
    by_vars = exprs(USUBJID),
    set_values_to = exprs(
      filled_in = 1
    )
  ) %>%
  derive_summary_records(
    filter = str_detect(PARAMCD, "GDS02[01][0-9]"),
    by_vars = exprs(USUBJID, AVISIT),
    analysis_var = AVAL,
    summary_fun = function(x) all(!is.na(x)),
    set_values_to = exprs(
      PARAMCD = "COMPLALL",
      PARAM = "Completed all questions?"
    )
  ) %>%
  mutate(
    AVALC = if_else(
      PARAMCD == "COMPLALL",
      if_else(
        AVAL == 1,
        "YES",
        "NO"
      ),
      AVALC
    )
  ) %>%
  filter(is.na(filled_in)) %>%
  select(-filled_in)

## ----echo=FALSE---------------------------------------------------------------
dataset_vignette(
  arrange(adgdssf, USUBJID, PARAMCD, ADY),
  display_vars = exprs(USUBJID, PARAMCD, PARAM, AVISIT, AVALC)
)

