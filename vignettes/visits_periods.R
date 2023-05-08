## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ---- warning=FALSE, message=FALSE--------------------------------------------
library(admiral)
library(tibble)
library(dplyr, warn.conflicts = FALSE)
library(lubridate)

## ---- echo=FALSE, warning=FALSE, message=FALSE--------------------------------
library(admiraldev)

## -----------------------------------------------------------------------------
windows <- tribble(
  ~AVISIT,    ~AWLO, ~AWHI, ~AVISITN, ~AWTARGET,
  "BASELINE",   -30,     1,        0,         1,
  "WEEK 1",       2,     7,        1,         5,
  "WEEK 2",       8,    15,        2,        11,
  "WEEK 3",      16,    22,        3,        19,
  "WEEK 4",      23,    30,        4,        26
)

## -----------------------------------------------------------------------------
adbds <- tribble(
  ~USUBJID, ~ADY,
  "1",       -33,
  "1",        -2,
  "1",         3,
  "1",        24,
  "2",        NA,
)

derive_vars_joined(
  adbds,
  dataset_add = windows,
  filter_join = AWLO <= ADY & ADY <= AWHI
)

## ----echo=FALSE---------------------------------------------------------------
phase_ref <- tribble(
  ~USUBJID, ~APHASEN, ~PHSDT,       ~PHEDT,       ~APHASE,
  "1",             1, "2021-01-04", "2021-02-06", "TREATMENT",
  "1",             2, "2021-02-07", "2021-03-07", "FUP",
  "2",             1, "2021-02-02", "2021-03-02", "TREATMENT"
) %>%
  mutate(
    STUDYID = "xyz",
    APHASEN = as.integer(APHASEN),
    across(matches("PH[ES]DT"), ymd)
  ) %>%
  select(STUDYID, everything())

phase_ref

## -----------------------------------------------------------------------------
adsl <- tibble(STUDYID = "xyz", USUBJID = c("1", "2"))

adsl <- derive_vars_period(
  adsl,
  dataset_ref = phase_ref,
  new_vars = exprs(PHwSDT = PHSDT, PHwEDT = PHEDT, APHASEw = APHASE)
) %>%
  select(STUDYID, USUBJID, PH1SDT, PH1EDT, PH2SDT, PH2EDT, APHASE1, APHASE2)

adsl

## -----------------------------------------------------------------------------
adae <- tribble(
  ~USUBJID, ~ASTDT,
  "1",      "2021-01-01",
  "1",      "2021-01-05",
  "1",      "2021-02-05",
  "1",      "2021-03-05",
  "1",      "2021-04-05",
  "2",      "2021-02-15",
) %>%
  mutate(ASTDT = ymd(ASTDT))

derive_vars_joined(
  adae,
  dataset_add = phase_ref,
  by_vars = exprs(USUBJID),
  filter_join = PHSDT <= ASTDT & ASTDT <= PHEDT
)

## -----------------------------------------------------------------------------
create_period_dataset(
  adsl,
  new_vars = exprs(PHSDT = PHwSDT, PHEDT = PHwEDT, APHASE = APHASEw)
)

## -----------------------------------------------------------------------------
# Add period variables to ADSL
period_ref <- tribble(
  ~USUBJID, ~APERIOD, ~APERSDT,     ~APEREDT,     ~TRTA,
  "1",             1, "2021-01-04", "2021-02-06", "DRUG A",
  "1",             2, "2021-02-07", "2021-03-07", "DRUG B",
  "2",             1, "2021-02-02", "2021-03-02", "DRUG B",
  "2",             2, "2021-03-03", "2021-04-01", "DRUG B"
) %>%
  mutate(
    STUDYID = "xyz",
    APERIOD = as.integer(APERIOD),
    across(ends_with("DT"), ymd)
  )

adsl <- derive_vars_period(
  adsl,
  dataset_ref = period_ref,
  new_vars = exprs(
    APxxSDT = APERSDT,
    APxxEDT = APEREDT,
    TRTxxA = TRTA
  )
) %>%
  select(
    STUDYID, USUBJID,
    TRT01A, TRT02A,
    AP01SDT, AP01EDT, AP02SDT, AP02EDT
  )

adsl

## -----------------------------------------------------------------------------
adae <- tribble(
  ~USUBJID, ~ASTDT,
  "1",      "2021-01-05",
  "1",      "2021-02-05",
  "1",      "2021-03-05",
  "1",      "2021-04-05",
  "2",      "2021-02-15",
  "2",      "2021-03-10",
) %>%
  mutate(
    ASTDT = ymd(ASTDT),
    STUDYID = "xyz"
  )

derive_vars_joined(
  adae,
  dataset_add = period_ref,
  by_vars = exprs(STUDYID, USUBJID),
  new_vars = exprs(APERIOD, TRTA),
  join_vars = exprs(APERSDT, APEREDT),
  filter_join = APERSDT <= ASTDT & ASTDT <= APEREDT
)

## -----------------------------------------------------------------------------
create_period_dataset(
  adsl,
  new_vars = exprs(APERSDT = APxxSDT, APEREDT = APxxEDT, TRTA = TRTxxA)
)

