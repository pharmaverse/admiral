## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

library(admiraldev)

## ---- warning=FALSE, message=FALSE--------------------------------------------
library(admiral)
library(lubridate)
library(tibble)
library(dplyr, warn.conflicts = FALSE)

## -----------------------------------------------------------------------------
impute_dtc_dtm(
  "2019-10",
  highest_imputation = "M",
  date_imputation = "01-01",
  time_imputation = "00:00:00"
)

## -----------------------------------------------------------------------------
impute_dtc_dtm(
  "2019-02",
  highest_imputation = "M",
  date_imputation = "02-31",
  time_imputation = "00:00:00"
)

## -----------------------------------------------------------------------------
impute_dtc_dtm(
  "2019-02",
  highest_imputation = "M",
  date_imputation = "last",
  time_imputation = "00:00:00"
)

## -----------------------------------------------------------------------------
dates <- c(
  "2019-02",
  "2019",
  "2019---01"
)
impute_dtc_dtm(
  dates,
  highest_imputation = "M",
  date_imputation = "mid",
  time_imputation = "00:00:00",
  preserve = FALSE
)
impute_dtc_dtm(
  dates,
  highest_imputation = "M",
  date_imputation = "mid",
  time_imputation = "00:00:00",
  preserve = TRUE
)

## -----------------------------------------------------------------------------
dates <- c(
  "2019-02",
  "2019",
  "2019---01"
)
impute_dtc_dtm(
  dates,
  highest_imputation = "M",
  date_imputation = "06-15",
  time_imputation = "00:00:00"
)

## -----------------------------------------------------------------------------
dates <- c(
  "2019-02-03T12:30:15",
  "2019-02-03T12:30",
  "2019-02-03",
  "2019-02",
  "2019"
)

# Do not impute
impute_dtc_dtm(
  dates,
  highest_imputation = "n"
)

# Impute seconds only
impute_dtc_dtm(
  dates,
  highest_imputation = "s"
)

# Impute time (hours, minutes, seconds) only
impute_dtc_dtm(
  dates,
  highest_imputation = "h"
)

# Impute days and time
impute_dtc_dtm(
  dates,
  highest_imputation = "D"
)

# Impute date (months and days) and time
impute_dtc_dtm(
  dates,
  highest_imputation = "M"
)

## -----------------------------------------------------------------------------
impute_dtc_dtm(
  "2019-02",
  highest_imputation = "M",
  date_imputation = "last",
  time_imputation = "last",
  max_dates = list(ymd("2019-01-14"), ymd("2019-02-25"))
)

## -----------------------------------------------------------------------------
# Impute year to first
impute_dtc_dtm(
  c("2019-02", NA),
  highest_imputation = "Y",
  min_dates = list(
    ymd("2019-01-14", NA),
    ymd("2019-02-25", "2020-01-01")
  )
)

# Impute year to last
impute_dtc_dtm(
  c("2019-02", NA),
  highest_imputation = "Y",
  date_imputation = "last",
  time_imputation = "last",
  max_dates = list(
    ymd("2019-01-14", NA),
    ymd("2019-02-25", "2020-01-01")
  )
)

## -----------------------------------------------------------------------------
ae <- tribble(
  ~AESTDTC,
  "2019-08-09T12:34:56",
  "2019-04-12",
  "2010-09",
  NA_character_
) %>%
  derive_vars_dtm(
    dtc = AESTDTC,
    new_vars_prefix = "AST",
    highest_imputation = "M",
    date_imputation = "first",
    time_imputation = "first"
  ) %>%
  derive_vars_dtm_to_dt(exprs(ASTDTM))

## ---- echo=FALSE--------------------------------------------------------------
dataset_vignette(ae)

## -----------------------------------------------------------------------------
ae <- tribble(
  ~AESTDTC,
  "2019-08-09T12:34:56",
  "2019-04-12",
  "2010-09",
  NA_character_
) %>%
  derive_vars_dt(
    dtc = AESTDTC,
    new_vars_prefix = "AST",
    highest_imputation = "M",
    date_imputation = "first"
  )

## ---- echo=FALSE--------------------------------------------------------------
dataset_vignette(ae)

## -----------------------------------------------------------------------------
ae <- tribble(
  ~AESTDTC,
  "2019-08-09T12:34:56",
  "2019-04-12",
  "2010-09",
  NA_character_
) %>%
  derive_vars_dtm(
    dtc = AESTDTC,
    new_vars_prefix = "AST",
    highest_imputation = "h",
    time_imputation = "first"
  )

## ---- echo=FALSE--------------------------------------------------------------
dataset_vignette(ae)

## -----------------------------------------------------------------------------
ae <- tribble(
  ~AESTDTC,              ~TRTSDTM,
  "2019-08-09T12:34:56", ymd_hms("2019-11-11T12:34:56"),
  "2019-10",             ymd_hms("2019-11-11T12:34:56"),
  "2019-11",             ymd_hms("2019-11-11T12:34:56"),
  "2019-12-04",          ymd_hms("2019-11-11T12:34:56")
) %>%
  derive_vars_dtm(
    dtc = AESTDTC,
    new_vars_prefix = "AST",
    highest_imputation = "M",
    date_imputation = "first",
    time_imputation = "first",
    min_dates = exprs(TRTSDTM)
  )

## ---- echo=FALSE--------------------------------------------------------------
dataset_vignette(ae)

## -----------------------------------------------------------------------------
ae <- tribble(
  ~AEENDTC,              ~DTHDT,            ~DCUTDT,
  "2019-08-09T12:34:56", ymd("2019-11-11"), ymd("2019-12-02"),
  "2019-11",             ymd("2019-11-11"), ymd("2019-12-02"),
  "2019-12",             NA,                ymd("2019-12-02"),
  "2019-12-04",          NA,                ymd("2019-12-02")
) %>%
  derive_vars_dtm(
    dtc = AEENDTC,
    new_vars_prefix = "AEN",
    highest_imputation = "M",
    date_imputation = "last",
    time_imputation = "last",
    max_dates = exprs(DTHDT, DCUTDT)
  )

## ---- echo=FALSE--------------------------------------------------------------
dataset_vignette(ae)

## -----------------------------------------------------------------------------
mh <- tribble(
  ~MHSTDTC,     ~TRTSDT,
  "2019-04",    ymd("2019-04-15"),
  "2019-04-01", ymd("2019-04-15"),
  "2019-05",    ymd("2019-04-15"),
  "2019-06-21", ymd("2019-04-15")
) %>%
  filter(
    convert_dtc_to_dt(
      MHSTDTC,
      highest_imputation = "M",
      date_imputation = "first"
    ) < TRTSDT
  )

## ---- echo=FALSE--------------------------------------------------------------
dataset_vignette(mh)

## -----------------------------------------------------------------------------
vs <- tribble(
  ~VSDTC,                ~VSTPT,
  "2019-08-09T12:34:56", NA,
  "2019-10-12",          "PRE-DOSE",
  "2019-11-10",          NA,
  "2019-12-04",          NA
) %>%
  slice_derivation(
    derivation = derive_vars_dtm,
    args = params(
      dtc = VSDTC,
      new_vars_prefix = "A"
    ),
    derivation_slice(
      filter = VSTPT == "PRE-DOSE",
      args = params(time_imputation = "first")
    ),
    derivation_slice(
      filter = TRUE,
      args = params(time_imputation = "last")
    )
  )

## ---- echo=FALSE--------------------------------------------------------------
dataset_vignette(vs)

