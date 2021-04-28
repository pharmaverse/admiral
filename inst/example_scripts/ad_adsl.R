# Name: ADSL
#
# Label: Subject Level Analysis Dataset
#
# Input: dm, ex
#

library(dplyr)
library(lubridate)
library(admiral)

# read in predecessor variables from dm
data("dm")
data("ex")

# derive treatment variables (TRT01P, TRT01A)
adsl <- dm %>%
  mutate(TRT01P = ARMCD, TRT01A = ACTARMCD) %>%

  # derive treatment start date (TRTSDTM, TRTSDT)
  derive_var_trtsdtm(dataset_ex = ex) %>%
  mutate(TRTSDT = date(TRTSDTM)) %>%

  # derive treatment end date (TRTEDTM, TRTEDT)
  derive_var_trtedtm(dataset_ex = ex) %>%
  mutate(TRTEDT = date(TRTEDTM)) %>%

  # derive treatment duration (TRTDURD)
  derive_var_trtdurd()

save(adsl, file = "data/adsl.rda", compress = TRUE)
