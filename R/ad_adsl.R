# Name: ADSL
#
# Label: Subject Level Analysis Dataset
#
# Input:
#   - dm
#   - ex
#

# read in predecessor variables from dm
data("dm")

# derive treatment variables (TRT01P, TRT01A)
adsl <- dm %>%
  mutate(TRT01P = ARMCD, TRT01A = ACTARMCD)

# derive treatment start date (TRTSDTM, TRTSDT)
data("ex")
adsl <- adsl %>%
  derive_var_trtsdtm(dataset_ex = ex) %>%
  mutate(TRTSDT = lubridate::date(TRTSDTM))

# derive treatment end date (TRTEDTM, TRTEDT)
adsl <- adsl %>%
  derive_var_trtedtm(dataset_ex = ex) %>%
  dplyr::mutate(TRTEDT = lubridate::date(TRTEDTM))

# derive treatment duration (TRTDURD)
 adsl <- adsl %>%
   derive_var_trtdurd()
