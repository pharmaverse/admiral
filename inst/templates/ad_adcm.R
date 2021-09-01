# Name: ADCM
#
# Label: Concomitant Medications Analysis Dataset
#
# Input: cm, adsl, suppcm, suppdm, ex
#

#devtools :: load_all()
#library(admiral)
#library(dplyr)
#library(lubridate)

# ---- Load source datasets ----

# Use e.g. haven::read_sas to read in .sas7bdat, or other suitable functions
#  as needed and assign to the variables below.

data("cm")
data("adsl")
data("ex")
adsl <- adsl
suppcm <- NULL
ex <- NULL

# ---- Derivations ----

adcm <- cm %>%
  # join supplementary qualifier variables
  #derive_suppqual_vars(suppcm) %>%

  # join adsl to cm
  left_join(adsl, by = c("STUDYID", "USUBJID")) %>%

  # derive analysis start time
  derive_vars_dtm(
    dtc = CMSTDTC,
    new_vars_prefix = "AST",
    date_imputation = "first",
    time_imputation = "first",
    min_dates = list(TRTSDT)
  ) %>%

  # derive analysis end time
  derive_vars_dtm(
    dtc = CMENDTC,
    new_vars_prefix = "AEN",
    date_imputation = "last",
    time_imputation = "last",
    max_dates = list(DTHDT, EOSDT)
  ) %>%

  # derive analysis end/start date
  mutate(
    ASTDT = date(ASTDTM),
    AENDT = date(AENDTM)
  ) %>%

  # derive analysis start relative day
  derive_var_astdy(
    reference_date = TRTSDT,
    date = ASTDT
  ) %>%

  # derive analysis end relative day
  derive_var_aendy(
    reference_date = TRTSDT,
    date = AENDT
  ) %>%

  # derive analysis duration (value and unit)
  derive_duration(
    new_var = ADURN,
    new_var_unit = ADURU,
    start_date = ASTDT,
    end_date = AENDT,
    in_unit = "days",
    out_unit = "days",
    add_one = TRUE,
    trunc_out = FALSE
  ) %>%
  # Derive Time Relative to Reference
  derive_var_atirel(flag_var = ASTTMF,
                    new_var = ATIREL
  ) %>%
  # derive On-Treatment flag
  derive_var_ontrtfl(
   date = ASTDT,
     ref_start_date = TRTSDT,
     ref_end_date = TRTEDT
   ) %>%
  # derive Pre-Treatment flag
  mutate(
    PREFL = ifelse(ASTDT < TRTSDT, "Y", "")
  ) %>%
  # derive Follow-Up flag
  mutate(
    FUPFL = ifelse(ASTDT > TRTSDT+28, "Y", "")
  ) %>%
  # derive Aphase and Aphasen Variable
  # Other timing variable can be derived similarly.
  mutate(
    APHASE = case_when(ONTRTFL == "Y"~"On-Treatment",
                       PREFL == "Y" ~ "Pre-Treatment",
                       FUPFL == "Y" ~ "Follow-Up"),
    APHASEN = case_when(ONTRTFL == "Y"~2,
                       PREFL == "Y" ~ 1,
                       FUPFL == "Y" ~ 3)
  ) %>%
  # Assign TRTP/TRTA
  mutate(
    TRTP = TRT01P,
    TRTA = TRT01A
  ) %>%
  # derive ANL01FL : This is variable is sponsor specific
  mutate(
    ANL01FL = ifelse(ONTRTFL == "Y","Y", "")
  ) %>%
  # derive occurrence flags
  derive_extreme_flag(
    new_var = AOCCIFL,
    by_vars = vars(USUBJID),
    order = vars(ASTDTM, CMSEQ),
    flag_filter = ANL01FL == "Y",
    mode = "last"
  )



# Note some of the roche specific variables are not included like ADCUTFL, ATC Variable, PERASTDY, APERIOD.I can add this variable based on
# the discussion.


# ---- Save output ----

saveRDS(adcm, file = "./ADCM.rds", compress = TRUE)


