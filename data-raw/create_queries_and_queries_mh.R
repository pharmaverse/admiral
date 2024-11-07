# PREFIX, e.g., SMQ01, CQ12
# GRPNAME, non NULL
# GRPID, could be NULL
# SCOPE, 'BROAD', 'NARROW', or NULL
# SRCVAR, e.g., AEDECOD, AELLT, ...
# TERMCHAR, non NULL

library(dplyr) # nolint

queries <- tibble::tribble(
  ~PREFIX, ~GRPNAME, ~GRPID, ~SCOPE,
  ~SCOPEN, ~SRCVAR, ~TERMCHAR, ~TERMNUM,
  "CQ01", "Dermatologic events", NA_integer_, NA_character_,
  NA_integer_, "AELLT", "APPLICATION SITE ERYTHEMA", NA_integer_,
  "CQ01", "Dermatologic events", NA_integer_, NA_character_,
  NA_integer_, "AELLT", "APPLICATION SITE PRURITUS", NA_integer_,
  "CQ01", "Dermatologic events", NA_integer_, NA_character_,
  NA_integer_, "AELLT", "ERYTHEMA", NA_integer_,
  "CQ01", "Dermatologic events", NA_integer_, NA_character_,
  NA_integer_, "AELLT", "LOCALIZED ERYTHEMA", NA_integer_,
  "CQ01", "Dermatologic events", NA_integer_, NA_character_,
  NA_integer_, "AELLT", "GENERALIZED PRURITUS", NA_integer_,
  "SMQ02", "Immune-Mediated Hypothyroidism", 20000160L, "BROAD",
  1L, "AEDECOD", "BIOPSY THYROID GLAND ABNORMAL", NA_integer_,
  "SMQ02", "Immune-Mediated Hypothyroidism", 20000160L, "BROAD",
  1L, "AEDECOD", "BLOOD THYROID STIMULATING HORMONE ABNORMAL", NA_integer_,
  "SMQ02", "Immune-Mediated Hypothyroidism", 20000160L, "NARROW",
  2L, "AEDECOD", "BIOPSY THYROID GLAND INCREASED", NA_integer_,
  "SMQ03", "Immune-Mediated Guillain-Barre Syndrome", 20000131L, "NARROW",
  2L, "AEDECOD", "GUILLAIN-BARRE SYNDROME", NA_integer_,
  "SMQ03", "Immune-Mediated Guillain-Barre Syndrome", 20000131L, "NARROW",
  2L, "AEDECOD", "MILLER FISHER SYNDROME", NA_integer_,
  "CQ04", "Immune-Mediated Adrenal Insufficiency",
  12150L, NA_character_, NA_integer_, "AEDECOD", "ADDISON'S DISEASE", NA_integer_,
  "CQ04", "Immune-Mediated Adrenal Insufficiency",
  12150L, NA_character_, NA_integer_, "AEDECOD", "ADRENAL ATROPHY", NA_integer_,
  "SMQ05", "Immune-Mediated Pneumonitis", 20000042L, "NARROW",
  2L, "AEDECOD", "ALVEOLAR PROTEINOSIS", NA_integer_,
  "SMQ05", "Immune-Mediated Pneumonitis", 20000042L, "NARROW",
  2L, "AEDECOD", "ALVEOLITIS", NA_integer_,
  "CQ06", "Immune-Mediated Colitis", 10009888L, NA_character_,
  NA_integer_, "AELLTCD", NA_character_, 1L
)

admiral:::assert_valid_queries(queries)

# Get previous dataset for comparison
queries_old <- admiral::queries

#  create queries.rda in data/
usethis::use_data(queries, overwrite = TRUE)

# Compare with previous version
diffdf::diffdf(
  base = queries_old,
  compare = queries,
  keys = c("PREFIX", "TERMCHAR", "TERMNUM")
)

# example to use for ADMH:
queries_mh <- queries %>%
  filter(SRCVAR %in% c("AELLT", "AEDECOD")) %>%
  mutate(SRCVAR = ifelse(SRCVAR == "AELLT", "MHLLT", "MHDECOD"))

admiral:::assert_valid_queries(queries)

# Get previous dataset for comparison
queries_mh_old <- admiral::queries_mh

usethis::use_data(queries_mh, overwrite = TRUE)

# Compare with previous version
diffdf::diffdf(
  base = queries_mh_old,
  compare = queries_mh,
  keys = c("PREFIX", "TERMCHAR", "TERMNUM")
)
