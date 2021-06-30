# VAR_PREFIX, e.g., SMQ01, CQ12
# QUERY_NAME, non NULL
# QUERY_ID, could be NULL
# QUERY_SCOPE, ‘BROAD’, ‘NARROW’, or NULL
# TERM_LEVEL, e.g., AEDECOD, AELLT, ...
# TERM_NAME, non NULL

queries <- tibble::tribble(
  ~VAR_PREFIX, ~QUERY_NAME, ~QUERY_ID, ~QUERY_SCOPE, ~QUERY_SCOPE_NUM, ~TERM_LEVEL, ~TERM_NAME,
  "CQ01", "Dermatologic events", NA_integer_, NA_character_,
  NA_integer_, "AELLT", "APPLICATION SITE ERYTHEMA",
  "CQ01", "Dermatologic events", NA_integer_, NA_character_,
  NA_integer_, "AELLT", "APPLICATION SITE PRURITUS",
  "CQ01", "Dermatologic events", NA_integer_, NA_character_,
  NA_integer_, "AELLT", "ERYTHEMA",
  "CQ01", "Dermatologic events", NA_integer_, NA_character_,
  NA_integer_, "AELLT", "LOCALIZED ERYTHEMA",
  "CQ01", "Dermatologic events", NA_integer_, NA_character_,
  NA_integer_, "AELLT", "GENERALIZED PRURITUS",
  "SMQ02", "Immune-Mediated Hypothyroidism", 20000160, "BROAD",
  1, "AEDECOD", "BIOPSY THYROID GLAND ABNORMAL",
  "SMQ02", "Immune-Mediated Hypothyroidism", 20000160, "BROAD",
  1, "AEDECOD", "BLOOD THYROID STIMULATING HORMONE ABNORMAL",
  "SMQ02", "Immune-Mediated Hypothyroidism", 20000160, "NARROW",
  1, "AEDECOD", "BIOPSY THYROID GLAND INCREASED",
  "SMQ03", "Immune-Mediated Guillain-Barre Syndrome", 20000131, "NARROW",
  2, "AEDECOD", "GUILLAIN-BARRE SYNDROME",
  "SMQ03", "Immune-Mediated Guillain-Barre Syndrome", 20000131, "NARROW",
  2, "AEDECOD", "MILLER FISHER SYNDROME",
  "CQ04", "Immune-Mediated Adrenal Insufficiency",
  12150, NA_character_, NA_integer_, "AEDECOD", "ADDISON'S DISEASE",
  "CQ04", "Immune-Mediated Adrenal Insufficiency",
  12150, NA_character_, NA_integer_, "AEDECOD", "ADRENAL ATROPHY",
  "SMQ05", "Immune-Mediated Pneumonitis", 20000042, "NARROW",
  2, "AEDECOD", "ALVEOLAR PROTEINOSIS",
  "SMQ05", "Immune-Mediated Pneumonitis", 20000042, "NARROW",
  2, "AEDECOD", "ALVEOLITIS",
  "CQ06", "Immune-Mediated Colitis", 10009888, NA_character_,
  NA_integer_, "AELLT", "COLITIS"
) %>% dplyr::mutate(
  TERM_ID = as.integer(as.factor(.data$TERM_NAME))
)

adae <- tibble::tribble(
  ~USUBJID, ~ASTDTM, ~AETERM, ~AESEQ, ~AEDECOD, ~AELLT,
  "01", "2020-06-02 23:59:59", "ERYTHEMA", 3,
  "Erythema", "Localized erythema",
  "02", "2020-06-05 23:59:59", "BASEDOW'S DISEASE", 5,
  "Basedow's disease", NA_character_,
  "02", "2020-06-05 23:59:59", "ALVEOLAR PROTEINOSIS", 1,
  "Alveolar proteinosis", NA_character_,
  "03", "2020-06-07 23:59:59", "SOME TERM", 2,
  "Some query", "Some term",
  "04", "2020-06-10 23:59:59", "APPLICATION SITE ERYTHEMA", 7,
  "APPLICATION SITE ERYTHEMA", "Application site erythema",
)

# try below:
derive_query_vars(adae, queries)
