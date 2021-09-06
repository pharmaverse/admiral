# VAR_PREFIX, e.g., SMQ01, CQ12
# QUERY_NAME, non NULL
# QUERY_ID, could be NULL
# QUERY_SCOPE, ‘BROAD’, ‘NARROW’, or NULL
# TERM_LEVEL, e.g., AEDECOD, AELLT, ...
# TERM_NAME, non NULL

queries <- tibble::tribble(
  ~VAR_PREFIX, ~QUERY_NAME, ~QUERY_ID, ~QUERY_SCOPE,
  ~QUERY_SCOPE_NUM, ~TERM_LEVEL, ~TERM_NAME, ~TERM_ID,
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
  1L, "AEDECOD", "BIOPSY THYROID GLAND INCREASED", NA_integer_,
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
derive_vars_query(adae, queries)
