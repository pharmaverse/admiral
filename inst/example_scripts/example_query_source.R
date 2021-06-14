# VAR_PREFIX, e.g., SMQ01, CQ12
# QUERY_NAME, non NULL
# QUERY_ID, could be NULL
# QUERY_SCOPE, ‘BROAD’, ‘NARROW’, or NULL
# TERM_LEVEL, e.g., AEDECOD, AELLT, ...
# TERM_NAME, non NULL

queries <- tibble::tribble(
  ~VAR_PREFIX, ~QUERY_NAME, ~QUERY_ID, ~QUERY_SCOPE, ~TERM_LEVEL, ~TERM_NAME,
  "CQ01", "Immune-Mediated Hepatitis (Diagnosis and Lab Abnormalities)", "20000008", "NARROW", "AEDECOD", "ALANINE AMINOTRANSFERASE ABNORMAL",
  "CQ01", "Immune-Mediated Hepatitis (Diagnosis and Lab Abnormalities)", "20000008", "NARROW", "AEDECOD", "ALANINE AMINOTRANSFERASE INCREASED",
  "CQ01", "Immune-Mediated Hepatitis (Diagnosis and Lab Abnormalities)", "20000008", "NARROW", "AEDECOD", "AMMONIA ABNORMAL",
  "SMQ02", "Immune-Mediated Hypothyroidism", "20000160", "BROAD", "AEDECOD", "BIOPSY THYROID GLAND ABNORMAL",
  "SMQ02", "Immune-Mediated Hypothyroidism", "20000160", "BROAD", "AEDECOD", "BLOOD THYROID STIMULATING HORMONE ABNORMAL",
  "SMQ03", "Immune-Mediated Hypothyroidism", "20000161", "NARROW", "AEDECOD", "BASEDOW'S DISEASE",
  "SMQ03", "Immune-Mediated Hypothyroidism", "20000161", "NARROW", "AEDECOD", "EXOPHTHALMOS",
  "CQ04", "Immune-Mediated Adrenal Insufficiency", "12150", NA_character_, "AEDECOD", "ADDISON'S DISEASE",
  "CQ04", "Immune-Mediated Adrenal Insufficiency", "12150", NA_character_, "AEDECOD", "ADRENAL ATROPHY",
  "SMQ05", "Immune-Mediated Pneumonitis", "20000042", "NARROW", "AEDECOD", "ALVEOLAR PROTEINOSIS",
  "SMQ05", "Immune-Mediated Pneumonitis", "20000042", "NARROW", "AEDECOD", "ALVEOLITIS",
  "CQ06", "Immune-Mediated Colitis", "10009888", NA_character_, "AELLT", "COLITIS"
)

# save(queries, file = "data/queries.rda")

adae <- tibble::tribble(
  ~USUBJID, ~ASTDTM, ~AETERM, ~AESEQ, ~AEDECOD, ~AELLT,
  "01", "2020-06-02 23:59:59", "ALANINE AMINOTRANSFERASE ABNORMAL", 3, "Alanine aminotransferase abnormal", NA_character_,
  "02", "2020-06-05 23:59:59", "BASEDOW'S DISEASE", 5, "Basedow's disease", NA_character_,
  "02", "2020-06-05 23:59:59", "ALVEOLAR PROTEINOSIS", 1, "Alveolar proteinosis", NA_character_,
  "03", "2020-06-07 23:59:59", "SOME TERM", 2, "Some query", "Some term"
)

# try below:
derive_query_vars(adae, queries)
