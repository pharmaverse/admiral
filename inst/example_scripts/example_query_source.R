# VAR_PREFIX, e.g., SMQ01, CQ12
# QUERY_NAME, non NULL
# QUERY_ID, could be NULL
# QUERY_SCOPE, ‘BROAD’, ‘NARROW’, or NULL
# TERM_LEVEL, e.g., AEDECOD, AELLT, ...
# TERM_NAME, non NULL

queries <- tibble::tribble(
  ~VAR_PREFIX, ~QUERY_NAME, ~QUERY_ID, ~QUERY_SCOPE, ~TERM_LEVEL, ~TERM_NAME,
  "CQ01", "Immune-Mediated Hepatitis (Diagnosis and Lab Abnormalities)", "20000008", "NARROW", "AEDECOD", "ALANINE AMINOTRANSFERASE ABNORMAL",
  "SMQ02", "Immune-Mediated Hypothyroidism", "20000160", "BROAD", "AEDECOD", "BIOPSY THYROID GLAND ABNORMAL",
  "SMQ03", "Immune-Mediated Hypothyroidism", "20000161", "NARROW", "AEDECOD", "BASEDOW'S DISEASE",
  "CQ04", "Immune-Mediated Adrenal Insufficiency", "12150", "", "AEDECOD", "ADDISON'S DISEASE",
  "SMQ05", "Immune-Mediated Pneumonitis", "20000042", "NARROW", "AEDECOD", "ALVEOLAR PROTEINOSIS"
)

# save(queries, file = "data/queries.rda")

adae <- tibble::tribble(
  ~USUBJID, ~ASTDTM, ~AETERM, ~AESEQ, ~AEDECOD,
  "01", "2020-06-02 23:59:59", "ALANINE AMINOTRANSFERASE ABNORMAL", 3, "Alanine aminotransferase abnormal",
  "02", "2020-06-05 23:59:59", "BASEDOW'S DISEASE", 5, "Basedow's disease",
  "03", "2020-06-05 23:59:59", "ALVEOLAR PROTEINOSIS", 1, "Alveolar proteinosis"
)

# try below:
derive_query_vars(adae, queries)
