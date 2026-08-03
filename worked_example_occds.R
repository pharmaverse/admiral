# =============================================================================
# Worked example: the OCCDS-specific summary  (issue #3160 POC)
#
# POC for the OCCDS section of summary_admiral_df_design_notes.md: what
# `summary()` should report when the dataset is an Occurrence Data Structure.
#
# It shows
#   1. a real ADAE    -- derived from pharmaversesdtm::ae the same way the
#                        admiral ADAE template does, with occurrence flags;
#                        the facts and the checks that pass
#   2. a broken ADAE  -- three classic OCCDS defects injected, so every check
#                        fires; none of them changes the record count, so the
#                        old summary would have shown nothing
#   3. a partial ADAE -- a mid-derivation dataset, to show that items are
#                        skipped silently when their variables are absent
#
# HOW TO RUN
#   From the package root:   Rscript worked_example_occds.R
# =============================================================================

suppressMessages({
  library(dplyr)
  library(cli)
})

devtools::load_all()

# =============================================================================
# 1. A real ADAE
# =============================================================================
cli::cli_h1("1. Real ADAE (derived from pharmaversesdtm::ae)")

# the standard template steps: merge the treatment dates from ADSL, impute the
# analysis dates, flag the treatment-emergent records, then flag the first
# occurrence per subject and per subject/body system among the emergent records
adae <- pharmaversesdtm::ae %>%
  convert_blanks_to_na() %>%
  derive_vars_merged(
    dataset_add = admiral_adsl,
    new_vars = exprs(TRTSDT, TRTEDT, TRT01A),
    by_vars = exprs(STUDYID, USUBJID)
  ) %>%
  derive_vars_dt(
    new_vars_prefix = "AST",
    dtc = AESTDTC,
    highest_imputation = "M"
  ) %>%
  derive_vars_dt(new_vars_prefix = "AEN", dtc = AEENDTC) %>%
  derive_var_trtemfl(
    start_date = ASTDT,
    end_date = AENDT,
    trt_start_date = TRTSDT,
    trt_end_date = TRTEDT
  ) %>%
  mutate(ASEQ = AESEQ) %>%
  restrict_derivation(
    derivation = derive_var_extreme_flag,
    args = params(
      by_vars = exprs(USUBJID),
      order = exprs(ASTDT, AESEQ),
      new_var = AOCCFL,
      mode = "first"
    ),
    filter = TRTEMFL == "Y"
  ) %>%
  restrict_derivation(
    derivation = derive_var_extreme_flag,
    args = params(
      by_vars = exprs(USUBJID, AEBODSYS),
      order = exprs(ASTDT, AESEQ),
      new_var = AOCCSFL,
      mode = "first"
    ),
    filter = TRTEMFL == "Y"
  ) %>%
  as_admiral_df(keys = c("USUBJID", "ASEQ"))

summary(adae)

cli::cli_bullets(c(
  "i" = "The term counts show the coding granularity at each level (body
         system, dictionary term, reported term) instead of the single
         {.var --DECOD} count reported before.",
  "i" = "The treatment-emergent line is the denominator shift: most AE outputs
         are restricted to {.var TRTEMFL == \"Y\"}, so this is the population
         the tables will actually see.",
  "i" = "The occurrence flags are unique per subject and level here because
         {.fun derive_var_extreme_flag} was restricted to the emergent records
         -- zero flags for a subject is legitimate, only a duplicate is a
         defect."
))

# =============================================================================
# 2. The same ADAE with defects injected
# =============================================================================
cli::cli_h1("2. Broken ADAE (defects injected)")

# each defect targets a row where it really is a defect: a second first-occurrence
# flag for a subject who already has one, a removed analysis start date, and an
# emergent record moved before treatment start
emergent <- which(adae$TRTEMFL == "Y" & !is.na(adae$ASTDT) & !is.na(adae$TRTSDT))
dup_row <- head(which(adae$TRTEMFL == "Y" & is.na(adae$AOCCFL)), 1)
na_row <- head(setdiff(emergent, dup_row), 1)
pre_row <- head(setdiff(emergent, c(dup_row, na_row)), 1)

broken <- adae %>%
  mutate(
    AOCCFL = if_else(row_number() == dup_row, "Y", AOCCFL),
    ASTDT = case_when(
      row_number() == na_row ~ as.Date(NA),
      row_number() == pre_row ~ TRTSDT - 10,
      TRUE ~ ASTDT
    )
  ) %>%
  as_admiral_df(keys = c("USUBJID", "ASEQ"))

summary(broken)

cli::cli_bullets(c(
  "!" = "None of the three defects changes the number of records, subjects, or
         terms -- the facts section of the two summaries is identical. Only
         the checks tell them apart.",
  "i" = "The duplicated {.var AOCCFL} is the classic silent one: the subject
         is counted twice in every {.val n (%)} of-subjects table downstream,
         and the record structure check cannot see it because the record count
         is unchanged."
))

# =============================================================================
# 3. A partially built ADAE
# =============================================================================
cli::cli_h1("3. Partial ADAE (mid-derivation)")

partial <- adae %>%
  select(STUDYID, USUBJID, ASEQ, AEBODSYS, AEDECOD, AETERM, AESEV, AESER) %>%
  as_admiral_df(keys = c("USUBJID", "ASEQ"))

summary(partial)

cli::cli_bullets(c(
  "i" = "No {.var TRTEMFL}, no dates, no occurrence flags: the corresponding
         lines and checks are dropped rather than reported as missing, so the
         summary stays useful while the dataset is being built.",
  "i" = "The term counts and the severity/seriousness breakdown only need the
         coded SDTM variables, which exist from the first derivation step."
))

cli::cli_h1("Done")
