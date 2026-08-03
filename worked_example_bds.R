# =============================================================================
# Worked example: the BDS-specific summary  (issue #3160 POC)
#
# POC for the BDS section of summary_admiral_df_design_notes.md: what
# `summary()` should report when the dataset is a Basic Data Structure.
#
# It shows
#   1. the real ADLB  -- the per-parameter table and the derived-record
#                        breakdown, plus the record-level checks
#   2. a broken ADLB  -- four classic BDS defects injected, so every check
#                        fires; none of them changes the record count, so the
#                        old summary would have shown nothing
#   3. a partial ADLB -- a mid-derivation dataset, to show that items are
#                        skipped silently when their variables are absent
#   4. a derivation   -- `derive_param_computed()` output summarized directly,
#                        the original motivation of #3160
#
# HOW TO RUN
#   From the package root:   Rscript worked_example_bds.R
# =============================================================================

suppressMessages({
  library(dplyr)
  library(cli)
})

devtools::load_all()

# =============================================================================
# 1. The real ADLB
# =============================================================================
cli::cli_h1("1. Real ADLB (admiral_adlb)")

adlb <- as_admiral_df(admiral_adlb)
summary(adlb)

cli::cli_bullets(c(
  "i" = "The per-parameter table replaces the flat {.var PARAMCD} list: an
         empty parameter, a missing-value spike, or an absurd min/max is
         visible as a row, which a list of codes can never show.",
  "i" = "The {.var DTYPE} line is exactly what the derivations added
         (last-observation, minimum, maximum records), separated from the
         observed data.",
  "i" = "The visit count includes the unscheduled and {.val POST-BASELINE}
         visits, which is why it is larger than the number of scheduled
         weeks."
))

# =============================================================================
# 2. The same ADLB with defects injected
# =============================================================================
cli::cli_h1("2. Broken ADLB (defects injected)")

# each defect targets a row where it really is a defect: a non-baseline ALT
# record turned into a second baseline, a change value with its baseline
# removed, one Week 2 record renumbered, and one ALT record relabelled
dup_base_row <- head(
  which(
    admiral_adlb$PARAMCD == "ALT" & is.na(admiral_adlb$ABLFL) &
      is.na(admiral_adlb$DTYPE) & admiral_adlb$USUBJID == admiral_adlb$USUBJID[1]
  ),
  1
)
chg_row <- head(which(!is.na(admiral_adlb$CHG)), 1)
visit_row <- head(which(admiral_adlb$AVISIT == "Week 2"), 1)
param_row <- head(which(admiral_adlb$PARAMCD == "ALT"), 1)

broken <- admiral_adlb %>%
  mutate(
    ABLFL = if_else(row_number() == dup_base_row, "Y", ABLFL),
    BASE = if_else(row_number() == chg_row, NA_real_, BASE),
    AVISITN = if_else(row_number() == visit_row, 99, AVISITN),
    PARAM = if_else(row_number() == param_row, "ALT (U/L)", PARAM)
  ) %>%
  as_admiral_df()

summary(broken)

cli::cli_bullets(c(
  "!" = "None of the four defects changes the number of records, subjects, or
         parameters -- the facts section of the two summaries is identical.
         Only the checks tell them apart.",
  "i" = "The duplicated baseline is the classic silent one: {.var BASE},
         {.var CHG}, and {.var PCHG} downstream become wrong without any
         derivation erroring."
))

# =============================================================================
# 3. A partially built ADLB
# =============================================================================
cli::cli_h1("3. Partial ADLB (mid-derivation)")

partial <- admiral_adlb %>%
  select(STUDYID, USUBJID, PARAMCD, PARAM, AVISIT, AVISITN, AVAL) %>%
  as_admiral_df()

summary(partial)

cli::cli_bullets(c(
  "i" = "No {.var DTYPE}, no {.var ABLFL}, no {.var BASE}/{.var CHG}: the
         corresponding line and checks are dropped rather than reported as
         missing, so the summary stays useful while the dataset is being
         built.",
  "i" = "The {.var AVISIT}/{.var AVISITN} and {.var PARAM} consistency checks
         still run -- they only need the identifier variables, which exist
         from the first derivation step."
))

# =============================================================================
# 4. Summarizing a derivation directly
# =============================================================================
cli::cli_h1("4. Output of derive_param_computed()")

ratio <- admiral_adlb %>%
  filter(is.na(DTYPE)) %>%
  derive_param_computed(
    by_vars = exprs(USUBJID, AVISIT, AVISITN),
    parameters = c("AST", "ALT"),
    set_values_to = exprs(
      AVAL = AVAL.AST / AVAL.ALT,
      PARAMCD = "ASTALT",
      PARAM = "AST/ALT ratio"
    )
  )

summary(ratio)

cli::cli_bullets(c(
  "i" = "This is the original #3160 use case: {.fun derive_param_computed}
         tags its output, so {.fun summary} works immediately after the
         derivation -- the new {.val ASTALT} row in the parameter table shows
         what was added, with its value range as a sanity check.",
  "i" = "The ratio rows carry no {.var PARAM} label duplication and no
         baseline flags, so the checks which can run still pass."
))

cli::cli_h1("Done")
