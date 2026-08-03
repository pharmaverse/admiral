# =============================================================================
# Worked example: the ADSL-specific summary  (issue #3160 POC)
#
# POC for the first item of summary_admiral_df_design_notes.md: what
# `summary()` should report when the dataset is subject level (ADSL).
#
# It shows
#   1. the real ADSL  -- the facts a programmer checks first, and a check which
#                        finds something genuine in the CDISC pilot data
#   2. a clean ADSL   -- the single confirmation line when the checks pass
#   3. a broken ADSL  -- two defects injected, so both checks fire
#   4. a partial ADSL -- a mid-derivation dataset, to show that every item is
#                        skipped silently when its variables are absent
#
# HOW TO RUN
#   From the package root:   Rscript worked_example_adsl.R
# =============================================================================

suppressMessages({
  library(dplyr)
  library(cli)
})

devtools::load_all()

# =============================================================================
# 1. The real ADSL
# =============================================================================
cli::cli_h1("1. Real ADSL (admiral_adsl)")

adsl <- as_admiral_df(admiral_adsl, keys = "USUBJID")
summary(adsl)

cli::cli_bullets(c(
  "i" = "The treatment and population lines are the denominators of nearly
         every downstream table, so they are worth seeing before anything is
         derived from this dataset.",
  "!" = "The 12 flagged subjects are not an artefact: in the CDISC pilot data
         they were planned for {.val Xanomeline High Dose} and actually took
         {.val Xanomeline Low Dose}. All 12 are treated and in the safety
         population, so this is a genuine dose deviation which changes which
         arm they are summarized under.",
  "i" = "Whether that is a finding or a known feature of the study is the
         programmer's call -- which is why the check reports rather than
         signals."
))

# =============================================================================
# 2. A clean ADSL
# =============================================================================
cli::cli_h1("2. Clean ADSL (dose deviations removed)")

clean <- admiral_adsl %>%
  filter(is.na(ARM) | is.na(ACTARM) | ARM == ACTARM) %>%
  as_admiral_df(keys = "USUBJID")

summary(clean)

cli::cli_bullets(c(
  "i" = "With nothing to report, the subject-level checks collapse to a single
         line -- a clean ADSL should not cost a screen of zeroes."
))

# =============================================================================
# 3. The same ADSL with defects injected
# =============================================================================
cli::cli_h1("3. Broken ADSL (defects injected)")

# target subjects whose planned arm is Placebo, so the injected value really
# is a mismatch (picking rows blindly can hit a subject already on that arm)
mis_rows <- head(which(clean$ARM == "Placebo"), 3)
saf_rows <- head(which(clean$SAFFL == "Y" & !is.na(clean$TRTSDT)), 2)

broken <- clean %>%
  mutate(
    # three subjects were dosed with something other than their planned arm
    ACTARM = if_else(
      row_number() %in% mis_rows, "Xanomeline High Dose", ACTARM
    ),
    # two safety-population subjects lost their treatment start date, so no
    # treatment-emergent record can be derived for them downstream
    TRTSDT = if_else(row_number() %in% saf_rows, as.Date(NA), TRTSDT)
  ) %>%
  as_admiral_df(keys = "USUBJID")

summary(broken)

cli::cli_bullets(c(
  "i" = "Both checks are diagnostics, not validations: they are reported and
         never signalled, so {.fun summary} stays safe to call anywhere.",
  "i" = "Neither defect changes the record count, the subject count, or the
         record structure -- the existing summary would have shown nothing."
))

# =============================================================================
# 4. A partially built ADSL
# =============================================================================
cli::cli_h1("4. Partial ADSL (mid-derivation)")

partial <- admiral_adsl %>%
  select(STUDYID, USUBJID, AGE, SEX, ARM) %>%
  as_admiral_df(keys = "USUBJID")

summary(partial)

cli::cli_bullets(c(
  "i" = "No treatment flags, no {.var TRTDURD}, no {.var ACTARM}: each item is
         dropped rather than reported as missing, so the summary is still
         readable while a dataset is being built.",
  "!" = "Note the checks line is absent too -- neither check could run, which
         is different from both passing."
))

cli::cli_h1("Done")
