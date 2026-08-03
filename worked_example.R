# =============================================================================
# Worked example: admiral_df summary() diagnostics  (issue #3160 POC)
#
# Demonstrates everything implemented so far:
#   1. as_admiral_df()          -- tag a dataset (optionally with keys)
#   2. get_admiral_df_type()    -- classify ADaM structure
#   3. summary.admiral_df()     -- quick diagnostic (subjects, obs, params, ...)
#   4. Data-driven grain check  -- minimal_unique_key() / infer_admiral_keys()
#         - minimal key (no padding when not needed)
#         - automatic expansion until unique
#         - surrogate-key guardrail (ASEQ never masks real duplicates)
#         - graded, sourced messaging: declared vs (inferred)
#   5. by_vars wiring           -- keys stamped by derive_param_computed()
#   6. metacore spec keys       -- get_admiral_keys() / set_admiral_keys(),
#                                  incl. the "admiral_keys" attribute they attach
#
# HOW TO RUN
#   From the package root:   Rscript worked_example.R
#
# NOTE: This script sources the REAL functions from R/admiral_df.R and only
# stubs the two admiral-internal helpers it relies on (get_admiral_option,
# vars2chr) plus the admiraldev assertions, so it runs without a full
# devtools::load_all() (handy where admiraldev < 1.5.0 blocks load_all).
# =============================================================================

suppressMessages({
  library(dplyr)
  library(stringr)
  library(cli)
  library(rlang)
  library(tibble)
})

devtools::load_all()

rule <- function(title) cli::cli_h1(title)

# =============================================================================
# 1. BDS -- inferred grain, single visit  => minimal key (no AVISIT padding)
# =============================================================================
rule("1. BDS, inferred grain, single visit (minimal key)")

advs <- tribble(
  ~STUDYID,  ~USUBJID,      ~PARAMCD, ~AVAL, ~AVISIT,
  "PILOT01", "01-701-1015", "DIABP",     51, "BASELINE",
  "PILOT01", "01-701-1015", "SYSBP",    121, "BASELINE",
  "PILOT01", "01-701-1028", "DIABP",     79, "BASELINE",
  "PILOT01", "01-701-1028", "SYSBP",    130, "BASELINE"
)
print(summary(as_admiral_df(advs)))
# -> Grain (USUBJID, PARAMCD) (inferred): unique   (AVISIT not needed here)

# =============================================================================
# 2. BDS -- inferred grain, multiple visits => auto-expansion to reach unique
# =============================================================================
rule("2. BDS, inferred grain, multiple visits (auto-expansion)")

advs_2visits <- tribble(
  ~STUDYID,  ~USUBJID,      ~PARAMCD, ~AVAL, ~AVISIT,
  "PILOT01", "01-701-1015", "DIABP",     51, "BASELINE",
  "PILOT01", "01-701-1015", "DIABP",     50, "WEEK 2",
  "PILOT01", "01-701-1015", "SYSBP",    121, "BASELINE",
  "PILOT01", "01-701-1015", "SYSBP",    119, "WEEK 2"
)
print(summary(as_admiral_df(advs_2visits)))
# -> Grain (USUBJID, PARAMCD, AVISIT) (inferred): unique

# =============================================================================
# 3. by_vars wiring -- keys declared by derive_param_computed()
#    (here we stamp them directly, exactly as the derivation now does internally:
#     keys = c(vars2chr(by_vars), "PARAMCD"))
# =============================================================================
rule("3. Declared grain via by_vars (as derive_param_computed stamps it)")

by_vars <- exprs(USUBJID, AVISIT)
map_records <- tribble(
  ~STUDYID,  ~USUBJID,      ~PARAMCD, ~AVAL, ~AVISIT,
  "PILOT01", "01-701-1015", "MAP",       74, "BASELINE",
  "PILOT01", "01-701-1028", "MAP",       96, "BASELINE"
)
advs_map <- bind_rows(advs, map_records)
advs_map <- as_admiral_df(advs_map, keys = c(vars2chr(by_vars), "PARAMCD"))
print(summary(advs_map))
# -> Grain (USUBJID, AVISIT, PARAMCD): unique     (no "(inferred)" tag)

# =============================================================================
# 4. Duplicate detection -- a derivation that added a record twice
# =============================================================================
rule("4. Declared grain with an unintended duplicate (real finding)")

advs_dup <- bind_rows(advs_map, map_records[1, ]) # add MAP for 1015 a 2nd time
advs_dup <- as_admiral_df(advs_dup, keys = c(vars2chr(by_vars), "PARAMCD"))
print(summary(advs_dup))
# -> Grain (...): 1 duplicate key combination  (red X)

# =============================================================================
# 5. Surrogate-key guardrail -- ASEQ present but must NOT mask a real duplicate
# =============================================================================
rule("5. Surrogate-key guardrail (ASEQ does not hide duplicates)")

advs_aseq <- tribble(
  ~STUDYID,  ~USUBJID,      ~PARAMCD, ~AVAL, ~AVISIT,     ~ASEQ,
  "PILOT01", "01-701-1015", "MAP",       74, "BASELINE",  1L,
  "PILOT01", "01-701-1015", "MAP",       75, "BASELINE",  2L # dup grain, new ASEQ
)
print(summary(as_admiral_df(advs_aseq))) # inferred path, ASEQ excluded
# -> Grain (USUBJID, PARAMCD, AVISIT) (inferred): 1 duplicate key combination

# =============================================================================
# 6. Other structures -- type detection tailors the output
# =============================================================================
rule("6a. ADSL (one row per subject)")
adsl <- tribble(
  ~STUDYID,  ~USUBJID,      ~TRT01P,     ~TRTSDT,
  "PILOT01", "01-701-1015", "Placebo",   as.Date("2024-01-10"),
  "PILOT01", "01-701-1028", "Xanomeline", as.Date("2024-01-11")
)
print(summary(as_admiral_df(adsl)))

rule("6b. TTE (events vs censored)")
adtte <- tribble(
  ~STUDYID,  ~USUBJID,      ~PARAMCD, ~CNSR, ~STARTDT,
  "PILOT01", "01-701-1015", "OS",         0, as.Date("2024-01-10"),
  "PILOT01", "01-701-1028", "OS",         1, as.Date("2024-01-11")
)
print(summary(as_admiral_df(adtte)))

rule("6c. OCCDS with ASEQ (grain = USUBJID + ASEQ, the genuine record key)")
adae <- tribble(
  ~STUDYID,  ~USUBJID,      ~AEDECOD,    ~TRTEMFL, ~ASEQ,
  "PILOT01", "01-701-1015", "HEADACHE",  "Y",      1L,
  "PILOT01", "01-701-1015", "NAUSEA",    "Y",      2L,
  "PILOT01", "01-701-1028", "HEADACHE",  "Y",      1L
)
print(summary(as_admiral_df(adae)))

rule("6d. OCCDS without ASEQ (grain cannot be checked -> warning)")
adae_no_seq <- select(adae, -ASEQ)
print(summary(as_admiral_df(adae_no_seq)))

# =============================================================================
# 7. metacore spec keys -- get_admiral_keys() / set_admiral_keys()
#    (skipped automatically if {metacore} or its example spec is unavailable)
# =============================================================================
rule("7. Grain from a metacore ADaM specification")

if (requireNamespace("metacore", quietly = TRUE)) {
  spec_path <- tryCatch(metacore::metacore_example("pilot_ADaM.rda"),
    error = function(e) NA_character_
  )
  if (!is.na(spec_path)) {
    load(spec_path) # loads an object named `metacore`
    spec <- metacore

    cli::cli_text("Keys declared in the ADSL spec:")
    print(get_admiral_keys(spec, "ADSL"))

    # before tagging: a plain working dataset carries no keys and is not an
    # `admiral_df`
    cli::cli_text("Before {.fun set_admiral_keys}:")
    cli::cli_bullets(c(
      " " = "class:          {.val {class(adsl)}}",
      " " = "'admiral_keys': {.val {attr(adsl, 'admiral_keys') %||% 'not set'}}"
    ))

    # attach the spec-defined keys to the working dataset
    adsl_spec <- set_admiral_keys(adsl, spec, "ADSL")

    # after tagging: the keys are stored in the `admiral_keys` attribute and the
    # dataset is tagged as an `admiral_df`
    # (`extra attributes` = everything beyond the usual names/row.names/class,
    # i.e. what `set_admiral_keys()` actually added)
    extra_attrs <- setdiff(
      names(attributes(adsl_spec)),
      c("names", "row.names", "class")
    )
    cli::cli_text("After {.fun set_admiral_keys}:")
    cli::cli_bullets(c(
      " " = "class:          {.val {class(adsl_spec)}}",
      " " = "'admiral_keys': {.val {attr(adsl_spec, 'admiral_keys')}}",
      " " = "extra attrs:    {.val {extra_attrs}}"
    ))

    # the attribute is exactly what the spec declares -- nothing is inferred
    stopifnot(identical(
      attr(adsl_spec, "admiral_keys"),
      get_admiral_keys(spec, "ADSL")
    ))

    # ...and `summary()` uses it: the grain is reported WITHOUT the "(inferred)"
    # tag, unlike sections 1, 2 and 5 above
    print(summary(adsl_spec))

    # CAVEAT: the attribute survives the row/column-preserving {dplyr} verbs
    # (filter(), mutate(), bind_rows()) but not reshaping ones such as
    # summarise(), so re-tag (or tag last) if the dataset is manipulated further
    keys_after <- function(data) attr(data, "admiral_keys") %||% "dropped"
    cli::cli_text("Attribute after further manipulation:")
    cli::cli_bullets(c(
      " " = "filter():    {.val {keys_after(filter(adsl_spec, TRT01P == 'Placebo'))}}",
      " " = "mutate():    {.val {keys_after(mutate(adsl_spec, X = 1))}}",
      " " = "summarise(): {.val {keys_after(summarise(adsl_spec, N = n()))}}"
    ))
  } else {
    cli::cli_alert_info("metacore example spec not found; skipping section 7.")
  }
} else {
  cli::cli_alert_info("{{metacore}} not installed; skipping section 7.")
}

cli::cli_h1("Done")
