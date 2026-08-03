# `summary(x, adsl =)` — comparing a dataset against ADSL, design notes

Brainstorm for the ADSL-linkage feature (2026-08-03, companion to issue #3160,
`summary_admiral_df_design_notes.md` design question 3, and
`keys_design_notes.md`). The idea: for non-ADSL types (`BDS`, `OCCDS`, `TTE`,
`other`), `summary()` accepts an **optional** `adsl =` argument and adds a
cross-dataset section. ADSL is the subject-level source of truth — population
membership, treatment assignment and dates, death, disposition — so a child
dataset can be wrong *relative to ADSL* in ways it can never reveal on its own.

Status: the **suggested first increment below is implemented** (2026-08-03,
`summarize_vs_adsl()` / `print_vs_adsl_summary()` and the
`summary(x, adsl =)` argument; see `admiral_df_showcase.qmd`). The remaining
catalogue items (TTE completeness, fatal-AE vs `DTHFL`, `STUDYID` agreement,
data-cut/consent windows, `population =`, the options channel) are still open.

## Why ADSL specifically

Three things only ADSL can provide:

1. **The full subject universe.** A child dataset only contains subjects with
   records; it cannot know who is *missing* from it, and it cannot supply the
   denominator for any "n (%) of subjects" statement.
2. **The authoritative subject-level values.** `TRTSDT`, `TRT01A`, `SAFFL`,
   `DTHDT` in a child dataset are *copies* (via `derive_vars_merged()`), and
   copies go stale when ADSL is refreshed after the child was built.
3. **Subject-level events the child must respect.** Death, end of study, data
   cut: records dated after them are defects the child dataset cannot detect
   alone.

## What can go wrong — the catalogue

Organized by the same facts-vs-red-flags principle as the per-type notes:
🚩 = worth surfacing only when non-zero, usually a genuine defect.

### A. Referential integrity (any type)

| Item | Notes |
|---|---|
| 🚩 **Orphan subjects** — in the dataset but not in ADSL | The single highest-value check. Causes: wrong ADSL version, screen failures left in, test/dummy subjects, `USUBJID` formatting drift (`"01-701-1015"` vs `"01701-1015"` — a formatting mismatch *presents* as orphans, so example values matter in the output). Always a defect. |
| **Coverage** — subjects in ADSL not in the dataset | Usually legitimate (a subject with no AEs is absent from ADAE), so this is a *fact*, not a flag: "225 of 254 ADSL subjects appear here (88.6%)". Restricting the denominator to the relevant population (`SAFFL == "Y"` for safety datasets) makes it more meaningful when the flag exists. |
| 🚩 TTE completeness | TTE is the exception where absence *is* a defect: every subject in the analysis population should contribute exactly one record per endpoint (`PARAMCD`). Report missing subjects per endpoint. |
| 🚩 `STUDYID` disagreement | Same `USUBJID`, different `STUDYID` between the two datasets — pooling/merge defect. Cheap to check when both have `STUDYID`. |

### B. Stale or wrong merged variables (any type)

The variables a child dataset shares with ADSL should be *identical* per
subject, because they were copied from it. Compare, per subject, every shared
subject-level variable from a curated list — e.g. `TRTSDT`, `TRTEDT`,
`TRT01P`, `TRT01A`, `SAFFL`, `ITTFL`, `DTHDT`, `EOSDT`, `AGE`, `SEX`, `RACE` —
and report:

| Item | Notes |
|---|---|
| 🚩 **Value disagreement per shared variable** | "TRTSDT differs from ADSL for 12 subjects" — the classic stale-merge signal: ADSL was re-derived after the child dataset was built, or the merge used the wrong keys (demographics disagreeing is near-proof of a bad join). This turns the invisible "built from last month's ADSL" into one line. |
| 🚩 Treatment mapping mismatch | Period-aware pairs: `TRTA`/`TRTP` in the child vs `TRT01A`/`TRT01P` in ADSL. A mismatch means wrong period handling or wrong merge. Only safe to check naively in single-period studies — multi-period needs `APERIOD`-aware mapping, so start with equality-when-single-period and skip otherwise. |
| Direction of authority | Report mismatches as "differs from ADSL", not "child is wrong" — the *cause* may be an outdated child **or** a newer ADSL; the wording should not presume. |

### C. Cross-dataset event logic (any type, strongest for OCCDS/TTE)

These need an ADSL variable the child typically does not carry:

| Item | Notes |
|---|---|
| 🚩 **Records after death** | `ASTDT`/`ADT` after ADSL `DTHDT`. Always a finding. |
| 🚩 Emergent without treatment | `TRTEMFL == "Y"` for a subject whose ADSL `TRTSDT` is missing (never treated). The existing OCCDS check uses the *merged* `TRTSDT`; with `adsl =` it works even when the child lacks the date columns, and against the authoritative values. |
| 🚩 Death inconsistency (OCCDS) | AE with fatal outcome (`AESDTH == "Y"` / `AEOUT == "FATAL"`) but ADSL `DTHFL != "Y"` — the two datasets disagree about whether the subject died. Classic cross-dataset QC finding. |
| 🚩 Records after data cut / end of study | `ADT`/`ASTDT` after `DCUTDT` (or `EOSDT`) when ADSL carries it. Advisory: some studies legitimately keep post-cut records, so this may belong behind a fact line rather than a red flag. |
| 🚩 Records before informed consent / randomization | `ADT < RFICDT` (or `RANDDT` for post-randomization datasets). Lower priority; medical-history-style datasets legitimately predate both. |
| 🚩 TTE event vs ADSL death | An OS event date differing from `DTHDT`, or `EVNTDESC` indicating death while `DTHFL != "Y"`. Endpoint-specific, so likely a later increment. |

### D. Denominators — the facts users actually want

Not defect checks but the reason programmers open ADSL alongside every review:

| Item | Notes |
|---|---|
| **Subjects line gains a percentage** | "Subjects (USUBJID): 225 (88.6% of 254 in ADSL)" — the existing line, upgraded in place. |
| **Per-arm coverage** | Subjects with records per `TRT01A`/`TRT01P` vs the arm N from ADSL: "Placebo 60/86, Xano HD 70/84, ...". An empty or thin arm is instantly visible, and it is the shape of every downstream incidence table. |
| **Per-parameter coverage (BDS)** | Add a `subjects/N` flavour to the per-parameter table — a parameter measured on half the population stands out only when the denominator is there. Possibly a printed column only when `adsl` is supplied. |
| **OCCDS incidence** | "Subjects with ≥1 treatment-emergent record: 220/254 (86.6%)" — the number every AE table leads with, and the child dataset alone cannot compute it. Answers the open item in the OCCDS notes ("subjects-with-occurrence as a percentage of the population"). |

## What the user should see

A `vs ADSL` block appended to the type-specific section, following the existing
conventions exactly:

- **Facts first**: the coverage line, per-arm coverage; percentage upgrades
  happen in place on existing lines rather than duplicating them.
- **Red flags via `print_summary_checks()`**: each failing check is one bullet
  with its count (orphans should include a few example `USUBJID`s — formatting
  drift is only diagnosable from examples); passing checks collapse into the
  one green line; a check whose variables are absent is not mentioned.
- **Details in the object, not on screen**: the offending `USUBJID`s and the
  per-variable mismatch table go into the returned `summary_admiral_df` object
  (e.g. `$vs_adsl$orphans`, `$vs_adsl$mismatches`) for drill-down, truncated on
  print — same pattern as `$bds$params`.
- **Absent argument = absent section**: no `adsl` supplied → no `vs ADSL`
  content at all, not "skipped" noise. Not run is different from passed.

## Design decisions to make

1. **Signature.** `summary(object, keys = NULL, adsl = NULL, ...)` with
   `check_dots_empty()` — the pattern established by `keys =`. Named formal,
   not `...`, so a typo (`adls =`) errors instead of being swallowed.
2. **Argument validation vs check philosophy.** The checks never signal, but
   *argument* validation may: if `adsl` is not a data frame, or is not unique
   by the subject key, the comparisons are ambiguous and an informative
   `cli_abort()` is fine (caller error, not data finding). If `adsl` lacks
   `USUBJID` entirely, abort likewise. Everything beyond that degrades
   gracefully: each item runs only when its variables exist on *both* sides.
3. **Join key.** `USUBJID` alone, consistent with the existing subject-count
   preference (other subject keys may be `NA` on derivation-added records).
   `STUDYID` agreement becomes its own check rather than part of the join key.
4. **`summary()` on an ADSL with `adsl =` supplied.** Two defensible options:
   warn-and-ignore, or actually compare (self-vs-frozen-copy is a real QC use:
   current ADSL vs previous data cut). Suggest warn-and-ignore for the first
   increment; the compare-two-ADSLs idea is a separate feature (`{diffdf}`
   territory).
5. **Population-aware denominators.** All-of-ADSL is the simple denominator;
   the *right* one is often a population (`SAFFL` for ADAE). First increment:
   report against all of ADSL, plus the `SAFFL` denominator when both sides
   allow it. A `population =` argument is a later refinement.
6. **Third channel.** Mirroring the keys design: a session-wide
   `set_admiral_options(adsl = )` (the `subject_keys` precedent) would give
   the zero-ceremony version — set once at the top of a script, every
   `summary()` call gains the section, per-call argument overrides. Also the
   natural home if `derive_vars_merged()` ever records provenance ("this
   dataset was merged from that ADSL"). Worth deciding *after* the argument
   proves its value.
7. **Cost.** Everything is one `USUBJID`-level join plus vectorized
   comparisons — single-pass, nothing quadratic, fine on a 2M-row ADLB.
8. **Scope boundary.** Same as ever: diagnostics, not validation, never a
   `{diffdf}` replacement, never errors on data.

## Suggested first increment

Highest signal, no new machinery beyond one join:

1. 🚩 Orphan subjects (with example `USUBJID`s in the object).
2. Coverage line + per-arm coverage (facts).
3. 🚩 Stale shared variables, starting with the treatment/date/flag list.
4. 🚩 Records after death (`DTHDT`).
5. OCCDS: incidence percentage; 🚩 emergent-without-treatment via ADSL.

TTE per-endpoint completeness and the death-consistency checks are the natural
second increment (they want the TTE per-`PARAMCD` rework first).
