# `summary.admiral_df()` — per-type content, design notes

Brainstorm for expanding what `summary()` reports about an admiral dataset,
depending on its ADaM type (`ADSL`, `BDS`, `OCCDS`, `TTE`, `other`) as detected
by `get_admiral_df_type()`. Companion to issue #3160.

Status: **notes for discussion, nothing implemented.** The current
implementation reports type, observations, subjects, the record structure
check, `PARAMCD`/`PARAM` values, `AVISIT` values, events/censored (TTE), and
the per-type ADSL, BDS, and OCCDS sections (see the dated implementation notes
below).

## Organizing principle: facts vs red flags

A clinical programmer running `summary()` after a derivation is asking two
questions:

1. **Did I create what I expected?** — counts, parameters, visits, date ranges.
2. **Did I break anything?** — duplicated baselines, orphaned subjects,
   negative survival times.

The second is where the value is, and it is what a generic `skimr`-style
summary cannot give you. Throughout the list below, 🚩 marks a **red flag**:
something that is worth surfacing only when it is non-zero, and which usually
indicates a genuine defect rather than a property of the data.

Red flags should stay advisory. If `summary()` ever errors, people stop calling
it, and it becomes useless as a browsing tool.

## Cross-cutting (any type)

| Item | Notes |
|---|---|
| Population / analysis flag counts | `SAFFL`, `ITTFL`, `EFFFL`, `FASFL`, `PPROTFL`, `RANDFL`. Almost every downstream TFL is subset by one of these. |
| Treatment distribution | N per `TRT01P` / `TRTA` / `TRTP`. Immediately shows an empty or misspelled arm. |
| Date coverage | min/max of the main date variable (`ADT`, `ASTDT`, `TRTSDT`), plus counts of imputation flags (`ADTF`, `ATMF`, `ASTDTF`). Imputed-date counts are a standing QC question. |
| Missingness | Only on the analysis-critical variables for the type (`AVAL`, `TRTSDT`, `PARAMCD`, `ASTDT`), not every column. |
| 🚩 Referential integrity with ADSL | Subjects here but not in ADSL, and optionally treated subjects in ADSL missing here. One of the highest-value checks in the whole list. Needs an `adsl` argument — see design questions. |
| 🚩 Records added vs input | If the object knows it came from a `derive_*` call: "N in, N out, N added". The single most useful line after a derivation. |

## ADSL

- N randomized / treated / completed, by arm.
- 🚩 `ARM` vs `ACTARM` mismatches (treatment misallocation). A real finding
  every time it is non-zero.
- 🚩 Missing `TRTSDT` among `SAFFL == "Y"`.
- Disposition: `EOSSTT` / `DCSREAS` counts.
- Treatment duration (`TRTDURD`): min / median / max.
- Demographics one-liner: age range, sex and race breakdown.
- Data cut date (`DCUTDT`) if present.

## BDS

- **Per-parameter table** instead of the current flat `PARAMCD` list: records,
  subjects, visits, missing `AVAL`, min / median / max of `AVAL` per parameter.
  This is the big one — most BDS defects are visible in a per-parameter row.
- 🚩 **Baseline flag integrity**: subject × parameter combinations with zero or
  more than one `ABLFL == "Y"`. Common, silent, and damaging downstream.
- 🚩 `CHG` / `PCHG` populated where `BASE` is missing.
- Derived record counts by `DTYPE` (LOCF, AVERAGE, MAXIMUM, ...) — shows
  exactly what the derivation added.
- Analysis flag counts (`ANL01FL`), and 🚩 whether the flagged subset is unique
  per subject / parameter / visit.
- Visit structure: number of visits, records per visit, unscheduled visits.
- 🚩 `AVISIT` ↔ `AVISITN` inconsistency (the same `AVISIT` with two different
  `AVISITN` values).
- 🚩 Multiple units (`AVALU`) or multiple `PARAM` labels for a single
  `PARAMCD`.
- Criterion flag counts (`CRIT1FL`) if present.

## OCCDS

- Subjects with at least one occurrence, and the percentage of the population.
- Distinct terms at each level: `AEBODSYS` (SOC) and `AEDECOD` (PT). Currently
  only one `--DECOD` level is reported.
- `TRTEMFL == "Y"` counts, and pre-treatment record counts.
- Severity / grade distribution (`AESEV`, `ATOXGR`), serious (`AESER`), related
  (`AEREL`), outcome.
- 🚩 **Occurrence flag integrity**: `AOCCFL`, `AOCC02FL` etc. should be exactly
  one `"Y"` per subject per level. The classic OCCDS defect, and the
  record-structure check does not catch it.
- 🚩 Missing `ASTDT`, or `ASTDT` before `TRTSDT` where `TRTEMFL == "Y"`.

**Implemented 2026-08-01** (`summarize_occds()` / `print_occds_summary()`,
mirroring the ADSL/BDS pattern): term counts per coding level
(`--BODSYS`/`--DECOD`/`--TERM`), treatment-emergent vs total records,
severity/grade distribution (`ASEV`/`AESEV`, falling back to
`ATOXGR`/`AETOXGR`), serious count (`AESER`), plus the three checks --
occurrence flag integrity (level derived from the flag's infix: `S` = per
subject and body system, `P` = per subject and dictionary term, otherwise per
subject; only *more than one* `"Y"` per group is a defect, since flags are
typically restricted to the emergent records), missing `ASTDT`, and
treatment-emergent records starting before `TRTSDT`. See
`worked_example_occds.R`. Still open: subjects-with-occurrence as a percentage
of the population (needs the `adsl =` argument), pre-treatment counts as their
own line, relatedness (`AEREL`, heterogeneous values) and outcome breakdowns.

## TTE

- Events / censored **per `PARAMCD`**. Each parameter is a separate endpoint,
  so the current single total is misleading with more than one endpoint.
- Censoring reason distribution (`CNSDTDSC`, `EVNTDESC`).
- Follow-up time: min / median / max `AVAL` with units, and median follow-up.
- 🚩 `AVAL <= 0`, or `ADT < STARTDT` (negative survival time).
- 🚩 Subjects appearing more than once per endpoint.

## Design questions

1. **Verbosity.** The full list would run well past one screen. Suggested: the
   default stays a short digest (counts plus red flags only), with the detail —
   per-parameter tables and so on — always computed into the returned object
   but truncated on print, or gated behind `summary(x, detail = TRUE)`.
2. **Cost.** These are cheap per-column operations, but on a 2M-row ADLB the
   per-parameter statistics add up. Everything should stay single-pass `dplyr`,
   with nothing quadratic.
3. **ADSL linkage.** Several of the best checks need ADSL. Options: an `adsl`
   argument, or a `get_admiral_option()` entry. Worth deciding early, because it
   changes the signature.
4. **Scope boundary.** Red flags are diagnostics, not validation. This should
   not grow into a replacement for `{diffdf}` or a study's QC scripts, and it
   must never error.
5. **Graceful degradation.** Every item must be skipped silently when the
   variables it needs are absent — datasets in progress are the normal case.

## Suggested next increment

High signal, no signature change, reuses the counting that already exists:

1. **BDS** — per-parameter table.
2. **TTE** — events / censored per parameter.
3. **OCCDS** — occurrence flag integrity.
4. **ADSL** — `ARM` vs `ACTARM` mismatches, and missing `TRTSDT` among
   `SAFFL == "Y"`.

Referential integrity with ADSL is the strongest remaining check, but it should
be separate work because of the argument change.

Open question before implementing: agree the shape of the returned object per
type first, since `print.summary_admiral_df()` and the tests both depend on it.

## Brainstorm 2026-08-01: more ADSL content

The first ADSL increment is implemented (treatment, populations, `TRTDURD`,
`EOSSTT`, demographics, `ARM`/`ACTARM` and `SAFFL`-without-`TRTSDT` checks).
Candidates for the next one, all present in `admiral_adsl`:

| Item | Notes |
|---|---|
| Deaths | `DTHFL == "Y"` count, `DTHCAUS` breakdown when populated. |
| Subject funnel | randomized (`RANDDT`/`RANDFL`) / treated (`TRTSDT`) / completed (`EOSSTT`) on one line; makes screen failures explicit instead of `NA: 52` in the `EOSSTT` line. |
| Disposition reason | `DCSREAS` next to `EOSSTT`. |
| Race/ethnicity | fold `RACE` into the demographics line. |
| Multi-period | detect `TRT02P` etc. and report the number of periods. |
| 🚩 `TRTEDT < TRTSDT` | negative treatment duration. |
| 🚩 `TRTDURD` recomputation | disagrees with `TRTEDT - TRTSDT + 1`. |
| 🚩 `TRTSDT < RANDDT` | treated before randomized. |
| 🚩 Death inconsistency | `DTHDT` present but `DTHFL != "Y"`, or death before `TRTSDT`. |
| 🚩 Flag domain | population flag values outside `Y`/`N`/`NA` (`"y"`, `"Yes"`, `""`). |

Each red flag is one entry in the `checks` list of `print_adsl_summary()`.
Suggested subset: deaths, funnel, `TRTEDT < TRTSDT`, flag domain — **implemented
2026-08-01**; the remaining items (`DCSREAS`, race, multi-period, the other
date checks) are still open.

## Brainstorm 2026-08-01: options via the method signature (the `...` question)

Do **not** consume options from raw `...`: a typo
(`summary(adsl, kyes = ...)`) is swallowed silently. Instead add named formals
to the method (the `summary.lm` pattern) and error on leftovers with
`rlang::check_dots_empty()`.

Candidate arguments, roughly in order of value:

1. `keys =` -- character vector *or* a `Metacore` object resolved via
   `get_admiral_keys()`. Precedence: `keys` arg > `admiral_keys` attribute >
   inference; report as `key_source = "supplied"`.
2. `metacore =` -- the whole spec: keys, plus mandatory variables from
   `ds_vars` not yet present, plus (BDS, later) expected vs actual `PARAMCD`
   from `value_spec`. If this exists, `keys =` may not need to accept Metacore
   objects.
3. `adsl =` -- enables the referential-integrity check for non-ADSL types
   (design question 3, answered without touching the generic).
4. `input =` -- the pre-derivation dataset, enabling "N in, N out, N added".
5. `detail =` -- digest vs full; but per design question 1 this likely belongs
   on `print()`, with everything always computed into the object.

Third channel: `set_admiral_options(metacore = spec)` for a session-wide spec
(the `subject_keys` precedent), with the per-call argument as override --
arguably the most admiral-idiomatic route.

Suggested next increment: `keys =` (character or Metacore) plus
`check_dots_empty()`, establishing the pattern the later arguments slot into.
Implemented; see `keys_design_notes.md` for the assessment of the argument
vs. attribute pathways (conclusion: keep both, precedence supplied >
declared > inferred).
