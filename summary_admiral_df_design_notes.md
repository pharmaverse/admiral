# `summary.admiral_df()` — per-type content, design notes

Brainstorm for expanding what `summary()` reports about an admiral dataset,
depending on its ADaM type (`ADSL`, `BDS`, `OCCDS`, `TTE`, `other`) as detected
by `get_admiral_df_type()`. Companion to issue #3160.

Status: **notes for discussion.** The current implementation reports type,
observations, subjects, the record structure check, `PARAMCD`/`PARAM` values,
`AVISIT` values, events/censored (TTE), and the per-type ADSL, BDS, and OCCDS
sections (see the dated implementation notes below).

**Read the 2026-08-03 decision below before acting on anything in this
document.** The red-flag catalogues in the per-type sections predate it and are
largely out of scope now; they are kept as the hand-over list for the ADaM
checks package rather than as a roadmap for this one.

## Decision 2026-08-03: core team review — checks move out of {admiral}

Presented the prototype to the core team. The feature itself was well received;
the objection was to **where the data-quality checks live**. Edoardo's position:
integrating them into {admiral} opens a can of worms — {admiral} is a
derivation package, and once it also adjudicates data quality it inherits an
open-ended obligation to every study's QC expectations. Better to outsource
them to a separate *ADaM checks* package.

Agreed direction: an **iterative, additive** summary. It adds a piece of
information only when the variables it needs are present, built up from the
variables actually found in the data, rather than attempting every check at
once. The action item was to remove the excessive checks over the following
week.

### The line that was drawn

Fifteen checks existed; **three** were kept. The kept ones share a property
that makes them {admiral}'s business rather than a checks package's:

> A check stays only if it is a **uniqueness or referential-integrity
> statement about keys and key-like flags** — the same kind of statement as the
> record structure check that motivated #3160, just applied to a subset of the
> records or across two datasets.

| Kept | Statement |
|---|---|
| `multiple_baselines` (BDS) | `ABLFL == "Y"` is unique within subject and parameter |
| `occ_flag_dups` (OCCDS) | an occurrence flag is unique within its level |
| `n_orphans` (vs ADSL) | every subject key points at an ADSL subject |

Everything else was a statement about *values* — domains, cross-variable
consistency, missingness, date ordering — and was removed:

| Removed | Section |
|---|---|
| `arm_mismatch`, `saffl_no_trtsdt`, `trt_end_before_start`, `flag_domain` | ADSL (now has no check section at all) |
| `chg_no_base`, `avisit_mismatch`, `param_inconsistent` | BDS |
| `missing_astdt`, `pre_trt_emergent` | OCCDS |
| `stale_total`/`stale_vars`, `n_after_death`, `n_emergent_untreated` | vs ADSL |

None of the *facts* were touched: treatment and population breakdowns, the
subject funnel, deaths, demographics, the BDS per-parameter table, `DTYPE`,
OCCDS terms/emergence/severity, and all the ADSL coverage denominators
(including OCCDS incidence) stay exactly as they were. The removal is visible
in the tests as `expect_null()` assertions on datasets that still contain the
injected defects — a defect present but deliberately unreported.

`print_summary_checks()` stays. It still partitions failed / passed / not-run,
and the three remaining checks still use it; it is also the slot any future
check drops into, so the machinery cost of the trim is zero.

### What this means for the rest of this document

The red-flag catalogues below are now a **hand-over list**, not a backlog. When
the ADaM checks package exists they are its specification, and the removed
implementations are recoverable from the git history of this branch. The
suggested-next-increment lists should be read as facts-only.

## Organizing principle: facts vs red flags

*(Superseded in part by the 2026-08-03 decision: red flags are still the right
way to think about the domain, but as of that decision {admiral} only carries
the key-integrity subset of them.)*

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
typically restricted to the emergent records; and, from 2026-08-07, per
`APERIOD` as well when the dataset has it -- see the {pharmaverseadam} sweep
below). See `worked_example_occds.R`.
The other two checks added at the same time — missing `ASTDT` and
treatment-emergent records starting before `TRTSDT` — were **removed
2026-08-03**; the flag integrity check is the only one left. Still open:
pre-treatment counts as their own line, relatedness (`AEREL`, heterogeneous
values) and outcome breakdowns. Subjects-with-occurrence as a percentage of the
population was answered by the `adsl =` argument (`incidence`).

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

*(Written before 2026-08-03. Items 1–3 are implemented; item 4 was implemented
and then removed by that decision. Kept for the record.)*

1. **BDS** — per-parameter table.
2. **TTE** — events / censored per parameter.
3. **OCCDS** — occurrence flag integrity.
4. ~~**ADSL** — `ARM` vs `ACTARM` mismatches, and missing `TRTSDT` among
   `SAFFL == "Y"`.~~ Removed 2026-08-03; belongs to the ADaM checks package.

Referential integrity with ADSL is the strongest remaining check, but it should
be separate work because of the argument change. *(Implemented 2026-08-03 as
`n_orphans`; it survived the trim because it is a key-integrity statement.)*

Open question before implementing: agree the shape of the returned object per
type first, since `print.summary_admiral_df()` and the tests both depend on it.

## Validation 2026-08-07: the {pharmaverseadam} sweep

`worked_example_pharmaverseadam.R` runs `summary()` over all 30 datasets in
{pharmaverseadam}, each non-ADSL one twice (with and without `adsl =`), for 58
calls in total. It is the first time the feature has been run on data nobody
wrote for it, and it is worth re-running after any change to type detection or
key inference.

Findings about **key inference** — two bugs and a design flaw, all now fixed —
are written up in `keys_design_notes.md`. Headline: reported duplicate records
fell from 15,729 to 1,110 and structure-check warnings from 6 to 0, with the
1,110 residual verified as genuine defects in the test data rather than
inference failures.

Findings about the **summary content** itself:

- Every type printed without error, and degradation held: no dataset produced a
  line it lacked the variables for.
- Zero orphan subjects across all 29 ADSL comparisons, once the ADSL is chosen
  by `USUBJID` overlap rather than by name (the therapeutic-area datasets do
  not share a subject universe, and the vaccine ones pair to `adsl_vaccine`).
- The coverage line earns its place: it makes visible that the
  therapeutic-area datasets cover thin slices of the shared ADSL
  (`advs_peds` 5 of 306, `adcoeq_metabolic` 5 of 306) — exactly what you want
  flagged before building a table from one of them.
- `adce_vaccine` was the only dataset where `occ_flag_dups` fired, and it was a
  **false positive** — see below. Fixed by adding `APERIOD` to the grouping;
  the check now passes on all four OCCDS datasets.
- Variable count added to the shared header (`Observations: n | Variables: n`)
  after running this — with 30 real datasets in front of you, the absence of
  the second dimension was conspicuous.

### The occurrence flag check assumed a single period

`adce_vaccine` reported a duplicate `AOCC01FL`: subject `ABC-1002` has two
`"Y"` records. They differ on `APERIOD` (1 and 2). It is a two-vaccination
study, and the flag marks the first reactogenicity event per subject *per
vaccination period* — one `"Y"` per subject and period in every case, which is
correct.

The check derived its grouping level from the flag's infix (`S` → body system,
`P` → dictionary term, anything else → subject), and `01` is a sponsor
numbering that says nothing about scope. There is no way to recover
period-scoping from the name, so it has to be read from the data. Fixed by
adding `APERIOD` to the grouping for *every* flag whenever it is present and
populated; a single-period study is unaffected because `APERIOD` is then
constant.

Note the same bug would have hit `AOCCSFL`/`AOCCPFL` in a multi-period study —
the flaw was not "numbered flags" but the assumption that flags are never
period-scoped, which is wrong for the whole vaccine therapeutic area.

The alternative considered and rejected: skip numbered flags entirely rather
than guess their level. That would silence the false positive without the
`APERIOD` handling, but it gives up on the most common flag form. The trade
favours grouping by period, because the check only reports *more than one*
`"Y"` as a defect — so the failure mode of over-grouping is a missed duplicate
(quiet), while under-grouping is a false alarm on correct data (loud, and
exactly what the 2026-08-03 review objected to).

### Sweep status after all fixes

Zero structure-check warnings, zero false-positive check failures, and three
datasets reporting genuine duplicate records (`adeg` 63, `advs` 39, `adpp`
1,008 — all verified by hand, see `keys_design_notes.md`). Re-run
`worked_example_pharmaverseadam.R` after any change to type detection, key
inference, or the checks.

## Post-2026-08-03 backlog

Facts only, in rough order of value:

1. **TTE** — the per-`PARAMCD` rework: events / censored per endpoint (the
   single total is misleading with more than one endpoint), censoring reason
   distribution, follow-up time. TTE is the one type with no per-type section
   at all, so this is the largest remaining gap.
2. **ADSL** — `DCSREAS` alongside `EOSSTT`, `RACE` folded into the demographics
   line, multi-period detection (`TRT02P` and beyond).
3. **OCCDS** — pre-treatment record counts as their own fact line;
   relatedness (`AEREL`) and outcome breakdowns.
4. **Cross-cutting** — date coverage (min/max of the main date variable) and
   imputation flag counts (`ADTF`, `ATMF`, `ASTDTF`), which are facts about
   what a derivation did, not judgements about whether it was right.
5. **`detail =`** — per design question 1, most likely on `print()` rather than
   `summary()`, with everything always computed into the object.

Deliberately not on this list: anything from the red-flag catalogues above.

## Brainstorm 2026-08-01: more ADSL content

The first ADSL increment is implemented (treatment, populations, `TRTDURD`,
`EOSSTT`, demographics, plus the `ARM`/`ACTARM` and `SAFFL`-without-`TRTSDT`
checks — **the checks were removed 2026-08-03**, the facts remain).
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
2026-08-01**. Superseded 2026-08-03: deaths and the funnel are facts and stayed;
`TRTEDT < TRTSDT` and the flag domain were removed with the rest of the ADSL
checks, and `print_adsl_summary()` no longer has a `checks` list. The remaining
fact items (`DCSREAS`, race, multi-period) are still open; the remaining date
checks are for the ADaM checks package.

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
