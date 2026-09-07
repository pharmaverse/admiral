# Rolling out the `admiral_df` / `summary()` feature — staged PR plan

The current branch (`3160_poc_liam`) lands the whole feature — metacore key
extraction, tagging, type detection, key inference, three per-type
summarizers, print methods, docs, worked examples — as one commit
("Everything and the kitchen sink"). That's ~130 files and 17k+ insertions,
which is too much for one review pass: a reviewer can't evaluate the ADSL
summary content and the OCCDS occurrence-flag heuristic and the `summary()`
dispatch decision all at once without one of them getting rubber-stamped.

This splits the branch into a dependency-ordered sequence of small,
independently mergeable PRs. Each stage is reviewable on its own, has its
own tests, and leaves `main` in a working state (nothing downstream depends
on unmerged stages until it exists). Companion to `summary_admiral_df_design_notes.md`,
`keys_design_notes.md`, `adsl_comparison_design_notes.md`, and
`summary_admiral_df_review_notes.md`.

## Guiding principles for splitting

- **Each PR should be reviewable in isolation** — a reviewer shouldn't need
  to hold the whole feature in their head to judge one PR.
- **Land plumbing before policy.** Generic infrastructure (tagging,
  attribute handling) before ADaM-specific heuristics (type detection, key
  inference) before presentation (print methods) before the highest-risk
  decision (hijacking `summary()`).
- **Every stage ships with its own tests** — don't defer test coverage to
  a later "add tests" PR.
- **Defer the riskiest, least-reversible decision as long as possible.**
  Attaching this to `summary()` (see `summary_admiral_df_review_notes.md`
  §3) is the one design choice that's hard to walk back once released, so
  it should be the last thing merged, once the reviewer has already seen
  and accepted everything it depends on.
- **Docs/vignettes/showcase ride with the functionality they document**, not
  as one giant doc dump at the end.

## Stage 1 — Extract dataset keys from a `metacore` spec (this is the starting point)

**Ships:** `get_admiral_keys()`, `set_admiral_keys()`, `as_admiral_df()`
(tagging only — just adds the class + optional `"admiral_keys"` attribute,
no `summary()` method yet).

**Why first:** it's fully self-contained. `metacore` is already a
`Suggests` dependency added in this branch; nothing else in the feature
(type detection, inference, summarizers) is needed for this to be useful
and correct on its own. It's also the piece with the clearest, narrowest
contract — "read `key_seq` from a `Metacore` object, store it on the
dataset" — so it's the easiest PR for a reviewer to fully verify.

**Reviewable surface:**
- `get_admiral_keys(metacore, dataset_name)` — thin wrapper over
  `metacore::get_keys()`, single-dataset-spec convenience, the
  `zero-keys-defined` warning path.
- `set_admiral_keys(dataset, metacore, dataset_name)` — applies the above
  and tags the dataset.
- `as_admiral_df(dataset, keys)` — the tagging primitive both of the above
  (and every later stage) build on.

**Tests:** the `get_admiral_keys()`/`set_admiral_keys()`/`as_admiral_df()`
subset of `test-admiral_df.R` (spec with keys, spec with no keys defined,
single- vs. multi-dataset spec, missing/invalid `metacore` object).

**Explicitly out of scope for this PR:** no `summary()` dispatch exists
yet, so tagging a dataset here has no visible effect beyond the attribute —
call that out in the PR description so reviewers don't go looking for
behavior that arrives in a later stage.

**Docs:** roxygen for the three functions, `NEWS.md` entry noting these are
new exported utilities with no behavioral effect yet (summary() lands
later).

## Stage 2 — `admiral_df` tag auditing

**Ships:** `check_admiral_df()` + `print.admiral_df_check()`.

**Why here:** it's the natural "does my tag survive this pipeline step"
companion to Stage 1, still has zero dependency on type detection or
summarization, and gives reviewers early visibility into the
attribute-fragility tradeoffs (documented at length in `admiral_df.R`'s
`set_admiral_keys()`/`check_admiral_df()` roxygen) *before* those tradeoffs
matter for a user-facing `summary()`.

**Tests:** all four `is_admiral_df`/`keys`/`stale_keys` combinations
documented in `check_admiral_df()`'s own `@details` (class+keys both
present, class only, class only with zero-length keys, keys only).

## Stage 3 — ADaM dataset type detection

**Ships:** `get_admiral_df_type()`, `is_adsl_structure()`.

**Why here:** pure, dependency-free classification logic (input: a data
frame; output: a string) with no dispatch or side effects yet — reviewable
purely as "do these heuristics classify a `PARAMCD`/`AVAL`/`--DECOD`-bearing
data frame correctly," independent of anything downstream.

**Reviewable surface:** the precedence table in `get_admiral_df_type()`'s
`@details`, and specifically the ambiguity risks already flagged in
`summary_admiral_df_review_notes.md` (OCCDS regex over-matching, ADSL
false-negative on a work-in-progress dataset) — worth resolving or
consciously deferring *here*, since every later stage inherits whatever
this stage decides.

**Tests:** one fixture per type, plus the currently-missing precedence
*interaction* test (a dataset simultaneously matching two types) and the
WIP-ADSL edge case from the review notes.

## Stage 4 — Record-structure inference fallback

**Ships:** `minimal_unique_key()`, `infer_admiral_keys()`.

**Why here:** builds directly on Stage 3's type detection and Stage 1's
key-attribute mechanism (this is the third rung of the `keys=` >
`set_admiral_keys()` > inferred precedence chain), but still has no
dependency on any `summarize_*`/print code — it can be fully tested by
asserting on its return vector directly.

**Reviewable surface:** this is the single highest-value-per-line-of-review
stage, since it's the piece `keys_design_notes.md` documents as having
caught real bugs against `{pharmaverseadam}`. Worth reviewing with that
document open, and worth landing the `worked_example_pharmaverseadam.R`
sweep as an automated (opt-in/`skip_on_cran`) regression test in this same
PR rather than leaving it a manual script — see
`summary_admiral_df_review_notes.md` §3 quick-win 4.

**Tests:** the existing per-type inference tests, plus the multi-period
`APERIOD` case, plus an all-`NA` `ASEQ` case (review notes finding #2).

## Stage 5 — `summary.admiral_df()` core (no per-type sections)

**Ships:** `summary.admiral_df()` and `print.summary_admiral_df()`, but
trimmed to only the type-independent facts: `type`, `n_obs`, `n_vars`,
`n_subjects`, and the record-structure/duplicate-key check. Per-type
sections (`$adsl`, `$bds`, `$occds`, `$vs_adsl`) are stubbed out /
deferred to later stages — this PR is purely about the dispatch mechanism
and the shared skeleton.

**Why here, and why last-but-one:** this is the PR that actually makes
`summary()` behave differently for a tagged data frame — the single
riskiest, least-reversible decision in the whole feature (see
`summary_admiral_df_review_notes.md` §3, "Reconsider hijacking `summary()`").
By this point the reviewer has already independently signed off on tagging,
type detection, and key inference, so this PR is scoped to exactly the
one open question that matters: *should this be `summary()`, or a
dedicated verb like `admiral_summary()`?* Isolating it here means that
question gets its own focused discussion/decision instead of being buried
under per-type formatting review.

**Tests:** `nrow() == 0` handling (review notes finding #6), the
`ds_name`-capture-via-`substitute()` behavior and its documented limits
(doesn't work through `%>%`), malformed-`adsl=`-argument abort path.

## Stage 6 — ADSL section

**Ships:** `summarize_adsl()`, `print_adsl_summary()`, `render_summary_table()`
(introduced here as the shared table renderer, since ADSL is the first
section that needs it — reused unchanged by Stages 7–8).

**Why here:** ADSL is the simplest per-type section (no uniqueness check,
per the 2026-08-03 scope decision — pure facts), so it's the best place to
get reviewer sign-off on the *shape* of a per-type section and the shared
table-rendering helper before the more heuristic-heavy BDS/OCCDS sections
build on the same pattern.

## Stage 7 — BDS section

**Ships:** `summarize_bds()`, `print_bds_summary()`.

**Why here, after ADSL:** introduces the first uniqueness check
(`multiple_baselines`) and the `max_params` truncation pattern. Land the
`APERIOD` grouping fix (review notes finding #1) as part of this PR, not
as a follow-up — the bug is about a check being introduced in this exact
PR, so there's no reason to ship it broken and fix it later.

## Stage 8 — OCCDS section

**Ships:** `summarize_occds()`, `print_occds_summary()`.

**Why last of the per-type sections:** the most heuristic-heavy of the
three (`count_dups`'s flag-infix inference, `--DECOD`/`--BODSYS` fallbacks)
— reviewing it after ADSL/BDS have established the pattern lets the
reviewer focus entirely on the OCCDS-specific logic. Include the
"ambiguous flag infix" fallback (review notes quick-win 5) here rather than
after the fact.

## Stage 9 — TTE handling + cross-dataset (`adsl=`) comparison

**Ships:** TTE-specific fields in `summary.admiral_df()` (`n_events`,
`n_censored`), `summarize_vs_adsl()`, `print_vs_adsl_summary()`.

**Why last:** both pieces are explicitly flagged as still-open design
questions in the design notes (TTE per-`PARAMCD` rework,
`set_admiral_options()`-based `adsl=` default). Landing them last means
they can absorb whatever conventions were settled in Stages 6–8 without
forcing an early decision on the two least-settled parts of the design.
Cap `orphans` for printing (review notes quick-win 6) as part of this PR.

## Stage 10 — Wire into `derive_param_computed()`, docs, examples

**Ships:** the `as_admiral_df()` call inside `derive_param_computed()`,
`admiral_df_showcase.qmd`, `worked_example*.R` scripts, and any remaining
vignette touch-ups.

**Why last:** these are consumers/demonstrations of the finished feature,
not part of it — reviewing them before Stages 1–9 have landed means
reviewing usage of an API that doesn't exist on `main` yet. Landing them
last also means the showcase/worked examples exercise the *actual* merged
behavior rather than a snapshot of the POC that may have changed during
review.

## Explicitly deferred (raise as follow-up issues, not blockers)

Several findings in `summary_admiral_df_review_notes.md` are real but
shouldn't gate this rollout — file them as separate issues once the
relevant stage lands:

- Consolidating the multiple full-table `distinct()`/`group_by()` passes
  into fewer (performance; needs a `bench::mark()` baseline first).
- Making the key-inference candidate list extensible via
  `set_admiral_options()`.
- A session-wide `adsl=` default via `set_admiral_options()`.
- Broadening `as_admiral_df()` tagging into other `derive_*` functions
  beyond `derive_param_computed()`.

## Summary table

| Stage | Ships | Depends on |
|---|---|---|
| 1 | `get_admiral_keys()`, `set_admiral_keys()`, `as_admiral_df()` (tagging) | — |
| 2 | `check_admiral_df()`, `print.admiral_df_check()` | 1 |
| 3 | `get_admiral_df_type()`, `is_adsl_structure()` | — |
| 4 | `minimal_unique_key()`, `infer_admiral_keys()` | 1, 3 |
| 5 | `summary.admiral_df()` core, `print.summary_admiral_df()` | 1, 3, 4 |
| 6 | `summarize_adsl()`, `print_adsl_summary()`, `render_summary_table()` | 5 |
| 7 | `summarize_bds()`, `print_bds_summary()` | 6 |
| 8 | `summarize_occds()`, `print_occds_summary()` | 6 |
| 9 | TTE fields, `summarize_vs_adsl()`, `print_vs_adsl_summary()` | 6–8 |
| 10 | `derive_param_computed()` wiring, showcase, worked examples, vignettes | 1–9 |
