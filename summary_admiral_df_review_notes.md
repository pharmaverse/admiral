# `summary.admiral_df()` — code review findings and suggested improvements

Review of the `admiral_df` diagnostic infrastructure (`R/admiral_df.R` and
supporting design docs) as it stands on `3160_poc_liam`. Companion to
`summary_admiral_df_design_notes.md`, `keys_design_notes.md`, and
`adsl_comparison_design_notes.md` — this document doesn't re-derive the
architecture, it audits the current implementation against it and proposes
concrete next steps.

Status: **notes for discussion**, written from a read-through of
`R/admiral_df.R`, the three `*_design_notes.md` files, `test-admiral_df.R`,
the `worked_example*.R` scripts, and `admiral_df_showcase.qmd`.

## 1. Architecture recap

- `as_admiral_df()` / `set_admiral_keys()` tag a data frame with the
  `admiral_df` S3 class (+ optional `"admiral_keys"` attribute), which makes
  `summary()` dispatch to `summary.admiral_df()`.
- `get_admiral_df_type()` classifies a dataset as `ADSL`/`BDS`/`OCCDS`/`TTE`/
  `other` from column-name heuristics (fixed precedence: `TTE` > `OCCDS` >
  `BDS` > `ADSL` > `other`).
- `infer_admiral_keys()` + `minimal_unique_key()` guess the record structure
  ("one record per...") when keys weren't declared, with precedence
  `keys=` argument > `"admiral_keys"` attribute (`set_admiral_keys()`/
  `metacore`) > inferred.
- `summarize_adsl()` / `summarize_bds()` / `summarize_occds()` /
  `summarize_vs_adsl()` each build a small list of facts plus exactly one
  uniqueness/referential-integrity check apiece (`multiple_baselines`,
  `occ_flag_dups`, `n_orphans`) — the scope agreed in the 2026-08-03 decision
  in `summary_admiral_df_design_notes.md`.
- `summary.admiral_df()` assembles everything into a `summary_admiral_df`
  object; `print.summary_admiral_df()` and friends render it.
- `check_admiral_df()` is a separate auditor for whether the tag/keys
  survived a pipeline.

Only one production derivation, `derive_param_computed()`, currently calls
`as_admiral_df()` — everywhere else the tag is applied manually.

## 2. Findings

### Correctness / bugs

| # | Finding | Where |
|---|---|---|
| 1 | `summarize_bds()$multiple_baselines` doesn't group by `APERIOD`, even though the *identical* false-positive bug was found and fixed for OCCDS's `occ_flag_dups`. Multi-period BDS studies (e.g. vaccine trials with repeat baselines per period) will false-positive on this today. | `admiral_df.R:817-824` vs. `943-951`; flagged but not generalized in `keys_design_notes.md:262` |
| 2 | `infer_admiral_keys()` prefers `ASEQ` for OCCDS record keys without checking it's populated. Other candidate picks (`trt_var`, `sev_var`) skip all-`NA` columns; `ASEQ` doesn't, so an empty sequence variable can silently "win" and produce spurious duplicate/unique results. | `admiral_df.R:354-358` vs. `665-670`, `917-922` |
| 3 | OCCDS type detection (`(DECOD\|TERM)$\|^AOCC.*FL$\|^TRTEMFL$`) matches any `*TERM`-suffixed column with no structural corroboration (e.g. `USUBJID` presence) — a non-OCCDS dataset can misclassify as OCCDS. | `admiral_df.R:74` |
| 4 | `is_adsl_structure()` can misfire on a legitimate work-in-progress ADSL (no `TRT01P`/`TRTSDT` yet, subject key not yet unique) — it degrades to the *wrong type* (`"other"`, whole ADSL section skipped) rather than a partial ADSL summary. | `admiral_df.R:98-105` |
| 5 | `summarize_occds()`'s `count_dups` infers occurrence-flag grouping level from a substring match on the flag name's infix (`"S"` vs `"P"`). A non-standard/ambiguous flag name (e.g. containing both letters, or a numbered sponsor-specific flag) is silently mis-scoped, with no "cannot determine level" fallback the way there is for a missing `--BODSYS`/`--DECOD` variable. | `admiral_df.R:952-974` vs. `959-965` |
| 6 | `nrow(dataset) == 0` is never exercised — an untested, plausible intermediate pipeline state across `get_admiral_df_type()`, `infer_admiral_keys()`, and `summary.admiral_df()`. | throughout |
| 7 | `n_orphans`/`orphans` stores every orphan `USUBJID` with no cap, unlike the `max_params` truncation already used for BDS params — a potential print/memory issue if a wrong-ADSL-version mistake produces thousands of orphans. | `admiral_df.R:1048-1050` vs. `1561` |

### Fragility / silent failure modes

- **Attribute loss is invisible.** The `"admiral_keys"` attribute is
  documented as dropped by `summarise()`, `group_by()`/`ungroup()`, pivots,
  non-first `bind_rows()`, and right-side joins — all common pipeline steps —
  with no warning at drop time (`admiral_df.R:212-224`). The only signal is
  the `(declared)` vs `(inferred)` annotation in printed output, easy to
  miss.
- **Class-dropped-but-attribute-kept is the worst case**: `summary()` then
  silently falls back to `summary.data.frame()` with zero indication the
  diagnostic is gone (`admiral_df.R:428-431`).
- Keys are stored as plain variable names, so `dplyr::rename()`/`select()`
  can silently produce stale keys; the `stale_keys` warning only fires at
  the next `summary()` call, so a wrong key can persist unnoticed through
  many pipeline steps.

### API design concerns

- **Hijacking `summary()`** is the biggest structural risk. `summary()` has
  a long-standing contract (numeric/statistical digest — cf.
  `summary.data.frame`, `summary.lm`). Any code that structurally inspects
  `summary(df)` rather than just printing it (generic QC scripts,
  `purrr::map(datasets, summary)`, etc.) gets a completely different object
  once a data frame is tagged `admiral_df`. The design notes debate *what*
  summary should contain but not *whether* `summary()` is the right generic
  to attach it to.
- **Narrow, inconsistent wiring**: only `derive_param_computed()` tags its
  output today, so the "declared keys" path — the reliable one — almost
  never fires in practice despite being the headline usage story in the
  design notes and showcase doc.
- **Naming idiom mismatch**: `as_admiral_df()`/`get_admiral_df_type()`/
  `infer_admiral_keys()` read as S3-scaffolding/internal-infrastructure
  names, a different idiom from admiral's verb-first public API
  (`derive_param_computed`, `get_summary_records`). `check_admiral_df()` in
  particular sits oddly next to `{admiraldev}`'s `assert_*`/`check_*`
  validators, which *error*; this one only reports.
- **Two independent, per-call precedence chains** (`keys=` vs. attribute vs.
  inference; `adsl=` vs. no session default) both currently resolve friction
  points with "call the function again," where the rest of admiral solves
  the equivalent problem via `set_admiral_options()` (e.g. `subject_keys`).
  Both `keys_design_notes.md:188-192` and
  `adsl_comparison_design_notes.md:141-147` flag this and defer it.
- **Broad exported surface for an evolving POC** — the returned-object shape
  is asserted on extensively in tests while the design notes describe the
  TTE section and `adsl=` defaulting as still open.

### Performance

- `minimal_unique_key()`'s `is_unique()` re-runs a full-table `distinct()`
  once per candidate added and once per candidate considered for removal —
  up to ~2× the length of the `extra` candidate list (13 variables) in full
  table passes, just to infer keys.
- `summary.admiral_df()` overall runs at least 5–7 independent full-table
  `distinct()`/`group_by()+summarise()` passes for a BDS dataset
  (`n_subjects`, key-uniqueness, `params`, `avisits`,
  `summarize_bds()`'s per-parameter grouping, `multiple_baselines`). None
  are quadratic, but the constant-factor cost is nontrivial for something
  meant to be a casual "run after every derivation" diagnostic on
  multi-million-row ADLB/ADPC-scale data. No `bench::mark()` regression test
  exists despite cost being flagged as an open concern in
  `summary_admiral_df_design_notes.md:184-186`.

### Test coverage gaps

- No test for the type-detection precedence *interaction* (verifying `TTE`
  actually wins over `OCCDS` when both signals are present) — only each
  type in isolation.
- No test for the `is_adsl_structure()` false-negative WIP-ADSL case.
- No test for `multiple_baselines` under a multi-period design (the exact
  gap in finding #1 above — OCCDS has this test, BDS doesn't).
- No test for large `n_orphans` (print truncation / object size).
- No `nrow() == 0` test.
- `worked_example_pharmaverseadam.R` — the only thing that has ever caught
  real inference bugs (see `keys_design_notes.md`) — is a manual script, not
  wired into `testthat`/CI, so a regression there currently can't fail a
  build.

## 3. Suggested improvements, prioritized

### Quick wins

1. Fix `multiple_baselines` in `summarize_bds()` to include `APERIOD` in its
   grouping when present and populated, mirroring the existing
   `occ_flag_dups` fix.
2. Guard `ASEQ` selection in `infer_admiral_keys()` against an all-`NA`
   column, consistent with `trt_var`/`sev_var` elsewhere in the file.
3. Add an explicit `nrow(dataset) == 0` fast path (and test) to
   `summary.admiral_df()`, `get_admiral_df_type()`, and
   `infer_admiral_keys()`.
4. Wire `worked_example_pharmaverseadam.R` into an opt-in
   (`skip_on_cran`-style) automated regression test.
5. Add a "cannot determine occurrence-flag level" fallback in
   `summarize_occds()`'s `count_dups` for ambiguous/non-standard flag
   infixes, rather than silently defaulting to a possibly-wrong grouping.
6. Cap/truncate `vs_adsl$orphans` for printing (mirroring the `max_params`
   pattern), keeping the full list in the returned object.
7. Tighten the OCCDS type-detection regex to require co-occurrence with
   `USUBJID` or another structural signal, not just a `*TERM`-suffixed
   column name.

### Bigger, structural

1. **Reconsider hijacking `summary()`.** Expose the diagnostic under a
   distinct name (e.g. `admiral_summary()`), keeping
   `summary.admiral_df()` as a thin, optional S3 convenience wrapper — this
   preserves the "quick diagnostic after a derivation" ergonomics without
   silently redefining `summary()`'s contract for any tagged tibble that
   leaks into a wider pipeline.
2. **Broaden tagging beyond `derive_param_computed()`**, either by wiring
   `as_admiral_df(keys=)` into other `derive_*` functions with well-defined
   output keys, or by providing a single wrapping helper that tags a
   derivation chain's result from its `by_vars` automatically.
3. **Make attribute loss visible.** Warn when stored metadata
   (`attr(object, "admiral_ds_name")`) is present but `"admiral_keys"` is
   `NULL` and the class survived, rather than silently falling back to
   inference with no signal that declared intent was lost.
4. **Make the key-inference candidate list extensible** via
   `set_admiral_options()` (mirroring the existing `subject_keys`
   precedent), since the `{pharmaverseadam}` sweep already proved real
   dataset shapes fall outside the hard-coded candidate list, and extending
   it today requires a package code change.
5. **Consolidate redundant full-table passes** in `summary.admiral_df()`
   into fewer shared aggregates, backed by a `bench::mark()`-based
   performance regression test.
6. **Stabilize the returned-object contract before wider release** — add an
   experimental lifecycle badge, or hold export until the TTE section and
   `set_admiral_options()` integration (candidate list, `adsl=` default) are
   settled.
7. **Add a session-wide `adsl=` default via `set_admiral_options()`**, as
   already proposed in `adsl_comparison_design_notes.md:141-147`, so the
   argument doesn't need repeating on every `summary()` call in an analysis
   script.
