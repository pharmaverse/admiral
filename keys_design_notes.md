# Key pathways in `summary.admiral_df()` — design assessment

Investigation of the two ways key variables reach the record-structure check
(2026-08-01, companion to issue #3160 and
`summary_admiral_df_design_notes.md`). Question: both the `keys =` argument and
the `admiral_keys` attribute exist — is that redundant, and which should stay?

**Unaffected by the 2026-08-03 core team decision.** All three key pathways
(inference, explicit user input, `{metacore}`) were presented at that session
and drew no objection; the discussion was about the data-quality checks, which
were trimmed to a core set — see `summary_admiral_df_design_notes.md`. If
anything the decision strengthens the case for this layer: the record structure
check is now the *paradigm* for what {admiral} will check, and everything below
determines what it checks against.

## What is in the code

`summary.admiral_df()` resolves keys in strict precedence order:

1. **Supplied** — the `keys =` argument (character vector, or a single-dataset
   `Metacore` object resolved via `get_admiral_keys()`). Reported as
   `key_source = "supplied"`, printed as `Structure (supplied): ...`.
   Transient: affects only that call, never touches the object.
2. **Declared** — the `admiral_keys` attribute. Three writers:
   `as_admiral_df(keys =)` (used inside `derive_param_computed()`, which tags
   its output with `by_vars` + `PARAMCD`), `set_admiral_keys()` (the metacore
   path, which also stores the dataset name), and a manual `attr<-`. Printed
   without annotation — it is the trusted default.
3. **Inferred** — `infer_admiral_keys()` guesses from type and data. Printed as
   `(inferred)`.

Both supplied and declared keys get the stale-key warning when they name absent
variables. `check_admiral_df()` audits the declared/inferred pathways only —
correct, since it audits the object, not future calls.

## Tradeoffs

The pathways differ on three axes: **who** sets the keys, **when**, and **how
long they last**.

**The attribute is the producer's channel.** The function that *creates* a
dataset knows its intended structure better than anyone downstream —
`derive_param_computed()` knows the output is one record per `by_vars` +
`PARAMCD`; the spec knows what ADVS's keys are supposed to be. Persisting that
knowledge is what makes the core #3160 story work: `summary(map)` right after a
derivation checks the *intended* structure with zero ceremony, and still works
three pipeline steps later. Weaknesses are the flip side of persistence: it is
invisible state (two visually identical tibbles can summarize differently), it
can be silently dropped (`summarise()`, `bind_rows()` when not first, pivots)
or go stale (`rename()`). The mitigations exist for exactly this: the stale-key
warning and `check_admiral_df()`.

**The argument is the consumer's channel.** It answers a different question:
not "is this dataset what its producer intended" but "is this dataset unique by
*the keys I say*". It is explicit at the call site (reproducible when reading
the script), cannot go stale, requires no mutation of the object for a one-off
question, and gives a direct metacore route without a tagging step. It is also
the escape hatch when the declared keys are wrong or a stricter hypothesis
should be tested. Weakness: the knowledge evaporates after the call, and
habitual use means repetition — forget it once and the summary silently falls
back to inference, which may validate a weaker key than intended.

## Which is better?

Neither subsumes the other: they are not two solutions to one problem but one
solution each to two different problems (persistent intent vs. call-time
question). The layering — *explicit argument overrides persistent metadata
overrides heuristic* — is a well-worn R idiom (`print(x, digits =)` over
`getOption("digits")`; contrast arguments over contrast attributes). Dropping
either has a real cost:

- Drop the argument → users must mutate an object just to ask a question, and
  the override for a wrong attribute disappears.
- Drop the attribute → the derivation tagging dies, and with it the
  zero-ceremony post-derivation check that motivated #3160. This one is
  load-bearing.

**Recommendation: keep both.** The genuine risk of dual pathways is
comprehension — "why did summary check *these* keys?" — and that is already
addressed: `key_source` is stored in the returned object, the print annotates
`(supplied)`/`(inferred)` while leaving the trusted declared case unmarked, and
`check_admiral_df()` explains the object's state. The precedence is the right
shape: most-explicit wins, heuristic last.

## One asymmetry to document

`summary(x, keys = spec)` is transient while `set_admiral_keys(x, spec)` is
persistent — same spec, different lifetime. That is coherent (it mirrors
argument-vs-options semantics), but it is the one place a user could plausibly
expect the call to have tagged the object when it did not. If that confusion
shows up in review, the fix is documentation (vignette sentence), not
collapsing the pathways.

## Evidence: the {pharmaverseadam} sweep (2026-08-07)

Until now pathway 3 (inference) had only ever been exercised against hand-built
fixtures written by the same person who wrote the inference. `worked_example_pharmaverseadam.R`
runs `summary()` over all 30 datasets in {pharmaverseadam} — 58 calls, since
every non-ADSL dataset is also summarized with `adsl =` — and it found two real
bugs plus one design flaw. The ADSL is chosen by `USUBJID` overlap rather than
by name, because the therapeutic-area datasets do not share a subject universe.

### Bug 1: OCCDS required the literal variable `ASEQ`

`infer_admiral_keys()` special-cases OCCDS because there is no analysis-value
structure to fall back on — two adverse events for one subject need not differ
in any analysis variable, so the sequence number *is* the record key. That part
is right. What was wrong was insisting on the name `ASEQ`: occurrence datasets
conventionally carry the SDTM domain sequence (`AESEQ`, `CMSEQ`, `MHSEQ`), and
`ASEQ` is the generic ADaM analysis sequence many sponsors never populate.
Three of the four OCCDS datasets had only the domain form, so the record
structure check — the feature that motivated #3160 — silently did not run on
the most commonly used occurrence datasets.

Fixed by preferring `ASEQ`, then falling back to a single `^[A-Z]{2}SEQ$`
match. The two-letter domain prefix is what keeps `SRCSEQ` out: that is
provenance added by `derive_vars_merged()`, not a record key. Several domain
sequences means the dataset merges two domains and the key is genuinely
ambiguous, so that case warns rather than guessing.

### Bug 2: the candidate list assumed findings-shaped data

The `extra` list (`AVISITN`, `AVISIT`, `ATPTN`, `ATPT`, `ADTM`, `ADT`,
`APERIOD`, `APERIODC`, `ASPID`) describes one dataset shape: findings keyed by
visit and timepoint. Three other shapes exist in {pharmaverseadam} and none was
covered:

| Dataset | Real key | What was missing |
|---|---|---|
| `adex` (exposure) | one record per dosing **interval** | `ASTDTM`/`ASTDT`/`AENDT` — keyed by interval start, not analysis date |
| `adppk` (population PK) | one record per **relative time** | `NFRLT`/`AFRLT` — no visit structure at all |
| `adpc` (PK concentration) | derived records alongside source | `DTYPE` — a derived record shares every analysis variable with its source |

All three fell back to `USUBJID + PARAMCD` and reported thousands of false
duplicates. `DTYPE` in particular affects any dataset with LOCF or averaged
records, which is most real BDS data.

### The design flaw: prefix-minimal is not subset-minimal

`minimal_unique_key()` walks the candidates in priority order and returns the
first unique **prefix**, so it carries every variable it passed over on the way
— whether or not that variable discriminates. Lengthening the candidate list to
fix bug 2 therefore made keys *worse*: `adpc` reached uniqueness at 14
variables, and the summary would have printed "one record per USUBJID,
PARAMCD, AVISITN, AVISIT, ATPTN, ATPT, NFRLT, AFRLT, ADTM, ADT, ASTDTM, ASTDT,
AENDT, DTYPE" — a claim about the intended structure that is simply false.

Fixed with a backward-elimination pass after uniqueness is reached: drop any
added variable whose removal preserves uniqueness, lowest priority first.
`adpc` collapses to `USUBJID + PARAMCD + AVISITN + ATPTN + DTYPE`. It also
shortened keys that were never broken — `adlb` dropped `AVISIT` as redundant
with `AVISITN`, `advs_metabolic` went from seven variables to four.

Note the interaction with an existing behaviour: pruning drops `AVISITN` from a
single-visit dataset, because there it genuinely does not discriminate. That is
the same "inferred key is narrower than the intended key" property inference
has always had, and the reason the declared and supplied pathways exist.

### Outcome

| | Before | After |
|---|---|---|
| Warnings (structure not checked) | 6 | 0 |
| Datasets reporting duplicates | 6 of 30 | 3 of 30 |
| Total duplicate records reported | 15,729 | 1,110 |

The 1,110 that remain were each checked by hand and are **true positives in the
test data**, not inference failures — recorded here so nobody re-diagnoses
them:

- `advs` (39) and `adeg` (63) — two SDTM `VISIT`s mapped onto one `AVISIT`,
  with `AVAL` differing in most of the pairs.
- `adpp` (1,008) — 336 groups of four records identical on every analysis
  variable, differing only by `SRCSEQ`. No semantic key separates them; this is
  a fanned-out merge in the data. `SRCSEQ` must not be added to the candidate
  list to make it disappear — that is precisely the surrogate-key masking
  `minimal_unique_key()` exists to prevent.

### What this says about the design

Inference is a heuristic whose quality is entirely a function of whether its
candidate list matches the dataset's shape, and there turned out to be three
shapes it had never heard of. Every one was found by running it on data nobody
wrote for it. That is an argument for keeping the declared and supplied
pathways prominent — not against inference, which now gets 27 of 30 real
datasets right with no configuration at all.

Open question it raises: should the candidate list be extensible (a
`set_admiral_options()` entry) rather than hard-coded? Sponsors with in-house
dataset shapes will hit exactly this. Against it: an option that changes what
the structure check means is a footgun, and the honest answer for an unusual
shape is `set_admiral_keys()`. Worth deciding only if someone asks.
