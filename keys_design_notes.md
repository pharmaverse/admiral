# Key pathways in `summary.admiral_df()` — design assessment

Investigation of the two ways key variables reach the record-structure check
(2026-08-01, companion to issue #3160 and
`summary_admiral_df_design_notes.md`). Question: both the `keys =` argument and
the `admiral_keys` attribute exist — is that redundant, and which should stay?

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
