# Demo: `admiral_df` class, dataset keys, and `summary()` override
#
# Showcases the functionality added on branch 3209-5-summary-override:
#   1. `as_admiral_df()`         - admiral dataset functions tag their output
#   2. `set_admiral_keys()`      - store key variables and a dataset name
#   3. `print.admiral_df()`      - header line with name, type, and keys
#   4. `get_admiral_df_type()`   - ADaM dataset type detection (internal)
#   5. `infer_admiral_keys()`    - record structure inferred from the data (internal)
#   6. `summary.admiral_df()`    - type-specific report + record structure check
#
# Run from the package root with the development version loaded.
# Not part of the package: do not commit (or add to .Rbuildignore).

devtools::load_all()
library(dplyr, warn.conflicts = FALSE)
library(tibble)

adsl <- pharmaverseadam::adsl
advs <- pharmaverseadam::advs
adae <- pharmaverseadam::adae
adtte <- pharmaverseadam::adtte_onco

# 1. Admiral functions tag their output -------------------------------------

# A plain tibble is not an admiral_df ...
class(advs)

# ... but the output of an admiral derivation is
advs_trt <- derive_vars_merged(
  select(advs, -TRT01P),
  dataset_add = adsl,
  by_vars = exprs(STUDYID, USUBJID),
  new_vars = exprs(TRT01P)
)
class(advs_trt)

# With no name or keys set it prints exactly as before (no header)
advs_trt %>% select(USUBJID, TRT01P, PARAMCD, AVISIT, AVAL) %>% head(3)

# Tagging by hand
class(as_admiral_df(adsl))

# 2. Setting keys -------------------------------------------------------------

advs_keyed <- set_admiral_keys(
  advs,
  keys = exprs(USUBJID, PARAMCD, AVISIT, ATPT),
  dataset_name = "ADVS"
)
attr(advs_keyed, "admiral_keys")
attr(advs_keyed, "admiral_ds_name")

# Keys from a spec arrive as a character vector: convert with chr2vars()
spec_keys <- c("USUBJID", "PARAMCD", "AVISIT", "ATPT")
advs_keyed <- set_admiral_keys(advs_keyed, keys = chr2vars(spec_keys))
attr(advs_keyed, "admiral_ds_name") # name kept from the previous call

# Passing the character vector directly gives a pointer to chr2vars()
try(set_admiral_keys(advs, keys = spec_keys))

# A typo in a key warns and leaves the previous keys untouched
advs_typo <- set_admiral_keys(advs_keyed, keys = exprs(USUBJID, PARMCD))
attr(advs_typo, "admiral_keys")

# Empty keys record "no key variables defined" (distinct from never set)
attr(set_admiral_keys(advs, keys = exprs()), "admiral_keys")

# 3. Print header ---------------------------------------------------------------

advs_keyed %>% select(USUBJID, PARAMCD, AVISIT, ATPT, AVAL)

set_admiral_keys(adsl, keys = exprs(USUBJID), dataset_name = "ADSL") %>%
  select(USUBJID, TRT01P, SAFFL)

set_admiral_keys(adsl, keys = exprs()) %>% select(USUBJID)

# Keys and class survive most dplyr verbs ...
advs_keyed %>%
  filter(PARAMCD == "SYSBP") %>%
  mutate(AVAL2 = AVAL * 2) %>%
  select(USUBJID, PARAMCD, AVISIT, ATPT, AVAL2)

# ... but group_by() drops the class (attributes kept) and summarise() drops both
advs_keyed %>% group_by(PARAMCD) %>% class()
advs_keyed %>% group_by(PARAMCD) %>% attr("admiral_keys")
advs_keyed %>% summarise(n = n()) %>% attributes() %>% names()

# 4. Dataset type detection (internal) ---------------------------------------

admiral:::get_admiral_df_type(adsl) # "ADSL"
admiral:::get_admiral_df_type(advs) # "BDS"
admiral:::get_admiral_df_type(adae) # "OCCDS"
admiral:::get_admiral_df_type(adtte) # "TTE"

# A lookup table sharing ADaM variable names is not mistaken for ADaM data
admiral:::get_admiral_df_type(distinct(adae, AEDECOD, AEBODSYS)) # "other"

# An empty dataset whose variables do not decide the type is unknown
admiral:::get_admiral_df_type(adsl[0, c("STUDYID", "USUBJID", "AGE")]) # NA

# 5. Inferring keys (internal) -------------------------------------------------

admiral:::infer_admiral_keys(adsl) # USUBJID (STUDYID is constant, so dropped)
admiral:::infer_admiral_keys(advs) # USUBJID, PARAMCD + only what is needed
admiral:::infer_admiral_keys(adae) # USUBJID, AESEQ / ASEQ
admiral:::infer_admiral_keys(adtte) # USUBJID, PARAMCD

# minimal_unique_key(): add candidates until unique, then drop redundant ones
toy <- tribble(
  ~USUBJID, ~PARAMCD, ~AVISIT,    ~ATPT,     ~AVAL,
  "1",      "SYSBP",  "BASELINE", "MORNING",   121,
  "1",      "SYSBP",  "BASELINE", "EVENING",   125,
  "1",      "SYSBP",  "WEEK 2",   "MORNING",   118,
  "1",      "SYSBP",  "WEEK 2",   "EVENING",   120
)
admiral:::minimal_unique_key(toy, c("USUBJID", "PARAMCD"), c("AVISIT", "ATPT"))
admiral:::minimal_unique_key(
  filter(toy, ATPT == "MORNING"),
  c("USUBJID", "PARAMCD"),
  c("ATPT", "AVISIT") # ATPT is tried first but dropped as redundant
)

# 6. summary() ---------------------------------------------------------------

## ADSL: treatment arms and population flags
adsl_df <- set_admiral_keys(adsl, keys = exprs(USUBJID), dataset_name = "ADSL")
summary(adsl_df)

## BDS: parameters and visits, keys inferred (name taken from the variable).
## Even the most specific inferred key leaves duplicates in pharmaverseadam's
## ADVS, which the summary flags
advs_df <- as_admiral_df(advs)
summary(advs_df)

## BDS: checked against declared keys. These keys ignore ADT and DTYPE, so
## unscheduled visits and derived AVERAGE records show up as duplicates
summary(advs_keyed)

## OCCDS: distinct terms and treatment-emergent records
adae_df <- set_admiral_keys(adae, keys = exprs(USUBJID, AESEQ), dataset_name = "ADAE")
summary(adae_df)

## TTE: events and censored per parameter
adtte_df <- as_admiral_df(adtte)
summary(adtte_df)

## Typical use: check a derivation added what was expected
advs_small <- tribble(
  ~USUBJID,      ~PARAMCD, ~AVAL, ~AVISIT,
  "01-701-1015", "DIABP",     51, "BASELINE",
  "01-701-1015", "SYSBP",    121, "BASELINE",
  "01-701-1015", "DIABP",     55, "WEEK 2",
  "01-701-1015", "SYSBP",    125, "WEEK 2",
  "01-701-1028", "DIABP",     79, "BASELINE",
  "01-701-1028", "SYSBP",    130, "BASELINE"
) %>%
  mutate(STUDYID = "PILOT01")

advs_map <- derive_param_computed(
  advs_small,
  by_vars = exprs(STUDYID, USUBJID, AVISIT),
  parameters = c("SYSBP", "DIABP"),
  set_values_to = exprs(
    AVAL = (AVAL.SYSBP + 2 * AVAL.DIABP) / 3,
    PARAMCD = "MAP"
  )
)
summary(advs_map) # MAP added for 2 subjects, 3 records; one record per key

## Catching duplicates: a fanned-out merge duplicates subject records.
## The declared keys keep it typed as ADSL, so the arm/flag details remain;
## without declared keys the records decide, and it would be "other"
adsl_dup <- bind_rows(adsl, filter(adsl, row_number() <= 5)) %>%
  set_admiral_keys(keys = exprs(USUBJID), dataset_name = "ADSL")
summary(adsl_dup)
admiral:::get_admiral_df_type(bind_rows(adsl, filter(adsl, row_number() <= 5))) # "other"

## Declared keys which no longer exist (e.g. after rename()) are reported, not checked
advs_keyed %>%
  rename(ANLVISIT = AVISIT) %>%
  summary()

## Declared empty keys
summary(set_admiral_keys(adsl, keys = exprs(), dataset_name = "ADSL"))

## Work in progress: only items whose variables exist are reported
adsl_wip <- as_admiral_df(select(adsl, STUDYID, USUBJID, AGE, SEX))
summary(adsl_wip)

## Grouped datasets are summarized ungrouped
summary(as_admiral_df(group_by(adsl, TRT01P)))

## Unknown type: empty dataset
summary(as_admiral_df(adsl[0, c("STUDYID", "USUBJID", "AGE")]))

## The summary object is a list, so its pieces are accessible programmatically
s <- summary(advs_keyed)
s$structure
s$details$params
