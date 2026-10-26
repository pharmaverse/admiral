# as_admiral_df ----
## Test 1: adds the admiral_df class to a data frame ----
test_that("as_admiral_df Test 1: adds the admiral_df class to a data frame", {
  input <- tibble::tribble(
    ~USUBJID, ~AVAL,
    "1",      10,
    "2",      20
  )

  actual <- as_admiral_df(input)

  expect_equal(class(actual), c("admiral_df", class(input)))
  # Test that the dataset is not changed
  expect_equal(unclass(actual), unclass(input))
})

## Test 2: is idempotent when the class is already present ----
test_that("as_admiral_df Test 2: is idempotent when the class is already present", {
  input <- tibble::tribble(
    ~USUBJID, ~AVAL,
    "1",      10
  )

  once <- as_admiral_df(input)
  twice <- as_admiral_df(once)

  expect_equal(class(once), class(twice))
  expect_equal(sum(class(twice) == "admiral_df"), 1L)
})

## Test 3: returns NULL unchanged ----
test_that("as_admiral_df Test 3: returns NULL unchanged", {
  expect_null(as_admiral_df(NULL))
})

# set_admiral_keys ----
## Test 4: stores the keys, the dataset name, and the admiral_df class ----
test_that("set_admiral_keys Test 4: stores the keys, the dataset name, and the admiral_df class", { # nolint
  input <- tibble::tribble(
    ~USUBJID, ~PARAMCD, ~AVISIT,    ~AVAL,
    "1",      "DIABP",  "BASELINE",    51,
    "1",      "SYSBP",  "BASELINE",   121
  )

  actual <- set_admiral_keys(
    input,
    keys = exprs(USUBJID, PARAMCD, AVISIT),
    dataset_name = "ADVS"
  )

  expect_equal(
    attr(actual, "admiral_keys"),
    c("USUBJID", "PARAMCD", "AVISIT")
  )
  expect_equal(attr(actual, "admiral_ds_name"), "ADVS")
  expect_s3_class(actual, "admiral_df")
  # Test that the data itself is not changed
  stripped <- actual
  attr(stripped, "admiral_keys") <- NULL
  attr(stripped, "admiral_ds_name") <- NULL
  class(stripped) <- class(input)
  expect_equal(stripped, input)
})

## Test 5: keys from a specification are passed via chr2vars() ----
test_that("set_admiral_keys Test 5: keys from a specification are passed via chr2vars()", {
  input <- tibble::tribble(
    ~USUBJID, ~PARAMCD, ~AVAL,
    "1",      "DIABP",     51
  )
  # keys extracted from a specification table carry a variable label
  spec_keys <- c("USUBJID", "PARAMCD")
  attr(spec_keys, "label") <- "Variable Name"

  from_exprs <- set_admiral_keys(input, keys = exprs(USUBJID, PARAMCD))
  from_spec <- set_admiral_keys(input, keys = chr2vars(spec_keys))

  expect_identical(attr(from_spec, "admiral_keys"), attr(from_exprs, "admiral_keys"))
})

## Test 6: the dataset name is stored, preserved, and replaced ----
test_that("set_admiral_keys Test 6: the dataset name is stored, preserved, and replaced", {
  input <- tibble::tribble(
    ~USUBJID, ~PARAMCD, ~AVAL,
    "1",      "DIABP",     51
  )

  # nothing is stored when there is no name to store
  unnamed <- set_admiral_keys(input, keys = exprs(USUBJID))
  expect_equal(attr(unnamed, "admiral_keys"), "USUBJID")
  expect_null(attr(unnamed, "admiral_ds_name"))

  # `dataset_name` is not repeated, so the name of the dataset is not lost by
  # a call which only revises its keys
  named <- set_admiral_keys(input, keys = exprs(USUBJID), dataset_name = "ADVS")
  rekeyed <- set_admiral_keys(named, keys = exprs(USUBJID, PARAMCD))
  expect_equal(attr(rekeyed, "admiral_keys"), c("USUBJID", "PARAMCD"))
  expect_equal(attr(rekeyed, "admiral_ds_name"), "ADVS")

  # both the keys and the name give way when a later call supplies them
  renamed <- set_admiral_keys(named, keys = exprs(USUBJID, PARAMCD), dataset_name = "ADSL")
  expect_equal(attr(renamed, "admiral_keys"), c("USUBJID", "PARAMCD"))
  expect_equal(attr(renamed, "admiral_ds_name"), "ADSL")
})

## Test 7: the keys are stored as a plain character vector ----
test_that("set_admiral_keys Test 7: the keys are stored as a plain character vector", {
  input <- tibble::tribble(
    ~USUBJID, ~AVAL,
    "1",      10
  )

  # `vars2chr()` returns a named vector, which must not end up on the stored
  # attribute
  actual <- set_admiral_keys(input, keys = exprs(USUBJID, AVAL))

  expect_identical(attr(actual, "admiral_keys"), c("USUBJID", "AVAL"))
})

## Test 8: empty keys are stored, keys not in the dataset are not ----
test_that("set_admiral_keys Test 8: empty keys are stored, keys not in the dataset are not", { # nolint
  input <- tibble::tribble(
    ~USUBJID, ~PARAMCD, ~AVAL,
    "1",      "DIABP",     51
  )

  # the empty attribute is a statement that the dataset has no key variables,
  # so it must be distinguishable from an attribute which was never set
  expect_silent(empty <- set_admiral_keys(input, keys = exprs()))
  expect_identical(attr(empty, "admiral_keys"), character(0))
  expect_false(is.null(attr(empty, "admiral_keys")))
  expect_s3_class(empty, "admiral_df")

  # a key which is not in the dataset is most likely a typo, so the keys are
  # not set and a warning says so; the name and the class are still applied
  expect_snapshot(
    typo <- set_admiral_keys(input, keys = exprs(USUBJID, PARMACD), dataset_name = "ADVS")
  )
  expect_null(attr(typo, "admiral_keys"))
  expect_identical(attr(typo, "admiral_ds_name"), "ADVS")
  expect_s3_class(typo, "admiral_df")

  # re-keying with a typo keeps the keys set previously, and the warning says
  # which they are
  keyed <- set_admiral_keys(input, keys = exprs(USUBJID, PARAMCD))
  expect_snapshot(
    rekeyed <- set_admiral_keys(keyed, keys = exprs(USUBJID, PARMACD, AVISIT))
  )
  expect_identical(attr(rekeyed, "admiral_keys"), c("USUBJID", "PARAMCD"))
})

## Test 9: the keys and the class survive a dplyr pipeline ----
test_that("set_admiral_keys Test 9: the keys and the class survive a dplyr pipeline", {
  input <- set_admiral_keys(
    tibble::tribble(
      ~USUBJID, ~PARAMCD, ~AVAL,
      "1",      "DIABP",     51,
      "2",      "DIABP",     79
    ),
    keys = exprs(USUBJID, PARAMCD),
    dataset_name = "ADVS"
  )
  adsl <- tibble::tribble(
    ~USUBJID, ~TRT01P,
    "1",      "Placebo",
    "2",      "Xanomeline"
  )

  # the usefulness of the attribute rests on it surviving the verbs a derivation
  # pipeline is built from -- which is a promise about `{dplyr}`, not about
  # admiral, so it is worth pinning against a dependency update
  actual <- input %>%
    mutate(BASE = AVAL) %>%
    filter(AVAL > 0) %>%
    arrange(USUBJID) %>%
    select(USUBJID, PARAMCD, AVAL, BASE) %>%
    left_join(adsl, by = "USUBJID")

  expect_equal(attr(actual, "admiral_keys"), c("USUBJID", "PARAMCD"))
  expect_equal(attr(actual, "admiral_ds_name"), "ADVS")
  expect_s3_class(actual, "admiral_df")
})

## Test 10: an error is issued for invalid arguments ----
test_that("set_admiral_keys Test 10: an error is issued for invalid arguments", {
  input <- tibble::tribble(
    ~USUBJID, ~AVAL,
    "1",      10
  )

  expect_error(set_admiral_keys(list(), keys = exprs(USUBJID)))
  expect_error(set_admiral_keys(input, keys = 1))
  expect_error(set_admiral_keys(input, keys = iris))
  expect_error(set_admiral_keys(input, keys = exprs(USUBJID + 1)))
  expect_error(set_admiral_keys(input, keys = exprs(USUBJID), dataset_name = c("A", "B")))

  # keys taken from a specification are a character vector, so the error
  # points to the conversion
  expect_snapshot(set_admiral_keys(input, keys = c("USUBJID", "PARAMCD")), error = TRUE)
})

# get_admiral_df_type ----
## Test 11: classifies the common ADaM structures ----
test_that("get_admiral_df_type Test 11: classifies the common ADaM structures", {
  bds <- tibble::tribble(
    ~STUDYID, ~USUBJID, ~PARAMCD, ~AVAL, ~AVISIT,
    "P",      "1",      "MAP",       90, "BASELINE"
  )
  expect_identical(get_admiral_df_type(bds), "BDS")

  # BDS is satisfied by either analysis value variable
  expect_identical(
    get_admiral_df_type(tibble::tibble(STUDYID = "P", USUBJID = "1", PARAMCD = "A", AVALC = "X")),
    "BDS"
  )
  # with neither, it is not a BDS dataset
  expect_identical(
    get_admiral_df_type(tibble::tibble(
      STUDYID = "P", USUBJID = c("1", "1"), PARAMCD = "A", FOO = 1
    )),
    "other"
  )

  # two records per subject, so that dropping `PARAMCD` below does not leave a
  # dataset which is subject-level by structure
  tte <- tibble::tribble(
    ~STUDYID, ~USUBJID, ~PARAMCD, ~AVAL, ~CNSR, ~STARTDT,
    "P",      "1",      "OS",       100,     0, as.Date("2020-01-01"),
    "P",      "1",      "PFS",       60,     1, as.Date("2020-01-01")
  )
  expect_identical(get_admiral_df_type(tte), "TTE")

  # the TTE signal needs all three of its variables: dropping any one term falls
  # back to BDS rather than staying TTE
  expect_identical(get_admiral_df_type(select(tte, -CNSR)), "BDS")
  expect_identical(get_admiral_df_type(select(tte, -STARTDT)), "BDS")
  expect_identical(get_admiral_df_type(select(tte, -PARAMCD)), "other")

  adsl <- tibble::tribble(
    ~STUDYID, ~USUBJID, ~TRT01P,
    "P",      "1",      "A",
    "P",      "2",      "B"
  )
  expect_identical(get_admiral_df_type(adsl), "ADSL")

  other <- tibble::tribble(
    ~FOO, ~BAR,
    1,    2
  )
  expect_identical(get_admiral_df_type(other), "other")

  # a summary of a BDS dataset shares its `PARAMCD` and `AVAL` but has no
  # subjects, so it is not a BDS dataset itself
  bds_summary <- tibble::tribble(
    ~TRT01A, ~PARAMCD, ~AVAL,
    "A",     "SYSBP",    125,
    "B",     "SYSBP",    125
  )
  expect_identical(get_admiral_df_type(bds_summary), "other")
})

## Test 12: each occurrence signal is recognized on its own ----
test_that("get_admiral_df_type Test 12: each occurrence signal is recognized on its own", {
  # two records for one subject, so that a dataset which fails every occurrence
  # test is not subject-level either and the fallback is visible as "other"
  occds <- function(...) {
    tibble::tibble(STUDYID = "P", USUBJID = c("1", "1"), ...)
  }

  # the occurrence flag family of the ADaM OCCDS implementation guide, plus the
  # sponsor-numbered form
  for (flag in c("AOCCFL", "AOCCIFL", "AOCCSFL", "AOCCPFL", "AOCCPIFL", "AOCC02FL")) {
    input <- occds()
    input[[flag]] <- c("Y", "N")
    expect_identical(get_admiral_df_type(input), "OCCDS", info = flag)
  }

  expect_identical(get_admiral_df_type(occds(TRTEMFL = c("Y", "N"))), "OCCDS")
  expect_identical(get_admiral_df_type(occds(AEDECOD = c("A", "B"))), "OCCDS")
  expect_identical(get_admiral_df_type(occds(MHTERM = c("A", "B"))), "OCCDS")
  expect_identical(get_admiral_df_type(occds(CMTERM = c("A", "B"))), "OCCDS")
})

## Test 13: an occurrence signal which is not corroborated does not classify ----
test_that("get_admiral_df_type Test 13: an occurrence signal which is not corroborated does not classify", { # nolint
  # `LONGTERM` ends in TERM but is not a `--TERM` variable, and the two-letter
  # SDTM domain prefix is what distinguishes the two
  not_occds <- tibble::tribble(
    ~STUDYID, ~USUBJID, ~LONGTERM, ~AVAL,
    "P",      "1",      "Y",          10,
    "P",      "1",      "N",          20
  )
  expect_identical(get_admiral_df_type(not_occds), "other")

  # the naming signal alone, with nothing identifying a subject, is not enough
  no_subject <- tibble::tribble(
    ~AEDECOD,   ~AESEV,
    "HEADACHE", "MILD"
  )
  expect_identical(get_admiral_df_type(no_subject), "other")

  # nor is a subset of the subject keys: `STUDYID` is constant within a study,
  # so a lookup table carrying it still identifies no subject
  study_only <- tibble::tribble(
    ~STUDYID, ~AEDECOD,
    "P",      "HEADACHE",
    "P",      "NAUSEA"
  )
  expect_identical(get_admiral_df_type(study_only), "other")

  # a lone severity variable is not an occurrence signal either
  severity_only <- tibble::tibble(
    STUDYID = "P", USUBJID = c("1", "1"), AESEV = c("MILD", "SEVERE")
  )
  expect_identical(get_admiral_df_type(severity_only), "other")
})

## Test 14: the type precedence resolves datasets matching more than one type ----
test_that("get_admiral_df_type Test 14: the type precedence resolves datasets matching more than one type", { # nolint
  # BDS over OCCDS: an occurrence dataset which has acquired a parameter is
  # reported per parameter, so `PARAMCD` excludes the OCCDS branch outright
  bds_and_occds <- tibble::tribble(
    ~STUDYID, ~USUBJID, ~PARAMCD, ~AVAL, ~AEDECOD,   ~TRTEMFL,
    "P",      "1",      "AEHLT",      1, "HEADACHE", "Y"
  )
  expect_identical(get_admiral_df_type(bds_and_occds), "BDS")

  # OCCDS over ADSL: one adverse event per subject is still an occurrence
  # dataset, even though it is trivially one record per subject and carries a
  # treatment variable merged from ADSL
  occds_and_adsl <- tibble::tribble(
    ~STUDYID, ~USUBJID, ~AEDECOD,   ~TRT01P,   ~TRTSDT,
    "P",      "1",      "HEADACHE", "Placebo", as.Date("2020-01-01"),
    "P",      "2",      "NAUSEA",   "Drug",    as.Date("2020-01-02")
  )
  expect_identical(get_admiral_df_type(occds_and_adsl), "OCCDS")
})

## Test 15: an empty dataset is typed from the variables it declares ----
test_that("get_admiral_df_type Test 15: an empty dataset is typed from the variables it declares", { # nolint
  empty_bds <- tibble::tibble(
    STUDYID = character(0), USUBJID = character(0), PARAMCD = character(0),
    AVAL = numeric(0)
  )
  expect_identical(get_admiral_df_type(empty_bds), "BDS")

  # ADSL is decided by the records, and an empty dataset cannot show whether it
  # has one record per subject -- the record and subject counts are trivially
  # equal -- so the type is unknown, whatever treatment variables it carries
  empty_adsl <- tibble::tibble(
    STUDYID = character(0), USUBJID = character(0), TRT01P = character(0)
  )
  expect_identical(get_admiral_df_type(empty_adsl), NA_character_)

  empty_shell <- tibble::tibble(
    STUDYID = character(0), USUBJID = character(0), AVAL = numeric(0)
  )
  expect_identical(get_admiral_df_type(empty_shell), NA_character_)

  # where the variables decide the type, the missing records do not matter:
  # without the subject keys it is not subject data at all
  empty_lookup <- tibble::tibble(AEDECOD = character(0), AETERM = character(0))
  expect_identical(get_admiral_df_type(empty_lookup), "other")

  # nothing is inferred for an empty dataset, whatever its type
  expect_identical(infer_admiral_keys(empty_adsl), character(0))
})

## Test 16: an error is issued for a non-data-frame dataset ----
test_that("get_admiral_df_type Test 16: an error is issued for a non-data-frame dataset", {
  expect_error(get_admiral_df_type(matrix(1:4, nrow = 2)))
  expect_error(get_admiral_df_type("ADSL"))
})

## Test 17: ADSL is decided by the records ----
test_that("get_admiral_df_type Test 17: ADSL is decided by the records", {
  one_per_subject <- tibble::tribble(
    ~STUDYID, ~USUBJID, ~AGE,
    "P",      "1",        63,
    "P",      "2",        71
  )
  expect_identical(get_admiral_df_type(one_per_subject), "ADSL")

  many_per_subject <- tibble::tribble(
    ~STUDYID, ~USUBJID, ~AVAL,
    "P",      "1",         10,
    "P",      "1",         20
  )
  expect_identical(get_admiral_df_type(many_per_subject), "other")

  # all the subject keys are needed: `STUDYID` is constant within a study, so a
  # single-study lookup table carrying it is one record per study, not subject
  study_only <- tibble::tribble(
    ~STUDYID, ~PARAM,
    "P",      "Code list A"
  )
  expect_identical(get_admiral_df_type(study_only), "other")
})

## Test 18: treatment variables are no evidence of ADSL ----
test_that("get_admiral_df_type Test 18: treatment variables are no evidence of ADSL", {
  # the shape every BDS/OCCDS derivation passes through: the ADSL variables have
  # been merged on, but `PARAMCD` has not been assigned yet. It carries every
  # treatment variable ADSL does while having many records per subject, so
  # nothing but the record structure distinguishes it
  mid_derivation <- tibble::tribble(
    ~STUDYID, ~USUBJID, ~VSTESTCD, ~VSSTRESN, ~TRT01P,   ~TRT01A,   ~TRTSDT,
    "P",      "1",      "SYSBP",         121, "Placebo", "Placebo", as.Date("2020-01-01"),
    "P",      "1",      "DIABP",          79, "Placebo", "Placebo", as.Date("2020-01-01"),
    "P",      "2",      "SYSBP",         130, "Drug",    "Drug",    as.Date("2020-01-02"),
    "P",      "2",      "DIABP",          85, "Drug",    "Drug",    as.Date("2020-01-02")
  )

  expect_identical(get_admiral_df_type(mid_derivation), "other")
  # the period-numbered form is equally no evidence
  expect_identical(get_admiral_df_type(mutate(mid_derivation, TRT02A = "Drug")), "other")

  # the finished ADSL, with the same variables and one record per subject, is
  expect_identical(
    get_admiral_df_type(distinct(mid_derivation, STUDYID, USUBJID, TRT01P, TRT01A, TRTSDT)),
    "ADSL"
  )
})

## Test 19: a grouped dataset is judged on its records overall ----
test_that("get_admiral_df_type Test 19: a grouped dataset is judged on its records overall", {
  many_per_subject <- tibble::tribble(
    ~STUDYID, ~USUBJID, ~AVISIT,    ~AVAL,
    "P",      "1",      "BASELINE",    10,
    "P",      "1",      "WEEK 2",      20
  )

  # `distinct()` on a grouped dataset silently adds the grouping variables, so
  # without `ungroup()` this would test one record per subject *per visit* and
  # report a findings dataset as subject-level
  expect_identical(get_admiral_df_type(group_by(many_per_subject, AVISIT)), "other")
})

## Test 20: the subject keys option is respected ----
test_that("get_admiral_df_type Test 20: the subject keys option is respected", {
  # one record per STUDYID + USUBJID (the default subject keys), but `SUBJID`
  # repeats across the two studies
  input <- tibble::tribble(
    ~STUDYID, ~USUBJID, ~SUBJID, ~AGE,
    "A",      "A-1",    "1",       63,
    "B",      "B-1",    "1",       71
  )

  expect_identical(get_admiral_df_type(input), "ADSL")

  subject_keys <- get_admiral_option("subject_keys")
  withr::defer(set_admiral_options(subject_keys = subject_keys))
  set_admiral_options(subject_keys = exprs(SUBJID))

  # what counts as one record per subject follows the option, so the same
  # dataset is no longer subject-level
  expect_identical(get_admiral_df_type(input), "other")
})

# minimal_unique_key ----
## Test 21: returns must_have when it is already unique ----
test_that("minimal_unique_key Test 21: returns must_have when it is already unique", {
  input <- tibble::tribble(
    ~USUBJID, ~PARAMCD, ~AVISIT,
    "1",      "DIABP",  "BASELINE",
    "2",      "DIABP",  "BASELINE"
  )

  expect_identical(
    minimal_unique_key(input, must_have = "USUBJID", optional = "AVISIT"),
    "USUBJID"
  )
})

## Test 22: drops candidates which do not discriminate ----
test_that("minimal_unique_key Test 22: drops candidates which do not discriminate", {
  # `AVISIT` is redundant with `AVISITN`, and `ATPT` is constant; both sit ahead
  # of `DTYPE` in the candidate order, so the walk carries them on its way to the
  # variable which actually separates the records
  input <- tibble::tribble(
    ~USUBJID, ~AVISITN, ~AVISIT,  ~ATPT, ~DTYPE,
    "1",             2, "WEEK 2", "PRE", NA,
    "1",             2, "WEEK 2", "PRE", "LOCF",
    "1",             4, "WEEK 4", "PRE", NA,
    "1",             4, "WEEK 4", "PRE", "LOCF"
  )

  expect_identical(
    minimal_unique_key(
      input,
      must_have = "USUBJID",
      optional = c("AVISITN", "AVISIT", "ATPT", "DTYPE")
    ),
    c("USUBJID", "AVISITN", "DTYPE")
  )
})

## Test 23: returns everything when uniqueness is never reached ----
test_that("minimal_unique_key Test 23: returns everything when uniqueness is never reached", {
  # wholly duplicated records: no combination of the candidates separates them,
  # so the caller sees duplicates against the returned key and can report them
  input <- tibble::tribble(
    ~USUBJID, ~PARAMCD, ~AVISIT,
    "1",      "DIABP",  "BASELINE",
    "1",      "DIABP",  "BASELINE"
  )

  expect_identical(
    minimal_unique_key(input, must_have = "USUBJID", optional = c("PARAMCD", "AVISIT")),
    c("USUBJID", "PARAMCD", "AVISIT")
  )
})

# infer_admiral_keys ----
## Test 24: the inferred key spans the common BDS shapes ----
test_that("infer_admiral_keys Test 24: the inferred key spans the common BDS shapes", {
  # exposure: one record per subject, parameter and dosing interval, keyed by
  # the interval start rather than by an analysis date or a visit
  adex <- tibble::tribble(
    ~STUDYID, ~USUBJID, ~PARAMCD, ~AVAL, ~ASTDT,                ~AENDT,
    "P",      "1",      "DOSE",      54, as.Date("2024-01-01"), as.Date("2024-01-31"),
    "P",      "1",      "DOSE",      81, as.Date("2024-02-01"), as.Date("2024-02-28")
  )
  expect_identical(infer_admiral_keys(adex), c("USUBJID", "PARAMCD", "ASTDT"))

  # population PK: keyed by relative time, with no visit structure at all
  adppk <- tibble::tribble(
    ~STUDYID, ~USUBJID, ~PARAMCD, ~NFRLT, ~AVAL,
    "P",      "1",      "CONC",        0,     0,
    "P",      "1",      "CONC",        1,    12,
    "P",      "1",      "CONC",        2,     8
  )
  expect_identical(infer_admiral_keys(adppk), c("USUBJID", "PARAMCD", "NFRLT"))

  # a derived record shares every analysis variable with the record it was
  # derived from, so only DTYPE separates the two within a visit
  locf <- tibble::tribble(
    ~STUDYID, ~USUBJID, ~PARAMCD, ~AVISITN, ~AVAL, ~DTYPE,
    "P",      "1",      "SYSBP",         2,   121, NA,
    "P",      "1",      "SYSBP",         2,   121, "LOCF",
    "P",      "1",      "SYSBP",         4,   130, NA,
    "P",      "1",      "SYSBP",         4,   130, "LOCF"
  )
  expect_identical(
    infer_admiral_keys(locf),
    c("USUBJID", "PARAMCD", "AVISITN", "DTYPE")
  )
})

## Test 25: multi-period designs are keyed by period ----
test_that("infer_admiral_keys Test 25: multi-period designs are keyed by period", {
  # the same visit recurs in each period, so nothing but APERIOD separates the
  # records -- a vaccine study shape
  input <- tibble::tribble(
    ~STUDYID, ~USUBJID, ~PARAMCD, ~AVISITN, ~APERIOD, ~AVAL,
    "P",      "1",      "SYSBP",         1,        1,   121,
    "P",      "1",      "SYSBP",         1,        2,   130,
    "P",      "1",      "SYSBP",         2,        1,   118,
    "P",      "1",      "SYSBP",         2,        2,   125
  )

  expect_identical(
    infer_admiral_keys(input),
    c("USUBJID", "PARAMCD", "AVISITN", "APERIOD")
  )
})

## Test 26: ADSL, TTE and unrecognized datasets ----
test_that("infer_admiral_keys Test 26: ADSL, TTE and unrecognized datasets", {
  adsl <- tibble::tribble(
    ~STUDYID, ~USUBJID, ~TRT01P,
    "P",      "1",      "Placebo",
    "P",      "2",      "Drug"
  )
  expect_identical(infer_admiral_keys(adsl), "USUBJID")

  adtte <- tibble::tribble(
    ~STUDYID, ~USUBJID, ~PARAMCD, ~AVAL, ~CNSR, ~STARTDT,
    "P",      "1",      "OS",       100,     0, as.Date("2020-01-01"),
    "P",      "1",      "PFS",       60,     1, as.Date("2020-01-01")
  )
  expect_identical(infer_admiral_keys(adtte), c("USUBJID", "PARAMCD"))

  # nothing is inferred for a dataset whose type could not be determined
  other <- tibble::tribble(
    ~FOO, ~BAR,
    1,    2,
    1,    3
  )
  expect_identical(infer_admiral_keys(other), character(0))
})

## Test 27: OCCDS is keyed by its sequence variable ----
test_that("infer_admiral_keys Test 27: OCCDS is keyed by its sequence variable", {
  occds <- function(...) {
    tibble::tibble(
      STUDYID = "P", USUBJID = c("1", "1"),
      AEDECOD = c("HEADACHE", "NAUSEA"), ...
    )
  }

  # `ASEQ` is preferred when present and populated
  expect_identical(
    infer_admiral_keys(occds(ASEQ = 1:2, AESEQ = 3:4)),
    c("USUBJID", "ASEQ")
  )

  # occurrence datasets most often carry only the domain sequence
  expect_identical(infer_admiral_keys(occds(AESEQ = 1:2)), c("USUBJID", "AESEQ"))
  expect_identical(infer_admiral_keys(occds(CMSEQ = 1:2)), c("USUBJID", "CMSEQ"))

  # `SRCSEQ` is provenance added by a merge, not a record key: three letters
  # before SEQ, so it must not be mistaken for a domain sequence
  expect_warning(
    result <- infer_admiral_keys(occds(SRCSEQ = 1:2)),
    regexp = "domain sequence"
  )
  expect_identical(result, character(0))

  # several domain sequences make the record key ambiguous rather than obvious
  expect_warning(
    ambiguous <- infer_admiral_keys(occds(AESEQ = 1:2, CESEQ = 1:2)),
    regexp = "ambiguous"
  )
  expect_identical(ambiguous, character(0))

  # an all-NA sequence has not been derived yet, so `ASEQ` gives way to the
  # populated domain sequence rather than keying the dataset on nothing
  expect_identical(
    infer_admiral_keys(occds(ASEQ = NA_integer_, AESEQ = 1:2)),
    c("USUBJID", "AESEQ")
  )
  expect_warning(
    empty_seq <- infer_admiral_keys(occds(ASEQ = NA_integer_)),
    regexp = "populated"
  )
  expect_identical(empty_seq, character(0))
})

## Test 28: nothing is inferred for a dataset with no records ----
test_that("infer_admiral_keys Test 28: nothing is inferred for a dataset with no records", {
  # every candidate key is trivially unique over no records, so inference would
  # report a structure the dataset has not demonstrated
  empty_bds <- tibble::tibble(
    STUDYID = character(0), USUBJID = character(0), PARAMCD = character(0),
    AVAL = numeric(0)
  )

  expect_identical(infer_admiral_keys(empty_bds), character(0))
})

## Test 29: the subject keys option is respected ----
test_that("infer_admiral_keys Test 29: the subject keys option is respected", {
  subject_keys <- get_admiral_option("subject_keys")
  withr::defer(set_admiral_options(subject_keys = subject_keys))

  # one study, so `STUDYID` distinguishes nothing and is left out of the
  # reported record structure
  one_study <- tibble::tribble(
    ~STUDYID, ~USUBJID, ~SUBJID, ~PARAMCD, ~AVAL,
    "A",      "A-1",    "1",     "SYSBP",    121,
    "A",      "A-2",    "2",     "SYSBP",    130
  )
  expect_identical(infer_admiral_keys(one_study), c("USUBJID", "PARAMCD"))

  # the subject part of the key follows the option, as the ADSL structure check
  # does
  set_admiral_options(subject_keys = exprs(SUBJID))
  expect_identical(infer_admiral_keys(one_study), c("SUBJID", "PARAMCD"))

  # `SUBJID` repeats across the two studies and neither variable identifies a
  # subject on its own, so both are doing real work here and both are kept
  pooled <- tibble::tribble(
    ~STUDYID, ~SUBJID, ~PARAMCD, ~AVAL,
    "A",      "1",     "SYSBP",    121,
    "A",      "2",     "SYSBP",    130,
    "B",      "1",     "SYSBP",    118,
    "B",      "2",     "SYSBP",    125
  )
  set_admiral_options(subject_keys = exprs(STUDYID, SUBJID))
  expect_identical(
    infer_admiral_keys(pooled),
    c("STUDYID", "SUBJID", "PARAMCD")
  )
})

## Test 30: a grouped dataset is keyed on its records overall ----
test_that("infer_admiral_keys Test 30: a grouped dataset is keyed on its records overall", {
  input <- tibble::tribble(
    ~STUDYID, ~USUBJID, ~PARAMCD, ~AVISITN, ~AVAL,
    "P",      "1",      "SYSBP",         2,   121,
    "P",      "1",      "SYSBP",         4,   130
  )

  # without `ungroup()` the uniqueness test would run within group, so the key
  # would stop at USUBJID + PARAMCD and the visit structure would be lost
  expect_identical(
    infer_admiral_keys(group_by(input, AVISITN)),
    c("USUBJID", "PARAMCD", "AVISITN")
  )
})

## Test 31: population PK needs both the nominal and the actual relative time ----
test_that("infer_admiral_keys Test 31: population PK needs both the nominal and the actual relative time", { # nolint
  # a draw taken off schedule shares its nominal time with the scheduled draw
  # and its actual time with a different nominal one, so neither variable keys
  # the dataset alone. Reproduces the record structure of `pharmaverseadam::adppk`.
  input <- tibble::tribble(
    ~STUDYID, ~USUBJID, ~PARAMCD, ~NFRLT, ~AFRLT, ~AVAL,
    "P",      "1",      "XAN",         3,      6, 25.77,
    "P",      "1",      "XAN",         6,      6,  1.86,
    "P",      "1",      "XAN",        36,     36,  0.31,
    "P",      "1",      "XAN",        36,     37,  0.27
  )

  expect_identical(
    infer_admiral_keys(input),
    c("USUBJID", "PARAMCD", "NFRLT", "AFRLT")
  )
})

## Test 32: visit, timepoint and DTYPE are all needed at once ----
test_that("infer_admiral_keys Test 32: visit, timepoint and DTYPE are all needed at once", { # nolint
  # visit, timepoint *and* `DTYPE` are all needed, while the relative time and
  # the analysis date carry no information the first two do not -- so uniqueness
  # is only reached at the end of the candidate list and the redundant
  # candidates picked up on the way have to be dropped again. Reproduces the
  # record structure of `pharmaverseadam::adpc`.
  input <- tibble::tribble(
    ~STUDYID, ~USUBJID, ~PARAMCD, ~AVISITN, ~ATPTN, ~NFRLT, ~ADT, ~DTYPE, ~AVAL,
    "P", "1", "XAN", 2, 24, 24, as.Date("2024-01-02"), NA, 0.0107,
    "P", "1", "XAN", 2, 36, 36, as.Date("2024-01-02"), NA, NA,
    "P", "1", "XAN", 2, 36, 36, as.Date("2024-01-02"), "HALFLLOQ", 0.005,
    "P", "1", "XAN", 3, 24, 24, as.Date("2024-01-09"), NA, 0.0098
  )

  expect_identical(
    infer_admiral_keys(input),
    c("USUBJID", "PARAMCD", "AVISITN", "ATPTN", "DTYPE")
  )
})

# print.admiral_df ----
## Test 33: the name, type and keys are printed in a header ----
test_that("print.admiral_df Test 33: the name, type and keys are printed in a header", {
  input <- tibble::tribble(
    ~STUDYID, ~USUBJID, ~PARAMCD, ~AVISIT,    ~AVAL,
    "P",      "1",      "DIABP",  "BASELINE",    51,
    "P",      "1",      "SYSBP",  "BASELINE",   121
  )

  expect_snapshot(
    print(set_admiral_keys(input, keys = exprs(USUBJID, PARAMCD, AVISIT), dataset_name = "ADVS"))
  )

  # either attribute is shown on its own
  expect_snapshot(print(set_admiral_keys(input, keys = exprs(USUBJID, PARAMCD, AVISIT))))

  # empty keys record that no keys are defined, so they are shown as such
  expect_snapshot(print(set_admiral_keys(input, keys = exprs(), dataset_name = "ADVS")))

  # a plain data frame prints with its own method below the header
  expect_snapshot(
    print(set_admiral_keys(as.data.frame(input), keys = exprs(USUBJID), dataset_name = "ADVS"))
  )

  # a type which cannot be determined is left out rather than shown as
  # "other": without the subject keys this is not subject data
  lookup <- tibble::tribble(
    ~PARAMCD, ~PARAM,
    "DIABP",  "Diastolic Blood Pressure (mmHg)"
  )
  expect_snapshot(print(set_admiral_keys(lookup, keys = exprs(PARAMCD), dataset_name = "PARAMS")))
})

## Test 34: a dataset with only the class prints as before ----
test_that("print.admiral_df Test 34: a dataset with only the class prints as before", {
  input <- tibble::tribble(
    ~USUBJID, ~AVAL,
    "1",      10
  )

  # every derivation tags its output with the class, so a header without
  # anything to show would change how every admiral dataset prints
  expect_identical(
    capture.output(print(as_admiral_df(input))),
    capture.output(print(input))
  )
  expect_invisible(print(as_admiral_df(input)))
})

# summary.admiral_df ----
## Test 35: the common facts and the declared record structure are reported ----
test_that("summary.admiral_df Test 35: the common facts and the declared record structure are reported", { # nolint
  input <- tibble::tribble(
    ~STUDYID, ~USUBJID, ~PARAMCD, ~AVISIT,    ~AVAL,
    "P",      "1",      "DIABP",  "BASELINE",    51,
    "P",      "1",      "SYSBP",  "BASELINE",   121,
    "P",      "2",      "DIABP",  "BASELINE",    79
  ) %>%
    set_admiral_keys(keys = exprs(USUBJID, PARAMCD, AVISIT), dataset_name = "ADVS")

  actual <- summary(input)

  expect_s3_class(actual, "summary_admiral_df")
  expect_identical(actual$name, "ADVS")
  expect_identical(actual$type, "BDS")
  expect_identical(actual$n_subjects, 2L)
  expect_identical(actual$n_obs, 3L)
  expect_identical(actual$n_vars, 5L)
  expect_identical(
    actual$structure,
    list(
      keys = c("USUBJID", "PARAMCD", "AVISIT"),
      source = "declared",
      missing_keys = character(0),
      n_duplicates = 0L
    )
  )
})

## Test 36: the keys are inferred if none were declared ----
test_that("summary.admiral_df Test 36: the keys are inferred if none were declared", {
  input <- tibble::tribble(
    ~STUDYID, ~USUBJID, ~PARAMCD, ~AVISIT,    ~AVAL,
    "P",      "1",      "SYSBP",  "BASELINE",   121,
    "P",      "1",      "SYSBP",  "WEEK 2",     118,
    "P",      "2",      "SYSBP",  "BASELINE",   130
  ) %>%
    as_admiral_df()

  actual <- summary(input)

  expect_identical(actual$name, "input")
  expect_identical(actual$structure$keys, c("USUBJID", "PARAMCD", "AVISIT"))
  expect_identical(actual$structure$source, "inferred")
  expect_identical(actual$structure$n_duplicates, 0L)
})

## Test 37: duplicate, missing and empty declared keys ----
test_that("summary.admiral_df Test 37: duplicate, missing and empty declared keys", {
  input <- tibble::tribble(
    ~STUDYID, ~USUBJID, ~PARAMCD, ~AVISIT,    ~AVAL,
    "P",      "1",      "SYSBP",  "BASELINE",   121,
    "P",      "1",      "SYSBP",  "WEEK 2",     118,
    "P",      "2",      "SYSBP",  "BASELINE",   130
  )

  duplicated <- summary(set_admiral_keys(input, keys = exprs(USUBJID, PARAMCD)))
  expect_identical(duplicated$structure$n_duplicates, 1L)

  # a renamed key cannot be checked, rather than being checked on part of the key
  renamed <- set_admiral_keys(input, keys = exprs(USUBJID, PARAMCD, AVISIT)) %>%
    dplyr::rename(VISIT = AVISIT)
  expect_identical(summary(renamed)$structure$missing_keys, "AVISIT")
  expect_identical(summary(renamed)$structure$n_duplicates, NA_integer_)

  no_keys <- summary(set_admiral_keys(input, keys = exprs()))
  expect_identical(no_keys$structure$keys, character(0))
  expect_identical(no_keys$structure$source, "declared")
})

## Test 38: the ADSL summary ----
test_that("summary.admiral_df Test 38: the ADSL summary", {
  input <- tibble::tribble(
    ~STUDYID, ~USUBJID, ~TRT01P,   ~SAFFL, ~ITTFL,
    "P",      "1",      "Placebo", "Y",    "Y",
    "P",      "2",      "Drug A",  "Y",    "Y",
    "P",      "3",      "Drug A",  "N",    "Y"
  ) %>%
    set_admiral_keys(keys = exprs(USUBJID), dataset_name = "ADSL")

  actual <- summary(input)

  expect_identical(actual$details$trt_var, "TRT01P")
  expect_identical(actual$details$trt$subjects, c(2L, 1L))
  expect_identical(actual$details$pops, c(SAFFL = 2L, ITTFL = 3L))
  expect_snapshot(print(actual))
})

## Test 39: the BDS summary ----
test_that("summary.admiral_df Test 39: the BDS summary", {
  input <- tibble::tribble(
    ~STUDYID, ~USUBJID, ~PARAMCD, ~AVISIT,    ~AVISITN, ~AVAL,
    "P",      "1",      "SYSBP",  "WEEK 2",          2,   118,
    "P",      "1",      "SYSBP",  "BASELINE",        0,   121,
    "P",      "1",      "DIABP",  "BASELINE",        0,    79,
    "P",      "2",      "SYSBP",  "BASELINE",        0,   130
  ) %>%
    as_admiral_df()

  actual <- summary(input)

  expect_identical(actual$details$params$PARAMCD, c("DIABP", "SYSBP"))
  expect_identical(actual$details$params$records, c(1L, 3L))
  expect_identical(actual$details$params$subjects, c(1L, 2L))
  # visits are in `AVISITN` order
  expect_identical(actual$details$visits, c("BASELINE", "WEEK 2"))
  expect_snapshot(print(actual))
})

## Test 40: the OCCDS summary ----
test_that("summary.admiral_df Test 40: the OCCDS summary", {
  input <- tibble::tribble(
    ~STUDYID, ~USUBJID, ~AESEQ, ~AEBODSYS,                ~AEDECOD,   ~TRTEMFL,
    "P",      "1",           1, "NERVOUS SYSTEM",         "HEADACHE", "Y",
    "P",      "1",           2, "GASTROINTESTINAL",       "NAUSEA",   NA,
    "P",      "2",           1, "NERVOUS SYSTEM",         "HEADACHE", "Y"
  ) %>%
    as_admiral_df()

  actual <- summary(input)

  expect_identical(actual$details$terms, c(AEBODSYS = 2L, AEDECOD = 2L))
  expect_identical(actual$details$trtem, c(emergent = 2L, total = 3L))
  expect_snapshot(print(actual))
})

## Test 41: the TTE summary ----
test_that("summary.admiral_df Test 41: the TTE summary", {
  input <- tibble::tribble(
    ~STUDYID, ~USUBJID, ~PARAMCD, ~STARTDT, ~CNSR, ~AVAL,
    "P", "1", "OS", lubridate::ymd("2020-01-01"), 0, 10,
    "P", "2", "OS", lubridate::ymd("2020-01-01"), 1, 12,
    "P", "1", "PFS", lubridate::ymd("2020-01-01"), 1, 3
  ) %>%
    as_admiral_df()

  actual <- summary(input)

  expect_identical(actual$details$events$PARAMCD, c("OS", "PFS"))
  expect_identical(actual$details$events$events, c(1L, 0L))
  expect_identical(actual$details$events$censored, c(1L, 1L))
  expect_snapshot(print(actual))
})

## Test 42: datasets with little to report do not error ----
test_that("summary.admiral_df Test 42: datasets with little to report do not error", {
  # not subject data: no subject count and no record structure
  lookup <- as_admiral_df(tibble::tibble(PARAMCD = c("DIABP", "SYSBP")))
  expect_snapshot(print(summary(lookup)))

  # no records: the type cannot be determined
  empty <- as_admiral_df(tibble::tibble(STUDYID = character(), USUBJID = character()))
  expect_identical(summary(empty)$type, NA_character_)
  expect_identical(summary(empty)$details, list())
  expect_snapshot(print(summary(empty)))

  # a work-in-progress ADSL without treatment or population flags
  adsl <- as_admiral_df(tibble::tibble(STUDYID = "P", USUBJID = c("1", "2"), TRT01P = NA))
  expect_identical(summary(adsl)$details, list())
})

## Test 43: print.summary_admiral_df returns its input invisibly ----
test_that("print.summary_admiral_df Test 43: returns its input invisibly", {
  input <- as_admiral_df(tibble::tibble(STUDYID = "P", USUBJID = "1"))
  actual <- summary(input)

  expect_output(res <- withVisible(print(actual)), "Summary of an admiral dataset")
  expect_false(res$visible)
  expect_identical(res$value, actual)
})
