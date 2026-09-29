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

## Test 4: returns NULL unchanged ----
test_that("as_admiral_df Test 4: returns NULL unchanged", {
  expect_null(as_admiral_df(NULL))
})

# set_admiral_keys ----
## Test 5: stores the keys, the dataset name, and the admiral_df class ----
test_that("set_admiral_keys Test 5: stores the keys, the dataset name, and the admiral_df class", { # nolint
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

## Test 6: keys are equally accepted as a character vector ----
test_that("set_admiral_keys Test 6: keys are equally accepted as a character vector", {
  input <- tibble::tribble(
    ~USUBJID, ~PARAMCD, ~AVAL,
    "1",      "DIABP",     51
  )

  from_exprs <- set_admiral_keys(input, keys = exprs(USUBJID, PARAMCD))
  from_chr <- set_admiral_keys(input, keys = c("USUBJID", "PARAMCD"))

  expect_equal(attr(from_chr, "admiral_keys"), c("USUBJID", "PARAMCD"))
  expect_equal(attr(from_exprs, "admiral_keys"), attr(from_chr, "admiral_keys"))
})

## Test 7: a dataset name is only stored when there is one to store ----
test_that("set_admiral_keys Test 7: a dataset name is only stored when there is one to store", { # nolint
  input <- tibble::tribble(
    ~USUBJID, ~AVAL,
    "1",      10
  )

  actual <- set_admiral_keys(input, keys = "USUBJID")

  expect_equal(attr(actual, "admiral_keys"), "USUBJID")
  expect_null(attr(actual, "admiral_ds_name"))
})

## Test 8: a stored dataset name survives re-keying ----
test_that("set_admiral_keys Test 8: a stored dataset name survives re-keying", {
  input <- tibble::tribble(
    ~USUBJID, ~PARAMCD, ~AVAL,
    "1",      "DIABP",     51
  )

  named <- set_admiral_keys(input, keys = "USUBJID", dataset_name = "ADVS")
  # `dataset_name` is not repeated, so the name of the dataset is not lost by
  # a call which only revises its keys
  rekeyed <- set_admiral_keys(named, keys = c("USUBJID", "PARAMCD"))

  expect_equal(attr(rekeyed, "admiral_keys"), c("USUBJID", "PARAMCD"))
  expect_equal(attr(rekeyed, "admiral_ds_name"), "ADVS")
})

## Test 9: keys replace those of a previous call ----
test_that("set_admiral_keys Test 9: keys replace those of a previous call", {
  input <- tibble::tribble(
    ~USUBJID, ~PARAMCD, ~AVAL,
    "1",      "DIABP",     51
  )

  once <- set_admiral_keys(input, keys = "USUBJID", dataset_name = "ADSL")
  twice <- set_admiral_keys(once, keys = c("USUBJID", "PARAMCD"), dataset_name = "ADVS")

  expect_equal(attr(twice, "admiral_keys"), c("USUBJID", "PARAMCD"))
  expect_equal(attr(twice, "admiral_ds_name"), "ADVS")
  expect_equal(sum(class(twice) == "admiral_df"), 1L)
})

## Test 10: attributes of the supplied keys are dropped ----
test_that("set_admiral_keys Test 10: attributes of the supplied keys are dropped", {
  input <- tibble::tribble(
    ~USUBJID, ~AVAL,
    "1",      10
  )
  # keys extracted from a specification table carry a variable label, which
  # must not end up on the stored attribute
  labelled_keys <- structure("USUBJID", label = "Unique Subject Identifier")

  actual <- set_admiral_keys(input, keys = labelled_keys)

  expect_identical(attr(actual, "admiral_keys"), "USUBJID")
})

## Test 11: an empty keys vector is stored quietly ----
test_that("set_admiral_keys Test 11: an empty keys vector is stored quietly", {
  input <- tibble::tribble(
    ~USUBJID, ~AVAL,
    "1",      10
  )

  expect_silent(actual <- set_admiral_keys(input, keys = character(0)))

  # the empty attribute is a statement that the dataset has no key variables,
  # so it must be distinguishable from an attribute which was never set
  expect_identical(attr(actual, "admiral_keys"), character(0))
  expect_false(is.null(attr(actual, "admiral_keys")))
  expect_s3_class(actual, "admiral_df")
})

## Test 12: keys which are not in the dataset are stored quietly ----
test_that("set_admiral_keys Test 12: keys which are not in the dataset are stored quietly", { # nolint
  input <- tibble::tribble(
    ~USUBJID, ~PARAMCD, ~AVAL,
    "1",      "DIABP",     51
  )

  # a partially derived dataset can declare the structure it is being built
  # towards; reporting keys which no longer match the dataset is the job of the
  # tooling which reads the attribute, once the derivation is finished
  expect_silent(
    actual <- set_admiral_keys(input, keys = c("USUBJID", "PARAMCD", "AVISIT"))
  )

  expect_equal(
    attr(actual, "admiral_keys"),
    c("USUBJID", "PARAMCD", "AVISIT")
  )
})

## Test 13: the keys and the class survive a dplyr pipeline ----
test_that("set_admiral_keys Test 13: the keys and the class survive a dplyr pipeline", {
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

## Test 14: an error is issued for invalid arguments ----
test_that("set_admiral_keys Test 14: an error is issued for invalid arguments", {
  input <- tibble::tribble(
    ~USUBJID, ~AVAL,
    "1",      10
  )

  expect_error(set_admiral_keys(list(), keys = "USUBJID"))
  expect_error(set_admiral_keys(input, keys = 1))
  expect_error(set_admiral_keys(input, keys = exprs(USUBJID + 1)))
  expect_error(set_admiral_keys(input, keys = "USUBJID", dataset_name = c("A", "B")))
})

# get_admiral_df_type ----
## Test 15: classifies the common ADaM structures ----
test_that("get_admiral_df_type Test 15: classifies the common ADaM structures", {
  bds <- tibble::tribble(
    ~USUBJID, ~PARAMCD, ~AVAL, ~AVISIT,
    "1",      "MAP",       90, "BASELINE"
  )
  tte <- tibble::tribble(
    ~USUBJID, ~PARAMCD, ~AVAL, ~CNSR, ~STARTDT,
    "1",      "OS",       100,     0, as.Date("2020-01-01")
  )
  occds <- tibble::tribble(
    ~USUBJID, ~AEDECOD,   ~TRTEMFL,
    "1",      "HEADACHE", "Y"
  )
  adsl <- tibble::tribble(
    ~STUDYID, ~USUBJID, ~TRT01P,
    "P",      "1",      "A",
    "P",      "2",      "B"
  )
  other <- tibble::tribble(
    ~FOO, ~BAR,
    1,    2
  )

  expect_identical(get_admiral_df_type(bds), "BDS")
  expect_identical(get_admiral_df_type(tte), "TTE")
  expect_identical(get_admiral_df_type(occds), "OCCDS")
  expect_identical(get_admiral_df_type(adsl), "ADSL")
  expect_identical(get_admiral_df_type(other), "other")
})

## Test 16: the type precedence resolves datasets matching more than one type ----
test_that("get_admiral_df_type Test 16: the type precedence resolves datasets matching more than one type", { # nolint
  # TTE over BDS: a time-to-event dataset is a BDS dataset by variable content
  tte_and_bds <- tibble::tribble(
    ~USUBJID, ~PARAMCD, ~AVAL, ~CNSR, ~STARTDT,
    "1",      "OS",       100,     0, as.Date("2020-01-01")
  )
  expect_identical(get_admiral_df_type(tte_and_bds), "TTE")

  # BDS over OCCDS: an occurrence dataset which has acquired a parameter is
  # reported per parameter, so `PARAMCD` excludes the OCCDS branch outright
  bds_and_occds <- tibble::tribble(
    ~USUBJID, ~PARAMCD, ~AVAL, ~AEDECOD,   ~TRTEMFL,
    "1",      "AEHLT",      1, "HEADACHE", "Y"
  )
  expect_identical(get_admiral_df_type(bds_and_occds), "BDS")

  # OCCDS over ADSL: one adverse event per subject is still an occurrence
  # dataset, even though it is trivially one record per subject and carries a
  # treatment variable merged from ADSL
  occds_and_adsl <- tibble::tribble(
    ~USUBJID, ~AEDECOD,   ~TRT01P,   ~TRTSDT,
    "1",      "HEADACHE", "Placebo", as.Date("2020-01-01"),
    "2",      "NAUSEA",   "Drug",    as.Date("2020-01-02")
  )
  expect_identical(get_admiral_df_type(occds_and_adsl), "OCCDS")
})

## Test 17: an occurrence dataset needs more than a TERM-suffixed variable ----
test_that("get_admiral_df_type Test 17: an occurrence dataset needs more than a TERM-suffixed variable", { # nolint
  # `LONGTERM` ends in TERM but is not a `--TERM` variable, and the two-letter
  # SDTM domain prefix is what distinguishes the two
  not_occds <- tibble::tribble(
    ~USUBJID, ~LONGTERM, ~AVAL,
    "1",      "Y",          10,
    "1",      "N",          20
  )
  expect_identical(get_admiral_df_type(not_occds), "other")

  # the naming signal alone, with nothing identifying a subject, is not enough
  no_subject <- tibble::tribble(
    ~AEDECOD,   ~AESEV,
    "HEADACHE", "MILD"
  )
  expect_identical(get_admiral_df_type(no_subject), "other")

  # a genuine `--TERM` variable on a dataset which identifies a subject is
  expect_identical(
    get_admiral_df_type(tibble::tribble(
      ~USUBJID, ~CMTERM,
      "1",      "ASPIRIN"
    )),
    "OCCDS"
  )
})

## Test 18: an empty dataset is typed from the variables it declares ----
test_that("get_admiral_df_type Test 18: an empty dataset is typed from the variables it declares", { # nolint
  empty_bds <- tibble::tibble(
    USUBJID = character(0), PARAMCD = character(0), AVAL = numeric(0)
  )
  expect_identical(get_admiral_df_type(empty_bds), "BDS")

  empty_adsl <- tibble::tibble(USUBJID = character(0), TRT01P = character(0))
  expect_identical(get_admiral_df_type(empty_adsl), "ADSL")

  # an empty dataset cannot demonstrate one record per subject -- the record and
  # subject counts are trivially equal -- so the subject key alone does not make
  # it subject-level
  empty_shell <- tibble::tibble(USUBJID = character(0), AVAL = numeric(0))
  expect_identical(get_admiral_df_type(empty_shell), "other")
})

# is_adsl_structure ----
## Test 19: the records decide the structure, not the treatment variables ----
test_that("is_adsl_structure Test 19: the records decide the structure, not the treatment variables", { # nolint
  one_per_subject <- tibble::tribble(
    ~STUDYID, ~USUBJID, ~AGE,
    "P",      "1",        63,
    "P",      "2",        71
  )
  expect_true(is_adsl_structure(one_per_subject))

  many_per_subject <- tibble::tribble(
    ~STUDYID, ~USUBJID, ~AVAL,
    "P",      "1",         10,
    "P",      "1",         20
  )
  expect_false(is_adsl_structure(many_per_subject))

  # the first step of a BDS/OCCDS derivation merges the treatment variables on
  # from ADSL, so carrying them is no evidence of a subject-level structure
  expect_false(is_adsl_structure(mutate(many_per_subject, TRT02A = "Drug")))
  expect_false(is_adsl_structure(mutate(many_per_subject, TRT01P = "Drug")))
  expect_false(is_adsl_structure(mutate(many_per_subject, TRTSDT = as.Date("2020-01-01"))))
})

## Test 20: a grouped dataset is judged on its records overall ----
test_that("is_adsl_structure Test 20: a grouped dataset is judged on its records overall", {
  many_per_subject <- tibble::tribble(
    ~STUDYID, ~USUBJID, ~AVISIT,    ~AVAL,
    "P",      "1",      "BASELINE",    10,
    "P",      "1",      "WEEK 2",      20
  )

  expect_false(is_adsl_structure(many_per_subject))

  # `distinct()` on a grouped dataset silently adds the grouping variables, so
  # without `ungroup()` this would test one record per subject *per visit* and
  # report a findings dataset as subject-level
  expect_false(is_adsl_structure(group_by(many_per_subject, AVISIT)))
  expect_identical(get_admiral_df_type(group_by(many_per_subject, AVISIT)), "other")
})

## Test 22: a partially derived ADSL is not recognized ----
test_that("is_adsl_structure Test 22: a partially derived ADSL is not recognized", {
  # a known limitation, pinned so that it stays a deliberate choice: until an
  # ADSL reaches one record per subject or declares a treatment variable,
  # nothing distinguishes it from any other record-level dataset
  wip_adsl <- tibble::tribble(
    ~STUDYID, ~USUBJID, ~AGE, ~SEX,
    "P",      "1",        63, "M",
    "P",      "1",        63, "M"
  )

  expect_false(is_adsl_structure(wip_adsl))
  expect_identical(get_admiral_df_type(wip_adsl), "other")
})

## Test 23: the subject keys option is respected ----
test_that("is_adsl_structure Test 23: the subject keys option is respected", {
  # one record per STUDYID + USUBJID (the default subject keys), but `SUBJID`
  # repeats across the two studies
  input <- tibble::tribble(
    ~STUDYID, ~USUBJID, ~SUBJID, ~AGE,
    "A",      "A-1",    "1",       63,
    "B",      "B-1",    "1",       71
  )

  expect_true(is_adsl_structure(input))

  subject_keys <- get_admiral_option("subject_keys")
  withr::defer(set_admiral_options(subject_keys = subject_keys))
  set_admiral_options(subject_keys = exprs(SUBJID))

  # what counts as one record per subject follows the option, so the same
  # dataset is no longer subject-level
  expect_false(is_adsl_structure(input))
})

# get_admiral_df_type ----
## Test 24: each occurrence signal is recognized on its own ----
test_that("get_admiral_df_type Test 24: each occurrence signal is recognized on its own", {
  # two records for one subject, so that a dataset which fails every occurrence
  # test is not subject-level either and the fallback is visible as "other"
  occds <- function(...) {
    tibble::tibble(USUBJID = c("1", "1"), ...)
  }

  # the occurrence flag family of the ADaM OCCDS implementation guide, plus the
  # sponsor-numbered form
  for (flag in c("AOCCFL", "AOCCIFL", "AOCCSFL", "AOCCPFL", "AOCCPIFL", "AOCC02FL")) {
    input <- occds()
    input[[flag]] <- c("Y", "N")
    expect_identical(get_admiral_df_type(input), "OCCDS", info = flag)
  }

  # and each of the other alternatives, in isolation rather than alongside one
  # another as in Test 15
  expect_identical(get_admiral_df_type(occds(TRTEMFL = c("Y", "N"))), "OCCDS")
  expect_identical(get_admiral_df_type(occds(AEDECOD = c("A", "B"))), "OCCDS")
  expect_identical(get_admiral_df_type(occds(MHTERM = c("A", "B"))), "OCCDS")

  # `SRCSEQ`-style provenance and a lone severity variable are not occurrence
  # signals
  expect_identical(get_admiral_df_type(occds(AESEV = c("MILD", "SEVERE"))), "other")
})

## Test 25: the TTE signal needs all three of its variables ----
test_that("get_admiral_df_type Test 25: the TTE signal needs all three of its variables", {
  # two records per subject, so that dropping `PARAMCD` does not leave a dataset
  # which is subject-level by structure
  tte <- tibble::tribble(
    ~USUBJID, ~PARAMCD, ~AVAL, ~CNSR, ~STARTDT,
    "1",      "OS",       100,     0, as.Date("2020-01-01"),
    "1",      "PFS",       60,     1, as.Date("2020-01-01")
  )
  expect_identical(get_admiral_df_type(tte), "TTE")

  # dropping any one term falls back to BDS rather than staying TTE
  expect_identical(get_admiral_df_type(select(tte, -CNSR)), "BDS")
  expect_identical(get_admiral_df_type(select(tte, -STARTDT)), "BDS")
  expect_identical(get_admiral_df_type(select(tte, -PARAMCD)), "other")

  # BDS is satisfied by either analysis value variable
  expect_identical(
    get_admiral_df_type(tibble::tibble(USUBJID = "1", PARAMCD = "A", AVALC = "X")),
    "BDS"
  )
  # with neither, it is not a BDS dataset
  expect_identical(
    get_admiral_df_type(tibble::tibble(USUBJID = c("1", "1"), PARAMCD = "A", FOO = 1)),
    "other"
  )
})

## Test 26: a findings dataset with treatment variables merged on is not ADSL ----
test_that("get_admiral_df_type Test 26: a findings dataset with treatment variables merged on is not ADSL", { # nolint
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

  expect_false(is_adsl_structure(mid_derivation))
  expect_identical(get_admiral_df_type(mid_derivation), "other")

  # the finished ADSL, with the same variables and one record per subject, is
  expect_identical(
    get_admiral_df_type(distinct(mid_derivation, STUDYID, USUBJID, TRT01P, TRT01A, TRTSDT)),
    "ADSL"
  )
})

## Test 27: an error is issued for a non-data-frame dataset ----
test_that("get_admiral_df_type Test 27: an error is issued for a non-data-frame dataset", {
  expect_error(get_admiral_df_type(matrix(1:4, nrow = 2)))
  expect_error(get_admiral_df_type("ADSL"))
})

# minimal_unique_key ----
## Test 28: returns must_have when it is already unique ----
test_that("minimal_unique_key Test 28: returns must_have when it is already unique", {
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

## Test 29: drops candidates which do not discriminate ----
test_that("minimal_unique_key Test 29: drops candidates which do not discriminate", {
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

## Test 30: returns everything when uniqueness is never reached ----
test_that("minimal_unique_key Test 30: returns everything when uniqueness is never reached", {
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
## Test 31: the inferred key spans the common BDS shapes ----
test_that("infer_admiral_keys Test 31: the inferred key spans the common BDS shapes", {
  # exposure: one record per subject, parameter and dosing interval, keyed by
  # the interval start rather than by an analysis date or a visit
  adex <- tibble::tribble(
    ~USUBJID, ~PARAMCD, ~AVAL, ~ASTDT,                ~AENDT,
    "1",      "DOSE",      54, as.Date("2024-01-01"), as.Date("2024-01-31"),
    "1",      "DOSE",      81, as.Date("2024-02-01"), as.Date("2024-02-28")
  )
  expect_identical(infer_admiral_keys(adex), c("USUBJID", "PARAMCD", "ASTDT"))

  # population PK: keyed by relative time, with no visit structure at all
  adppk <- tibble::tribble(
    ~USUBJID, ~PARAMCD, ~NFRLT, ~AVAL,
    "1",      "CONC",        0,     0,
    "1",      "CONC",        1,    12,
    "1",      "CONC",        2,     8
  )
  expect_identical(infer_admiral_keys(adppk), c("USUBJID", "PARAMCD", "NFRLT"))

  # a derived record shares every analysis variable with the record it was
  # derived from, so only DTYPE separates the two within a visit
  locf <- tibble::tribble(
    ~USUBJID, ~PARAMCD, ~AVISITN, ~AVAL, ~DTYPE,
    "1",      "SYSBP",         2,   121, NA,
    "1",      "SYSBP",         2,   121, "LOCF",
    "1",      "SYSBP",         4,   130, NA,
    "1",      "SYSBP",         4,   130, "LOCF"
  )
  expect_identical(
    infer_admiral_keys(locf),
    c("USUBJID", "PARAMCD", "AVISITN", "DTYPE")
  )
})

## Test 32: multi-period designs are keyed by period ----
test_that("infer_admiral_keys Test 32: multi-period designs are keyed by period", {
  # the same visit recurs in each period, so nothing but APERIOD separates the
  # records -- a vaccine study shape
  input <- tibble::tribble(
    ~USUBJID, ~PARAMCD, ~AVISITN, ~APERIOD, ~AVAL,
    "1",      "SYSBP",         1,        1,   121,
    "1",      "SYSBP",         1,        2,   130,
    "1",      "SYSBP",         2,        1,   118,
    "1",      "SYSBP",         2,        2,   125
  )

  expect_identical(
    infer_admiral_keys(input),
    c("USUBJID", "PARAMCD", "AVISITN", "APERIOD")
  )
})

## Test 33: ADSL, TTE and unrecognized datasets ----
test_that("infer_admiral_keys Test 33: ADSL, TTE and unrecognized datasets", {
  adsl <- tibble::tribble(
    ~STUDYID, ~USUBJID, ~TRT01P,
    "P",      "1",      "Placebo",
    "P",      "2",      "Drug"
  )
  expect_identical(infer_admiral_keys(adsl), "USUBJID")

  adtte <- tibble::tribble(
    ~USUBJID, ~PARAMCD, ~AVAL, ~CNSR, ~STARTDT,
    "1",      "OS",       100,     0, as.Date("2020-01-01"),
    "1",      "PFS",       60,     1, as.Date("2020-01-01")
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

## Test 34: OCCDS is keyed by its sequence variable ----
test_that("infer_admiral_keys Test 34: OCCDS is keyed by its sequence variable", {
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
})

## Test 35: an unusable OCCDS sequence is reported rather than guessed ----
test_that("infer_admiral_keys Test 35: an unusable OCCDS sequence is reported rather than guessed", { # nolint
  occds <- function(...) {
    tibble::tibble(
      STUDYID = "P", USUBJID = c("1", "1"),
      AEDECOD = c("HEADACHE", "NAUSEA"), ...
    )
  }

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

## Test 36: nothing is inferred for a dataset with no records ----
test_that("infer_admiral_keys Test 36: nothing is inferred for a dataset with no records", {
  # every candidate key is trivially unique over no records, so inference would
  # report a structure the dataset has not demonstrated
  empty_bds <- tibble::tibble(
    USUBJID = character(0), PARAMCD = character(0), AVAL = numeric(0)
  )

  expect_identical(infer_admiral_keys(empty_bds), character(0))
})

## Test 37: a grouped dataset is keyed on its records overall ----
test_that("infer_admiral_keys Test 37: a grouped dataset is keyed on its records overall", {
  input <- tibble::tribble(
    ~USUBJID, ~PARAMCD, ~AVISITN, ~AVAL,
    "1",      "SYSBP",         2,   121,
    "1",      "SYSBP",         4,   130
  )

  # without `ungroup()` the uniqueness test would run within group, so the key
  # would stop at USUBJID + PARAMCD and the visit structure would be lost
  expect_identical(
    infer_admiral_keys(group_by(input, AVISITN)),
    c("USUBJID", "PARAMCD", "AVISITN")
  )
})

## Test 38: inference holds up against the {pharmaverseadam} datasets ----
test_that("infer_admiral_keys Test 38: inference holds up against the {pharmaverseadam} datasets", { # nolint
  skip_on_cran()
  skip_if_not_installed("pharmaverseadam")

  ds_names <- utils::data(package = "pharmaverseadam")$results[, "Item"]
  load_ds <- function(nm) {
    env <- new.env()
    suppressWarnings(utils::data(list = nm, package = "pharmaverseadam", envir = env))
    env[[nm]]
  }

  keys <- lapply(setNames(ds_names, ds_names), function(nm) {
    infer_admiral_keys(load_ds(nm))
  })

  # every one of these is a real dataset written by somebody who had never heard
  # of this feature, so a shape the candidate list does not cover shows up here
  # as a dataset whose record structure cannot be checked at all
  expect_identical(names(keys)[lengths(keys) == 0], character(0))

  # the shapes the candidate list was extended to cover, each of which was once
  # inferred as USUBJID + PARAMCD and reported thousands of false duplicates
  expect_identical(keys$adex, c("USUBJID", "PARAMCD", "ASTDTM")) # dosing interval
  expect_identical(keys$adppk, c("USUBJID", "PARAMCD", "NFRLT", "AFRLT")) # relative time
  expect_identical(
    keys$adpc, # derived records alongside their source
    c("USUBJID", "PARAMCD", "AVISITN", "ATPTN", "DTYPE")
  )
  # occurrence datasets carrying only the SDTM domain sequence
  expect_identical(keys$adae, c("USUBJID", "AESEQ"))
  expect_identical(keys$adcm, c("USUBJID", "CMSEQ"))
  expect_identical(keys$admh, c("USUBJID", "MHSEQ"))

  # the remaining duplicates were each checked by hand and are true positives in
  # the test data, not inference failures -- see `keys_design_notes.md`. Pinned
  # so that an inference regression cannot hide among them.
  dups <- vapply(ds_names, function(nm) {
    ds <- load_ds(nm)
    nrow(ds) - nrow(distinct(ds, !!!syms(keys[[nm]])))
  }, integer(1))
  expect_identical(
    dups[dups > 0],
    c(adeg = 63L, adpp = 1008L, advs = 39L)
  )
})
