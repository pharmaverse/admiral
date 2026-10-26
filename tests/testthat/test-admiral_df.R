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

# print.admiral_df ----
## Test 21: the name, type and keys are printed in a header ----
test_that("print.admiral_df Test 21: the name, type and keys are printed in a header", {
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

## Test 22: a dataset with only the class prints as before ----
test_that("print.admiral_df Test 22: a dataset with only the class prints as before", {
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
