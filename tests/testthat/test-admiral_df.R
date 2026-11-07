# new_admiral_df ----
## Test 1: adds the admiral_df class to a data frame ----
test_that("new_admiral_df Test 1: adds the admiral_df class to a data frame", {
  input <- tibble::tribble(
    ~USUBJID, ~AVAL,
    "1",      10,
    "2",      20
  )

  actual <- new_admiral_df(input)

  expect_equal(class(actual), c("admiral_df", class(input)))
  # Test that the dataset is not changed
  expect_equal(unclass(actual), unclass(input))
})

## Test 2: is idempotent when the class is already present ----
test_that("new_admiral_df Test 2: is idempotent when the class is already present", {
  input <- tibble::tribble(
    ~USUBJID, ~AVAL,
    "1",      10
  )

  once <- new_admiral_df(input)
  twice <- new_admiral_df(once)

  expect_equal(class(once), class(twice))
  expect_equal(sum(class(twice) == "admiral_df"), 1L)
})

## Test 3: returns NULL unchanged ----
test_that("new_admiral_df Test 3: returns NULL unchanged", {
  expect_null(new_admiral_df(NULL))
})

## Test 4: keeps the metadata already set ----
test_that("new_admiral_df Test 4: keeps the metadata already set", {
  input <- admiral_df(
    tibble::tribble(
      ~USUBJID, ~AVAL,
      "1",      10
    ),
    dataset_keys = exprs(USUBJID),
    dataset_name = "ADSL",
    dataset_class = "ADSL"
  )

  # every derivation tags its output, which must not lose the metadata the user
  # set on its input
  actual <- new_admiral_df(input)

  expect_identical(attr(actual, "admiral_dataset_keys"), "USUBJID")
  expect_identical(attr(actual, "admiral_dataset_name"), "ADSL")
  expect_identical(attr(actual, "admiral_dataset_class"), "ADSL")
})

# admiral_df ----
## Test 5: stores the keys, name, ADaM class, and the admiral_df class ----
test_that("admiral_df Test 5: stores the keys, name, ADaM class, and the admiral_df class", {
  input <- tibble::tribble(
    ~USUBJID, ~PARAMCD, ~AVISIT,    ~AVAL,
    "1",      "DIABP",  "BASELINE",    51,
    "1",      "SYSBP",  "BASELINE",   121
  )

  actual <- admiral_df(
    input,
    dataset_keys = exprs(USUBJID, PARAMCD, AVISIT),
    dataset_name = "ADVS",
    dataset_class = "BDS"
  )

  expect_equal(
    attr(actual, "admiral_dataset_keys"),
    c("USUBJID", "PARAMCD", "AVISIT")
  )
  expect_equal(attr(actual, "admiral_dataset_name"), "ADVS")
  expect_equal(attr(actual, "admiral_dataset_class"), "BDS")
  expect_s3_class(actual, "admiral_df")
  # Test that the data itself is not changed
  stripped <- actual
  attr(stripped, "admiral_dataset_keys") <- NULL
  attr(stripped, "admiral_dataset_name") <- NULL
  attr(stripped, "admiral_dataset_class") <- NULL
  class(stripped) <- class(input)
  expect_equal(stripped, input)
})

## Test 6: keys from a specification are passed via chr2vars() ----
test_that("admiral_df Test 6: keys from a specification are passed via chr2vars()", {
  input <- tibble::tribble(
    ~USUBJID, ~PARAMCD, ~AVAL,
    "1",      "DIABP",     51
  )
  # keys extracted from a specification table carry a variable label
  spec_keys <- c("USUBJID", "PARAMCD")
  attr(spec_keys, "label") <- "Variable Name"

  from_exprs <- admiral_df(input, dataset_keys = exprs(USUBJID, PARAMCD))
  from_spec <- admiral_df(input, dataset_keys = chr2vars(spec_keys))

  expect_identical(
    attr(from_spec, "admiral_dataset_keys"),
    attr(from_exprs, "admiral_dataset_keys")
  )
})

## Test 7: the metadata is stored, preserved, and replaced ----
test_that("admiral_df Test 7: the metadata is stored, preserved, and replaced", {
  input <- tibble::tribble(
    ~USUBJID, ~PARAMCD, ~AVAL,
    "1",      "DIABP",     51
  )

  # nothing is stored when there is nothing to store
  bare <- admiral_df(input)
  expect_s3_class(bare, "admiral_df")
  expect_null(attr(bare, "admiral_dataset_keys"))
  expect_null(attr(bare, "admiral_dataset_name"))
  expect_null(attr(bare, "admiral_dataset_class"))

  # an argument which is not repeated is not lost by a call which only revises
  # another
  named <- admiral_df(
    input,
    dataset_keys = exprs(USUBJID),
    dataset_name = "ADVS",
    dataset_class = "BDS"
  )
  rekeyed <- admiral_df(named, dataset_keys = exprs(USUBJID, PARAMCD))
  expect_equal(attr(rekeyed, "admiral_dataset_keys"), c("USUBJID", "PARAMCD"))
  expect_equal(attr(rekeyed, "admiral_dataset_name"), "ADVS")
  expect_equal(attr(rekeyed, "admiral_dataset_class"), "BDS")

  retyped <- admiral_df(named, dataset_class = "TTE")
  expect_equal(attr(retyped, "admiral_dataset_keys"), "USUBJID")
  expect_equal(attr(retyped, "admiral_dataset_class"), "TTE")

  # all give way when a later call supplies them
  renamed <- admiral_df(
    named,
    dataset_keys = exprs(USUBJID, PARAMCD),
    dataset_name = "ADSL",
    dataset_class = "ADSL"
  )
  expect_equal(attr(renamed, "admiral_dataset_keys"), c("USUBJID", "PARAMCD"))
  expect_equal(attr(renamed, "admiral_dataset_name"), "ADSL")
  expect_equal(attr(renamed, "admiral_dataset_class"), "ADSL")
})

## Test 8: the keys are stored as a plain character vector ----
test_that("admiral_df Test 8: the keys are stored as a plain character vector", {
  input <- tibble::tribble(
    ~USUBJID, ~AVAL,
    "1",      10
  )

  # `vars2chr()` returns a named vector, which must not end up on the stored
  # attribute
  actual <- admiral_df(input, dataset_keys = exprs(USUBJID, AVAL))

  expect_identical(attr(actual, "admiral_dataset_keys"), c("USUBJID", "AVAL"))
})

## Test 9: empty keys are stored, keys not in the dataset are not ----
test_that("admiral_df Test 9: empty keys are stored, keys not in the dataset are not", {
  input <- tibble::tribble(
    ~USUBJID, ~PARAMCD, ~AVAL,
    "1",      "DIABP",     51
  )

  # the empty attribute is a statement that the dataset has no key variables,
  # so it must be distinguishable from an attribute which was never set
  expect_silent(empty <- admiral_df(input, dataset_keys = exprs()))
  expect_identical(attr(empty, "admiral_dataset_keys"), character(0))
  expect_false(is.null(attr(empty, "admiral_dataset_keys")))
  expect_s3_class(empty, "admiral_df")

  # a key which is not in the dataset is most likely a typo, so the keys are
  # not set and a warning says so; the name and the class are still applied
  expect_snapshot(
    typo <- admiral_df(input, dataset_keys = exprs(USUBJID, PARMACD), dataset_name = "ADVS")
  )
  expect_null(attr(typo, "admiral_dataset_keys"))
  expect_identical(attr(typo, "admiral_dataset_name"), "ADVS")
  expect_s3_class(typo, "admiral_df")

  # re-keying with a typo keeps the keys set previously, and the warning says
  # which they are
  keyed <- admiral_df(input, dataset_keys = exprs(USUBJID, PARAMCD))
  expect_snapshot(
    rekeyed <- admiral_df(keyed, dataset_keys = exprs(USUBJID, PARMACD, AVISIT))
  )
  expect_identical(attr(rekeyed, "admiral_dataset_keys"), c("USUBJID", "PARAMCD"))
})

## Test 10: the metadata and the class survive a dplyr pipeline ----
test_that("admiral_df Test 10: the metadata and the class survive a dplyr pipeline", {
  input <- admiral_df(
    tibble::tribble(
      ~USUBJID, ~PARAMCD, ~AVAL,
      "1",      "DIABP",     51,
      "2",      "DIABP",     79
    ),
    dataset_keys = exprs(USUBJID, PARAMCD),
    dataset_name = "ADVS",
    dataset_class = "BDS"
  )
  adsl <- tibble::tribble(
    ~USUBJID, ~TRT01P,
    "1",      "Placebo",
    "2",      "Xanomeline"
  )

  # the usefulness of the attributes rests on them surviving the verbs a
  # derivation pipeline is built from -- which is a promise about `{dplyr}`,
  # not about admiral, so it is worth pinning against a dependency update
  actual <- input %>%
    mutate(BASE = AVAL) %>%
    filter(AVAL > 0) %>%
    arrange(USUBJID) %>%
    select(USUBJID, PARAMCD, AVAL, BASE) %>%
    left_join(adsl, by = "USUBJID")

  expect_equal(attr(actual, "admiral_dataset_keys"), c("USUBJID", "PARAMCD"))
  expect_equal(attr(actual, "admiral_dataset_name"), "ADVS")
  expect_equal(attr(actual, "admiral_dataset_class"), "BDS")
  expect_s3_class(actual, "admiral_df")
})

## Test 11: an error is issued for invalid arguments ----
test_that("admiral_df Test 11: an error is issued for invalid arguments", {
  input <- tibble::tribble(
    ~USUBJID, ~AVAL,
    "1",      10
  )

  expect_error(admiral_df(list(), dataset_keys = exprs(USUBJID)), class = "assert_data_frame")
  expect_error(admiral_df(input, dataset_keys = 1), class = "assert_vars")
  expect_error(admiral_df(input, dataset_keys = iris), class = "assert_vars")
  expect_error(admiral_df(input, dataset_keys = exprs(USUBJID + 1)), class = "assert_vars")
  expect_error(
    admiral_df(input, dataset_keys = exprs(USUBJID), dataset_name = c("A", "B")),
    class = "assert_character_scalar"
  )
  expect_error(
    admiral_df(input, dataset_class = c("BDS", "TTE")),
    class = "assert_character_scalar"
  )

  # the class is one of the ADaM dataset classes
  expect_snapshot(admiral_df(input, dataset_class = "ADAE"), error = TRUE)
})

## Test 12: device classes are accepted ----
test_that("admiral_df Test 12: device classes are accepted", {
  input <- tibble::tribble(
    ~STUDYID, ~SPDEVID, ~PARAMCD, ~AVAL,
    "P",      "D-1",    "BATTV",    3.1,
    "P",      "D-2",    "BATTV",    2.9
  )

  for (dataset_class in c("ADDL", "MDBDS", "MDTTE", "MDOCCDS")) {
    expect_silent(actual <- admiral_df(input, dataset_class = dataset_class))
    expect_identical(attr(actual, "admiral_dataset_class"), dataset_class)
  }
})

# infer_dataset_class ----
## Test 13: classifies the common ADaM structures ----
test_that("infer_dataset_class Test 13: classifies the common ADaM structures", {
  bds <- tibble::tribble(
    ~STUDYID, ~USUBJID, ~PARAMCD, ~AVAL, ~AVISIT,
    "P",      "1",      "MAP",       90, "BASELINE"
  )
  expect_identical(infer_dataset_class(bds), "BDS")

  # BDS is satisfied by either analysis value variable
  expect_identical(
    infer_dataset_class(
      tibble::tibble(STUDYID = "P", USUBJID = "1", PARAMCD = "A", AVALC = "X")
    ),
    "BDS"
  )
  # with neither, it is not a BDS dataset
  expect_identical(
    infer_dataset_class(tibble::tibble(
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
  expect_identical(infer_dataset_class(tte), "TTE")

  # the TTE signal needs all three of its variables: dropping any one term falls
  # back to BDS rather than staying TTE
  expect_identical(infer_dataset_class(select(tte, -CNSR)), "BDS")
  expect_identical(infer_dataset_class(select(tte, -STARTDT)), "BDS")
  expect_identical(infer_dataset_class(select(tte, -PARAMCD)), "other")

  adsl <- tibble::tribble(
    ~STUDYID, ~USUBJID, ~TRT01P,
    "P",      "1",      "A",
    "P",      "2",      "B"
  )
  expect_identical(infer_dataset_class(adsl), "ADSL")

  other <- tibble::tribble(
    ~FOO, ~BAR,
    1,    2
  )
  expect_identical(infer_dataset_class(other), "other")

  # a summary of a BDS dataset shares its `PARAMCD` and `AVAL` but has no
  # subjects, so it is not a BDS dataset itself
  bds_summary <- tibble::tribble(
    ~TRT01A, ~PARAMCD, ~AVAL,
    "A",     "SYSBP",    125,
    "B",     "SYSBP",    125
  )
  expect_identical(infer_dataset_class(bds_summary), "other")
})

## Test 14: each occurrence signal is recognized on its own ----
test_that("infer_dataset_class Test 14: each occurrence signal is recognized on its own", {
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
    expect_identical(infer_dataset_class(input), "OCCDS", info = flag)
  }

  expect_identical(infer_dataset_class(occds(TRTEMFL = c("Y", "N"))), "OCCDS")
  expect_identical(infer_dataset_class(occds(AEDECOD = c("A", "B"))), "OCCDS")
  expect_identical(infer_dataset_class(occds(MHTERM = c("A", "B"))), "OCCDS")
  expect_identical(infer_dataset_class(occds(CMTERM = c("A", "B"))), "OCCDS")
})

## Test 15: an occurrence signal which is not corroborated does not classify ----
test_that("infer_dataset_class Test 15: an occurrence signal which is not corroborated does not classify", { # nolint
  # `LONGTERM` ends in TERM but is not a `--TERM` variable, and the two-letter
  # SDTM domain prefix is what distinguishes the two
  not_occds <- tibble::tribble(
    ~STUDYID, ~USUBJID, ~LONGTERM, ~AVAL,
    "P",      "1",      "Y",          10,
    "P",      "1",      "N",          20
  )
  expect_identical(infer_dataset_class(not_occds), "other")

  # the naming signal alone, with nothing identifying a subject, is not enough
  no_subject <- tibble::tribble(
    ~AEDECOD,   ~AESEV,
    "HEADACHE", "MILD"
  )
  expect_identical(infer_dataset_class(no_subject), "other")

  # nor is a subset of the subject keys: `STUDYID` is constant within a study,
  # so a lookup table carrying it still identifies no subject
  study_only <- tibble::tribble(
    ~STUDYID, ~AEDECOD,
    "P",      "HEADACHE",
    "P",      "NAUSEA"
  )
  expect_identical(infer_dataset_class(study_only), "other")

  # a lone severity variable is not an occurrence signal either
  severity_only <- tibble::tibble(
    STUDYID = "P", USUBJID = c("1", "1"), AESEV = c("MILD", "SEVERE")
  )
  expect_identical(infer_dataset_class(severity_only), "other")
})

## Test 16: the class precedence resolves datasets matching more than one class ----
test_that("infer_dataset_class Test 16: the class precedence resolves datasets matching more than one class", { # nolint
  # BDS over OCCDS: an occurrence dataset which has acquired a parameter is
  # reported per parameter, so `PARAMCD` excludes the OCCDS branch outright
  bds_and_occds <- tibble::tribble(
    ~STUDYID, ~USUBJID, ~PARAMCD, ~AVAL, ~AEDECOD,   ~TRTEMFL,
    "P",      "1",      "AEHLT",      1, "HEADACHE", "Y"
  )
  expect_identical(infer_dataset_class(bds_and_occds), "BDS")

  # OCCDS over ADSL: one adverse event per subject is still an occurrence
  # dataset, even though it is trivially one record per subject and carries a
  # treatment variable merged from ADSL
  occds_and_adsl <- tibble::tribble(
    ~STUDYID, ~USUBJID, ~AEDECOD,   ~TRT01P,   ~TRTSDT,
    "P",      "1",      "HEADACHE", "Placebo", as.Date("2020-01-01"),
    "P",      "2",      "NAUSEA",   "Drug",    as.Date("2020-01-02")
  )
  expect_identical(infer_dataset_class(occds_and_adsl), "OCCDS")
})

## Test 17: an empty dataset is classified from the variables it declares ----
test_that("infer_dataset_class Test 17: an empty dataset is classified from the variables it declares", { # nolint
  empty_bds <- tibble::tibble(
    STUDYID = character(0), USUBJID = character(0), PARAMCD = character(0),
    AVAL = numeric(0)
  )
  expect_identical(infer_dataset_class(empty_bds), "BDS")

  # ADSL is decided by the records, and an empty dataset cannot show whether it
  # has one record per subject -- the record and subject counts are trivially
  # equal -- so the class is unknown, whatever treatment variables it carries
  empty_adsl <- tibble::tibble(
    STUDYID = character(0), USUBJID = character(0), TRT01P = character(0)
  )
  expect_identical(infer_dataset_class(empty_adsl), NA_character_)

  empty_shell <- tibble::tibble(
    STUDYID = character(0), USUBJID = character(0), AVAL = numeric(0)
  )
  expect_identical(infer_dataset_class(empty_shell), NA_character_)

  # the same holds with `PARAMCD` but no analysis value: with records this is
  # "other", without any it is unknown
  empty_paramcd <- tibble::tibble(
    STUDYID = character(0), USUBJID = character(0), PARAMCD = character(0)
  )
  expect_identical(infer_dataset_class(empty_paramcd), NA_character_)

  # where the variables decide the class, the missing records do not matter:
  # without the subject keys it is not subject data at all
  empty_lookup <- tibble::tibble(AEDECOD = character(0), AETERM = character(0))
  expect_identical(infer_dataset_class(empty_lookup), "other")
})

## Test 18: an error is issued for a non-data-frame dataset ----
test_that("infer_dataset_class Test 18: an error is issued for a non-data-frame dataset", {
  expect_error(infer_dataset_class(matrix(1:4, nrow = 2)), class = "assert_data_frame")
  expect_error(infer_dataset_class("ADSL"), class = "assert_data_frame")
})

## Test 19: ADSL is decided by the records ----
test_that("infer_dataset_class Test 19: ADSL is decided by the records", {
  one_per_subject <- tibble::tribble(
    ~STUDYID, ~USUBJID, ~AGE,
    "P",      "1",        63,
    "P",      "2",        71
  )
  expect_identical(infer_dataset_class(one_per_subject), "ADSL")

  many_per_subject <- tibble::tribble(
    ~STUDYID, ~USUBJID, ~AVAL,
    "P",      "1",         10,
    "P",      "1",         20
  )
  expect_identical(infer_dataset_class(many_per_subject), "other")

  # all the subject keys are needed: `STUDYID` is constant within a study, so a
  # single-study lookup table carrying it is one record per study, not subject
  study_only <- tibble::tribble(
    ~STUDYID, ~PARAM,
    "P",      "Code list A"
  )
  expect_identical(infer_dataset_class(study_only), "other")
})

## Test 20: treatment variables are no evidence of ADSL ----
test_that("infer_dataset_class Test 20: treatment variables are no evidence of ADSL", {
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

  expect_identical(infer_dataset_class(mid_derivation), "other")
  # the period-numbered form is equally no evidence
  expect_identical(infer_dataset_class(mutate(mid_derivation, TRT02A = "Drug")), "other")

  # the finished ADSL, with the same variables and one record per subject, is
  expect_identical(
    infer_dataset_class(distinct(mid_derivation, STUDYID, USUBJID, TRT01P, TRT01A, TRTSDT)),
    "ADSL"
  )
})

## Test 21: a grouped dataset is judged on its records overall ----
test_that("infer_dataset_class Test 21: a grouped dataset is judged on its records overall", {
  many_per_subject <- tibble::tribble(
    ~STUDYID, ~USUBJID, ~AVISIT,    ~AVAL,
    "P",      "1",      "BASELINE",    10,
    "P",      "1",      "WEEK 2",      20
  )

  # `distinct()` on a grouped dataset silently adds the grouping variables, so
  # without `ungroup()` this would test one record per subject *per visit* and
  # report a findings dataset as subject-level
  expect_identical(infer_dataset_class(group_by(many_per_subject, AVISIT)), "other")
})

## Test 22: the subject keys option is respected ----
test_that("infer_dataset_class Test 22: the subject keys option is respected", {
  # one record per STUDYID + USUBJID (the default subject keys), but `SUBJID`
  # repeats across the two studies
  input <- tibble::tribble(
    ~STUDYID, ~USUBJID, ~SUBJID, ~AGE,
    "A",      "A-1",    "1",       63,
    "B",      "B-1",    "1",       71
  )

  expect_identical(infer_dataset_class(input), "ADSL")

  subject_keys <- get_admiral_option("subject_keys")
  withr::defer(set_admiral_options(subject_keys = subject_keys))
  set_admiral_options(subject_keys = exprs(SUBJID))

  # what counts as one record per subject follows the option, so the same
  # dataset is no longer subject-level
  expect_identical(infer_dataset_class(input), "other")
})

## Test 23: device datasets are inferred as the device classes ----
test_that("infer_dataset_class Test 23: device datasets are inferred as the device classes", {
  # device data: `SPDEVID`, and no subject keys
  addl <- tibble::tribble(
    ~STUDYID, ~SPDEVID, ~DEVTYPE,
    "P",      "D-1",    "INSULIN PUMP",
    "P",      "D-2",    "INSULIN PUMP"
  )
  expect_identical(infer_dataset_class(addl), "ADDL")

  mdbds <- tibble::tribble(
    ~STUDYID, ~SPDEVID, ~PARAMCD, ~AVISITN, ~AVAL,
    "P",      "D-1",    "BATTV",         0,   3.1,
    "P",      "D-1",    "BATTV",         4,   2.9
  )
  expect_identical(infer_dataset_class(mdbds), "MDBDS")

  mdoccds <- tibble::tribble(
    ~STUDYID, ~SPDEVID, ~DEDECOD,    ~ASEQ,
    "P",      "D-1",    "OCCLUSION",     1,
    "P",      "D-1",    "OCCLUSION",     2
  )
  expect_identical(infer_dataset_class(mdoccds), "MDOCCDS")

  # TTE is tested before BDS for device data too
  mdtte <- tibble::tribble(
    ~STUDYID, ~SPDEVID, ~PARAMCD, ~AVAL, ~CNSR, ~STARTDT,
    "P",      "D-1",    "TTFAIL",   120,     0, lubridate::ymd("2022-01-01")
  )
  expect_identical(infer_dataset_class(mdtte), "MDTTE")

  # more than one record per device, and no parameter or occurrence signal
  device_other <- tibble::tribble(
    ~STUDYID, ~SPDEVID, ~LOTNUM,
    "P",      "D-1",    "L-1",
    "P",      "D-1",    "L-2"
  )
  expect_identical(infer_dataset_class(device_other), "other")

  # `USUBJID` is conditional in a device dataset, so the subject keys do not
  # decide the class: populated for every record, missing for some, or absent,
  # a dataset with `SPDEVID` for every record is device data
  mdbds_subject <- tibble::tribble(
    ~STUDYID, ~USUBJID, ~SPDEVID, ~PARAMCD, ~AVAL,
    "P",      "1",      "D-1",    "BATTV",    3.1,
    "P",      "2",      "D-2",    "BATTV",    2.9
  )
  expect_identical(infer_dataset_class(mdbds_subject), "MDBDS")
  expect_identical(
    infer_dataset_class(mutate(mdbds_subject, USUBJID = c("1", NA))),
    "MDBDS"
  )
  mdoccds_subject <- tibble::tribble(
    ~STUDYID, ~USUBJID, ~SPDEVID, ~DEDECOD,
    "P",      "1",      "D-1",    "OCCLUSION",
    "P",      "1",      "D-1",    "LEAK"
  )
  expect_identical(infer_dataset_class(mdoccds_subject), "MDOCCDS")

  # one record per device is ADDL, not ADSL, whether each subject has one
  # device or several
  addl_one_per_subject <- tibble::tribble(
    ~STUDYID, ~USUBJID, ~SPDEVID, ~DEVTYPE,
    "P",      "1",      "D-1",    "INSULIN PUMP",
    "P",      "2",      "D-2",    "INSULIN PUMP"
  )
  expect_identical(infer_dataset_class(addl_one_per_subject), "ADDL")
  addl_many_per_subject <- tibble::tribble(
    ~STUDYID, ~USUBJID, ~SPDEVID, ~DEVTYPE,
    "P",      "1",      "D-1",    "INSULIN PUMP",
    "P",      "1",      "D-2",    "CGM SENSOR"
  )
  expect_identical(infer_dataset_class(addl_many_per_subject), "ADDL")

  # a subject dataset records `SPDEVID` only for the records involving a
  # device, e.g. an adverse event dataset, so it is subject data -- whether the
  # missing values are `NA` or blank, as in data read from SAS
  adae <- tibble::tribble(
    ~STUDYID, ~USUBJID, ~SPDEVID, ~AEDECOD,
    "P",      "1",      "D-1",    "HEADACHE",
    "P",      "2",      NA,       "NAUSEA"
  )
  expect_identical(infer_dataset_class(adae), "OCCDS")
  expect_identical(infer_dataset_class(mutate(adae, SPDEVID = c("D-1", ""))), "OCCDS")

  # without the subject keys, a partly populated `SPDEVID` is neither class
  expect_identical(infer_dataset_class(select(adae, -USUBJID)), "other")

  # with no records, the variables decide as far as they can
  empty_mdbds <- tibble::tibble(
    STUDYID = character(0), SPDEVID = character(0), PARAMCD = character(0),
    AVAL = numeric(0)
  )
  expect_identical(infer_dataset_class(empty_mdbds), "MDBDS")
  empty_addl <- tibble::tibble(STUDYID = character(0), SPDEVID = character(0))
  expect_identical(infer_dataset_class(empty_addl), NA_character_)
})

# get_dataset_class ----
## Test 24: the declared class takes precedence over the inferred one ----
test_that("get_dataset_class Test 24: the declared class takes precedence over the inferred one", {
  input <- tibble::tribble(
    ~STUDYID, ~USUBJID, ~SPDEVID, ~AEDECOD,
    "P",      "1",      "D-1",    "DEVICE SITE PAIN"
  )

  # without a declared class, the class is inferred: `SPDEVID` is populated for
  # every record, so this is inferred as device data
  expect_identical(get_dataset_class(input), "MDOCCDS")

  # a declared class is returned as is, even where the inferred class differs --
  # here an adverse event dataset restricted to device-related events
  expect_identical(
    get_dataset_class(admiral_df(input, dataset_class = "OCCDS")),
    "OCCDS"
  )
})

# get_dataset_keys ----
## Test 25: the stored keys are returned as symbols ----
test_that("get_dataset_keys Test 25: the stored keys are returned as symbols", {
  input <- tibble::tribble(
    ~USUBJID, ~PARAMCD, ~AVAL,
    "1",      "DIABP",     51
  )

  # keys which have not been set are not inferred
  expect_null(get_dataset_keys(input))
  expect_null(get_dataset_keys(admiral_df(input, dataset_name = "ADVS")))

  # the keys come back in the form `admiral_df()` takes them, so they can be
  # passed straight to a `by_vars` argument
  keyed <- admiral_df(input, dataset_keys = exprs(USUBJID, PARAMCD))
  expect_equal(get_dataset_keys(keyed), exprs(USUBJID, PARAMCD), ignore_attr = TRUE)
  expect_identical(
    get_dataset_keys(admiral_df(input, dataset_keys = get_dataset_keys(keyed))),
    get_dataset_keys(keyed)
  )

  # empty keys are returned as an empty list, distinct from `NULL`
  expect_identical(get_dataset_keys(admiral_df(input, dataset_keys = exprs())), list())

  expect_error(get_dataset_keys("ADVS"), class = "assert_data_frame")
})

# get_dataset_name ----
## Test 26: the stored name is returned ----
test_that("get_dataset_name Test 26: the stored name is returned", {
  input <- tibble::tribble(
    ~USUBJID, ~AVAL,
    "1",      10
  )

  expect_null(get_dataset_name(input))
  expect_identical(get_dataset_name(admiral_df(input, dataset_name = "ADVS")), "ADVS")

  expect_error(get_dataset_name("ADVS"), class = "assert_data_frame")
})

# print.admiral_df ----
## Test 27: the name, class and keys are printed in a header ----
test_that("print.admiral_df Test 27: the name, class and keys are printed in a header", {
  input <- tibble::tribble(
    ~STUDYID, ~USUBJID, ~PARAMCD, ~AVISIT,    ~AVAL,
    "P",      "1",      "DIABP",  "BASELINE",    51,
    "P",      "1",      "SYSBP",  "BASELINE",   121
  )

  expect_snapshot(
    print(admiral_df(
      input,
      dataset_keys = exprs(USUBJID, PARAMCD, AVISIT),
      dataset_name = "ADVS",
      dataset_class = "BDS"
    ))
  )

  # each attribute is shown on its own; without a declared class, the inferred
  # class is shown and marked as such
  expect_snapshot(print(admiral_df(input, dataset_keys = exprs(USUBJID, PARAMCD, AVISIT))))
  expect_snapshot(print(admiral_df(input, dataset_class = "BDS")))

  # a declared class is shown instead of the inferred class, even where they differ
  expect_snapshot(print(admiral_df(input, dataset_keys = exprs(USUBJID), dataset_class = "TTE")))

  # an inferred class which is not an ADaM dataset class is left out rather than
  # shown as "other": without the subject keys this is not subject data
  lookup <- tibble::tribble(
    ~PARAMCD, ~PARAM,
    "DIABP",  "Diastolic Blood Pressure (mmHg)"
  )
  expect_snapshot(
    print(admiral_df(lookup, dataset_keys = exprs(PARAMCD), dataset_name = "PARAMS"))
  )

  # empty keys record that no keys are defined, so they are shown as such
  expect_snapshot(print(admiral_df(input, dataset_keys = exprs(), dataset_name = "ADVS")))

  # a plain data frame prints with its own method below the header
  expect_snapshot(
    print(admiral_df(as.data.frame(input), dataset_keys = exprs(USUBJID), dataset_name = "ADVS"))
  )
})

## Test 28: a dataset with only the class prints as before ----
test_that("print.admiral_df Test 28: a dataset with only the class prints as before", {
  input <- tibble::tribble(
    ~USUBJID, ~AVAL,
    "1",      10
  )

  # every derivation tags its output with the class, so a header without
  # anything to show would change how every admiral dataset prints
  expect_identical(
    capture.output(print(new_admiral_df(input))),
    capture.output(print(input))
  )
  expect_invisible(print(new_admiral_df(input)))
})
