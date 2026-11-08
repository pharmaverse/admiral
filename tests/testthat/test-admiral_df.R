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

# summary.admiral_df ----
## Test 29: the common facts and the declared record structure are reported ----
test_that("summary.admiral_df Test 29: the common facts and the declared record structure are reported", { # nolint
  input <- tibble::tribble(
    ~STUDYID, ~USUBJID, ~PARAMCD, ~AVISIT,    ~AVAL,
    "P",      "1",      "DIABP",  "BASELINE",    51,
    "P",      "1",      "SYSBP",  "BASELINE",   121,
    "P",      "2",      "DIABP",  "BASELINE",    79
  ) %>%
    admiral_df(dataset_keys = exprs(USUBJID, PARAMCD, AVISIT), dataset_name = "ADVS")

  actual <- summary(input)

  expect_s3_class(actual, "summary_admiral_df")
  expect_identical(actual$name, "ADVS")
  expect_identical(actual$dataset_class, "BDS")
  expect_identical(actual$n_subjects, 2L)
  expect_identical(actual$n_obs, 3L)
  expect_identical(actual$n_vars, 5L)
  expect_identical(
    actual$structure,
    list(
      keys = c("USUBJID", "PARAMCD", "AVISIT"),
      missing_keys = character(0),
      n_duplicates = 0L
    )
  )
})

## Test 30: the structure is not checked if no keys were set ----
test_that("summary.admiral_df Test 30: the structure is not checked if no keys were set", {
  input <- tibble::tribble(
    ~STUDYID, ~USUBJID, ~PARAMCD, ~AVISIT,    ~AVAL,
    "P",      "1",      "SYSBP",  "BASELINE",   121,
    "P",      "1",      "SYSBP",  "WEEK 2",     118,
    "P",      "2",      "SYSBP",  "BASELINE",   130
  ) %>%
    admiral_df()

  actual <- summary(input)

  expect_identical(actual$name, "input")
  # keys are never inferred
  expect_null(actual$structure$keys)
  expect_identical(actual$structure$n_duplicates, NA_integer_)
  expect_snapshot(print(actual))
})

## Test 31: duplicate, missing and empty declared keys ----
test_that("summary.admiral_df Test 31: duplicate, missing and empty declared keys", {
  input <- tibble::tribble(
    ~STUDYID, ~USUBJID, ~PARAMCD, ~AVISIT,    ~AVAL,
    "P",      "1",      "SYSBP",  "BASELINE",   121,
    "P",      "1",      "SYSBP",  "WEEK 2",     118,
    "P",      "2",      "SYSBP",  "BASELINE",   130
  )

  duplicated <- summary(admiral_df(input, dataset_keys = exprs(USUBJID, PARAMCD)))
  expect_identical(duplicated$structure$n_duplicates, 1L)

  # a renamed key cannot be checked, rather than being checked on part of the key
  renamed <- admiral_df(input, dataset_keys = exprs(USUBJID, PARAMCD, AVISIT)) %>%
    dplyr::rename(VISIT = AVISIT)
  expect_identical(summary(renamed)$structure$missing_keys, "AVISIT")
  expect_identical(summary(renamed)$structure$n_duplicates, NA_integer_)

  no_keys <- summary(admiral_df(input, dataset_keys = exprs()))
  expect_identical(no_keys$structure$keys, character(0))

  # each outcome of the check is printed as its own kind of bullet
  expect_snapshot(print(duplicated))
  expect_snapshot(print(summary(renamed)))
  expect_snapshot(print(no_keys))
})

## Test 32: the ADSL summary ----
test_that("summary.admiral_df Test 32: the ADSL summary", {
  input <- tibble::tribble(
    ~STUDYID, ~USUBJID, ~TRT01P,   ~SAFFL, ~ITTFL,
    "P",      "1",      "Placebo", "Y",    "Y",
    "P",      "2",      "Drug A",  "Y",    "Y",
    "P",      "3",      "Drug A",  "N",    "Y"
  ) %>%
    # the class is declared, as duplicated subject records would otherwise make
    # the inferred class "other"
    admiral_df(dataset_keys = exprs(USUBJID), dataset_name = "ADSL", dataset_class = "ADSL")

  actual <- summary(input)

  expect_identical(actual$details$trt_var, "TRT01P")
  expect_identical(actual$details$trt$subjects, c(2L, 1L))
  expect_identical(actual$details$pops, c(SAFFL = 2L, ITTFL = 3L))
  expect_snapshot(print(actual))

  # duplicated subject records are reported without losing the ADSL details
  duplicated <- summary(dplyr::bind_rows(input, input[1, ]))
  expect_identical(duplicated$dataset_class, "ADSL")
  expect_identical(duplicated$structure$n_duplicates, 1L)
  expect_identical(duplicated$details$trt$subjects, c(2L, 1L))
  # the population flags count subjects, not the duplicated records
  expect_identical(duplicated$details$pops, c(SAFFL = 2L, ITTFL = 3L))
})

## Test 33: the BDS summary ----
test_that("summary.admiral_df Test 33: the BDS summary", {
  input <- tibble::tribble(
    ~STUDYID, ~USUBJID, ~PARAMCD, ~AVISIT,    ~AVISITN, ~AVAL,
    "P",      "1",      "SYSBP",  "WEEK 2",          2,   118,
    "P",      "1",      "SYSBP",  "BASELINE",        0,   121,
    "P",      "1",      "DIABP",  "BASELINE",        0,    79,
    "P",      "2",      "SYSBP",  "BASELINE",        0,   130
  ) %>%
    admiral_df()

  actual <- summary(input)

  expect_identical(actual$details$params$PARAMCD, c("DIABP", "SYSBP"))
  expect_identical(actual$details$params$records, c(1L, 3L))
  expect_identical(actual$details$params$subjects, c(1L, 2L))
  # visits are in `AVISITN` order
  expect_identical(actual$details$visits, c("BASELINE", "WEEK 2"))
  expect_snapshot(print(actual))

  # a long parameter table is truncated on printing, but kept in the object
  many_params <- tibble::tibble(
    STUDYID = "P", USUBJID = "1", PARAMCD = sprintf("P%02d", 1:12), AVAL = 1
  ) %>%
    admiral_df()
  many_summary <- summary(many_params)
  expect_identical(nrow(many_summary$details$params), 12L)
  expect_snapshot(print(many_summary))
})

## Test 34: the OCCDS summary ----
test_that("summary.admiral_df Test 34: the OCCDS summary", {
  input <- tibble::tribble(
    ~STUDYID, ~USUBJID, ~AESEQ, ~AEBODSYS,                ~AEDECOD,   ~TRTEMFL,
    "P",      "1",           1, "NERVOUS SYSTEM",         "HEADACHE", "Y",
    "P",      "1",           2, "GASTROINTESTINAL",       "NAUSEA",   NA,
    "P",      "2",           1, "NERVOUS SYSTEM",         "HEADACHE", "Y"
  ) %>%
    admiral_df()

  actual <- summary(input)

  expect_identical(actual$details$terms, c(AEBODSYS = 2L, AEDECOD = 2L))
  expect_identical(actual$details$trtem, c(emergent = 2L, total = 3L))
  expect_snapshot(print(actual))
})

## Test 35: the TTE summary ----
test_that("summary.admiral_df Test 35: the TTE summary", {
  input <- tibble::tribble(
    ~STUDYID, ~USUBJID, ~PARAMCD, ~STARTDT, ~CNSR, ~AVAL,
    "P", "1", "OS", lubridate::ymd("2020-01-01"), 0, 10,
    "P", "2", "OS", lubridate::ymd("2020-01-01"), 1, 12,
    "P", "1", "PFS", lubridate::ymd("2020-01-01"), 1, 3
  ) %>%
    admiral_df()

  actual <- summary(input)

  expect_identical(actual$details$events$PARAMCD, c("OS", "PFS"))
  expect_identical(actual$details$events$events, c(1L, 0L))
  expect_identical(actual$details$events$censored, c(1L, 1L))
  expect_snapshot(print(actual))
})

## Test 36: datasets with little to report do not error ----
test_that("summary.admiral_df Test 36: datasets with little to report do not error", {
  # not subject data: no subject count and no record structure
  lookup <- admiral_df(tibble::tibble(PARAMCD = c("DIABP", "SYSBP")))
  expect_snapshot(print(summary(lookup)))

  # no records: the class cannot be determined
  empty <- admiral_df(tibble::tibble(STUDYID = character(), USUBJID = character()))
  expect_identical(summary(empty)$dataset_class, NA_character_)
  expect_identical(summary(empty)$details, list())
  expect_snapshot(print(summary(empty)))

  # a work-in-progress ADSL without treatment or population flags
  adsl <- admiral_df(tibble::tibble(STUDYID = "P", USUBJID = c("1", "2"), TRT01P = NA))
  expect_identical(summary(adsl)$details, list())
})

## Test 37: print.summary_admiral_df returns its input invisibly ----
test_that("print.summary_admiral_df Test 37: returns its input invisibly", {
  input <- admiral_df(tibble::tibble(STUDYID = "P", USUBJID = "1"))
  actual <- summary(input)

  # cli prints through messages
  msgs <- capture.output(res <- withVisible(print(actual)), type = "message")
  expect_match(msgs[1], "ADSL summary")
  expect_false(res$visible)
  expect_identical(res$value, actual)
})

# print_summary_table ----
## Test 38: tables are aligned, counts are not rounded, NA is shown as `.` ----
test_that("print_summary_table Test 38: tables are aligned, counts are not rounded, NA is shown as `.`", { # nolint
  tbl <- tibble::tibble(
    PARAMCD = c("DIABP", "SYSBP", NA),
    records = c(11509L, 7L, 1L),
    mean = c(80.123456, NA, 1),
    subjects = c(254, 3, NA)
  )

  expect_snapshot(print_summary_table(tbl))
  expect_snapshot(print_summary_table(tbl, max_rows = 1, what = "parameter"))
})

## Test 39: the ADSL subject flow, outcome and demographics ----
test_that("summary.admiral_df Test 39: the ADSL subject flow, outcome and demographics", { # nolint
  input <- tibble::tibble(
    STUDYID = "P",
    USUBJID = c("1", "2", "3", "4"),
    RANDDT = lubridate::ymd(c("2020-01-01", "2020-01-01", "2020-01-01", NA)),
    TRTSDT = lubridate::ymd(c("2020-01-02", "2020-01-03", NA, NA)),
    TRTDURD = c(100, 20, NA, NA),
    EOSSTT = c("COMPLETED", "DISCONTINUED", "DISCONTINUED", NA),
    DTHFL = c("N", "Y", "N", NA),
    DTHCAUS = c(NA, "ADVERSE EVENT", NA, NA),
    AGE = c(63, 71, 55, 40),
    SEX = c("F", "M", "F", "M")
  ) %>%
    # the class is declared, as duplicated subject records would otherwise make
    # the inferred class "other"
    admiral_df(dataset_keys = exprs(USUBJID), dataset_name = "ADSL", dataset_class = "ADSL")

  actual <- summary(input)$details

  expect_identical(actual$flow, c(randomized = 3L, treated = 2L, completed = 1L))
  expect_identical(actual$trtdur, c(median = 60, min = 20, max = 100))
  expect_identical(actual$eosstt, c(COMPLETED = 1L, DISCONTINUED = 2L, "NA" = 1L))
  # the causes are counted among the deaths only
  expect_identical(actual$n_deaths, 1L)
  expect_identical(actual$dthcaus, c("ADVERSE EVENT" = 1L))
  expect_identical(actual$age, c(median = 59, min = 40, max = 71))
  expect_identical(actual$sex, c(F = 2L, M = 2L))
  expect_snapshot(print(summary(input)))

  # duplicated subject records do not inflate the counts or skew the
  # distributions
  duplicated <- summary(dplyr::bind_rows(input, input[2, ], input[2, ]))$details
  expect_identical(duplicated[names(actual)], actual)

  # the flow falls back to `RANDFL` and `TRTSDTM`; items whose variables are
  # missing or not derived yet are left out
  fallback <- tibble::tibble(
    STUDYID = "P", USUBJID = c("1", "2"), RANDFL = c("Y", NA),
    TRTSDTM = lubridate::ymd_hms("2020-01-01 10:00:00", NA), DTHFL = NA, AGE = NA_real_
  ) %>%
    admiral_df()
  fallback_details <- summary(fallback)$details
  expect_identical(fallback_details$flow, c(randomized = 1L, treated = 1L))
  expect_null(fallback_details$n_deaths)
  expect_null(fallback_details$age)
  expect_null(fallback_details$sex)
})

## Test 40: the BDS per-parameter statistics and derived records ----
test_that("summary.admiral_df Test 40: the BDS per-parameter statistics and derived records", { # nolint
  input <- tibble::tribble(
    ~STUDYID, ~USUBJID, ~PARAMCD, ~PARAM,                  ~AVISIT,    ~AVAL, ~DTYPE,
    "P",      "1",      "SYSBP",  "Systolic BP (mmHg)",    "BASELINE",   120, NA,
    "P",      "1",      "SYSBP",  "Systolic BP (mmHg)",    "WEEK 2",      NA, NA,
    "P",      "1",      "SYSBP",  "Systolic BP (mmHg)",    "WEEK 2",     125, "LOCF",
    "P",      "2",      "SYSBP",  "Systolic BP (mm Hg)",   "BASELINE",  130.5, NA,
    "P",      "2",      "DIABP",  "Diastolic BP (mmHg)",   NA,            NA, NA
  ) %>%
    admiral_df()

  actual <- summary(input)$details

  expect_identical(
    actual$params,
    tibble::tibble(
      PARAMCD = c("DIABP", "SYSBP"),
      # every `PARAM` of a `PARAMCD` is shown
      PARAM = c("Diastolic BP (mmHg)", "Systolic BP (mmHg) / Systolic BP (mm Hg)"),
      records = c(1L, 4L),
      subjects = c(1L, 2L),
      visits = c(0L, 2L),
      missing = c(1L, 1L),
      min = c(NA, 120),
      median = c(NA, 125),
      max = c(NA, 130.5)
    )
  )
  expect_identical(actual$dtype, c(LOCF = 1L))
  # `PARAM` is not printed, the fractional maximum is
  expect_snapshot(print(summary(input)))

  # without `AVAL`, `AVISIT`, `PARAM`, or `DTYPE` only the counts are reported
  minimal <- summary(admiral_df(tibble::tibble(
    STUDYID = "P", USUBJID = "1", PARAMCD = "A", AVALC = "X"
  )))$details
  expect_named(minimal$params, c("PARAMCD", "records", "subjects"))
  expect_null(minimal$dtype)
})

## Test 41: a declared class does not require the variables of that class ----
test_that("summary.admiral_df Test 41: a declared class does not require the variables of that class", { # nolint
  # a declared class is not checked against the variables, so the details of
  # the class must not assume them
  input <- tibble::tibble(FOO = 1:2)
  for (dataset_class in admiral_dataset_classes) {
    actual <- summary(admiral_df(input, dataset_class = dataset_class))
    expect_identical(actual$dataset_class, dataset_class)
    expect_identical(actual$details, list())
    expect_no_error(capture.output(print(actual), type = "message"))
  }

  # TTE needs `CNSR` as well as `PARAMCD`
  tte <- admiral_df(tibble::tibble(PARAMCD = "OS", AVAL = 1), dataset_class = "TTE")
  expect_null(summary(tte)$details$events)

  # without a subject variable none of the ADSL items, which all count
  # subjects, can be reported
  adsl <- admiral_df(
    tibble::tibble(TRT01P = c("A", "B"), SAFFL = c("Y", "N"), AGE = c(40, 50)),
    dataset_class = "ADSL"
  )
  expect_identical(summary(adsl)$details, list())
})

## Test 42: more than one baseline per subject and parameter is reported ----
test_that("summary.admiral_df Test 42: more than one baseline per subject and parameter is reported", { # nolint
  input <- tibble::tribble(
    ~STUDYID, ~USUBJID, ~PARAMCD, ~AVISIT,     ~ABLFL, ~AVAL,
    "P",      "1",      "SYSBP",  "SCREENING", "Y",      120,
    "P",      "1",      "SYSBP",  "BASELINE",  "Y",      121,
    "P",      "1",      "DIABP",  "BASELINE",  "Y",       80,
    "P",      "1",      "DIABP",  "WEEK 2",    NA,        82,
    "P",      "2",      "SYSBP",  "WEEK 2",    NA,       130
  ) %>%
    admiral_df()

  actual <- summary(input)$details
  expect_identical(actual$baseline_by, c("USUBJID", "PARAMCD"))
  # each combination is counted once, however many baselines it has; subject 2
  # has no baseline, which is not a defect
  expect_identical(actual$multiple_baselines, 1L)
  expect_identical(summary(input[-1, ])$details$multiple_baselines, 0L)
  expect_snapshot(print(summary(input)))
  expect_snapshot(print(summary(input[-1, ])))

  # a baseline per baseline type, timepoint, or period is not a duplicate
  for (by_var in c("BASETYPE", "ATPT", "APERIOD")) {
    by_input <- dplyr::mutate(input, !!by_var := c("A", "B", "A", "A", "A"))
    by_actual <- summary(by_input)$details
    expect_identical(by_actual$baseline_by, c("USUBJID", "PARAMCD", by_var))
    expect_identical(by_actual$multiple_baselines, 0L)
  }

  # without `ABLFL` the check is not run, and not mentioned
  no_ablfl <- summary(dplyr::select(input, -ABLFL))
  expect_null(no_ablfl$details$multiple_baselines)
  expect_no_match(
    paste(capture.output(print(no_ablfl), type = "message"), collapse = "\n"),
    "Checks|baseline"
  )
})

# print_summary_checks ----
## Test 43: failed checks are listed, passed ones combined, absent ones left out ----
test_that("print_summary_checks Test 43: failed checks are listed, passed ones combined, absent ones left out", { # nolint
  checks <- list(
    a = list(failed = "{d$a} a failure{?s}", passed = "no a"),
    b = list(failed = "{d$b} b failure{?s}", passed = "no b"),
    c = list(failed = "{d$c} c failure{?s}", passed = "no c"),
    d = list(failed = "{d$d} d failure{?s}", passed = "no d")
  )
  d <- list(a = 2L, b = 0L, c = 0L)

  expect_snapshot(print_summary_checks(d, checks))
  # nothing is printed if no check ran
  expect_silent(print_summary_checks(list(), checks))
})

## Test 44: the OCCDS severity and seriousness ----
test_that("summary.admiral_df Test 44: the OCCDS severity and seriousness", {
  input <- tibble::tribble(
    ~STUDYID, ~USUBJID, ~AEDECOD,   ~ASEV,      ~AESEV,     ~AETOXGR, ~AESER,
    "P",      "1",      "HEADACHE", "MILD",     "MILD",     "1",      "N",
    "P",      "1",      "NAUSEA",   "SEVERE",   "MODERATE", "3",      "Y",
    "P",      "2",      "HEADACHE", NA,         NA,         NA,       NA
  ) %>%
    admiral_df()

  # the analysis severity is preferred, and a missing value is counted
  actual <- summary(input)$details
  expect_identical(actual$sev_var, "ASEV")
  expect_identical(actual$sev, c(MILD = 1L, SEVERE = 1L, "NA" = 1L))
  expect_identical(actual$ser_var, "AESER")
  expect_identical(actual$serious, c(serious = 1L, total = 3L))
  expect_snapshot(print(summary(input)))

  # a severity which is not derived yet gives way to the next, then to the
  # toxicity grade
  expect_identical(
    summary(dplyr::mutate(input, ASEV = NA))$details$sev_var,
    "AESEV"
  )
  expect_identical(
    summary(dplyr::select(input, -ASEV, -AESEV))$details$sev_var,
    "AETOXGR"
  )

  # neither is reported without the variables
  none <- summary(dplyr::select(input, -ASEV, -AESEV, -AETOXGR, -AESER))$details
  expect_null(none$sev)
  expect_null(none$serious)
})

# check_occurrence_flags ----
## Test 45: each flag is checked within the group its name defines ----
test_that("check_occurrence_flags Test 45: each flag is checked within the group its name defines", { # nolint
  input <- tibble::tribble(
    ~USUBJID, ~AEBODSYS, ~AEDECOD,   ~AOCCFL, ~AOCCSFL, ~AOCCPFL, ~AOCCPIFL, ~AOCC01FL,
    "1",      "NERV",    "HEADACHE", "Y",     "Y",      "Y",      "Y",       "Y",
    "1",      "NERV",    "DIZZY",    NA,      NA,       "Y",      "Y",       "Y",
    "1",      "GI",      "NAUSEA",   NA,      "Y",      "Y",      "Y",       NA,
    "2",      "NERV",    "HEADACHE", "Y",     "Y",      "Y",      NA,        NA
  )

  actual <- check_occurrence_flags(input, "USUBJID")
  expect_identical(
    actual$occ_flag_by,
    list(
      AOCCFL = "USUBJID",
      AOCCSFL = c("USUBJID", "AEBODSYS"),
      AOCCPFL = c("USUBJID", "AEDECOD"),
      AOCCPIFL = c("USUBJID", "AEDECOD")
    )
  )
  expect_identical(
    actual$occ_flag_dups,
    c(AOCCFL = 0L, AOCCSFL = 0L, AOCCPFL = 0L, AOCCPIFL = 0L)
  )
  # the level of a sponsor-defined flag cannot be told from its name
  expect_identical(actual$occ_flags_unchecked, "AOCC01FL")

  # a second "Y" in the same group, each group counted once
  dup <- dplyr::bind_rows(input, input[1:2, ])
  expect_identical(
    check_occurrence_flags(dup, "USUBJID")$occ_flag_dups,
    c(AOCCFL = 1L, AOCCSFL = 1L, AOCCPFL = 2L, AOCCPIFL = 2L)
  )

  # the first occurrence is flagged per period where there are periods
  per_period <- dplyr::mutate(dup, APERIOD = c(1, 1, 1, 1, 2, 2))
  per_period_actual <- check_occurrence_flags(per_period, "USUBJID")
  expect_identical(per_period_actual$occ_flag_by$AOCCFL, c("USUBJID", "APERIOD"))
  expect_true(all(per_period_actual$occ_flag_dups == 0L))

  # without the level variable the flag cannot be checked
  no_bodsys <- check_occurrence_flags(dplyr::select(input, -AEBODSYS), "USUBJID")
  expect_identical(no_bodsys$occ_flags_unchecked, c("AOCCSFL", "AOCC01FL"))

  # a single grouping variable works on a plain data frame too
  expect_identical(
    check_occurrence_flags(as.data.frame(dup[, c("USUBJID", "AOCCFL")]), "USUBJID"),
    list(occ_flag_by = list(AOCCFL = "USUBJID"), occ_flag_dups = c(AOCCFL = 1L))
  )

  # nothing to check without flags or subjects
  expect_identical(check_occurrence_flags(input[, 1:3], "USUBJID"), list())
  expect_identical(check_occurrence_flags(input, character(0)), list())

  expect_snapshot(print(summary(admiral_df(dplyr::mutate(dup, STUDYID = "P")))))
})

## Test 46: the dataset is compared with the supplied ADSL ----
test_that("summary.admiral_df Test 46: the dataset is compared with the supplied ADSL", {
  adsl <- tibble::tribble(
    ~STUDYID, ~USUBJID, ~SAFFL, ~TRT01P,   ~TRT01A,
    "P",      "1",      "Y",    "Placebo", "Placebo",
    "P",      "2",      "Y",    "Drug A",  "Drug A",
    "P",      "3",      "N",    "Drug A",  NA,
    "P",      "4",      "Y",    "Drug A",  "Drug A"
  )
  adae <- tibble::tribble(
    ~STUDYID, ~USUBJID, ~AEDECOD,   ~TRTEMFL,
    "P",      "1",      "HEADACHE", "Y",
    "P",      "1",      "NAUSEA",   "Y",
    "P",      "2",      "NAUSEA",   NA,
    "P",      "9",      "NAUSEA",   "Y",
    # a record which belongs to no subject is not compared
    NA,       NA,       "NAUSEA",   "Y"
  ) %>%
    admiral_df()

  actual <- summary(adae, adsl = adsl)$vs_adsl
  expect_identical(actual$n_adsl, 4L)
  expect_identical(actual$n_common, 2L)
  expect_identical(actual$n_orphans, 1L)
  expect_identical(actual$orphans, tibble::tibble(USUBJID = "9"))
  expect_identical(actual$saffl, c(covered = 2L, total = 3L))
  # the actual treatment ahead of the planned one; subjects without one are
  # counted too
  expect_identical(actual$arm_var, "TRT01A")
  expect_identical(
    actual$arm,
    tibble::tibble(
      TRT01A = c("Drug A", "Placebo", NA),
      subjects = c(1L, 1L, 0L),
      total = c(2L, 1L, 1L)
    )
  )
  # treatment-emergent subjects out of the safety population
  expect_identical(actual$incidence, c(subjects = 1L, total = 3L))
  expect_identical(actual$incidence_denom, "in the safety population")
  expect_snapshot(print(summary(adae, adsl = adsl)))

  # without a safety population, out of all ADSL subjects
  no_saffl <- summary(adae, adsl = dplyr::select(adsl, -SAFFL))$vs_adsl
  expect_null(no_saffl$saffl)
  expect_identical(no_saffl$incidence, c(subjects = 1L, total = 4L))
  expect_identical(no_saffl$incidence_denom, "in ADSL")

  # the incidence is for occurrence data only
  advs <- admiral_df(
    tibble::tibble(STUDYID = "P", USUBJID = "1", PARAMCD = "SYSBP", AVAL = 120)
  )
  expect_null(summary(advs, adsl = adsl)$vs_adsl$incidence)
  expect_snapshot(print(summary(advs, adsl = adsl)))

  # only the first and last few of the subjects not in ADSL are printed
  many_orphans <- admiral_df(tibble::tibble(
    STUDYID = "P", USUBJID = as.character(11:20), AEDECOD = "NAUSEA"
  ))
  expect_identical(summary(many_orphans, adsl = adsl)$vs_adsl$n_orphans, 10L)
  expect_snapshot(print(summary(many_orphans, adsl = adsl)))
})

## Test 47: an ADSL which cannot be compared with is an error ----
test_that("summary.admiral_df Test 47: an ADSL which cannot be compared with is an error", {
  adsl <- tibble::tibble(STUDYID = "P", USUBJID = c("1", "2"))
  adae <- admiral_df(tibble::tibble(STUDYID = "P", USUBJID = "1", AEDECOD = "NAUSEA"))

  expect_snapshot(summary(adae, adsl = adsl[c(1, 1, 2), ]), error = TRUE)
  expect_error(
    summary(adae, adsl = dplyr::select(adsl, -USUBJID)),
    class = "assert_data_frame"
  )
  expect_error(summary(adae, adsl = "ADSL"), class = "assert_data_frame")

  # a dataset which cannot be compared is only a warning, as the rest of the
  # summary is valid
  expect_snapshot(actual <- summary(admiral_df(adsl), adsl = adsl))
  expect_null(actual$vs_adsl)
  lookup <- admiral_df(tibble::tibble(PARAMCD = "SYSBP"))
  expect_warning(
    no_subjects <- summary(lookup, adsl = adsl),
    regexp = "no variables to match subjects by"
  )
  expect_null(no_subjects$vs_adsl)

  # a grouped ADSL is compared as a whole
  expect_identical(
    summary(adae, adsl = dplyr::group_by(adsl, STUDYID))$vs_adsl$n_common,
    1L
  )
})

## Test 48: a record without a subject is not counted as a subject ----
test_that("summary.admiral_df Test 48: a record without a subject is not counted as a subject", { # nolint
  # e.g. records added by a derivation which does not populate `USUBJID`; a
  # blank counts as missing, as in data read from SAS
  advs <- tibble::tibble(
    STUDYID = "P",
    USUBJID = c("1", "2", NA, ""),
    PARAMCD = "SYSBP",
    AVAL = c(120, 130, 125, 128)
  ) %>%
    admiral_df()

  actual <- summary(advs)
  expect_identical(actual$n_subjects, 2L)
  expect_identical(actual$details$params$subjects, 2L)

  adsl <- tibble::tibble(
    STUDYID = "P",
    USUBJID = c("1", "2", NA),
    SAFFL = "Y",
    TRT01P = "Placebo",
    AGE = c(40, 50, 99)
  ) %>%
    admiral_df(dataset_class = "ADSL")

  adsl_details <- summary(adsl)$details
  expect_identical(adsl_details$trt$subjects, 2L)
  expect_identical(adsl_details$pops, c(SAFFL = 2L))
  expect_identical(adsl_details$age, c(median = 45, min = 40, max = 50))

  # the subject keys are used without `USUBJID`, numeric ones too
  subject_keys <- get_admiral_option("subject_keys")
  withr::defer(set_admiral_options(subject_keys = subject_keys))
  set_admiral_options(subject_keys = exprs(STUDYID, SUBJID))
  by_subjid <- admiral_df(tibble::tibble(STUDYID = "P", SUBJID = c(1, NA, 2), AGE = 1:3))
  expect_identical(summary(by_subjid)$n_subjects, 2L)
})

## Test 49: the checks leave out records without a subject ----
test_that("summary.admiral_df Test 49: the checks leave out records without a subject", {
  # two baselines without a subject belong to no subject, so they are not a
  # second baseline of anyone
  advs <- tibble::tibble(
    STUDYID = "P",
    USUBJID = c("1", NA, NA, "", ""),
    PARAMCD = "SYSBP",
    ABLFL = "Y",
    AVAL = c(120, 125, 126, 127, 128)
  ) %>%
    admiral_df()
  expect_identical(summary(advs)$details$multiple_baselines, 0L)

  adae <- tibble::tibble(
    STUDYID = "P",
    USUBJID = c("1", NA, NA, "", ""),
    AEDECOD = "NAUSEA",
    AOCCFL = "Y"
  )
  expect_identical(
    check_occurrence_flags(adae, "USUBJID")$occ_flag_dups,
    c(AOCCFL = 0L)
  )

  # duplicates of a real subject are still reported
  expect_identical(
    check_occurrence_flags(dplyr::bind_rows(adae, adae[1, ]), "USUBJID")$occ_flag_dups,
    c(AOCCFL = 1L)
  )
})
