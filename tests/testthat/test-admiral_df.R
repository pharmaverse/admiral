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
  input <- set_admiral_attrs(
    tibble::tribble(
      ~USUBJID, ~AVAL,
      "1",      10
    ),
    keys = exprs(USUBJID),
    dataset_name = "ADSL",
    dataset_type = "ADSL"
  )

  # every derivation tags its output, which must not lose the metadata the user
  # set on its input
  actual <- new_admiral_df(input)

  expect_identical(attr(actual, "admiral_keys"), "USUBJID")
  expect_identical(attr(actual, "admiral_ds_name"), "ADSL")
  expect_identical(attr(actual, "admiral_ds_type"), "ADSL")
})

# set_admiral_attrs ----
## Test 5: stores the keys, name, type, and the admiral_df class ----
test_that("set_admiral_attrs Test 5: stores the keys, name, type, and the admiral_df class", {
  input <- tibble::tribble(
    ~USUBJID, ~PARAMCD, ~AVISIT,    ~AVAL,
    "1",      "DIABP",  "BASELINE",    51,
    "1",      "SYSBP",  "BASELINE",   121
  )

  actual <- set_admiral_attrs(
    input,
    keys = exprs(USUBJID, PARAMCD, AVISIT),
    dataset_name = "ADVS",
    dataset_type = "BDS"
  )

  expect_equal(
    attr(actual, "admiral_keys"),
    c("USUBJID", "PARAMCD", "AVISIT")
  )
  expect_equal(attr(actual, "admiral_ds_name"), "ADVS")
  expect_equal(attr(actual, "admiral_ds_type"), "BDS")
  expect_s3_class(actual, "admiral_df")
  # Test that the data itself is not changed
  stripped <- actual
  attr(stripped, "admiral_keys") <- NULL
  attr(stripped, "admiral_ds_name") <- NULL
  attr(stripped, "admiral_ds_type") <- NULL
  class(stripped) <- class(input)
  expect_equal(stripped, input)
})

## Test 6: keys from a specification are passed via chr2vars() ----
test_that("set_admiral_attrs Test 6: keys from a specification are passed via chr2vars()", {
  input <- tibble::tribble(
    ~USUBJID, ~PARAMCD, ~AVAL,
    "1",      "DIABP",     51
  )
  # keys extracted from a specification table carry a variable label
  spec_keys <- c("USUBJID", "PARAMCD")
  attr(spec_keys, "label") <- "Variable Name"

  from_exprs <- set_admiral_attrs(input, keys = exprs(USUBJID, PARAMCD))
  from_spec <- set_admiral_attrs(input, keys = chr2vars(spec_keys))

  expect_identical(attr(from_spec, "admiral_keys"), attr(from_exprs, "admiral_keys"))
})

## Test 7: the metadata is stored, preserved, and replaced ----
test_that("set_admiral_attrs Test 7: the metadata is stored, preserved, and replaced", {
  input <- tibble::tribble(
    ~USUBJID, ~PARAMCD, ~AVAL,
    "1",      "DIABP",     51
  )

  # nothing is stored when there is nothing to store
  bare <- set_admiral_attrs(input)
  expect_s3_class(bare, "admiral_df")
  expect_null(attr(bare, "admiral_keys"))
  expect_null(attr(bare, "admiral_ds_name"))
  expect_null(attr(bare, "admiral_ds_type"))

  # an argument which is not repeated is not lost by a call which only revises
  # another
  named <- set_admiral_attrs(
    input,
    keys = exprs(USUBJID),
    dataset_name = "ADVS",
    dataset_type = "BDS"
  )
  rekeyed <- set_admiral_attrs(named, keys = exprs(USUBJID, PARAMCD))
  expect_equal(attr(rekeyed, "admiral_keys"), c("USUBJID", "PARAMCD"))
  expect_equal(attr(rekeyed, "admiral_ds_name"), "ADVS")
  expect_equal(attr(rekeyed, "admiral_ds_type"), "BDS")

  retyped <- set_admiral_attrs(named, dataset_type = "TTE")
  expect_equal(attr(retyped, "admiral_keys"), "USUBJID")
  expect_equal(attr(retyped, "admiral_ds_type"), "TTE")

  # all give way when a later call supplies them
  renamed <- set_admiral_attrs(
    named,
    keys = exprs(USUBJID, PARAMCD),
    dataset_name = "ADSL",
    dataset_type = "ADSL"
  )
  expect_equal(attr(renamed, "admiral_keys"), c("USUBJID", "PARAMCD"))
  expect_equal(attr(renamed, "admiral_ds_name"), "ADSL")
  expect_equal(attr(renamed, "admiral_ds_type"), "ADSL")
})

## Test 8: the keys are stored as a plain character vector ----
test_that("set_admiral_attrs Test 8: the keys are stored as a plain character vector", {
  input <- tibble::tribble(
    ~USUBJID, ~AVAL,
    "1",      10
  )

  # `vars2chr()` returns a named vector, which must not end up on the stored
  # attribute
  actual <- set_admiral_attrs(input, keys = exprs(USUBJID, AVAL))

  expect_identical(attr(actual, "admiral_keys"), c("USUBJID", "AVAL"))
})

## Test 9: empty keys are stored, keys not in the dataset are not ----
test_that("set_admiral_attrs Test 9: empty keys are stored, keys not in the dataset are not", {
  input <- tibble::tribble(
    ~USUBJID, ~PARAMCD, ~AVAL,
    "1",      "DIABP",     51
  )

  # the empty attribute is a statement that the dataset has no key variables,
  # so it must be distinguishable from an attribute which was never set
  expect_silent(empty <- set_admiral_attrs(input, keys = exprs()))
  expect_identical(attr(empty, "admiral_keys"), character(0))
  expect_false(is.null(attr(empty, "admiral_keys")))
  expect_s3_class(empty, "admiral_df")

  # a key which is not in the dataset is most likely a typo, so the keys are
  # not set and a warning says so; the name and the class are still applied
  expect_snapshot(
    typo <- set_admiral_attrs(input, keys = exprs(USUBJID, PARMACD), dataset_name = "ADVS")
  )
  expect_null(attr(typo, "admiral_keys"))
  expect_identical(attr(typo, "admiral_ds_name"), "ADVS")
  expect_s3_class(typo, "admiral_df")

  # re-keying with a typo keeps the keys set previously, and the warning says
  # which they are
  keyed <- set_admiral_attrs(input, keys = exprs(USUBJID, PARAMCD))
  expect_snapshot(
    rekeyed <- set_admiral_attrs(keyed, keys = exprs(USUBJID, PARMACD, AVISIT))
  )
  expect_identical(attr(rekeyed, "admiral_keys"), c("USUBJID", "PARAMCD"))
})

## Test 10: the metadata and the class survive a dplyr pipeline ----
test_that("set_admiral_attrs Test 10: the metadata and the class survive a dplyr pipeline", {
  input <- set_admiral_attrs(
    tibble::tribble(
      ~USUBJID, ~PARAMCD, ~AVAL,
      "1",      "DIABP",     51,
      "2",      "DIABP",     79
    ),
    keys = exprs(USUBJID, PARAMCD),
    dataset_name = "ADVS",
    dataset_type = "BDS"
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

  expect_equal(attr(actual, "admiral_keys"), c("USUBJID", "PARAMCD"))
  expect_equal(attr(actual, "admiral_ds_name"), "ADVS")
  expect_equal(attr(actual, "admiral_ds_type"), "BDS")
  expect_s3_class(actual, "admiral_df")
})

## Test 11: an error is issued for invalid arguments ----
test_that("set_admiral_attrs Test 11: an error is issued for invalid arguments", {
  input <- tibble::tribble(
    ~USUBJID, ~AVAL,
    "1",      10
  )

  expect_error(set_admiral_attrs(list(), keys = exprs(USUBJID)))
  expect_error(set_admiral_attrs(input, keys = 1))
  expect_error(set_admiral_attrs(input, keys = iris))
  expect_error(set_admiral_attrs(input, keys = exprs(USUBJID + 1)))
  expect_error(set_admiral_attrs(input, keys = exprs(USUBJID), dataset_name = c("A", "B")))
  expect_error(set_admiral_attrs(input, dataset_type = c("BDS", "TTE")))

  # keys taken from a specification are a character vector, so the error
  # points to the conversion
  expect_snapshot(set_admiral_attrs(input, keys = c("USUBJID", "PARAMCD")), error = TRUE)

  # the type is one of the ADaM dataset types
  expect_snapshot(set_admiral_attrs(input, dataset_type = "ADAE"), error = TRUE)
})

## Test 12: device types are accepted, and checked for SPDEVID ----
test_that("set_admiral_attrs Test 12: device types are accepted, and checked for SPDEVID", {
  input <- tibble::tribble(
    ~STUDYID, ~SPDEVID, ~PARAMCD, ~AVAL,
    "P",      "D-1",    "BATTV",    3.1,
    "P",      "D-2",    "BATTV",    2.9
  )

  for (type in c("ADDL", "MDBDS", "MDOCCDS")) {
    expect_silent(actual <- set_admiral_attrs(input, dataset_type = type))
    expect_identical(attr(actual, "admiral_ds_type"), type)
  }

  # a device type for a dataset without the device identifier is most likely
  # the wrong type, so a warning says so; the type is still stored, as the
  # declaration is the user's to make
  no_device <- select(input, -SPDEVID)
  expect_snapshot(actual <- set_admiral_attrs(no_device, dataset_type = "MDBDS"))
  expect_identical(attr(actual, "admiral_ds_type"), "MDBDS")

  # a non-device type is not checked for the device identifier
  expect_silent(set_admiral_attrs(no_device, dataset_type = "BDS"))
})

# print.admiral_df ----
## Test 13: the name, type and keys are printed in a header ----
test_that("print.admiral_df Test 13: the name, type and keys are printed in a header", {
  input <- tibble::tribble(
    ~STUDYID, ~USUBJID, ~PARAMCD, ~AVISIT,    ~AVAL,
    "P",      "1",      "DIABP",  "BASELINE",    51,
    "P",      "1",      "SYSBP",  "BASELINE",   121
  )

  expect_snapshot(
    print(set_admiral_attrs(
      input,
      keys = exprs(USUBJID, PARAMCD, AVISIT),
      dataset_name = "ADVS",
      dataset_type = "BDS"
    ))
  )

  # each attribute is shown on its own
  expect_snapshot(print(set_admiral_attrs(input, keys = exprs(USUBJID, PARAMCD, AVISIT))))
  expect_snapshot(print(set_admiral_attrs(input, dataset_type = "BDS")))

  # empty keys record that no keys are defined, so they are shown as such
  expect_snapshot(print(set_admiral_attrs(input, keys = exprs(), dataset_name = "ADVS")))

  # a plain data frame prints with its own method below the header
  expect_snapshot(
    print(set_admiral_attrs(as.data.frame(input), keys = exprs(USUBJID), dataset_name = "ADVS"))
  )
})

## Test 14: a dataset with only the class prints as before ----
test_that("print.admiral_df Test 14: a dataset with only the class prints as before", {
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
