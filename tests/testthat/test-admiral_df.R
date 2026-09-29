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
