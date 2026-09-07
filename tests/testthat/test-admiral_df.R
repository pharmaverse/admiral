# admiral test guidelines loaded

# helper: load the `{metacore}` pilot ADaM example spec, optionally selecting
# a single dataset. The example spec itself triggers informational metacore
# warnings about empty columns; they are not the subject of these tests.
load_pilot_adam_spec <- function() {
  spec_env <- new.env()
  load(metacore::metacore_example("pilot_ADaM.rda"), envir = spec_env)
  spec_env[[ls(spec_env)[1]]]
}

load_pilot_adam_dataset_spec <- function(dataset_name) {
  suppressWarnings(
    suppressMessages(metacore::select_dataset(load_pilot_adam_spec(), dataset_name))
  )
}

# as_admiral_df ----
## Test 1: adds the admiral_df class while preserving tibble classes ----
test_that("as_admiral_df Test 1: adds the admiral_df class while preserving tibble classes", {
  input <- tibble::tibble(USUBJID = c("1", "2"), AVAL = c(1, 2))

  result <- as_admiral_df(input)

  expect_s3_class(result, "admiral_df")
  expect_s3_class(result, "tbl_df")
  expect_identical(as.data.frame(result), as.data.frame(input))
})

## Test 2: is idempotent and returns NULL unchanged ----
test_that("as_admiral_df Test 2: is idempotent and returns NULL unchanged", {
  input <- as_admiral_df(tibble::tibble(USUBJID = "1"))

  expect_identical(class(as_admiral_df(input)), class(input))
  expect_null(as_admiral_df(NULL))
})

## Test 3: stores the keys attribute when supplied ----
test_that("as_admiral_df Test 3: stores the keys attribute when supplied", {
  input <- tibble::tibble(USUBJID = c("1", "2"), PARAMCD = c("A", "B"))

  result <- as_admiral_df(input, keys = c("USUBJID", "PARAMCD"))

  expect_identical(attr(result, "admiral_keys"), c("USUBJID", "PARAMCD"))

  # a NULL `keys` leaves an existing attribute untouched
  untouched <- as_admiral_df(result)
  expect_identical(attr(untouched, "admiral_keys"), c("USUBJID", "PARAMCD"))
})

# get_admiral_keys ----
## Test 4: extracts the keys defined for a single dataset ----
test_that("get_admiral_keys Test 4: extracts the keys defined for a single dataset", {
  skip_if_not_installed("metacore")

  adsl_spec <- load_pilot_adam_dataset_spec("ADSL")
  adae_spec <- load_pilot_adam_dataset_spec("ADAE")

  adsl_keys <- get_admiral_keys(adsl_spec)
  adae_keys <- get_admiral_keys(adae_spec)

  expect_identical(adsl_keys, "USUBJID")
  expect_identical(adae_keys, c("USUBJID", "AETERM", "ASTDT", "AESEQ"))
})

## Test 5: dataset_name is required when the spec has more than one dataset ----
test_that("get_admiral_keys Test 5: dataset_name is required when the spec has more than one dataset", { # nolint
  skip_if_not_installed("metacore")

  spec <- load_pilot_adam_spec()

  expect_error(get_admiral_keys(spec), regexp = "dataset_name")
  expect_identical(get_admiral_keys(spec, "ADSL"), get_admiral_keys(load_pilot_adam_dataset_spec("ADSL")))
})

## Test 6: warns and returns character(0) when no keys are defined ----
test_that("get_admiral_keys Test 6: warns and returns character(0) when no keys are defined", {
  skip_if_not_installed("metacore")

  spec <- load_pilot_adam_dataset_spec("ADSL")
  # simulate a spec with no keys defined for the dataset without mutating the
  # read-only `ds_vars` table directly
  local_mocked_bindings(
    get_keys = function(...) tibble::tibble(variable = character(0)),
    .package = "metacore"
  )

  expect_warning(
    result <- get_admiral_keys(spec),
    regexp = "No key variables"
  )
  expect_identical(result, character(0))
})

## Test 7: errors for a missing/invalid metacore object ----
test_that("get_admiral_keys Test 7: errors for a missing/invalid metacore object", {
  expect_error(get_admiral_keys(list()), regexp = "Metacore")
  expect_error(get_admiral_keys(NULL), regexp = "Metacore")
})

# set_admiral_keys ----
## Test 8: tags the dataset with keys, dataset name, and the admiral_df class ----
test_that("set_admiral_keys Test 8: tags the dataset with keys, dataset name, and the admiral_df class", { # nolint
  skip_if_not_installed("metacore")

  adsl_spec <- load_pilot_adam_dataset_spec("ADSL")
  input <- tibble::tibble(USUBJID = c("1", "2"))

  result <- set_admiral_keys(input, adsl_spec)

  expect_s3_class(result, "admiral_df")
  expect_identical(attr(result, "admiral_keys"), get_admiral_keys(adsl_spec))
  expect_identical(attr(result, "admiral_ds_name"), "ADSL")
})

## Test 9: dataset_name is passed through to get_admiral_keys() ----
test_that("set_admiral_keys Test 9: dataset_name is passed through to get_admiral_keys()", {
  skip_if_not_installed("metacore")

  spec <- load_pilot_adam_spec()
  input <- tibble::tibble(USUBJID = c("1", "2"))

  result <- set_admiral_keys(input, spec, "ADSL")

  expect_identical(attr(result, "admiral_keys"), get_admiral_keys(spec, "ADSL"))
  expect_identical(attr(result, "admiral_ds_name"), "ADSL")
})

## Test 10: errors for a non-data-frame dataset ----
test_that("set_admiral_keys Test 10: errors for a non-data-frame dataset", {
  skip_if_not_installed("metacore")

  adsl_spec <- load_pilot_adam_dataset_spec("ADSL")

  expect_error(set_admiral_keys(list(), adsl_spec), regexp = "data.frame")
})
