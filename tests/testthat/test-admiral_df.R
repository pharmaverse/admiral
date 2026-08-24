# admiral test guidelines loaded

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

# derive_param_computed ----
## Test 3: output is tagged with the admiral_df class ----
test_that("derive_param_computed Test 3: output is tagged with the admiral_df class", {
  input <- tibble::tribble(
    ~USUBJID,      ~PARAMCD, ~AVAL, ~VISIT,
    "01-701-1015", "DIABP",     51, "BASELINE",
    "01-701-1015", "SYSBP",    121, "BASELINE",
    "01-701-1028", "DIABP",     79, "BASELINE",
    "01-701-1028", "SYSBP",    130, "BASELINE"
  )

  result <- derive_param_computed(
    input,
    by_vars = exprs(USUBJID, VISIT),
    parameters = c("SYSBP", "DIABP"),
    set_values_to = exprs(
      AVAL = (AVAL.SYSBP + 2 * AVAL.DIABP) / 3,
      PARAMCD = "MAP"
    )
  )

  expect_s3_class(result, "admiral_df")
  expect_s3_class(result, "tbl_df")
})

# get_admiral_df_type ----
## Test 4: classifies the common ADaM structures ----
test_that("get_admiral_df_type Test 4: classifies the common ADaM structures", {
  bds <- tibble::tibble(
    USUBJID = "1", PARAMCD = "MAP", AVAL = 90, AVISIT = "BASELINE"
  )
  tte <- tibble::tibble(
    USUBJID = "1", PARAMCD = "OS", AVAL = 100, CNSR = 0,
    STARTDT = as.Date("2020-01-01")
  )
  occds <- tibble::tibble(
    USUBJID = "1", AEDECOD = "HEADACHE", TRTEMFL = "Y"
  )
  adsl <- tibble::tibble(
    STUDYID = "P", USUBJID = c("1", "2"), TRT01P = c("A", "B")
  )
  other <- tibble::tibble(FOO = 1, BAR = 2)

  expect_identical(get_admiral_df_type(bds), "BDS")
  expect_identical(get_admiral_df_type(tte), "TTE")
  expect_identical(get_admiral_df_type(occds), "OCCDS")
  expect_identical(get_admiral_df_type(adsl), "ADSL")
  expect_identical(get_admiral_df_type(other), "other")
})

# summary.admiral_df ----
## Test 5: returns the expected diagnostics for a BDS dataset ----
test_that("summary.admiral_df Test 5: returns the expected diagnostics for a BDS dataset", {
  input <- as_admiral_df(tibble::tribble(
    ~STUDYID,  ~USUBJID, ~PARAMCD, ~AVAL, ~AVISIT,
    "PILOT01", "1",      "DIABP",     51, "BASELINE",
    "PILOT01", "1",      "SYSBP",    121, "BASELINE",
    "PILOT01", "1",      "MAP",       74, "BASELINE",
    "PILOT01", "2",      "DIABP",     79, "BASELINE",
    "PILOT01", "2",      "SYSBP",    130, "BASELINE",
    "PILOT01", "2",      "MAP",       96, "WEEK 2"
  ))

  result <- summary(input)

  expect_s3_class(result, "summary_admiral_df")
  expect_identical(result$type, "BDS")
  expect_identical(result$n_obs, 6L)
  expect_identical(result$n_subjects, 2L)
  expect_identical(result$params$PARAMCD, c("DIABP", "MAP", "SYSBP"))
  expect_identical(result$avisits, c("BASELINE", "WEEK 2"))
})

## Test 6: counts subjects by USUBJID when subject keys are partly NA ----
test_that("summary.admiral_df Test 6: counts subjects by USUBJID when subject keys are partly NA", {
  # records added by derivations may leave STUDYID as NA; subject count must
  # still be driven by USUBJID
  input <- as_admiral_df(tibble::tribble(
    ~STUDYID,        ~USUBJID, ~PARAMCD, ~AVAL,
    "PILOT01",       "1",      "SYSBP",    121,
    NA_character_,   "1",      "MAP",       74,
    "PILOT01",       "2",      "SYSBP",    130
  ))

  expect_identical(summary(input)$n_subjects, 2L)
})

## Test 7: formatted output is stable ----
test_that("summary.admiral_df Test 7: formatted output is stable", {
  input <- as_admiral_df(tibble::tribble(
    ~STUDYID,  ~USUBJID, ~PARAMCD, ~AVAL, ~AVISIT,
    "PILOT01", "1",      "DIABP",     51, "BASELINE",
    "PILOT01", "1",      "SYSBP",    121, "BASELINE",
    "PILOT01", "2",      "DIABP",     79, "BASELINE"
  ))

  expect_snapshot(print(summary(input)))
})
