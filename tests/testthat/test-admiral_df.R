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

## Test 7: warns if declared keys are no longer in the dataset ----
test_that("summary.admiral_df Test 7: warns if declared keys are no longer in the dataset", {
  input <- as_admiral_df(
    tibble::tribble(
      ~STUDYID,  ~USUBJID, ~PARAMCD, ~AVAL,
      "PILOT01", "1",      "SYSBP",    121,
      "PILOT01", "2",      "SYSBP",    130
    ),
    keys = c("USUBJID", "PARAMCD", "AVISIT")
  )

  # `AVISIT` is not in the dataset, the structure must be checked on the rest
  expect_warning(
    result <- summary(input),
    regexp = "AVISIT"
  )
  expect_identical(result$keys, c("USUBJID", "PARAMCD"))
  expect_identical(result$key_source, "declared")

  # keys surviving `rename()` are stale, i.e. no key is left
  renamed <- dplyr::rename(input, SUBJID = USUBJID)
  expect_warning(
    result_renamed <- summary(dplyr::select(renamed, -PARAMCD)),
    regexp = "USUBJID"
  )
  expect_null(result_renamed$keys)

  # inferred keys must not trigger the warning
  expect_no_warning(summary(as_admiral_df(tibble::tribble(
    ~STUDYID,  ~USUBJID, ~PARAMCD, ~AVAL,
    "PILOT01", "1",      "SYSBP",    121,
    "PILOT01", "2",      "SYSBP",    130
  ))))
})

## Test 8: formatted output is stable ----
test_that("summary.admiral_df Test 8: formatted output is stable", {
  input <- as_admiral_df(tibble::tribble(
    ~STUDYID,  ~USUBJID, ~PARAMCD, ~AVAL, ~AVISIT,
    "PILOT01", "1",      "DIABP",     51, "BASELINE",
    "PILOT01", "1",      "SYSBP",    121, "BASELINE",
    "PILOT01", "2",      "DIABP",     79, "BASELINE"
  ))

  expect_snapshot(print(summary(input)))
})

# check_admiral_df ----
## Test 9: reports a fully tagged dataset ----
test_that("check_admiral_df Test 9: reports a fully tagged dataset", {
  input <- as_admiral_df(
    tibble::tribble(
      ~STUDYID,  ~USUBJID, ~PARAMCD, ~AVAL,
      "PILOT01", "1",      "SYSBP",    121,
      "PILOT01", "2",      "SYSBP",    130
    ),
    keys = c("USUBJID", "PARAMCD")
  )

  result <- check_admiral_df(input)

  expect_s3_class(result, "admiral_df_check")
  expect_true(result$is_admiral_df)
  expect_identical(result$keys, c("USUBJID", "PARAMCD"))
  expect_identical(result$key_source, "declared")
  expect_identical(result$stale_keys, character(0))
  expect_true(result$structure_checked)
})

## Test 10: distinguishes absent from empty keys ----
test_that("check_admiral_df Test 10: distinguishes absent from empty keys", {
  input <- as_admiral_df(tibble::tribble(
    ~STUDYID,  ~USUBJID, ~PARAMCD, ~AVAL,
    "PILOT01", "1",      "SYSBP",    121,
    "PILOT01", "2",      "SYSBP",    130
  ))

  # no attribute at all: the keys are inferred, so the structure is checked
  absent <- check_admiral_df(input)
  expect_null(absent$declared_keys)
  expect_identical(absent$key_source, "inferred")
  expect_identical(absent$keys, c("USUBJID", "PARAMCD"))
  expect_true(absent$structure_checked)

  # attribute set but empty (e.g. a spec without key variables): the inference
  # fallback is skipped and no key is left, so nothing is checked at all
  empty <- input
  attr(empty, "admiral_keys") <- character(0)
  empty <- check_admiral_df(empty)
  expect_identical(empty$declared_keys, character(0))
  expect_identical(empty$key_source, "none")
  expect_false(empty$structure_checked)
})

## Test 11: reports stale keys and an orphaned attribute ----
test_that("check_admiral_df Test 11: reports stale keys and an orphaned attribute", {
  input <- as_admiral_df(
    tibble::tribble(
      ~STUDYID,  ~USUBJID, ~PARAMCD, ~AVAL,
      "PILOT01", "1",      "SYSBP",    121,
      "PILOT01", "2",      "SYSBP",    130
    ),
    keys = c("USUBJID", "PARAMCD")
  )

  # partially stale: the structure is checked on the remaining key only
  partial <- check_admiral_df(dplyr::rename(input, SUBJID = USUBJID))
  expect_identical(partial$stale_keys, "USUBJID")
  expect_identical(partial$keys, "PARAMCD")
  expect_identical(partial$key_source, "declared")
  expect_true(partial$structure_checked)

  # fully stale: no declared key is left, and nothing is inferred instead
  full <- check_admiral_df(
    dplyr::select(dplyr::rename(input, SUBJID = USUBJID), -PARAMCD)
  )
  expect_identical(full$stale_keys, c("USUBJID", "PARAMCD"))
  expect_identical(full$keys, character(0))
  expect_false(full$structure_checked)

  # class dropped but attribute kept: `summary()` would not dispatch to
  # `summary.admiral_df()`, so nothing is checked despite valid keys
  orphaned <- check_admiral_df(tibble::as_tibble(input))
  expect_false(orphaned$is_admiral_df)
  expect_identical(orphaned$declared_keys, c("USUBJID", "PARAMCD"))
  expect_false(orphaned$structure_checked)
})

## Test 12: does not warn when no keys can be inferred ----
test_that("check_admiral_df Test 12: does not warn when no keys can be inferred", {
  # `infer_admiral_keys()` warns for OCCDS without `ASEQ`; the reporter must
  # stay silent and report the empty result instead
  input <- as_admiral_df(tibble::tribble(
    ~STUDYID,  ~USUBJID, ~AEDECOD,
    "PILOT01", "1",      "HEADACHE"
  ))

  expect_no_warning(result <- check_admiral_df(input))
  expect_identical(result$keys, character(0))
  expect_identical(result$key_source, "none")
  expect_false(result$structure_checked)
})

# print.admiral_df_check ----
## Test 13: formatted output is stable ----
test_that("print.admiral_df_check Test 13: formatted output is stable", {
  input <- as_admiral_df(
    tibble::tribble(
      ~STUDYID,  ~USUBJID, ~PARAMCD, ~AVAL,
      "PILOT01", "1",      "SYSBP",    121,
      "PILOT01", "2",      "SYSBP",    130
    ),
    keys = c("USUBJID", "PARAMCD")
  )

  expect_snapshot(print(check_admiral_df(input)))
  expect_snapshot(print(check_admiral_df(dplyr::rename(input, SUBJID = USUBJID))))
  expect_snapshot(print(check_admiral_df(tibble::as_tibble(input))))
})

# summarize_adsl ----
## Test 14: treats all-NA variables as not yet derived ----
test_that("summarize_adsl Test 14: treats all-NA variables as not yet derived", {
  input <- as_admiral_df(tibble::tribble(
    ~STUDYID,  ~USUBJID, ~TRT01P,       ~TRT01A,   ~EOSSTT,       ~SEX,          ~AGE,
    "PILOT01", "1",      NA_character_, "Placebo", NA_character_, NA_character_, NA_real_,
    "PILOT01", "2",      NA_character_, "Drug A",  NA_character_, NA_character_, NA_real_
  ))

  result <- summary(input)$adsl

  # an all-NA `TRT01P` is skipped in favour of the populated `TRT01A`
  expect_identical(result$trt_var, "TRT01A")
  expect_identical(result$trt$n, c(1L, 1L))
  expect_null(result$eosstt)
  expect_null(result$sex)
  expect_null(result$age)

  # no populated treatment candidate at all: no treatment element
  all_na <- as_admiral_df(tibble::tribble(
    ~STUDYID,  ~USUBJID, ~TRT01P,       ~TRTSDT,
    "PILOT01", "1",      NA_character_, as.Date("2024-01-01")
  ))
  expect_null(summary(all_na)$adsl$trt_var)
})

## Test 15: the subject-level check lines name what ran ----
test_that("print.summary_admiral_df Test 15: the subject-level check lines name what ran", {
  clean <- as_admiral_df(tibble::tribble(
    ~STUDYID,  ~USUBJID, ~ARM,      ~ACTARM,   ~SAFFL, ~TRTSDT,
    "PILOT01", "1",      "Placebo", "Placebo", "Y",    as.Date("2024-01-01"),
    "PILOT01", "2",      "Drug A",  "Drug A",  "Y",    as.Date("2024-01-02")
  ))
  # both checks pass: one confirmation line naming each check
  expect_snapshot(print(summary(clean)))

  broken <- as_admiral_df(tibble::tribble(
    ~STUDYID,  ~USUBJID, ~ARM,      ~ACTARM,  ~SAFFL, ~TRTSDT,
    "PILOT01", "1",      "Placebo", "Drug A", "Y",    as.Date("2024-01-01"),
    "PILOT01", "2",      "Drug A",  "Drug A", "Y",    as.Date(NA)
  ))
  # both checks fail: one detailed bullet each, no confirmation line
  expect_snapshot(print(summary(broken)))

  # one fails, one passes: the passing check is still reported by name
  mixed <- as_admiral_df(tibble::tribble(
    ~STUDYID,  ~USUBJID, ~ARM,      ~ACTARM,  ~SAFFL, ~TRTSDT,
    "PILOT01", "1",      "Placebo", "Drug A", "Y",    as.Date("2024-01-01"),
    "PILOT01", "2",      "Drug A",  "Drug A", "Y",    as.Date("2024-01-02")
  ))
  expect_snapshot(print(summary(mixed)))
})

## Test 16: supplied keys override the attribute and inference ----
test_that("summary.admiral_df Test 16: supplied keys override the attribute and inference", {
  input <- as_admiral_df(
    tibble::tribble(
      ~STUDYID,  ~USUBJID, ~PARAMCD, ~AVISIT,    ~AVAL,
      "PILOT01", "1",      "SYSBP",  "BASELINE",   121,
      "PILOT01", "1",      "SYSBP",  "WEEK 2",     130
    ),
    keys = c("USUBJID", "PARAMCD", "AVISIT")
  )

  result <- summary(input, keys = c("USUBJID", "PARAMCD"))
  expect_identical(result$key_source, "supplied")
  expect_identical(result$keys, c("USUBJID", "PARAMCD"))
  expect_identical(result$n_duplicate_keys, 1L)

  # a zero-length vector skips the record structure check entirely
  expect_null(summary(input, keys = character(0))$keys)

  # supplied keys which are not in the dataset warn like stale declared keys
  expect_warning(
    result_stale <- summary(input, keys = c("USUBJID", "PARAMCD", "ATPT")),
    regexp = "supplied"
  )
  expect_identical(result_stale$keys, c("USUBJID", "PARAMCD"))

  # the dots must stay empty so a typoed option fails instead of being swallowed
  expect_error(
    summary(input, kyes = "USUBJID"),
    class = "rlib_error_dots_nonempty"
  )
})

## Test 17: keys can be read from a metacore specification ----
test_that("summary.admiral_df Test 17: keys can be read from a metacore specification", {
  skip_if_not_installed("metacore")

  spec_env <- new.env()
  load(metacore::metacore_example("pilot_ADaM.rda"), envir = spec_env)
  spec <- spec_env[[ls(spec_env)[1]]]
  # the example spec itself triggers informational metacore warnings about
  # empty columns; they are not the subject of this test
  suppressWarnings(
    suppressMessages(adsl_spec <- metacore::select_dataset(spec, "ADSL"))
  )
  spec_keys <- get_admiral_keys(adsl_spec)

  # a dataset containing exactly the spec keys, unique per row
  input <- tibble::as_tibble(setNames(
    lapply(spec_keys, function(k) c("a", "b")),
    spec_keys
  ))
  input$TRTSDT <- as.Date("2024-01-01")
  input <- as_admiral_df(input)

  result <- summary(input, keys = adsl_spec)
  expect_identical(result$key_source, "supplied")
  expect_identical(result$keys, spec_keys)
  expect_identical(result$n_duplicate_keys, 0L)

  # a specification describing more than one dataset needs the dataset selected
  expect_error(summary(input, keys = spec), regexp = "dataset_name")
})

## Test 18: reports subject flow, deaths, and the date/flag-domain checks ----
test_that("summarize_adsl Test 18: reports subject flow, deaths, and the date/flag-domain checks", {
  input <- as_admiral_df(tibble::tibble(
    STUDYID = "PILOT01",
    USUBJID = c("1", "2", "3", "4"),
    RANDDT = as.Date(c("2024-01-01", "2024-01-02", "2024-01-03", NA)),
    TRTSDT = as.Date(c("2024-01-05", "2024-01-06", "2024-01-07", NA)),
    # subject 2 ends treatment before starting it
    TRTEDT = as.Date(c("2024-02-01", "2024-01-01", "2024-02-03", NA)),
    EOSSTT = c("COMPLETED", "DISCONTINUED", "COMPLETED", NA),
    DTHFL = c(NA, "Y", NA, NA),
    DTHCAUS = c(NA, "ADVERSE EVENT", NA, NA),
    # subject 3 has a flag value outside the Y/N/NA domain
    SAFFL = c("Y", "Y", "Yes", NA)
  ))

  result <- summary(input)$adsl

  expect_identical(
    result$funnel,
    c(randomized = 3L, treated = 3L, completed = 2L)
  )
  expect_identical(result$n_dth, 1L)
  expect_identical(result$dthcaus$DTHCAUS, "ADVERSE EVENT")
  expect_identical(result$trt_end_before_start, 1L)
  expect_identical(result$flag_domain, 1L)
  expect_identical(result$flag_domain_vars, "SAFFL")
})

# summarize_bds ----
## Test 19: builds the per-parameter table and DTYPE breakdown ----
test_that("summarize_bds Test 19: builds the per-parameter table and DTYPE breakdown", {
  input <- as_admiral_df(tibble::tribble(
    ~USUBJID, ~PARAMCD, ~PARAM,         ~AVISIT,    ~AVISITN, ~AVAL, ~DTYPE,
    "1",      "DIABP",  "Diastolic BP", "BASELINE",        0,    60, NA,
    "1",      "DIABP",  "Diastolic BP", "WEEK 2",          2,    65, NA,
    "1",      "DIABP",  "Diastolic BP", "WEEK 4",          4,    NA, "LOCF",
    "2",      "DIABP",  "Diastolic BP", "BASELINE",        0,    70, NA,
    "1",      "SYSBP",  "Systolic BP",  "BASELINE",        0,   120, NA
  ))

  result <- summary(input)$bds

  expect_identical(
    result$params,
    tibble::tribble(
      ~PARAMCD, ~records, ~subjects, ~visits, ~missing, ~min, ~median, ~max,
      "DIABP",        4L,        2L,      3L,       1L,   60,      65,   70,
      "SYSBP",        1L,        1L,      1L,       0L,  120,     120,  120
    )
  )
  expect_identical(result$dtype$DTYPE, "LOCF")
  expect_identical(result$dtype$n, 1L)

  # consistent AVISIT/AVISITN and PARAM: the checks ran and passed
  expect_identical(result$avisit_mismatch, 0L)
  expect_identical(result$param_inconsistent, 0L)
  # no ABLFL, CHG, or BASE: those checks did not run
  expect_null(result$multiple_baselines)
  expect_null(result$chg_no_base)
})

## Test 20: the BDS checks find injected defects ----
test_that("summarize_bds Test 20: the BDS checks find injected defects", {
  input <- as_admiral_df(tibble::tribble(
    ~USUBJID, ~PARAMCD, ~PARAM,        ~AVISIT,      ~AVISITN, ~AVAL, ~ABLFL, ~BASE, ~CHG,
    "1",      "SYSBP",  "Systolic BP", "BASELINE",          0,   120, "Y",      120,   NA,
    # a second baseline for the same subject/parameter, and a second AVISIT
    # for AVISITN == 0
    "1",      "SYSBP",  "Systolic BP", "BASELINE 2",        0,   118, "Y",      120,   NA,
    # CHG without BASE, and a second PARAM label for the same PARAMCD
    "1",      "SYSBP",  "Sys BP",      "WEEK 2",            2,   130, NA,        NA,   10,
    "2",      "SYSBP",  "Systolic BP", "BASELINE",          0,   115, "Y",      115,   NA
  ))

  result <- summary(input)$bds

  expect_identical(result$multiple_baselines, 1L)
  expect_identical(result$chg_no_base, 1L)
  expect_identical(result$avisit_mismatch, 1L)
  expect_identical(result$param_inconsistent, 1L)
})

## Test 22: the dataset name is captured for the heading ----
test_that("summary.admiral_df Test 22: the dataset name is captured for the heading", {
  advs <- as_admiral_df(tibble::tribble(
    ~STUDYID,  ~USUBJID, ~PARAMCD, ~AVAL,
    "PILOT01", "1",      "SYSBP",    121
  ))

  # the variable name the summary was called on
  expect_identical(summary(advs)$ds_name, "advs")

  # a stored name (as set by `set_admiral_keys()`) wins over the variable name
  attr(advs, "admiral_ds_name") <- "ADVS"
  expect_identical(summary(advs)$ds_name, "ADVS")

  # no name is better than a wrong one: the magrittr placeholder and complex
  # expressions are not used
  advs2 <- as_admiral_df(tibble::tribble(
    ~STUDYID,  ~USUBJID, ~PARAMCD, ~AVAL,
    "PILOT01", "1",      "SYSBP",    121
  ))
  expect_null((advs2 %>% summary())$ds_name)
  expect_null(summary(as_admiral_df(tibble::as_tibble(advs2)))$ds_name)
})

## Test 21: BDS formatted output is stable ----
test_that("print.summary_admiral_df Test 21: BDS formatted output is stable", {
  input <- as_admiral_df(tibble::tribble(
    ~USUBJID, ~PARAMCD, ~PARAM,         ~AVISIT,    ~AVISITN, ~AVAL, ~DTYPE,
    "1",      "DIABP",  "Diastolic BP", "BASELINE",        0,    60, NA,
    "1",      "DIABP",  "Diastolic BP", "WEEK 2",          2,    65, NA,
    "1",      "DIABP",  "Diastolic BP", "WEEK 4",          4,    NA, "LOCF",
    "2",      "DIABP",  "Diastolic BP", "BASELINE",        0,    70, NA,
    "1",      "SYSBP",  "Systolic BP",  "BASELINE",        0,   120, NA
  ))

  expect_snapshot(print(summary(input)))
})

# summarize_occds ----
## Test 23: reports terms, treatment-emergent, severity, and serious facts ----
test_that("summarize_occds Test 23: reports terms, treatment-emergent, severity, and serious facts", { # nolint
  input <- as_admiral_df(tibble::tribble(
    ~USUBJID, ~ASEQ, ~AEBODSYS, ~AEDECOD,   ~AETERM,    ~TRTEMFL,      ~AESEV,     ~AESER, ~ASTDT,                 ~TRTSDT, # nolint
    "1",         1L, "NERV",    "HEADACHE", "Headache", "Y",           "MILD",     "N",    as.Date("2024-01-10"), as.Date("2024-01-01"), # nolint
    "1",         2L, "GI",      "NAUSEA",   "Nausea",   "Y",           "MODERATE", "Y",    as.Date("2024-01-12"), as.Date("2024-01-01"), # nolint
    "2",         1L, "NERV",    "HEADACHE", "headache", NA_character_, "MILD",     "N",    as.Date("2023-12-30"), as.Date("2024-01-01") # nolint
  ))

  result <- summary(input)$occds

  expect_identical(
    result$terms,
    c(AEBODSYS = 2L, AEDECOD = 2L, AETERM = 3L)
  )
  expect_identical(result$trtem, c(emergent = 2L, total = 3L))
  expect_identical(result$sev_var, "AESEV")
  expect_identical(result$sev$AESEV, c("MILD", "MODERATE"))
  expect_identical(result$sev$n, c(2L, 1L))
  expect_identical(result$n_serious, 1L)

  # the date checks ran and passed; the pre-treatment record is not flagged
  # treatment emergent, so it is not a defect
  expect_identical(result$missing_astdt, 0L)
  expect_identical(result$pre_trt_emergent, 0L)
  # no occurrence flags: that check did not run
  expect_null(result$occ_flag_dups)
})

## Test 24: the OCCDS checks find injected defects ----
test_that("summarize_occds Test 24: the OCCDS checks find injected defects", {
  input <- as_admiral_df(tibble::tribble(
    ~USUBJID, ~ASEQ, ~AEBODSYS, ~AEDECOD,    ~TRTEMFL, ~AOCCFL,       ~AOCCSFL,      ~ASTDT,                 ~TRTSDT, # nolint
    # a treatment-emergent record starting before treatment
    "1",         1L, "NERV",    "HEADACHE",  "Y",      "Y",           "Y",           as.Date("2023-12-30"), as.Date("2024-01-01"), # nolint
    # a second AOCCFL "Y" for the subject, a second AOCCSFL "Y" for the same
    # body system, and a missing start date
    "1",         2L, "NERV",    "DIZZINESS", "Y",      "Y",           "Y",           as.Date(NA),           as.Date("2024-01-01"), # nolint
    "2",         1L, "GI",      "NAUSEA",    "Y",      "Y",           "Y",           as.Date("2024-02-05"), as.Date("2024-02-01") # nolint
  ))

  result <- summary(input)$occds

  expect_identical(result$occ_flag_dups, 2L)
  expect_identical(result$occ_flag_vars, c("AOCCFL", "AOCCSFL"))
  expect_identical(result$missing_astdt, 1L)
  expect_identical(result$pre_trt_emergent, 1L)

  # a flag whose level variable is absent cannot be judged and is skipped:
  # AOCCPFL needs a --DECOD variable
  no_decod <- as_admiral_df(tibble::tribble(
    ~USUBJID, ~ASEQ, ~AETERM,    ~AOCCPFL,
    "1",         1L, "Headache", "Y",
    "1",         2L, "Nausea",   "Y"
  ))
  expect_null(summary(no_decod)$occds$occ_flag_dups)
})

# summarize_vs_adsl ----
## Test 26: reports coverage and denominators from ADSL ----
test_that("summarize_vs_adsl Test 26: reports coverage and denominators from ADSL", {
  adsl <- tibble::tribble(
    ~USUBJID, ~SAFFL, ~TRT01A,       ~TRTSDT,               ~DTHDT,
    "1",      "Y",    "Placebo",     as.Date("2024-01-01"), as.Date(NA),
    "2",      "Y",    "Active",      as.Date("2024-01-02"), as.Date("2024-03-01"),
    "3",      "N",    "Placebo",     as.Date(NA),           as.Date(NA)
  )
  input <- as_admiral_df(tibble::tribble(
    ~USUBJID, ~PARAMCD, ~AVAL, ~ADT,                  ~TRTSDT,
    "1",      "SYSBP",    121, as.Date("2024-01-10"), as.Date("2024-01-01"),
    "2",      "SYSBP",    130, as.Date("2024-02-01"), as.Date("2024-01-02")
  ))

  result <- summary(input, adsl = adsl)$vs_adsl

  expect_identical(result$n_adsl, 3L)
  expect_identical(result$n_common, 2L)
  expect_identical(result$saffl, c(covered = 2L, total = 2L))
  expect_identical(result$arm_var, "TRT01A")
  expect_identical(result$arm_coverage$TRT01A, c("Active", "Placebo"))
  expect_identical(result$arm_coverage$subjects, c(1L, 1L))
  expect_identical(result$arm_coverage$total, c(1L, 2L))

  # subject 3 having no records is legitimate absence, not a defect; the
  # shared TRTSDT matches and no record is dated after death
  expect_identical(result$n_orphans, 0L)
  expect_identical(result$stale_total, 0L)
  expect_identical(result$stale_vars$variable, character(0))
  expect_identical(result$after_death_var, "ADT")
  expect_identical(result$n_after_death, 0L)
})

## Test 27: the ADSL comparison checks find injected defects ----
test_that("summarize_vs_adsl Test 27: the ADSL comparison checks find injected defects", {
  adsl <- tibble::tribble(
    ~USUBJID, ~TRTSDT,               ~DTHDT,
    "1",      as.Date("2024-01-01"), as.Date("2024-02-01"),
    "2",      as.Date("2024-01-02"), as.Date(NA)
  )
  input <- as_admiral_df(tibble::tribble(
    ~USUBJID, ~PARAMCD, ~AVAL, ~ADT,                  ~TRTSDT,
    # TRTSDT disagrees with ADSL (stale merge) and ADT is after death
    "1",      "SYSBP",    121, as.Date("2024-02-10"), as.Date("2024-01-05"),
    "2",      "SYSBP",    130, as.Date("2024-01-20"), as.Date("2024-01-02"),
    # a subject with no ADSL entry at all
    "9",      "SYSBP",    115, as.Date("2024-01-15"), as.Date(NA)
  ))

  result <- summary(input, adsl = adsl)$vs_adsl

  expect_identical(result$n_orphans, 1L)
  expect_identical(result$orphans, "9")
  # the orphan's NA TRTSDT is not additionally counted as a mismatch
  expect_identical(result$stale_total, 1L)
  expect_identical(result$stale_vars$variable, "TRTSDT")
  expect_identical(result$stale_vars$subjects, 1L)
  expect_identical(result$n_after_death, 1L)
})

## Test 28: OCCDS gains the incidence and emergent-untreated items ----
test_that("summarize_vs_adsl Test 28: OCCDS gains the incidence and emergent-untreated items", { # nolint
  adsl <- tibble::tribble(
    ~USUBJID, ~SAFFL, ~TRTSDT,
    "1",      "Y",    as.Date("2024-01-01"),
    "2",      "Y",    as.Date(NA),
    "3",      "Y",    as.Date("2024-01-03"),
    "4",      "N",    as.Date(NA)
  )
  input <- as_admiral_df(tibble::tribble(
    ~USUBJID, ~ASEQ, ~AEDECOD,   ~TRTEMFL,
    "1",         1L, "HEADACHE", "Y",
    "1",         2L, "NAUSEA",   NA_character_,
    # emergent records for a subject ADSL says was never treated
    "2",         1L, "DIZZINESS", "Y"
  ))

  result <- summary(input, adsl = adsl)$vs_adsl

  expect_identical(result$incidence, c(subjects = 2L, total = 3L))
  expect_identical(result$incidence_denom, "safety population")
  expect_identical(result$n_emergent_untreated, 1L)

  # without SAFFL the denominator falls back to all of ADSL
  no_saffl <- summary(input, adsl = select(adsl, -SAFFL))$vs_adsl
  expect_identical(no_saffl$incidence, c(subjects = 2L, total = 4L))
  expect_identical(no_saffl$incidence_denom, "ADSL")
})

## Test 29: a malformed adsl errors and an ADSL-type object warns ----
test_that("summarize_vs_adsl Test 29: a malformed adsl errors and an ADSL-type object warns", { # nolint
  input <- as_admiral_df(tibble::tribble(
    ~USUBJID, ~PARAMCD, ~AVAL,
    "1",      "SYSBP",    121
  ))

  # duplicate subjects make every comparison ambiguous: caller error
  dup_adsl <- tibble::tibble(USUBJID = c("1", "1"), SAFFL = c("Y", "Y"))
  expect_error(
    summary(input, adsl = dup_adsl),
    regexp = "one record per"
  )

  # missing USUBJID in adsl: caller error
  expect_error(
    summary(input, adsl = tibble::tibble(SUBJID = "1")),
    class = "assert_data_frame"
  )

  # an ADSL compared against an ADSL is ignored with a warning
  adsl_input <- as_admiral_df(tibble::tibble(
    USUBJID = c("1", "2"), TRT01P = c("A", "B")
  ))
  expect_warning(
    result <- summary(adsl_input, adsl = tibble::tibble(USUBJID = "1")),
    regexp = "ignored"
  )
  expect_null(result$vs_adsl)
})

## Test 30: ADSL comparison formatted output is stable ----
test_that("print.summary_admiral_df Test 30: ADSL comparison formatted output is stable", { # nolint
  adsl <- tibble::tribble(
    ~USUBJID, ~SAFFL, ~TRT01A,   ~TRTSDT,               ~DTHDT,
    "1",      "Y",    "Placebo", as.Date("2024-01-01"), as.Date(NA),
    "2",      "Y",    "Active",  as.Date("2024-01-02"), as.Date("2024-02-01"),
    "3",      "N",    "Placebo", as.Date(NA),           as.Date(NA)
  )
  input <- as_admiral_df(tibble::tribble(
    ~USUBJID, ~PARAMCD, ~AVAL, ~ADT,                  ~TRTSDT,
    # TRTSDT stale vs ADSL
    "1",      "SYSBP",    121, as.Date("2024-01-10"), as.Date("2024-01-05"),
    # dated after death
    "2",      "SYSBP",    130, as.Date("2024-02-10"), as.Date("2024-01-02"),
    # orphan
    "9",      "SYSBP",    115, as.Date("2024-01-15"), as.Date(NA)
  ))

  expect_snapshot(print(summary(input, adsl = adsl)))
})

## Test 25: OCCDS formatted output is stable ----
test_that("print.summary_admiral_df Test 25: OCCDS formatted output is stable", {
  input <- as_admiral_df(tibble::tribble(
    ~USUBJID, ~ASEQ, ~AEBODSYS,          ~AEDECOD,    ~TRTEMFL,      ~AESEV,     ~AESER, ~AOCCFL,       ~ASTDT,                 ~TRTSDT, # nolint
    "1",         1L, "NERVOUS SYSTEM",   "HEADACHE",  "Y",           "MILD",     "N",    "Y",           as.Date("2024-01-10"), as.Date("2024-01-01"), # nolint
    "1",         2L, "NERVOUS SYSTEM",   "DIZZINESS", "Y",           "MODERATE", "N",    NA_character_, as.Date("2024-01-15"), as.Date("2024-01-01"), # nolint
    "2",         1L, "GASTROINTESTINAL", "NAUSEA",    "Y",           "SEVERE",   "Y",    "Y",           as.Date("2024-02-05"), as.Date("2024-02-01"), # nolint
    "2",         2L, "GASTROINTESTINAL", "NAUSEA",    NA_character_, "MILD",     "N",    NA_character_, as.Date("2024-01-25"), as.Date("2024-02-01") # nolint
  ))

  expect_snapshot(print(summary(input)))
})
