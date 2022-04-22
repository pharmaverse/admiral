# slice_derivation ----
## slice_derivation Test 1: slice derivation ----
test_that("slice_derivation Test 1: slice derivation", {
  advs <- tibble::tribble(
    ~USUBJID, ~VSDTC,       ~VSTPT,             ~VSSEQ,
    "1",      "2020-04-16", NA_character_,      1,
    "1",      "2020-04-16", "BEFORE TREATMENT", 2
  )

  actual <- slice_derivation(
    advs,
    derivation = derive_vars_dtm,
    args = params(dtc = VSDTC,
                  new_vars_prefix = "A"),
    derivation_slice(
      filter = str_detect(VSTPT, "PRE|BEFORE"),
      args = params(time_imputation = "first")
    ),
    derivation_slice(args = params(time_imputation = "last"))
  )

  expected <- mutate(advs,
                     ADTM = c(ymd_hms("2020-04-16 23:59:59"), ymd_hms("2020-04-16 00:00:00")),
                     ATMF = "H")

  expect_dfs_equal(base = expected,
                   compare = actual,
                   keys = c("USUBJID", "VSSEQ"))
})

## slice_derivation Test 2: non matching observations ----
test_that("slice_derivation Test 2: non matching observations", {
  advs <- tibble::tribble(
    ~USUBJID, ~VSDTC,       ~VSTPT,             ~VSSEQ,
    "1",      "2020-04-16", NA_character_,      1,
    "1",      "2020-04-16", "BEFORE TREATMENT", 2
  )

  actual <- slice_derivation(
    advs,
    derivation = derive_vars_dtm,
    args = params(dtc = VSDTC,
                  new_vars_prefix = "A"),
    derivation_slice(
      filter = str_detect(VSTPT, "PRE|BEFORE"),
      args = params(time_imputation = "first")
    )
  )

  expected <- mutate(advs,
                     ADTM = c(ymd_hms(NA), ymd_hms("2020-04-16 00:00:00")),
                     ATMF = c(NA_character_, "H")) %>%
    mutate(ADTM = as_iso_dtm(ADTM))

  expect_dfs_equal(base = expected,
                   compare = actual,
                   keys = c("USUBJID", "VSSEQ"))
})

# print.derivation_slice ----
## print.derivation_slice Test 1: `derivation_slice` objects are printed as intended ----
test_that("print.derivation_slice Test1: `derivation_slice` objects are printed as intended", {
  slice <-
    derivation_slice(filter = AVISITN > 0,
                     args = params(new_var = CHG))
  expected_print_output <- c(
    "<derivation_slice> object",
    "filter: AVISITN > 0 ",
    "args:",
    "$new_var",
    "CHG",
    "",
    "attr(,\"class\")",
    "[1] \"params\" \"list\"  "
  )
  expect_identical(capture.output(print(slice)), expected_print_output)
})
