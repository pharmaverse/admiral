## Test 1: wihtout by_vars ----
test_that("derive_var_merged_ef_msrc Test 1: wihtout by_vars", {
  expected <- tibble::tribble(
    ~USUBJID, ~CANCTRFL,
    "1",      "Y",
    "2",      "Y",
    "3",      "Y",
    "4",      NA_character_
  )

  adsl <- select(expected, -CANCTRFL)

  cm <- tibble::tribble(
    ~USUBJID, ~CMCAT,        ~CMSEQ,
    "1",      "ANTI-CANCER",      1,
    "1",      "GENERAL",          2,
    "2",      "GENERAL",          1,
    "3",      "ANTI-CANCER",      1
  )

  pr <- tibble::tribble(
    ~USUBJID, ~PRSEQ,
    "2",      1,
    "3",      1
  )

  actual <- derive_var_merged_ef_msrc(
    adsl,
    flag_events = list(
      flag_event(
        dataset_name = "cm",
        condition = CMCAT == "ANTI-CANCER"
      ),
      flag_event(
        dataset_name = "pr"
      )
    ),
    source_datasets = list(cm = cm, pr = pr),
    by_vars = exprs(USUBJID),
    new_var = CANCTRFL
  )

  expect_dfs_equal(
    base = expected,
    compare = actual,
    keys = "USUBJID"
  )
})

## Test 2: with by_vars ----
test_that("derive_var_merged_ef_msrc Test 2: with by_vars", {
  expected <- tibble::tribble(
    ~USUBJID, ~EXLNKID, ~EXADJ,         ~DOSADJFL,
    "1",      "1",      "AE",           "Y",
    "1",      "2",      NA_character_,  "N",
    "1",      "3",      NA_character_,  "Y",
    "2",      "1",      NA_character_,  "N",
    "3",      "1",      NA_character_,  "Y"
  )

  adex <- select(expected, -DOSADJFL)

  ec <- tibble::tribble(
    ~USUBJID, ~ECLNKID, ~ECADJ,
    "1",      "3",      "AE",
    "3",      "1",      NA_character_
  )

  fa <- tibble::tribble(
    ~USUBJID, ~FALNKID, ~FATESTCD, ~FAOBJ,            ~FASTRESC,
    "3",      "1",      "OCCUR",   "DOSE ADJUSTMENT", "Y"
  )

  actual <- derive_var_merged_ef_msrc(
    adex,
    flag_events = list(
      flag_event(
        dataset_name = "ex",
        condition = !is.na(EXADJ)
      ),
      flag_event(
        dataset_name = "ec",
        condition = !is.na(ECADJ),
        by_vars = exprs(USUBJID, EXLNKID = ECLNKID)
      ),
      flag_event(
        dataset_name = "fa",
        condition = FATESTCD == "OCCUR" & FAOBJ == "DOSE ADJUSTMENT" & FASTRESC == "Y",
        by_vars = exprs(USUBJID, EXLNKID = FALNKID)
      )
    ),
    source_datasets = list(ex = adex, ec = ec, fa = fa),
    by_vars = exprs(USUBJID, EXLNKID),
    new_var = DOSADJFL,
    false_value = "N"
  )

  expect_dfs_equal(
    base = expected,
    compare = actual,
    keys = c("USUBJID", "EXLNKID")
  )
})
