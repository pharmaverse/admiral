# Extracted from test-admiral_df.R:550

# test -------------------------------------------------------------------------
addl <- tibble::tribble(
    ~STUDYID, ~SPDEVID, ~DEVTYPE,
    "P",      "D-1",    "INSULIN PUMP",
    "P",      "D-2",    "INSULIN PUMP"
  )
expect_identical(guess_admiral_df_class(addl), "ADDL")
mdbds <- tibble::tribble(
    ~STUDYID, ~SPDEVID, ~PARAMCD, ~AVISITN, ~AVAL,
    "P",      "D-1",    "BATTV",         0,   3.1,
    "P",      "D-1",    "BATTV",         4,   2.9
  )
expect_identical(guess_admiral_df_class(mdbds), "MDBDS")
mdoccds <- tibble::tribble(
    ~STUDYID, ~SPDEVID, ~DEDECOD,    ~ASEQ,
    "P",      "D-1",    "OCCLUSION",     1,
    "P",      "D-1",    "OCCLUSION",     2
  )
expect_identical(guess_admiral_df_class(mdoccds), "MDOCCDS")
mdtte <- tibble::tribble(
    ~STUDYID, ~SPDEVID, ~PARAMCD, ~AVAL, ~CNSR, ~STARTDT,
    "P",      "D-1",    "TTFAIL",   120,     0, lubridate::ymd("2022-01-01")
  )
expect_identical(guess_admiral_df_class(mdtte), "MDBDS")
