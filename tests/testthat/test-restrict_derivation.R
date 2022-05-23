# restrict_derivation ----
## restrict_derivation Test 1: restrict derivation with parameters ----
test_that("restrict_derivation Test 1: restrict derivation with parameters", {
  adlb <- tibble::tribble(
    ~USUBJID, ~AVISITN, ~AVAL, ~ABLFL,
    "1",            -1,   113, NA_character_,
    "1",             0,   113, "Y",
    "1",             3,   117, NA_character_,
    "2",             0,    95, "Y",
    "3",             0,   111, "Y",
    "3",             1,   101, NA_character_,
    "3",             2,   123, NA_character_
  )

  actual <- restrict_derivation(
    adlb,
    derivation = derive_var_base,
    args = params(by_vars = vars(USUBJID)),
    filter = AVISITN >= 0
  )

  expected <- mutate(adlb, BASE = c(NA, 113, 113, 95, 111, 111, 111))

  expect_dfs_equal(
    base = expected,
    compare = actual,
    keys = c("USUBJID", "AVISITN")
  )
})

## restrict_derivation Test 2: restrict derivation without parameters ----
test_that("restrict_derivation Test 2: restrict derivation without parameters", {
  adlb <- tibble::tribble(
    ~USUBJID, ~AVISITN, ~AVAL, ~ABLFL,        ~BASE,
    "1",            -1,   113, NA_character_,    NA,
    "1",             0,   113, "Y",             113,
    "1",             3,   117, NA_character_,   113,
    "2",             0,    95, "Y",              95,
    "3",             1,   101, NA_character_,    NA,
    "3",             2,   123, NA_character_,    NA
  )

  actual <- restrict_derivation(
    adlb,
    derivation = derive_var_chg,
    filter = AVISITN > 0
  )

  expected <- mutate(adlb, CHG = c(NA, NA, 4, NA, NA, NA))

  expect_dfs_equal(
    base = expected,
    compare = actual,
    keys = c("USUBJID", "AVISITN")
  )
})
