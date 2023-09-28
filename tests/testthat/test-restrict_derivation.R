# restrict_derivation ----
## Test 1: restrict derivation with parameters ----
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
    args = params(by_vars = exprs(USUBJID)),
    filter = AVISITN >= 0
  )

  expected <- mutate(adlb, BASE = c(NA, 113, 113, 95, 111, 111, 111))

  expect_dfs_equal(
    base = expected,
    compare = actual,
    keys = c("USUBJID", "AVISITN")
  )
})

## Test 2: restrict derivation without parameters ----
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

## Test 3: access functions from the parent environment ----
test_that("restrict_derivation Test 3: access functions from the parent environment", {
  my_derivation <- function(dataset, new_var) {
    mutate(
      dataset,
      !!enexpr(new_var) := 42
    )
  }

  my_data <- tibble::tribble(
    ~PARAMCD,
    "A",
    "B"
  )

  expect_silent({
    restrict_derivation(
      my_data,
      derivation = my_derivation,
      args = params(
        new_var = X
      ),
      filter = PARAMCD == "A"
    )
  })
})

## Test 4: allow dplyr functions ----
test_that("restrict_derivation Test 4: allow dplyr functions", {
  adlb <- tibble::tribble(
    ~USUBJID, ~AVISITN, ~AVAL,
    "1",            -1,   113,
    "1",             0,   113,
    "1",             3,   117,
    "2",             0,    95,
    "3",             0,   111,
    "3",             1,   101,
    "3",             2,   123
  )

  actual <- restrict_derivation(
    adlb,
    derivation = mutate,
    args = params(AVAL = AVAL + 1),
    filter = USUBJID == "1"
  )

  expected <- adlb %>%
    mutate(AVAL = ifelse(USUBJID == "1", AVAL + 1, AVAL))

  expect_dfs_equal(
    base = expected,
    compare = actual,
    keys = c("USUBJID", "AVISITN")
  )
})
