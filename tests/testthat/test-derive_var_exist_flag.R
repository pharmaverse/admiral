library("tibble")
library("dplyr")
library("testthat")
library("rlang")
library("admiral")
library("admiraldev")
adsl <- tibble::tribble(
  ~USUBJID, ~SEX, ~COUNTRY,
  "ST42-1", "F",  "AUT",
  "ST42-2", "M",  "MWI",
  "ST42-3", "M",  "NOR",
  "ST42-4", "F",  "UGA"
) %>% mutate(STUDYID = "ST42")

advs <- tibble::tribble(
  ~USUBJID, ~PARAMCD, ~AVISIT,    ~AVAL,
  "ST42-1", "WEIGHT", "BASELINE", 66,
  "ST42-1", "WEIGHT", "Week 2",   68,
  "ST42-2", "WEIGHT", "BASELINE", 88,
  "ST42-3", "WEIGHT", "Week 2",   55,
  "ST42-3", "WEIGHT", "Week 4",   50
) %>% mutate(STUDYID = "ST42")

# derive_var_exist_flag ----
## Test 1: generate existence flag ----
test_that("derive_var_exist_flag Test 1: generate existence flag", {
  
  actual <- derive_var_exist_flag(
    advs,
    new_var = VSEVALFL,
    condition = advs$AVISIT == "BASELINE"
  )
    
  expected <- mutate(advs,VSEVALFL = c(1, 0, 1, 0, 0))

  expect_dfs_equal(
    base = expected,
    compare = actual,
    keys = c("USUBJID","AVISIT")
  )
})
