library("tibble")
adsl <- tibble::tribble(
  ~USUBJID, ~SEX, ~COUNTRY,
  "ST42-1", "F",  "AUT",
  "ST42-2", "M",  "MWI",
  "ST42-3", "M",  "NOR",
  "ST42-4", "F",  "UGA"
) %>% mutate(STUDYID = "ST42")

adsl1 <- tibble::tribble(
  ~ID, ~SEX, ~COUNTRY,
  "ST42-1", "F", "AUT",
  "ST42-2", "M", "MWI",
  "ST42-3", "M", "NOR",
  "ST42-4", "F", "UGA"
) %>% mutate(STUDYID = "ST42")

adsl2 <- tibble::tribble(
  ~ID, ~SEX, ~COUNTRY,
  "ST42-1", "F", "AUT",
  "ST42-1", "F", "NOR",
  "ST42-2", "M", "MWI",
  "ST42-3", "M", "NOR",
  "ST42-4", "F", "UGA"
) %>% mutate(STUDYID = "ST42")


advs <- tibble::tribble(
  ~USUBJID, ~PARAMCD, ~AVISIT,    ~AVAL,
  "ST42-1", "WEIGHT", "BASELINE", 66,
  "ST42-1", "WEIGHT", "Week 2",   68,
  "ST42-2", "WEIGHT", "BASELINE", 88,
  "ST42-3", "WEIGHT", "Week 2",   55,
  "ST42-3", "WEIGHT", "Week 4",   50
) %>% mutate(STUDYID = "ST42")

advs1 <- tibble::tribble(
  ~ID, ~PARAMCD, ~AVISIT, ~AVAL,
  "ST42-1", "WEIGHT", "BASELINE", 66,
  "ST42-1", "WEIGHT", "Week 2", 68,
  "ST42-2", "WEIGHT", "BASELINE", 88,
  "ST42-3", "WEIGHT", "Week 2", 55,
  "ST42-3", "WEIGHT", "Week 4", 50
) %>% mutate(STUDYID = "ST42")


ex <- tibble::tribble(
  ~USUBJID, ~EXSTDTC,
  "ST42-1", "2020-12-07",
  "ST42-1", "2020-12-14",
  "ST42-2", "2021-01-12T12:00:00",
  "ST42-2", "2021-01-26T13:21",
  "ST42-3", "2021-03-02"
) %>% mutate(STUDYID = "ST42")

vs <- tibble::tribble(
  ~USUBJID, ~VSTESTCD, ~VSTEST, ~VSORRES, ~VSSEQ,
  "ST42-1", "DIABP", "Diastolic Blood Pressure", 64, 1,
  "ST42-1", "DIABP", "Diastolic Blood Pressure", 83, 2,
  "ST42-1", "WEIGHT", "Weight", 120, 3,
  "ST42-2", "WEIGHT", "Weight", 110, 1,
  "ST42-2", "HEIGHT", "Height", 58, 2
) %>% mutate(STUDYID = "ST42")
# derive_var_merged_exist_flag ----
## Test 14: merge existence flag ----
test_that("derive_var_exist_flag Test 14: merge existence flag", {
    
  actual <- derive_var_exist_flag(
    advs,
    new_var = VSEVALFL,
    condition = advs$AVISIT == "BASELINE"
  )
    
  expected <-
    data.frame(VSEVALFL = c("Y", "Y", NA_character_, NA_character_))

  print(actual)
  print(expected)
  expect_dfs_equal(
    base = expected,
    compare = actual,
    keys = "USUBJID"
  )
})

## Test 15: by_vars with rename ----
test_that("derive_var_merged_exist_flag Test 15: by_vars with rename", {
  expected <-
    c("Y", "Y", NA_character_, NA_character_)
    print(expected)
  actual <- derive_var_exist_flag(
    adsl,
    new_var = VSEVALFL,
    condition = AVISIT == "BASELINE",
  )

  
    
    print(expected)

  expect_dfs_equal(
    base = expected,
    compare = actual,
    keys = "USUBJID"
  )
})
