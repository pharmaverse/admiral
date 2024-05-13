advs <- tibble::tribble(
  ~USUBJID, ~PARAMCD, ~AVISIT,    ~AVAL,
  "ST42-1", "WEIGHT", "BASELINE", 66,
  "ST42-1", "WEIGHT", "Week 2",   68,
  "ST42-2", "WEIGHT", "BASELINE", 88,
  "ST42-3", "WEIGHT", "Week 2",   55,
  "ST42-3", "WEIGHT", "Week 4",   50
) %>% mutate(STUDYID = "ST42")

# get_flagged_records ----
## Test 1: generate existence flag ----
test_that("get_flagged_records Test 1: generate existence flag", {
  actual <- get_flagged_records(
    dataset = advs,
    new_var = VSEVALFL,
    condition = AVISIT == "BASELINE"
  )

  expected <- mutate(advs, VSEVALFL = c(1, 0, 1, 0, 0))

  expect_dfs_equal(
    base = expected,
    compare = actual,
    keys = c("USUBJID", "AVISIT")
  )
})
