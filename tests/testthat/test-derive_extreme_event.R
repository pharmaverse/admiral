# derive_extreme_records ----
## Test 1: add the worst event for each group ----
test_that("derive_extreme_records Test 1: add the worst event for each group", {
  input <- tibble::tribble(
    ~USUBJID, ~PARAMCD,       ~AVALC,
    "1",      "NO SLEEP",     "N",
    "1",      "WAKE UP",      "N",
    "1",      "FALL ASLEEP",  "N",
    "2",      "NO SLEEP",     "N",
    "2",      "WAKE UP",      "Y",
    "2",      "FALL ASLEEP",  "N",
    "3",      "NO SLEEP",     NA_character_,
    "3",      "WAKE UP",      "N"
  )

  expected_output <- bind_rows(
    input,
    tibble::tribble(
      ~USUBJID, ~PARAMCD, ~AVALC,                            ~AVAL,
      "1",      "WSP",    "No sleeping problems",                4,
      "2",      "WSP",    "Waking up more than three times",     2,
      "3",      "WSP",    "Missing",                            99
    )
  )

  actual_output <- derive_extreme_event(
    input,
    by_vars = exprs(USUBJID),
    events = list(
      event(
        condition = PARAMCD == "NO SLEEP" & AVALC == "Y",
        set_values_to = exprs(AVALC = "No sleep", AVAL = 1)
      ),
      event(
        condition = PARAMCD == "WAKE UP" & AVALC == "Y",
        set_values_to = exprs(AVALC = "Waking up more than three times", AVAL = 2)
      ),
      event(
        condition = PARAMCD == "FALL ASLEEP" & AVALC == "Y",
        set_values_to = exprs(AVALC = "More than 30 mins to fall asleep", AVAL = 3)
      ),
      event(
        condition = all(AVALC == "N"),
        set_values_to = exprs(
          AVALC = "No sleeping problems", AVAL = 4
        )
      ),
      event(
        condition = TRUE,
        set_values_to = exprs(AVALC = "Missing", AVAL = 99)
      )
    ),
    order = exprs(AVAL),
    mode = "first",
    set_values_to = exprs(
      PARAMCD = "WSP"
    ),
    check_type = "none"
  )

  expect_dfs_equal(
    base = expected_output,
    compare = actual_output,
    keys = c("USUBJID", "PARAMCD")
  )
})
