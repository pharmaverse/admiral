adqs <- tibble::tribble(
  ~STUDYID, ~USUBJID, ~PARAMCD, ~AVISIT, ~AVISITN, ~AVAL,
  "ADMIRAL01", "01-701-1015", "ITEM1", "WEEK 10", 100, 5,
  "ADMIRAL01", "01-701-1015", "ITEM2", "WEEK 10", 100, 1,
  "ADMIRAL01", "01-701-1015", "ITEM3", "WEEK 10", 100, 2,
  "ADMIRAL01", "01-701-1015", "ITEM4", "WEEK 10", 100, 4,
  "ADMIRAL01", "01-701-1015", "ITEM5", "WEEK 10", 100, 4,
  "ADMIRAL01", "01-701-1015", "ITEM1", "WEEK 20", 200, 1,
  "ADMIRAL01", "01-701-1015", "ITEM2", "WEEK 20", 200, 3,
  "ADMIRAL01", "01-701-1015", "ITEM3", "WEEK 20", 200, 5,
  "ADMIRAL01", "01-701-1015", "ITEM4", "WEEK 20", 200, 5,
  "ADMIRAL01", "01-701-1015", "ITEM5", "WEEK 20", 200, 2,
  "ADMIRAL01", "01-701-1281", "ITEM1", "WEEK 10", 100, 4,
  "ADMIRAL01", "01-701-1281", "ITEM2", "WEEK 10", 100, 1,
  "ADMIRAL01", "01-701-1281", "ITEM3", "WEEK 10", 100, 2,
  "ADMIRAL01", "01-701-1281", "ITEM4", "WEEK 10", 100, 3,
  "ADMIRAL01", "01-701-1281", "ITEM5", "WEEK 10", 100, 1,
  "ADMIRAL01", "01-701-1281", "ITEM1", "WEEK 20", 200, 2,
  "ADMIRAL01", "01-701-1281", "ITEM2", "WEEK 20", 200, 4,
  "ADMIRAL01", "01-701-1281", "ITEM3", "WEEK 20", 200, 5,
  "ADMIRAL01", "01-701-1281", "ITEM4", "WEEK 20", 200, 2,
  "ADMIRAL01", "01-701-1281", "ITEM5", "WEEK 20", 200, 3
)

## Test 1: compute_scale() works as expected ----
test_that("compute_scale Test 1: compute_scale works as expected", {

  input <- adqs %>%
    filter(USUBJID == "01-701-1015" & AVISITN == 100) %>%
    pull()

  expected_output <- mean((input - 1) * 25, na.rm = TRUE)


  expect_equal(compute_scale(input, c(1, 5), c(0, 100), min_n = 4), expected_output)

})

## Test 2: compute_scale() works as expected within derive_summary_records() ----
test_that("compute_scale Test 2: compute_scale works as expected", {

  expected_output <- bind_rows(adqs,
                              adqs %>%
                                filter(PARAMCD %in% c("ITEM1", "ITEM2", "ITEM3", "ITEM4")) %>%
                                  group_by(STUDYID, USUBJID, AVISIT, AVISITN) %>%
                                  summarise(AVAL = mean(AVAL, na.rm = TRUE)) %>%
                                  mutate(
                                    PARAMCD = "ITEMAVG",
                                    AVAL = (AVAL - 1) * 25
                                  )
                              )

  expect_equal(derive_summary_records(
                adqs,
                by_vars = vars(STUDYID, USUBJID, AVISIT, AVISITN),
                filter = (PARAMCD %in% c("ITEM1", "ITEM2", "ITEM3", "ITEM4")),
                analysis_var = AVAL,
                summary_fun = function(x) compute_scale(x, c(1, 5), c(0, 100), min_n = 3),
                set_values_to = vars(PARAMCD = "ITEMAVG")
              ),
  expected_output)

})

# add test for when AVAL values are outside of source_range (look expect_failure() from testthat - does this test whether an error is thrown?)
