# Input test dataset
adqs <- tibble::tribble(
  ~STUDYID, ~USUBJID, ~PARAMCD, ~AVISIT, ~AVISITN, ~AVAL,
  "ADMIRAL01", "01-701-1015", "ITEM1", "WEEK 10", 100, 5,
  "ADMIRAL01", "01-701-1015", "ITEM2", "WEEK 10", 100, 1,
  "ADMIRAL01", "01-701-1015", "ITEM3", "WEEK 10", 100, 2,
  "ADMIRAL01", "01-701-1015", "ITEM4", "WEEK 10", 100, NA,
  "ADMIRAL01", "01-701-1015", "ITEM1", "WEEK 20", 200, 1,
  "ADMIRAL01", "01-701-1015", "ITEM2", "WEEK 20", 200, 3,
  "ADMIRAL01", "01-701-1015", "ITEM3", "WEEK 20", 200, 5,
  "ADMIRAL01", "01-701-1015", "ITEM4", "WEEK 20", 200, 5,
  "ADMIRAL01", "01-701-1015", "ITEM5", "WEEK 20", 200, 1,
  "ADMIRAL01", "01-701-1015", "ITEM6", "WEEK 20", 200, 3,
  "ADMIRAL01", "01-701-1015", "ITEM7", "WEEK 20", 200, 3,
  "ADMIRAL01", "01-701-1281", "ITEM1", "WEEK 10", 100, 4,
)

## Test 1: compute_scale() works as expected ----
test_that("compute_scale Test 1: compute_scale works as expected", {
  input <- adqs %>%
    filter(USUBJID == "01-701-1015" & AVISITN == 100) %>%
    pull()

  expected_output <- mean((input - 1) * 25, na.rm = TRUE)

  expect_equal(compute_scale(input,
                             c(1, 5),
                             c(0, 100),
                             min_n = 3),
               expected_output)
})

## Test 2: scale is flipped if flip_direction == TRUE ----
test_that("compute_scale Test 2: scale is flipped if flip_direction == TRUE", {
  input <- adqs %>%
    filter(USUBJID == "01-701-1015" &
             AVISITN == 200 &
             PARAMCD %in% c("ITEM5", "ITEM6", "ITEM7")) %>%
    pull()

  expected_output <- 100 - (mean((input - 1) * 50, na.rm = TRUE))

  expect_equal(compute_scale(input,
                             c(1, 3),
                             c(0, 100),
                             flip_direction = TRUE,
                             min_n = 3),
               expected_output)
})

## Test 3: result is missing if min_n is not met ----
test_that("compute_scale Test 3: result is missing if min_n is not met", {
  input <- adqs %>%
    filter(USUBJID == "01-701-1281" & AVISITN == 100) %>%
    pull()

  expected_output <- NA

  expect_equal(compute_scale(input,
                             c(1, 5),
                             c(0, 100),
                             min_n = 2),
               expected_output)
})

## Test 4: compute_scale() works as expected within derive_summary_records() ----
test_that("compute_scale Test 4: compute_scale() works as expected within
          derive_summary_records()", {
  expected_output <- bind_rows(
    adqs,
    adqs %>%
      filter(PARAMCD %in% c("ITEM1", "ITEM2", "ITEM3")) %>%
      group_by(STUDYID, USUBJID, AVISIT, AVISITN) %>%
      summarise(n = n(), AVAL = mean(AVAL, na.rm = TRUE)) %>%
      mutate(
        PARAMCD = "ITEMAVG",
        AVAL = if_else(n >= 3, 100 - ((AVAL - 1) * 25), NA)
      ) %>%
      select(-n)
  )

  expect_equal(
    derive_summary_records(
      adqs,
      by_vars = exprs(STUDYID, USUBJID, AVISIT, AVISITN),
      filter = (PARAMCD %in% c("ITEM1", "ITEM2", "ITEM3")),
      analysis_var = AVAL,
      summary_fun = function(x)
        compute_scale(x, c(1, 5), c(0, 100), flip_direction = TRUE, min_n = 3),
      set_values_to = exprs(PARAMCD = "ITEMAVG")
    ),
    expected_output
  )
})
