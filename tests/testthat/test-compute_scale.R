## Test 1: compute_scale works as expected ----
test_that("compute_scale Test 1: compute_scale works as expected", {
  input <- c(5, 1, 2, NA)

  expected_output <- mean((input - 1) * 25, na.rm = TRUE)

  expect_equal(
    compute_scale(input,
      c(1, 5),
      c(0, 100),
      min_n = 3
    ),
    expected_output
  )
})

## Test 2: scale is flipped if flip_direction == TRUE ----
test_that("compute_scale Test 2: scale is flipped if flip_direction == TRUE", {
  input <- c(1, 3, 3)

  expected_output <- 100 - (mean((input - 1) * 50, na.rm = TRUE))

  expect_equal(
    compute_scale(input,
      c(1, 3),
      c(0, 100),
      flip_direction = TRUE,
      min_n = 3
    ),
    expected_output
  )
})

## Test 3: result is missing if min_n is not met ----
test_that("compute_scale Test 3: result is missing if min_n is not met", {
  input <- c(4)

  expected_output <- NA

  expect_equal(
    compute_scale(input,
      c(1, 5),
      c(0, 100),
      min_n = 2
    ),
    expected_output
  )
})

## Test 4: no transformation if source and target range aren't specified ----
test_that("compute_scale Test 4: no transformation if source and target range aren't specified", {
  input <- c(1, 3, 5, 5, 1, 3, 3)

  expected_output <- mean(input, na.rm = TRUE)

  expect_equal(
    compute_scale(input,
      min_n = 4
    ),
    expected_output
  )
})

## Test 5: works as expected within derive_summary_records() ----
test_that("compute_scale Test 5: works as expected within derive_summary_records()", {
  input <- tibble::tribble(
    ~STUDYID, ~USUBJID, ~PARAMCD, ~AVISIT, ~AVISITN, ~AVAL,
    "ADMIRAL01", "01-701-1015", "ITEM1", "WEEK 10", 100, 5,
    "ADMIRAL01", "01-701-1015", "ITEM2", "WEEK 10", 100, 1,
    "ADMIRAL01", "01-701-1015", "ITEM3", "WEEK 10", 100, 2,
    "ADMIRAL01", "01-701-1015", "ITEM4", "WEEK 10", 100, NA,
    "ADMIRAL01", "01-701-1015", "ITEM1", "WEEK 20", 200, 1,
    "ADMIRAL01", "01-701-1015", "ITEM2", "WEEK 20", 200, 3,
    "ADMIRAL01", "01-701-1015", "ITEM3", "WEEK 20", 200, 5,
    "ADMIRAL01", "01-701-1015", "ITEM4", "WEEK 20", 200, 5,
    "ADMIRAL01", "01-701-1281", "ITEM1", "WEEK 10", 100, 4,
  )

  expected_output <- bind_rows(
    input,
    input %>%
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
      input,
      dataset_add = input,
      by_vars = exprs(STUDYID, USUBJID, AVISIT, AVISITN),
      filter_add = (PARAMCD %in% c("ITEM1", "ITEM2", "ITEM3")),
      set_values_to = exprs(
        AVAL = compute_scale(AVAL, c(1, 5), c(0, 100), flip_direction = TRUE, min_n = 3),
        PARAMCD = "ITEMAVG"
      )
    ),
    expected_output
  )
})

## Test 6: error if source_range is supplied, but not target_range ----
test_that("compute_scale Test 6: error if source_range is supplied, but not target_range", {
  input <- c(1, 3, 5, 5, 1, 3, 3)

  expect_snapshot(
    error = TRUE,
    compute_scale(input,
      source_range = c(1, 5),
      min_n = 2
    )
  )
})

## Test 7: error if target_range is supplied, but not source_range ----
test_that("compute_scale Test 7: error if target_range is supplied, but not source_range", {
  input <- c(1, 3, 5, 5, 1, 3, 3)

  expect_snapshot(
    error = TRUE,
    compute_scale(input,
      target_range = c(0, 100),
      min_n = 2
    )
  )
})
