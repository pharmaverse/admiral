dm <- dplyr::tribble(
  ~USUBJID,
  "01-701-1015",
  "01-701-1016",
)

# get_new_tmp_var ----
## Test 1: creating temporary variables works ----
test_that("get_new_tmp_var Test 1: creating temporary variables works", {
  expect_identical(get_new_tmp_var(dm), sym("tmp_var_1"))
})

## Test 2: errors if name does not start with tmp ----
test_that("get_new_tmp_var Test 2: errors if name does not start with tmp", {
  expect_snapshot(
    error = TRUE,
    get_new_tmp_var(dm, prefix = "my_var")
  )
})

## Test 3: the temporary variable counter is increased correctly ----
test_that("get_new_tmp_var Test 3: the temporary variable counter is increased correctly", {
  dm$tmp_var_1 <- NA
  dm$tmp_var_2 <- NA
  expect_identical(get_new_tmp_var(dm), sym("tmp_var_3"))
})

# remove_tmp_vars ----
## Test 4: no variables are removed when no tmp vars are present ----
test_that("remove_tmp_vars Test 4: no variables are removed when no tmp vars are present", {
  expect_identical(dm, remove_tmp_vars(dm))
})

## Test 5: removing temporary variables works ----
test_that("remove_tmp_vars Test 5: removing temporary variables works", {
  tmp_var <- get_new_tmp_var(dm)
  dm <- mutate(dm, !!tmp_var := NA)
  do_something <- function(dataset) {
    tmp_var_1 <- get_new_tmp_var(dm)
    tmp_var_2 <- get_new_tmp_var(dm)
    dm <- mutate(dm, !!tmp_var_1 := NA, !!tmp_var_2 := NA)
    remove_tmp_vars(dm)
  }
  expect_identical(colnames(dm), colnames(do_something(dm)))
})

## Test 6: removing temp variables works with the pipe operator ----
test_that("remove_tmp_vars Test 6: removing temp variables works with the pipe operator", {
  tmp_var <- get_new_tmp_var(dm)
  dm <- mutate(dm, !!tmp_var := NA)
  do_something_with_pipe <- function(dataset) {
    tmp_var_1 <- get_new_tmp_var(dm)
    tmp_var_2 <- get_new_tmp_var(dm)
    dm %>%
      mutate(!!tmp_var_1 := NA, !!tmp_var_2 := NA) %>%
      remove_tmp_vars()
  }
  expect_identical(colnames(dm), colnames(do_something_with_pipe(dm)))
})


## Test 7: running get_new_tmp_var on NULL dataset creates generic variable ----
test_that("running get_new_tmp_var on NULL dataset creates generic variable", {
  df <- NULL
  expect_identical(get_new_tmp_var(df), sym("tmp_var_1"))
})
