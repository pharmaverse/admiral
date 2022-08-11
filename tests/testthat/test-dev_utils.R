library(admiraltest)
data(admiral_dm)
dm <- select(admiral_dm, USUBJID)

test_that("creating temporary variables works", {
  expect_identical(get_new_tmp_var(dm), sym("tmp_var_1"))
})

test_that("the temporary variable counter is increased correctly", {
  dm$tmp_var_1 <- NA
  dm$tmp_var_2 <- NA
  expect_identical(get_new_tmp_var(dm), sym("tmp_var_3"))
})

test_that("no variables are removed when no tmp vars are present", {
  expect_identical(dm, remove_tmp_vars(dm))
})

test_that("removing temporary variables works", {
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

test_that("removing temporary variables works when using the pipe operator", {
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
