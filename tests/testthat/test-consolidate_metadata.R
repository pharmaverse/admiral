## Test 1: consolidation works ----
test_that("consolidate_metadata Test 1: consolidation works", {
  glob <- tibble::tribble(
    ~id, ~val,
    1,   "glob_val_1",
    2,   "glob_val_2"
  )
  proj <- tibble::tribble(
    ~id, ~val,
    2,   "proj_val2"
  )
  stud <- tibble::tribble(
    ~id, ~val,
    3,   "stud_val_3"
  )

  expected <- tibble::tribble(
    ~id, ~val,         ~SOURCE,
    1,   "glob_val_1", "global",
    2,   "proj_val2",  "project",
    3,   "stud_val_3", "study"
  )

  expect_dfs_equal(
    base = expected,
    comp = consolidate_metadata(
      datasets = list(
        global = glob,
        project = proj,
        study = stud
      ),
      key_vars = exprs(id)
    ),
    keys = c("id")
  )
})

## Test 2: error if key vars are not unique ----
test_that("consolidate_metadata Test 2: error if key vars are not unique", {
  glob <- tibble::tribble(
    ~id, ~val,
    1,   "glob_val_1a",
    1,   "glob_val_1b",
    2,   "glob_val_2"
  )
  stud <- tibble::tribble(
    ~id, ~val,
    3,   "stud_val_3"
  )

  expect_error(
    consolidate_metadata(
      datasets = list(
        global = glob,
        study = stud
      ),
      key_vars = exprs(id)
    ),
    "Dataset contains duplicate records with respect to"
  )
})

## Test 3: warn if variables differ ----
test_that("consolidate_metadata Test 4: warn if variables differ", {
  glob <- tibble::tribble(
    ~id, ~val,
    1,   "glob_val_1",
    2,   "glob_val_2"
  )
  stud <- tibble::tribble(
    ~var, ~id, ~val,
    "abc",  3, "stud_val_3"
  )

  expect_snapshot(
    consolidate_metadata(
      datasets = list(
        global = glob,
        study = stud
      ),
      key_vars = exprs(id)
    )
  )
})
