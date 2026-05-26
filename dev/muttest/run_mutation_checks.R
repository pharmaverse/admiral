#!/usr/bin/env Rscript

if (!requireNamespace("muttest", quietly = TRUE)) {
  stop(
    paste(
      "Package `muttest` is required.",
      "Install it first with install.packages('muttest')."
    ),
    call. = FALSE
  )
}

source_files <- c(
  "R/compute_age_years.R",
  "R/compute_scale.R",
  "R/filter_exist.R"
)

test_files <- c(
  "tests/testthat/test-compute_age_years.R",
  "tests/testthat/test-compute_scale.R",
  "tests/testthat/test-filter_exist.R"
)

missing_files <- c(
  source_files[!file.exists(source_files)],
  test_files[!file.exists(test_files)]
)

if (length(missing_files) > 0L) {
  stop(
    paste(
      "Missing expected files:",
      paste(missing_files, collapse = ", ")
    ),
    call. = FALSE
  )
}

plan <- muttest::muttest_plan(
  source_files = source_files,
  mutators = c(
    muttest::comparison_operators(),
    muttest::logical_operators(),
    muttest::condition_mutations()
  )
)

message(
  "Running muttest against:\n",
  paste0(" - ", source_files, collapse = "\n"),
  "\nwith matching tests in tests/testthat."
)

result <- muttest::muttest(
  plan = plan,
  path = "tests/testthat",
  test_strategy = muttest::FileTestStrategy$new(load_package = "source"),
  workers = 1
)

print(result)
