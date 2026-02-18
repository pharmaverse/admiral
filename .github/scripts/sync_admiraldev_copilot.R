#!/usr/bin/env Rscript

#' Sync admiraldev Documentation for AI Assistants (v2.1 - Simplified)
#'
#' Creates AI assistant context files for admiral development:
#' - Root AGENT.md (universal AI assistant context)
#' - tests/testthat/AGENT.md (testing-specific context)
#' - .github/copilot-instructions.md (GitHub Copilot specific, legacy)

# Load required packages
suppressPackageStartupMessages({
  library(glue)
})

# Configuration
project_root <- if (file.exists(".git")) "." else ".."
agent_md_root <- file.path(project_root, "AGENT.md")
agent_md_tests <- file.path(project_root, "tests", "testthat", "AGENT.md")
copilot_instructions <- file.path(project_root, ".github", "copilot-instructions.md")
github_raw_base <- "https://raw.githubusercontent.com/pharmaverse/admiraldev/main/vignettes"

# Create directories if needed
dir.create(dirname(copilot_instructions), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(agent_md_tests), recursive = TRUE, showWarnings = FALSE)

cat("Creating AI assistant context files for admiral development...\n")

# Download admiraldev content (with fallbacks)
download_admiraldev_content <- function() {
  programming_strategy <- NULL
  unit_test_guidance <- NULL

  # Try to download programming strategy
  tryCatch(
    {
      programming_strategy <<- readLines(
        file.path(github_raw_base, "programming_strategy.Rmd"),
        warn = FALSE
      )
      cat("✓ Downloaded programming_strategy.Rmd\n")
    },
    error = function(e) {
      cat("⚠ Could not download programming_strategy.Rmd\n")
    }
  )

  # Try to download unit test guidance
  tryCatch(
    {
      unit_test_guidance <<- readLines(
        file.path(github_raw_base, "unit_test_guidance.Rmd"),
        warn = FALSE
      )
      cat("✓ Downloaded unit_test_guidance.Rmd\n")
    },
    error = function(e) {
      cat("⚠ Could not download unit_test_guidance.Rmd\n")
    }
  )

  list(
    programming_strategy = programming_strategy,
    unit_test_guidance = unit_test_guidance
  )
}

# Create root AGENT.md with actual admiraldev content
create_root_agent_md <- function(admiraldev_content) {
  # Create the file content as separate parts to avoid glue complexity
  header_lines <- c(
    "# Admiral Development Guidelines for AI Assistants",
    "",
    "This file provides context for AI coding assistants (GitHub Copilot, Gemini, Claude, etc.) about admiral ecosystem standards and best practices.",
    "",
    paste("**Auto-generated:**", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
    "**Source:** admiraldev package vignettes",
    "**Update script:** `source('.github/scripts/sync_admiraldev_copilot.R')`",
    "",
    "## Purpose",
    "",
    "Help AI assistants provide admiral-compliant code suggestions:",
    "- ✅ Consistent function naming and patterns",
    "- ✅ Proper argument handling and validation",
    "- ✅ CDISC/ADaM context and conventions",
    "- ✅ Testing and documentation standards",
    "",
    "---",
    ""
  )

  # Programming strategy content
  programming_lines <- c()
  if (!is.null(admiraldev_content$programming_strategy)) {
    programming_lines <- c(
      "# Admiral Programming Strategy",
      "",
      "**Source:** admiraldev programming_strategy.Rmd vignette",
      "",
      admiraldev_content$programming_strategy,
      "",
      "---",
      ""
    )
  } else {
    # Fallback content
    programming_lines <- c(
      "# Admiral Programming Strategy",
      "",
      "**Note:** Could not download latest admiraldev content. Using essential guidelines.",
      "**Full documentation:** https://pharmaverse.github.io/admiraldev/articles/programming_strategy.html",
      "",
      "## Function Design Principles",
      "",
      "- **Modularity:** All code follows a modular approach with clearly separated steps",
      "- **Avoid Copy and Paste:** Similar code should be put into separate functions",
      "- **Checks:** Meaningful error messages with clear reference to failing input",
      "- **Flexibility:** Functions should be as flexible as possible without reducing usability",
      "",
      "## Function Naming Convention",
      "",
      "Function names should start with a verb and use snake_case:",
      "",
      "- `derive_*` - Add rows/columns to datasets",
      "- `derive_var_*` - Add single variable",
      "- `derive_vars_*` - Add multiple variables",
      "- `compute_*` - Vector operations, return vectors",
      "- `assert_*` / `warn_*` - Input validation",
      "- `filter_*` - Filter observations",
      "",
      "---",
      ""
    )
  }

  # Footer content
  footer_lines <- c(
    "## Testing Guidelines",
    "",
    "For admiral unit testing guidelines, see `tests/testthat/AGENT.md`",
    "",
    "## Full Documentation",
    "",
    "- **Programming Strategy:** https://pharmaverse.github.io/admiraldev/articles/programming_strategy.html",
    "- **Unit Test Guidance:** https://pharmaverse.github.io/admiraldev/articles/unit_test_guidance.html",
    "- **Admiral Website:** https://pharmaverse.github.io/admiral/",
    "",
    "---",
    "",
    "*This file helps AI assistants understand admiral patterns for better code suggestions.*"
  )

  # Combine all content
  all_lines <- c(header_lines, programming_lines, footer_lines)
  writeLines(all_lines, agent_md_root)
  cat(glue("✓ Created: {agent_md_root}\n"))
}

# Create tests/testthat/AGENT.md with actual admiraldev content
create_tests_agent_md <- function(admiraldev_content) {
  # Create header
  header_lines <- c(
    "# Admiral Unit Testing Guidelines for AI Assistants",
    "",
    "Context for AI assistants when working with admiral unit tests in `tests/testthat/`.",
    "",
    paste("**Auto-generated:**", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
    "**Source:** admiraldev unit testing guidance vignette",
    "",
    "---",
    ""
  )

  # Testing content
  testing_lines <- c()
  if (!is.null(admiraldev_content$unit_test_guidance)) {
    testing_lines <- c(
      "# Admiral Unit Testing Guidance",
      "",
      "**Source:** admiraldev unit_test_guidance.Rmd vignette",
      "",
      admiraldev_content$unit_test_guidance,
      "",
      "---",
      ""
    )
  } else {
    # Fallback content
    testing_lines <- c(
      "# Admiral Unit Testing Guidelines",
      "",
      "**Note:** Could not download latest admiraldev content. Using essential guidelines.",
      "**Full documentation:** https://pharmaverse.github.io/admiraldev/articles/unit_test_guidance.html",
      "",
      "## Test Structure and Organization",
      "",
      "### File Naming and Organization",
      "- Match function names: `test-derive_var_example.R`",
      "- One test file per main function (can include helpers)",
      "- Use descriptive test names explaining the scenario",
      "",
      "### Admiral Test Requirements",
      "",
      "1. **Happy path** - basic functionality works",
      "2. **Error conditions** - invalid inputs fail appropriately",
      "3. **Edge cases** - empty data, boundary conditions",
      "4. **Argument validation** - all parameters properly validated",
      "5. **Output structure** - correct columns, types, ungrouped",
      "",
      "### Test Data Best Practices",
      "",
      "- Use `tribble()` for test data creation",
      "- Test that output is ungrouped with `expect_false(is.grouped_df(result))`",
      "- Test meaningful error messages",
      "- Aim for 80%+ test coverage",
      "",
      "---",
      ""
    )
  }

  # Footer
  footer_lines <- c(
    "*For complete guidance: https://pharmaverse.github.io/admiraldev/articles/unit_test_guidance.html*"
  )

  # Combine all content
  all_lines <- c(header_lines, testing_lines, footer_lines)
  writeLines(all_lines, agent_md_tests)
  cat(glue("✓ Created: {agent_md_tests}\n"))
}

# Create optimized .github/copilot-instructions.md for GitHub Copilot
create_copilot_instructions <- function(admiraldev_content) {
  # Simple static content optimized for Copilot
  content_lines <- c(
    "# GitHub Copilot Instructions - admiral Development",
    "",
    paste("**Auto-generated:**", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
    "**Optimized for:** GitHub Copilot code completion",
    "**Complete guidelines:** See `AGENT.md` files for full context",
    "",
    "⚠️ **DO NOT EDIT MANUALLY** - Run `source('.github/scripts/sync_admiraldev_copilot.R')` to update",
    "",
    "## Admiral Code Completion Patterns",
    "",
    "### Function Names (verb_object_detail)",
    "- `derive_var_base()` - Single variable derivation",
    "- `derive_vars_merged()` - Multiple variables from merge",
    "- `derive_param_tte()` - Parameter derivation",
    "- `compute_age_years()` - Vector computation",
    "- `assert_data_frame()` - Input validation",
    "",
    "### Argument Patterns",
    "Standard admiral function signature:",
    "- `dataset` - Always first argument",
    "- `by_vars = exprs(...)` - Grouping variables",
    "- `order = exprs(...)` - Sorting expressions",
    "- `new_var = VAR` - New variable (symbol)",
    "- `filter = PARAM == \"VALUE\"` - Filtering expression",
    "",
    "### Common Expressions",
    "- `subject_keys = exprs(STUDYID, USUBJID)`",
    "- `by_vars = exprs(USUBJID, PARAMCD)`",
    "- `order = exprs(AVISITN, desc(ADY))`",
    "- `filter = PARAMCD == \"TEMP\"`",
    "",
    "### Test Patterns",
    "- Use `tribble()` for test data",
    "- `expect_true(\"DERIVED\" %in% names(result))`",
    "- `expect_false(is.grouped_df(result))`",
    "- `expect_error(..., \"Required variable.*missing\")`",
    "",
    "---",
    "",
    "*This file is optimized for GitHub Copilot code completion. For complete admiral development guidelines, see `AGENT.md` and `tests/testthat/AGENT.md`.*"
  )

  writeLines(content_lines, copilot_instructions)
  cat(glue("✓ Created: {copilot_instructions}\n"))
}

# Main execution
admiraldev_content <- download_admiraldev_content()

cat("\nCreating AI assistant context files:\n")
create_root_agent_md(admiraldev_content)
create_tests_agent_md(admiraldev_content)
create_copilot_instructions(admiraldev_content)

cat(glue("\n🎉 Success! Created AI assistant context files:
📄 {agent_md_root} (universal AI assistants)
📄 {agent_md_tests} (testing context)
📄 {copilot_instructions} (Copilot legacy)

These files will help AI assistants provide better admiral-compliant suggestions.
"))
