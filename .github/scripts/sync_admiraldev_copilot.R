#!/usr/bin/env Rscript

#' Sync admiraldev Documentation for AI Assistants (v2.1 - Simplified)
#'
#' Creates AI assistant context files for admiral development:
#' - Root AGENTS.md (universal AI assistant context)
#' - tests/testthat/AGENTS.md (testing-specific context)
#' - .github/copilot-instructions.md (GitHub Copilot specific, legacy)

# Load required packages
suppressPackageStartupMessages({
  library(glue)
})

# Configuration
project_root <- if (file.exists(".git")) "." else ".."
agent_md_root <- file.path(project_root, "AGENTS.md")
agent_md_tests <- file.path(project_root, "tests", "testthat", "AGENTS.md")
copilot_instructions <- file.path(project_root, ".github", "copilot-instructions.md")
github_raw_base <- "https://raw.githubusercontent.com/pharmaverse/admiraldev/main/vignettes"

# Create directories if needed
dir.create(dirname(copilot_instructions), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(agent_md_tests), recursive = TRUE, showWarnings = FALSE)

cat("Creating AI assistant context files for admiral development...\n")
cat(glue("GitHub raw base URL: {github_raw_base}\n"))

# Test if we can reach the URLs
cat("Testing network connectivity...\n")
test_url <- file.path(github_raw_base, "programming_strategy.Rmd")
cat(glue("Testing URL: {test_url}\n"))

# Download admiraldev content (with multiple methods and better error handling)
download_admiraldev_content <- function() {
  programming_strategy <- NULL
  unit_test_guidance <- NULL

  # Test basic connectivity first
  cat("Testing network connectivity to GitHub...\n")
  can_connect <- FALSE
  tryCatch(
    {
      # Simple connectivity test
      con <- url("https://github.com", open = "rt")
      close(con)
      can_connect <- TRUE
      cat("✓ GitHub connectivity confirmed\n")
    },
    error = function(e) {
      cat("⚠ GitHub connectivity issue\n")
    }
  )

  if (!can_connect) {
    cat("Network connectivity issues - using fallback content\n")
    return(list(programming_strategy = NULL, unit_test_guidance = NULL))
  }

  # Function to try multiple download methods with timeouts
  try_download <- function(url, description) {
    content <- NULL

    # Method 1: readLines with timeout
    tryCatch(
      {
        con <- url(url, open = "rt")
        on.exit(try(close(con), silent = TRUE))
        content <- readLines(con, warn = FALSE, timeout = 30)
        close(con)
        on.exit() # Cancel the on.exit

        if (length(content) > 10) { # Basic sanity check
          cat(glue("✓ Downloaded {description} (readLines, {length(content)} lines)\n"))
          return(content)
        }
      },
      error = function(e) {
        cat(glue("⚠ Method 1 failed for {description}: {substr(e$message, 1, 50)}...\n"))
      }
    )

    # Method 2: download.file with timeout
    tryCatch(
      {
        temp_file <- tempfile(fileext = ".Rmd")
        result <- download.file(url, temp_file, quiet = TRUE, method = "auto", timeout = 30)

        if (result == 0 && file.exists(temp_file) && file.size(temp_file) > 100) {
          content <- readLines(temp_file, warn = FALSE)
          unlink(temp_file)
          if (length(content) > 10) {
            cat(glue("✓ Downloaded {description} (download.file, {length(content)} lines)\n"))
            return(content)
          }
        }
        unlink(temp_file)
      },
      error = function(e) {
        cat(glue("⚠ Method 2 failed for {description}: {substr(e$message, 1, 50)}...\n"))
      }
    )

    # Method 3: curl with timeout if available
    if (Sys.which("curl") != "") {
      tryCatch(
        {
          temp_file <- tempfile(fileext = ".Rmd")
          result <- system2("curl", c("-s", "-L", "--connect-timeout", "10", "--max-time", "30", "-o", temp_file, url),
            stdout = FALSE, stderr = FALSE
          )

          if (result == 0 && file.exists(temp_file) && file.size(temp_file) > 100) {
            content <- readLines(temp_file, warn = FALSE)
            unlink(temp_file)
            if (length(content) > 10) {
              cat(glue("✓ Downloaded {description} (curl, {length(content)} lines)\n"))
              return(content)
            }
          }
          unlink(temp_file)
        },
        error = function(e) {
          cat(glue("⚠ Method 3 (curl) failed for {description}: {substr(e$message, 1, 50)}...\n"))
        }
      )
    }

    cat(glue("✗ All download methods failed for {description}\n"))
    return(NULL)
  }

  # Try to download programming strategy
  programming_strategy <- try_download(
    file.path(github_raw_base, "programming_strategy.Rmd"),
    "programming_strategy.Rmd"
  )

  # Try to download unit test guidance
  unit_test_guidance <- try_download(
    file.path(github_raw_base, "unit_test_guidance.Rmd"),
    "unit_test_guidance.Rmd"
  )

  # Summary
  downloaded_count <- sum(!is.null(programming_strategy), !is.null(unit_test_guidance))
  cat(glue("Successfully downloaded {downloaded_count}/2 admiraldev vignettes\n"))

  list(
    programming_strategy = programming_strategy,
    unit_test_guidance = unit_test_guidance
  )
}

# Fix relative links in downloaded vignette content to use full admiraldev URLs
fix_relative_links <- function(lines) {
  if (is.null(lines)) {
    return(lines)
  }
  admiraldev_articles <- "https://pharmaverse.github.io/admiraldev/articles"
  admiral_reference <- "https://pharmaverse.github.io/admiral/reference"

  # Replace relative links to admiraldev articles (e.g. programming_strategy.html#anchor)
  lines <- gsub(
    "\\(([a-z_]+\\.html)(#[^)]*)?\\)",
    paste0("(", admiraldev_articles, "/\\1\\2)"),
    lines
  )
  # Replace ../reference/index.html links
  lines <- gsub(
    "\\(\\.\\./(reference/[^)]+)\\)",
    paste0("(", "https://pharmaverse.github.io/admiral/", "\\1)"),
    lines
  )
  lines
}

# Create root AGENTS.md with actual admiraldev content
create_root_agent_md <- function(admiraldev_content) {
  # Create the file content as separate parts to avoid glue complexity
  header_lines <- c(
    "# Admiral Development Guidelines for AI Assistants",
    "",
    "This file provides context for AI coding assistants (GitHub Copilot, Gemini, Claude, etc.) about admiral ecosystem standards and best practices.",
    "",
    "**Auto-generated:** See commit history for last update date",
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
      fix_relative_links(admiraldev_content$programming_strategy),
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
    "For admiral unit testing guidelines, see `tests/testthat/AGENTS.md`",
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

# Create tests/testthat/AGENTS.md with actual admiraldev content
create_tests_agent_md <- function(admiraldev_content) {
  # Create header
  header_lines <- c(
    "# Admiral Unit Testing Guidelines for AI Assistants",
    "",
    "Context for AI assistants when working with admiral unit tests in `tests/testthat/`.",
    "",
    "**Auto-generated:** See commit history for last update date",
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
      fix_relative_links(admiraldev_content$unit_test_guidance),
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
    "**Auto-generated:** See commit history for last update date",
    "**Optimized for:** GitHub Copilot code completion",
    "**Complete guidelines:** See `AGENTS.md` files for full context",
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
    "*This file is optimized for GitHub Copilot code completion. For complete admiral development guidelines, see `AGENTS.md` and `tests/testthat/AGENTS.md`.*"
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
