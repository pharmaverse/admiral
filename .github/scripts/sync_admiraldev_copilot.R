#!/usr/bin/env Rscript

#' Sync admiraldev Documentation for AI Assistants (v2.0)
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
  tryCatch({
    programming_strategy <<- readLines(
      file.path(github_raw_base, "programming_strategy.Rmd"), 
      warn = FALSE
    )
    cat("✓ Downloaded programming_strategy.Rmd\n")
  }, error = function(e) {
    cat("⚠ Could not download programming_strategy.Rmd\n")
  })
  
  # Try to download unit test guidance  
  tryCatch({
    unit_test_guidance <<- readLines(
      file.path(github_raw_base, "unit_test_guidance.Rmd"),
      warn = FALSE
    )
    cat("✓ Downloaded unit_test_guidance.Rmd\n")
  }, error = function(e) {
    cat("⚠ Could not download unit_test_guidance.Rmd\n")
  })
  
  list(
    programming_strategy = programming_strategy,
    unit_test_guidance = unit_test_guidance
  )
}

# Create root AGENT.md with actual admiraldev content
create_root_agent_md <- function(admiraldev_content) {
  
  # Create header
  header <- glue("# Admiral Development Guidelines for AI Assistants

This file provides context for AI coding assistants (GitHub Copilot, Gemini, Claude, etc.) 
about admiral ecosystem standards and best practices.

**Auto-generated:** {format(Sys.time(), '%Y-%m-%d %H:%M:%S')}  
**Source:** admiraldev package vignettes  
**Update script:** `source('.github/scripts/sync_admiraldev_copilot.R')`

## Purpose

Help AI assistants provide admiral-compliant code suggestions:
- ✅ Consistent function naming and patterns  
- ✅ Proper argument handling and validation
- ✅ CDISC/ADaM context and conventions
- ✅ Testing and documentation standards

---

")

  # Add programming strategy content if available
  programming_section <- ""
  if (!is.null(admiraldev_content$programming_strategy)) {
    programming_section <- glue("
# Admiral Programming Strategy

**Source:** admiraldev programming_strategy.Rmd vignette

{paste(admiraldev_content$programming_strategy, collapse = '\n')}

---

")
  } else {
    # Fallback content if download failed
    programming_section <- glue("
# Admiral Programming Strategy

**Note:** Could not download latest admiraldev content. Using essential guidelines.  
**Full documentation:** https://pharmaverse.github.io/admiraldev/articles/programming_strategy.html

## Function Design Principles

- **Modularity:** All code follows a modular approach with clearly separated steps
- **Avoid Copy and Paste:** Similar code should be put into separate functions  
- **Checks:** Meaningful error messages with clear reference to failing input
- **Flexibility:** Functions should be as flexible as possible without reducing usability

## Function Naming Convention

Function names should start with a verb and use snake_case:

| Prefix | Purpose |
|--------|---------|
| `derive_*` | Add rows/columns to datasets |
| `derive_var_*` | Add single variable |
| `derive_vars_*` | Add multiple variables |
| `compute_*` | Vector operations, return vectors |
| `assert_*` / `warn_*` | Input validation |
| `filter_*` | Filter observations |

## Function Arguments

**Standard order:** `dataset`, `by_vars`, `order`, `new_var`, `filter`, additional args

**Key patterns:**
```r
# Variables as symbols, not strings
new_var = AVAL                    # ✓ Correct  
new_var = \"AVAL\"                # ✗ Wrong

# Multiple variables in exprs()
by_vars = exprs(USUBJID, PARAMCD)
order = exprs(AVISITN, desc(AVAL))
```

---

")
  }

  # Add testing reference
  testing_section <- glue("
## Testing Guidelines

For admiral unit testing guidelines, see `tests/testthat/AGENT.md`

Key testing patterns:
- Use `tribble()` for test data creation
- Test happy path and error conditions  
- Validate that output is ungrouped
- Test argument validation with meaningful error messages

---

")

  # Add footer
  footer <- glue("
## Full Documentation

- **Programming Strategy:** https://pharmaverse.github.io/admiraldev/articles/programming_strategy.html
- **Unit Test Guidance:** https://pharmaverse.github.io/admiraldev/articles/unit_test_guidance.html  
- **Admiral Website:** https://pharmaverse.github.io/admiral/

---

*This file helps AI assistants understand admiral patterns for better code suggestions.*
")

  # Combine all sections
  full_content <- paste0(header, programming_section, testing_section, footer)
  
  writeLines(full_content, agent_md_root)
  cat(glue("✓ Created: {agent_md_root}\n"))
}

# Create tests/testthat/AGENT.md with actual admiraldev content
create_tests_agent_md <- function(admiraldev_content) {
  
  # Create header
  header <- glue("# Admiral Unit Testing Guidelines for AI Assistants

Context for AI assistants when working with admiral unit tests in `tests/testthat/`.

**Auto-generated:** {format(Sys.time(), '%Y-%m-%d %H:%M:%S')}  
**Source:** admiraldev unit testing guidance vignette

---

")

  # Add unit testing content if available
  testing_section <- ""
  if (!is.null(admiraldev_content$unit_test_guidance)) {
    testing_section <- glue("
# Admiral Unit Testing Guidance

**Source:** admiraldev unit_test_guidance.Rmd vignette

{paste(admiraldev_content$unit_test_guidance, collapse = '\n')}

---

")
  } else {
    # Fallback content if download failed  
    testing_section <- glue("
# Admiral Unit Testing Guidelines

**Note:** Could not download latest admiraldev content. Using essential guidelines.  
**Full documentation:** https://pharmaverse.github.io/admiraldev/articles/unit_test_guidance.html

## Test Structure and Organization

### File Naming and Organization
- Match function names: `test-derive_var_example.R`
- One test file per main function (can include helpers)
- Use descriptive test names explaining the scenario

### Test Framework Pattern
```r
# Standard admiral test structure
test_that(\"function_name() works with basic input\", {{
  # Setup test data using tribble()
  input <- tribble(
    ~USUBJID, ~PARAMCD, ~AVAL,
    \"S001\",   \"TEMP\",   37.1,
    \"S002\",   \"TEMP\",   36.8
  )
  
  # Call function
  result <- my_function(input, new_var = DERIVED)
  
  # Test expectations
  expect_equal(result$DERIVED, expected_values)
  expect_true(\"DERIVED\" %in% names(result))
  expect_equal(nrow(result), nrow(input))
}})
```

## Admiral Test Requirements

### Must Test Categories
1. **Happy path** - basic functionality works
2. **Error conditions** - invalid inputs fail appropriately
3. **Edge cases** - empty data, boundary conditions  
4. **Argument validation** - all parameters properly validated
5. **Output structure** - correct columns, types, ungrouped

### Standard Test Patterns

**Test argument validation:**
```r
test_that(\"function errors with invalid arguments\", {{
  data <- tribble(~USUBJID, ~AVAL, \"S001\", 10)
  
  # Symbol vs string validation
  expect_error(
    my_function(data, new_var = \"INVALID\"),
    \"new_var.*must be a symbol\"
  )
  
  # Required variables  
  expect_error(
    my_function(select(data, -AVAL)),
    \"Required variable.*AVAL.*missing\"
  )
}})
```

**Test grouped data handling:**
```r  
test_that(\"function errors with grouped data\", {{
  grouped_data <- tribble(
    ~USUBJID, ~AVAL,
    \"S001\", 10
  ) %>% group_by(USUBJID)
  
  expect_error(
    my_function(grouped_data),
    \"grouped.*not supported\"
  )
}})
```

**Test output structure:**
```r
test_that(\"function returns correct output structure\", {{
  input <- tribble(~USUBJID, ~AVAL, \"S001\", 10)
  result <- my_function(input, new_var = DERIVED)
  
  # Output is ungrouped
  expect_false(is.grouped_df(result))
  
  # New variable created
  expect_true(\"DERIVED\" %in% names(result))
  
  # Original variables preserved
  expect_true(all(names(input) %in% names(result)))
}})
```

## Test Data Best Practices

### Use tribble() for Clarity
```r
# ✓ Preferred - easy to read and modify
test_data <- tribble(
  ~USUBJID, ~PARAMCD, ~AVAL, ~AVISITN,
  \"S001\",   \"TEMP\",   37.1,  1,
  \"S001\",   \"TEMP\",   37.5,  2,
  \"S002\",   \"PULSE\",  72,    1
)

# ✗ Avoid - harder to scan in test code  
test_data <- data.frame(
  USUBJID = c(\"S001\", \"S001\", \"S002\"),
  # ... etc
)
```

### Library Loading
```r
# Standard test file setup
library(dplyr, warn.conflicts = FALSE)
library(admiral)  # Load package being tested

# Use functions directly in tests, not package::function
# ✓ my_function(data)
# ✗ admiral::my_function(data)
```

## Error Testing Patterns

### Meaningful Error Messages
```r
# Test specific error content users will see
expect_error(
  my_function(invalid_input),
  \"AVAL.*required.*missing\"  # Key words for user understanding
)

# Test error classes if using custom errors
expect_error(
  my_function(invalid_input), 
  class = \"admiral_error\"
)
```

## Coverage Guidelines

- **Aim for 80%+ coverage** on all exported functions
- **Focus on logic branches** - different code paths
- **Test edge cases** - empty data, single rows, boundary values
- **Mock external dependencies** if any (rare in admiral)

---

")
  }

  # Add footer
  footer <- glue("
*For complete guidance: https://pharmaverse.github.io/admiraldev/articles/unit_test_guidance.html*
")

  # Combine all sections
  full_content <- paste0(header, testing_section, footer)
  
  writeLines(full_content, agent_md_tests)
  cat(glue("✓ Created: {agent_md_tests}\n"))
}

# Create optimized .github/copilot-instructions.md for GitHub Copilot
create_copilot_instructions <- function(admiraldev_content) {
  content <- glue("# GitHub Copilot Instructions - admiral Development

**Auto-generated:** {format(Sys.time(), '%Y-%m-%d %H:%M:%S')}  
**Optimized for:** GitHub Copilot code completion  
**Complete guidelines:** See `AGENT.md` files for full context

⚠️ **DO NOT EDIT MANUALLY** - Run `source('.github/scripts/sync_admiraldev_copilot.R')` to update

## Admiral Code Completion Patterns

### Function Names (verb_object_detail)
```r
# Derivation functions
derive_var_base()           # Single variable derivation
derive_vars_merged()        # Multiple variables from merge
derive_param_tte()          # Parameter derivation

# Computation functions  
compute_age_years()         # Vector computation
compute_duration()          # Duration calculation

# Validation functions
assert_data_frame()         # Input validation
assert_vars()              # Variable validation
```

### Argument Patterns
```r
# Standard admiral function signature
my_function(
  dataset,                  # Always first
  by_vars = exprs(...),    # Grouping
  order = exprs(...),      # Sorting  
  new_var = VAR,           # New variable (symbol)
  filter = PARAM == \"VAL\", # Filtering
  set_values_to = exprs(...)# Value setting
)
```

### Common Expressions
```r
# Subject identification
subject_keys = exprs(STUDYID, USUBJID)

# Grouping patterns
by_vars = exprs(USUBJID, PARAMCD)
by_vars = exprs(USUBJID, PARAMCD, AVISIT)

# Ordering patterns  
order = exprs(AVISITN, ADTM)
order = exprs(PARAMCD, desc(ADY))

# Filtering patterns
filter = PARAMCD == \"TEMP\"
filter = AVISITN == 1
filter = !is.na(AVAL)
```

### Test Patterns
```r
# Test data with tribble
test_data <- tribble(
  ~USUBJID, ~PARAMCD, ~AVAL, ~AVISITN,
  \"S001\",   \"TEMP\",   37.1,  1,
  \"S001\",   \"TEMP\",   37.5,  2
)

# Standard test structure  
test_that(\"function_name() works with basic input\", {{
  result <- my_function(test_data, new_var = DERIVED)
  expect_true(\"DERIVED\" %in% names(result))
  expect_equal(nrow(result), nrow(test_data))
}})

# Error testing
test_that(\"function_name() errors with invalid input\", {{
  expect_error(
    my_function(invalid_data),
    \"Required variable.*missing\"
  )
}})
```

### Documentation Template
```r
#' Derive Example Variable
#'
#' @param dataset Input dataset
#' @param new_var New variable name
#'
#' @returns Dataset with new_var added
#'
#' @family der_gen  
#' @keywords der_gen
#'
#' @export
#'
#' @examples
#' library(dplyr, warn.conflicts = FALSE)
#' 
#' data <- tribble(~USUBJID, ~AVAL, \"S001\", 10)
#' derive_var_example(data, new_var = DERIVED)
```

---

*This file is optimized for GitHub Copilot code completion. For complete admiral development guidelines, see `AGENT.md` and `tests/testthat/AGENT.md`.*
")

  writeLines(content, copilot_instructions)
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
