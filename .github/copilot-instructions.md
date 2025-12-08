# GitHub Copilot Instructions for admiral

## Repository Overview

`{admiral}` (ADaM in R Asset Library) is an open-source R package that provides a modular toolbox for creating CDISC-compliant Analysis Data Model (ADaM) datasets. This package is part of the pharmaverse ecosystem and follows strict pharmaceutical industry standards for clinical trial data analysis.

**Key Points:**
- This is a core package in the admiral family with multiple extension packages for specific therapeutic areas
- All code must be production-ready and suitable for regulatory submission
- The package emphasizes readability, modularity, and comprehensive documentation
- Development follows the tidyverse style guide and admiral-specific conventions

## Essential Resources

Before making any code changes, familiarize yourself with these critical documents from `{admiraldev}`:

1. **[Programming Strategy](https://pharmaverse.github.io/admiraldev/articles/programming_strategy.html)** - MANDATORY reading for all code contributions. Covers:
   - Function design principles
   - Input/output conventions
   - Error handling and assertions
   - Naming conventions
   - Deprecation guidance
   - Function categorization and keywords

2. **[Unit Test Guidance](https://pharmaverse.github.io/admiraldev/articles/unit_test_guidance.html)** - Required for all new functions:
   - Writing comprehensive unit tests
   - Test structure and naming conventions
   - Testing edge cases and error conditions
   - Using testthat framework

3. **[Writing Vignettes](https://pharmaverse.github.io/admiraldev/articles/writing_vignettes.html)** - For documentation updates:
   - Vignette structure and style
   - Code examples and best practices
   - Integration with existing documentation

4. **[Pull Request Review Guidance](https://pharmaverse.github.io/admiraldev/articles/pr_review_guidance.html)** - Before submitting PRs

## Package Structure

```
admiral/
├── R/                       # Source code (all .R files)
├── tests/testthat/          # Unit tests (test-*.R files)
├── man/                     # Auto-generated documentation (DO NOT EDIT MANUALLY)
├── vignettes/               # Long-form documentation (.Rmd files)
├── inst/templates/          # ADaM dataset templates (ad_*.R files)
├── data/                    # Package data
├── DESCRIPTION              # Package metadata
├── NAMESPACE                # Auto-generated exports (DO NOT EDIT MANUALLY)
└── NEWS.md                  # Changelog (UPDATE for all changes)
```

## Code Style and Conventions

### General Principles (Admiral Manifesto)
1. **Usability** - Documentation is paramount; functions must be easy to understand and use
2. **Simplicity** - Functions have clear, focused purposes; avoid multi-purpose functions
3. **Findability** - Consistent naming conventions; avoid redundant functions
4. **Readability** - Follow Programming Strategy; code should be understandable to reviewers

### Coding Standards

**Style:**
- Follow the [tidyverse style guide](https://style.tidyverse.org/)
- Use tidyverse packages (dplyr, tidyr, purrr, stringr, etc.) over base R equivalents
- Prefer `%>%` pipe operator for readability
- Maximum line length: 100 characters (enforced by lintr)
- Use snake_case for function and variable names

**Function Design:**
- All functions must have roxygen2 documentation with `@export`, `@param`, `@return`, `@examples`
- Include `@keywords` and `@family` tags (see Programming Strategy for categorization)
- Start functions with comprehensive input validation using `assert_*()` functions from admiraldev
- Functions should not modify global state or have unexpected side effects
- Use descriptive parameter names; avoid abbreviations unless standard in CDISC

**Common Patterns:**
```r
# Example function structure
#' Brief Title
#'
#' Detailed description
#'
#' @param dataset Input dataset (`data.frame`)
#' @param new_var Name of new variable to create
#'
#' @return Input dataset with additional column
#'
#' @family utils_fmt
#' @keywords utils_fmt
#'
#' @export
#'
#' @examples
#' library(dplyr)
#' # Example code here
my_function <- function(dataset, new_var) {
  # Input validation
  assert_data_frame(dataset)
  assert_character_scalar(new_var)
  
  # Main logic using tidyverse
  dataset %>%
    mutate(!!sym(new_var) := computed_value)
}
```

### File Organization

**R/ directory:**
- One function per file OR related functions grouped logically
- File names should match primary function name (e.g., `derive_vars_dt.R`)
- Utility functions go in specific utility files (e.g., `user_utils.R`, `dev_utils.R`)

**tests/testthat/ directory:**
- One test file per source file: `test-<source_file_name>.R`
- Use descriptive test names: `"function_name Test N: description"`
- Test structure:
  ```r
  test_that("function_name Test 1: expected behavior description", {
    # Arrange
    input <- ...
    expected <- ...
    
    # Act & Assert
    expect_equal(function_name(input), expected)
  })
  ```

## Development Workflow

### Before Starting
1. Create a GitHub issue describing the change
2. Create a feature branch from `main` following naming convention: `<issue_number>_<brief_description>`
3. Familiarize yourself with the Programming Strategy (linked above)

### Making Changes

**For New Functions:**
1. Write the function in appropriate R/ file following programming strategy
2. Add comprehensive roxygen2 documentation
3. Add unit tests in tests/testthat/
4. Add examples demonstrating usage
5. Run `devtools::document()` to update NAMESPACE and man/ files
6. Update NEWS.md with a bullet point describing the change
7. Ensure NAMESPACE exports are in alphabetical order
8. If relevant, update vignettes and/or inst/templates/

**For Bug Fixes:**
1. Add a test that reproduces the bug (should fail)
2. Fix the bug
3. Verify the test now passes
4. Update NEWS.md
5. Update documentation if behavior changed

### Building and Testing

**Essential Commands:**
```r
# Install dependencies
install.packages("devtools")
install.packages("admiraldev")

# Load package for development
devtools::load_all()

# Update documentation (run after changing roxygen comments)
devtools::document()

# Run tests
devtools::test()

# Run specific test file
devtools::test_file("tests/testthat/test-<filename>.R")

# Check package (R CMD check)
devtools::check()

# Build pkgdown site locally
pkgdown::build_site()

# Code style check
lintr::lint_package()

# Style code automatically
styler::style_pkg()
```

**Before Committing:**
- [ ] Run `devtools::document()` if you modified roxygen comments
- [ ] Run `devtools::test()` - all tests must pass
- [ ] Run `devtools::check()` - should return 0 errors, 0 warnings, 0 notes
- [ ] Run `lintr::lint_package()` - fix any linting issues
- [ ] Run `styler::style_pkg()` if needed
- [ ] Update NEWS.md with your changes
- [ ] Build and review pkgdown site if you changed documentation: `pkgdown::build_site()`

## Documentation Requirements

### Roxygen2 Headers
Every exported function MUST include:
- `@param` for each parameter with description and expected type
- `@return` describing what the function returns
- `@export` tag for exported functions
- `@keywords` for categorization (see Programming Strategy)
- `@family` for grouping related functions
- `@examples` with working code examples that demonstrate usage
- Consider adding `@details` for complex functions

### Examples
- Must be executable code that runs without error
- Use existing test data from `{pharmaversesdtm}` package when possible
- Include `library()` calls for required packages (admiral, dplyr, etc.)
- Show realistic use cases
- Keep examples concise but comprehensive

### Vignettes
- Written in R Markdown (.Rmd)
- Follow structure outlined in Writing Vignettes guide
- Include working code examples with output
- Reference relevant functions using markdown links
- Update existing vignettes if your changes affect them

## Testing Guidelines

### Test Coverage
- **Every function must have unit tests**
- Test normal expected behavior
- Test edge cases (empty inputs, NAs, boundary values)
- Test error conditions (invalid inputs should trigger informative errors)
- Aim for >80% code coverage

### Test Structure
```r
# Pattern for tests
# function_name ----
## Test 1: description ----
test_that("function_name Test 1: description", {
  input <- tibble::tribble(
    ~VAR1, ~VAR2,
    "A",   1,
    "B",   2
  )
  
  expected <- tibble::tribble(
    ~VAR1, ~VAR2, ~NEW_VAR,
    "A",   1,     "computed_A",
    "B",   2,     "computed_B"
  )
  
  expect_equal(function_name(input), expected)
})

## Test 2: error handling ----
test_that("function_name Test 2: errors on invalid input", {
  expect_error(
    function_name(invalid_input),
    "expected error message pattern"
  )
})
```

### Snapshot Testing
- For functions with complex output, consider using `testthat::expect_snapshot()`
- Snapshots stored in `tests/testthat/_snaps/`

## Common Patterns and Conventions

### Working with Data
```r
# Use tidyverse for data manipulation
dataset %>%
  filter(!is.na(PARAM)) %>%
  mutate(NEW_VAR = case_when(...)) %>%
  arrange(USUBJID, ADT)

# Use admiral's derive_* functions for ADaM derivations
dataset <- dataset %>%
  derive_vars_dt(
    dtc = EXSTDTC,
    new_vars_prefix = "EXST"
  )
```

### Handling Missing Values
- SAS blank strings (`""`) should be converted to `NA` using `convert_blanks_to_na()`
- This is standard practice in all admiral templates
- Use `na.rm = TRUE` thoughtfully; document behavior with NAs

### Assertions
Use admiraldev assertion functions extensively:
```r
# Common assertions
assert_data_frame(dataset)
assert_character_scalar(var_name)
assert_vars(by_vars)
assert_date_vector(date_var)
assert_logical_scalar(flag, optional = TRUE)
```

### Variable Naming
- Follow CDISC standards for ADaM variable names (uppercase)
- Use descriptive names for internal variables (snake_case)
- Prefix parameters that take variable names with `*_var` or `*_vars`

## NAMESPACE Management

- **NEVER edit NAMESPACE manually** - it's auto-generated by roxygen2
- Exports MUST be in alphabetical order
- Run `devtools::document()` to regenerate
- Add `@export` tag in roxygen comments to export a function
- If you see NAMESPACE conflicts, check that your roxygen tags are correct

## NEWS.md Updates

Always update NEWS.md for your changes:

```md
# admiral 1.x.x.xxxx

## New Features
- New function `new_function()` to do X (#issue_number)

## Updates of Existing Functions
- `existing_function()` now handles Y differently (#issue_number)

## Bug Fixes  
- Fixed bug in `buggy_function()` where Z occurred (#issue_number)

## Documentation
- Updated vignette X to include Y (#issue_number)
```

## Pull Request Checklist

Before submitting a PR, verify:
- [ ] All tests pass (`devtools::test()`)
- [ ] R CMD check passes with 0 errors, 0 warnings, 0 notes (`devtools::check()`)
- [ ] Code follows tidyverse style guide
- [ ] All new functions have comprehensive roxygen documentation
- [ ] All new functions have unit tests
- [ ] NAMESPACE is properly updated (via `devtools::document()`)
- [ ] NEWS.md is updated
- [ ] If applicable: vignettes are updated
- [ ] If applicable: inst/templates/ are updated
- [ ] If applicable: deprecated functions follow deprecation guidance
- [ ] pkgdown site builds successfully (`pkgdown::build_site()`)
- [ ] No lintr warnings (`lintr::lint_package()`)

## Special Considerations

### Deprecation
When removing or replacing functions:
- Follow the [deprecation guidance](https://pharmaverse.github.io/admiraldev/articles/programming_strategy.html#deprecation)
- Mark as deprecated for at least one release cycle
- Use `lifecycle` package functions for deprecation warnings
- Update documentation to point users to replacement

### Dependencies
- Minimize new dependencies
- Use packages already in DESCRIPTION when possible
- Justify any new dependencies in the PR
- Prefer tidyverse packages when adding new dependencies

### Performance
- Consider performance for large datasets (typical clinical trials have 1000s of rows)
- Avoid loops when vectorization is possible
- Use dplyr for efficient data manipulation

### ADaM Standards
- Follow CDISC ADaM Implementation Guide
- Variable names, labels, and metadata must follow ADaM conventions
- Consult existing admiral functions for patterns

## Getting Help

- **Slack**: [pharmaverse Slack workspace](https://pharmaverse.slack.com/) - #admiral channel
- **GitHub Issues**: [admiral issues](https://github.com/pharmaverse/admiral/issues)
- **Documentation**: [admiral website](https://pharmaverse.github.io/admiral/)
- **Programming Strategy**: Reference link at top of this document

## Quick Reference

**Most Important Principle**: Readability and documentation are paramount. Code will be reviewed by health authorities and users need to understand and trust the implementations. When in doubt, prioritize clarity over cleverness.

**Key Files to Review for Patterns:**
- `R/user_utils.R` - Core utility functions with S3 methods
- `R/derive_vars_dt.R` - Date derivations (common pattern)
- `tests/testthat/test-user_utils.R` - Testing patterns
- `inst/templates/ad_adsl.R` - Complete ADaM dataset template
- `vignettes/adsl.Rmd` - Vignette example

**Remember**: This package is used in regulatory submissions. Quality, correctness, and documentation are more important than speed of implementation.
