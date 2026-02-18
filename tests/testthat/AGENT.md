# Admiral Unit Testing Guidelines for AI Assistants

Context for AI assistants when working with admiral unit tests in `tests/testthat/`.

**Auto-generated:** 2026-02-18 18:26:30  
**Source:** admiraldev unit testing guidance vignette

---
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
test_that("function_name() works with basic input", {
  # Setup test data using tribble()
  input <- tribble(
    ~USUBJID, ~PARAMCD, ~AVAL,
    "S001",   "TEMP",   37.1,
    "S002",   "TEMP",   36.8
  )
  
  # Call function
  result <- my_function(input, new_var = DERIVED)
  
  # Test expectations
  expect_equal(result$DERIVED, expected_values)
  expect_true("DERIVED" %in% names(result))
  expect_equal(nrow(result), nrow(input))
})
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
test_that("function errors with invalid arguments", {
  data <- tribble(~USUBJID, ~AVAL, "S001", 10)
  
  # Symbol vs string validation
  expect_error(
    my_function(data, new_var = "INVALID"),
    "new_var.*must be a symbol"
  )
  
  # Required variables  
  expect_error(
    my_function(select(data, -AVAL)),
    "Required variable.*AVAL.*missing"
  )
})
```

**Test grouped data handling:**
```r  
test_that("function errors with grouped data", {
  grouped_data <- tribble(
    ~USUBJID, ~AVAL,
    "S001", 10
  ) %>% group_by(USUBJID)
  
  expect_error(
    my_function(grouped_data),
    "grouped.*not supported"
  )
})
```

**Test output structure:**
```r
test_that("function returns correct output structure", {
  input <- tribble(~USUBJID, ~AVAL, "S001", 10)
  result <- my_function(input, new_var = DERIVED)
  
  # Output is ungrouped
  expect_false(is.grouped_df(result))
  
  # New variable created
  expect_true("DERIVED" %in% names(result))
  
  # Original variables preserved
  expect_true(all(names(input) %in% names(result)))
})
```

## Test Data Best Practices

### Use tribble() for Clarity
```r
# ✓ Preferred - easy to read and modify
test_data <- tribble(
  ~USUBJID, ~PARAMCD, ~AVAL, ~AVISITN,
  "S001",   "TEMP",   37.1,  1,
  "S001",   "TEMP",   37.5,  2,
  "S002",   "PULSE",  72,    1
)

# ✗ Avoid - harder to scan in test code  
test_data <- data.frame(
  USUBJID = c("S001", "S001", "S002"),
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
  "AVAL.*required.*missing"  # Key words for user understanding
)

# Test error classes if using custom errors
expect_error(
  my_function(invalid_input), 
  class = "admiral_error"
)
```

## Coverage Guidelines

- **Aim for 80%+ coverage** on all exported functions
- **Focus on logic branches** - different code paths
- **Test edge cases** - empty data, single rows, boundary values
- **Mock external dependencies** if any (rare in admiral)

---
*For complete guidance: https://pharmaverse.github.io/admiraldev/articles/unit_test_guidance.html*
