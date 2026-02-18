# GitHub Copilot Instructions - admiral Development

**Auto-generated:** 2026-02-18 18:26:30  
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
  filter = PARAM == "VAL", # Filtering
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
filter = PARAMCD == "TEMP"
filter = AVISITN == 1
filter = !is.na(AVAL)
```

### Test Patterns
```r
# Test data with tribble
test_data <- tribble(
  ~USUBJID, ~PARAMCD, ~AVAL, ~AVISITN,
  "S001",   "TEMP",   37.1,  1,
  "S001",   "TEMP",   37.5,  2
)

# Standard test structure  
test_that("function_name() works with basic input", {
  result <- my_function(test_data, new_var = DERIVED)
  expect_true("DERIVED" %in% names(result))
  expect_equal(nrow(result), nrow(test_data))
})

# Error testing
test_that("function_name() errors with invalid input", {
  expect_error(
    my_function(invalid_data),
    "Required variable.*missing"
  )
})
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
#' data <- tribble(~USUBJID, ~AVAL, "S001", 10)
#' derive_var_example(data, new_var = DERIVED)
```

---

*This file is optimized for GitHub Copilot code completion. For complete admiral development guidelines, see `AGENT.md` and `tests/testthat/AGENT.md`.*
