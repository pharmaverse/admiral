# admiraldev 0.3.0

## New Features
  - New function `process_set_values_to()` for creating the variables specified
  by the `set_value_to` argument and catching errors (#70)
  
## Updates of Existing Functions
  - Fix datatable styling for documentation (#197)
  - The `assert_character_vector()` function gained a `named` argument to check
  that all elements of the vector are named. (#70)
  - The `assert_list_of()` function gained a `named` argument to check that all
  elements of the list are named. (#203)
  - The `quote_fun` argument of `enumerate()` was extended such that `NULL` can
  be specified to request no quoting of the elements. (#203)
  - The `assert_list_of()` function was enhanced such that it also considers the
  type of the element, e.g., to check if a value is a list of symbols. (#208)

## Breaking Changes
- The default value of the `optional` argument in `assert_date_vector()`,
`assert_list_of()`, and `assert_s3_class()` was changed from `TRUE` to `FALSE`
to make the default behavior consistent. (#87)
- admiral functions no longer expect list of quosures created by `vars()` but
list of expressions created by `exprs()`. Thus the following functions and
arguments were deprecated:
    - `quo_c()` and `replace_symbol_in_quo()`
    - the `quosures` argument in `get_source_vars()`,
    `replace_values_by_names()`, and `vars2chr()`

## Documentation
  - New section in programming strategy regarding comments (#71)
  - Removed requirement to add `@author` tags to code scripts from programming 
  strategy, as we will only be tracking authors in the DESCRIPTION file. Authors
  have been removed from function documentation in line with this update. 
  (#206, #210)

# admiraldev 0.2.0

## New Features
  - Developer addin for formatting tests to admiral programming standards (#73)
  - New functions `replace_symbol_in_quo()` and `add_suffix_to_vars()` (#106)
  - New function `assert_atomic_vector()` (#98)
  - New keyword/family `create_aux` for functions creating auxiliary datasets (#126)
  - New function `assert_date_vector()` (#129)
  - New function `assert_same_type()` (#176)
  - Remove dependency on `{assertthat}` (#149)
  - Test coverage for `admiraldev` have increased from 45% to approximately 100% (#94, #95, #96, #98, #101, #103)
  - _Environment_ objects were consolidated into a single `admiraldev_environment` object under `R/admiraldev_environment.R`. (#179)

## Updates of Existing Functions
  - `expect_names` argument added to `assert_vars()` to check if all variables are named (#117)
  - Remove `dplyr` function exports and migration of user facing function `negate_vars()` to admiral  (#83)
  
## Breaking Changes
  - No longer compatible with admiral (<0.9)
  
## Documentation
  - New vignette for our package release strategy (#79) 
  - Updated multiple roxygen headers (#116, #133, #134, #141, #145, #172)
  - Description on how admiral options work for certain function inputs, i.e `subject_keys` (#133)

## Various
  - PR Checklist Template updated (#172)
  - New authors/contributors (#158)

# admiraldev 0.1.0

## New Features

  - Developer specific functions brought over from `{admiral}`
  - Developer specific vignettes brought over from `{admiral}`
  - New `{admiraldev}` website created

## Updates of Existing Functions
  - NA
  
## Breaking Changes
  - NA

## Documentation
  - NA

## Various
  - NA

