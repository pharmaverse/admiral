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

