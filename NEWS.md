# admiral 0.6.0

## New Features

- `derive_vars_dy()` derives the analysis day from one or more `--DT(M)` variables
(#700)

## Updates of Existing Functions

- The `derive_last_dose()` function has been split into a general function 
`derive_vars_last_dose()` and three wrapper functions `derive_var_last_dose_amt()`, 
`derive_var_last_dose_date()`, and `derive_var_last_dose_grp()` (#385)

- `derive_var_ontrtfl()` now has a `new_var` parameter to support the derivation of `ONTRxxFL` and `ONTRTwFL` variables (#721)

- `derive_vars_dtm()`, `derive_var_disposition` and `derive_var_lstalvdt` now have `preserve` argument.  A user can preserve partial dates when doing date imputation, e.g. `2019---07` would become `2019-06-07` by setting `preserve` to `TRUE` when doing date_imputation (#592)

- `derive_vars_dtm()` now has `ignore_seconds_flag` argument so users can suppress `"S"` flag if seconds are not recorded in the data (#589)

## Breaking Changes

- `derive_agegr_ema()`, `derive_agegr_fda()`, `derive_disposition_dt()`,
`derive_disposition_status()`,`derive_extreme_flag()`, `derive_worst_flag()`,
`derive_obs_number()`, `derive_disposition_reason()` have been deprecated and
renamed in favor of `derive_var_agegr_ema()`, `derive_var_agegr_fda()`,
`derive_var_disposition_dt()`, `derive_var_disposition_status()`,
`derive_var_extreme_flag()`, `derive_var_worst_flag()`,
`derive_var_last_dose()`, `derive_var_obs_number()`, and
`derive_vars_disposition_reason()` respectively (#738)

- `derive_var_basec()` and `derive_baseline()` have been deprecated in favor of the extended `derive_var_base()` function (#695)

- `derive_params_exposure()` has been deprecated and renamed as `derive_param_exposure()` (#722)

- The `derive_last_dose()` function has been deprecated in favor of
`derive_var_last_dose_date()` (#385)

- The behavior of all functions providing the `date_imputation` parameter, e.g., `derive_vars_dtm()` and
`derive_vars_dt()`, has changed for `date_imputation = "mid"`. Before the date was imputed as June 15th
if both month and day were missing. Now it is imputed as June 30th. For the old behavior please specify
`date_imputation = "06-15"`. Please note the behavior has not changed if only the day is missing. In
this case the day is imputed as `15` (#592)

- `derive_var_ontrtfl()` now has a `new_var` parameter to support the derivation of `ONTRxxFL` and `ONTRTwFL` variables (#721)

- The following functions and parameters, which were deprecated in previous {admiral} versions, were removed (#513):
  
  - `derive_aage()`, `derive_duration()`, `derive_query_vars()`, and
  `derive_suppqual_vars()` function
  - `fns` and `filter_rows` parameters in `derive_summary_records()`
  - `date_var` and `traceabilty_vars` parameters in `dthcaus_source()`
  - `flag_filter` parameter in `derive_extreme_flag()`
  - `flag_filter` parameter in `derive_var_extreme_flag()`
  - `date_var` parameter in `lstalvdt_source()`
  - `date` parameter in `derive_var_ontrtfl()`
  
- `derive_var_agegr_fda()` has been updated to use ranges <18, 18-64, >=65 (#829)

## Documentation

- README and site homepage has been updated with important new section around expectations of {admiral}, as well as other useful references such as links to conference talks (#868 & #802)

- New vignette [Development Process](../articles/development_process.html) and improvements made to contribution vignettes (#765 & #758)

- Updated [Pull Request Review Guidance](../articles/pr_review_guidance.html) on using `task-list-completed` workflow (#817)

## Various

- GitHub repo moved to pharmaverse org and associated broken site links fixed (#803 & #820)

- Examples have been added for `format_reason_default`, `format_eoxxstt_default`, `extend_source_datasets` and `filter_date_sources` (#745)

# admiral 0.5.0

- The first truly open source release licensed under Apache 2.0 (#680)

- New vignette [Contributing to admiral](../articles/contribution_model.html) (#679)

- New vignette [Unit Test Guidance](../articles/unit_test_guidance.html) (#679)

- Broken links in README have been fixed (#564)

# admiral 0.4.0

## New Features

### General

- `derive_vars_dtm_to_tm()` enables the easy conversion of datetime to time variables (#551)

- `derive_var_age_years()` derives age in years from a variable providing the
age in different units (#569)

### BDS

- `derive_param_tte()` derives time-to-event-parameters (#546)

- For common time-to-event endpoints [event and censoring source
objects](../reference/index.html#section-pre-defined-time-to-event-sources) are
provided (#612)

### Developer

- `assert_list_element()` checks if an element of a list of lists/classes
fulfills a condition

- `assert_one_to_one()` checks if there is a one to one mapping between two
lists of variables

- `negate_vars()` negates a list of variables to remove them from a dataset with
`select()`

## Updates of Existing Functions

- Unit checks in `derive_param_*()` functions are no longer case sensitive (#631)

- `derive_agegr_ema()` and `derive_agegr_fda()` gain a `age_unit` parameter used
to specify the unit of the input age (#569)

## Breaking Changes

- All SDTM datasets have been moved to the {admiral.test} package (#639)

- The `min_dates` and `max_dates` parameters of `derive_vars_dt()` and `derive_vars_dtm()` no longer expect a `list()` but `vars()` as input (#405)

## Bug Fixes

- `derive_vars_dtm()` no longer shifts the time of the input `--DTC` variable (#436)

## Documentation

- New vignette [Creating a BDS Time-to-Event ADaM](../articles/bds_tte.html) (#549)

- New vignette [Queries Dataset Documentation](../articles/queries_dataset.html) (#561)

- New vignette [Writing Vignettes](../articles/writing_vignettes.html) (#334)

- New vignette [Pull Request Review Guidance](../articles/pr_review_guidance.html) (#554)

- A section on handling missing values when working with {admiral} has been added to the "Get Started" vignette (#577)

- Package installation instructions have been added to the README (#558)

- The documentation of `derive_vars_dtm()` falsely stated that the `flag_imputation` parameter should be either `TRUE` or `FALSE`. It now correctly states that the possible values are `"time"`, `"date"` or `"auto"` (#539)


# admiral 0.3.0

## New Features

### General

- `convert_blanks_to_na()` can be used to convert SAS blanks, i.e. `""`, into proper R `NA` values (#482)

- `call_derivation()` enables users to call the same function multiple times with some parameters being fixed across iterations and others varying (#403)

- `derive_vars_dtm_to_dt()` enables the easy conversion of datetime to date variables (#376)

- `derive_var_ontrtfl()` can now handle events with a start and end date rather than just a single assessment date (#395)

- `derive_worst_flag()` enables the flagging of worst records (#300)

### BDS

- `derive_derived_param()` can be used to derive a new parameter from existing parameters in a BDS dataset (#325)

- `derive_param_bmi()`, `derive_param_bsa()` and `derive_param_map()` enables the derivation of the body mass index, body surface area and mean arterial pressure parameters respectively (#368)

- `derive_param_qtc()` enables the derivation of corrected QT intervals according to the formula of Bazett, Fridericia or Sagie (#325)

- `derive_param_rr()` enables the derivation of the RR interval (#325)

- `derive_params_exposure()` enables the derivation of summary exposure parameters (#400)

- `derive_param_doseint()` enables the derivation of dose intensity (#179)

### OCCDS

- `derive_var_atirel()` enables the derivation of the "Analysis Time Relative to Reference" (#397)

- `derive_vars_atc()` can be used to add ATC variables from FACM to ADCM (#396)

## Updates of Existing Functions

- `derive_var_anrind()` now checks whether the `AVAL` variable is present in the input dataset (#486)

- All derivation functions check whether the input dataset is grouped and throw an error if it is (#408)

- `use_ad_template()` has been refactored to no longer make use of the {usethis} package which is no longer a dependency of {admiral} (#433)

- A performance issue in `derive_vars_dt()` has been resolved (#384)

## Breaking Changes

- The `drop_values_from` parameter has been removed from `derive_summary_records()` (#425)

- The format of the `date_imputation` parameter of `derive_vars_dt()` and `derive_vars_dtm()` has been changed from "dd-mm" to "mm-dd". Thus, "01-12" now refers to January 12th rather than December 1st (#492)

- Several functions have been renamed. The old names are now deprecated. They can still be used but a warning will be issued (#507)
    - `derive_aage()` -> `derive_vars_aage()`
    - `derive_duration()` -> `derive_vars_duration()`
    - `derive_query_vars()` -> `derive_vars_query()`
    - `derive_suppqual_vars()` -> `derive_vars_suppqual()`
    
- The `date_var` parameter of `lstalvdt_source()` has been renamed to `date`

- The `filter_rows` parameter of `derive_summary_records()` has been renamed to `filter`. The `fns` parameter has been deprecated in favor of `analysis_var` and `summary_fun` (#491)

- The `date_var` and `traceabilty_vars` parameters of `dthcaus_source()` have been renamed to `date` and `traceability_vars`, respectively (#493)

- The `flag_filter` parameter of `derive_extreme_flag()` has been renamed to `filter` (#487)

## Bug Fixes

- `derive_agegr_fda()` used to return `NA` for ages less than or equal 18. It now returns `<=18` (#441)

## Documentation

- New vignette on "Date and Time Imputation" has been created (#198)

- A "Guidance for git Usage" has been created (#266)

- "OCCDS" has been added as a new section on the reference page on the package website (#485)

- The Programming Strategy has been updated (#495)

- A search feature has been added to the package website (#438)

- New template scripts for ADEX (#181), ADCM (#268) and ADEG (#258) have been created

- New vignette for programming ADEX has been created (#372)

- A section on how to create query variables (e.g. SMQs in ADAE) has been added to the Occurrence datasets vignette (#370)

- The BDS vignette has been updated to incorporate examples of ADVS and ADEG specific functions (#371)


# admiral 0.2.1

- Fixed a critical bug in `use_ad_template()` that prevented the function from being usable at all (#326)


# admiral 0.2.0

## New Features

### General

- Function argument checks have been completely re-written to provide clearer error messages to users (#263, #288)

- SDTM `SUPP--` datasets can be merged onto their parent domain using `derive_suppqual_vars()` (#145)

- In case a derivation detects duplicate records after applying a `filter`, the dataset of duplicate records is made available to users via `get_duplicates_dataset()` (#202)

- `derive_vars_dt()` and `derive_vars_dtm()` gain a `min_dates` and `max_dates` parameter which can be used to ensure that the imputed date(time) is not before the `min_dates` nor after the `max_dates`, e.g. avoid that `AENDT` is after the data cut date or `ASTDT` is before the first treatment date (#158)

- `use_ad_template()` can be used to open a template script for an ADaM dataset; all available templates can be displayed using `list_all_templates()` (#110)

### ADSL

- EMA and FDA defined age groupings can be derived using `derive_agegr_ema()` and `derive_agegr_fda()`, respectively

- Disposition Status can be derived using `derive_disposition_status()` (#92)

- Disposition Reason can be derived using `derive_disposition_reason()`

- Disposition Dates can be derived using `derive_disposition_dt()` (#91)

- Date Last Known Alive can be derived using `derive_var_lstalvdt()` (#94)

- Cause of Death can be derived using `derive_var_dthcaus()` (#93)

### BDS

- Summary records for BDS datasets, e.g. with `DTYPE == "AVERAGE"`, can be derived using `derive_summary_records()` (#177)

### OCCDS

- Last Dose Date(time) can be derived using `derive_last_dose()`

## Breaking Changes

- `derive_merged_vars()` has been removed from {admiral} in favor of smaller special purpose functions, e.g. `derive_disposition_status()` (#167)

- Function arguments no longer accept expressions created with `expr()` or `exprs()` as inputs; instead filter expressions can be passed "as is" and multiple variables have to be wrapped inside `vars()` (#187)

  **Old:**

  ```r
  derive_extreme_flag(
    vs,
    new_var = LASTFL,
    by_vars = exprs(USUBJID, VSTESTCD, VISIT),
    order = exprs(VSTPTNUM),
    mode = "last",
    flag_filter = expr(VISIT != "BASELINE")
  )
  ```

  **New:**

  ```r
  derive_extreme_flag(
    vs,
    new_var = LASTFL,
    by_vars = vars(USUBJID, VSTESTCD, VISIT),
    order = vars(VSTPTNUM),
    mode = "last",
    flag_filter = VISIT != "BASELINE"
  )
  ```
  
- `read_dap_m3()` and `initialize()` have been migrated to {admiral.roche} (#272)

- The `start_date` and `end_date` parameters of `derive_var_ady()`, `derive_var_aendy()` and `derive_var_astdy()` have been renamed to `reference_date` and `date`, respectively (#121)


## Bug Fixes

- `derive_var_basetype()` no longer drops records which do not match any condition defined in the `basetype` argument (#226)

- Join warnings like "Column `USUBJID` has different attributes on LHS and RHS of join when using left_join()" are no longer displayed (#271)

- For datetimes with time imputed to "00:00:00" the time part is now displayed (#206)


## Documentation

- [Frequently Asked Questions](../articles/faq.html)

- [Creating ADSL](../articles/adsl.html)

- [Creating a BDS Finding ADaM](../articles/bds_finding.html)

- [Creating an OCCDS ADaM](../articles/occds.html)
