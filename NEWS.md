# admiral 0.3.0

- Several functions have been renamed. The old names are now deprecated. They can still be used but a warning will be issued (#507). 
    - `derive_aage()` -> `derive_vars_aage()`
    - `derive_duration()` -> `derive_vars_duration()`
    - `derive_query_vars()` -> `derive_vars_query()`
    - `derive_suppqual_vars()` -> `derive_vars_suppqual()`

- The Programming Strategy has been updated (#495)

- The `date_var` and `traceabilty_vars` parameters of `dthcaus_source()` have been renamed to `date` and `traceability_vars`, respectively. The old names are now deprecated. They can still be used but a warning will be issued (#493).

- The format of the `date_imputation` parameter of `derive_vars_dt()` and `derive_vars_dtm()` has been changed from "dd-mm" to "mm-dd". Thus, "01-12" now refers to January 12th rather than December 1st (#492).

- The `filter_rows` parameter of `derive_summary_records()` has been renamed to `filter`. The `fns` parameter has been deprecated in favor of `analysis_var` and `summary_fun` (#491).

- The `flag_filter` parameter of `derive_extreme_flag()` has been renamed to `filter` (#487).

- `derive_var_anrind()` checks whether the `AVAL` variable is present in the input dataset (#486).

- "OCCDS" has been added as a new section on the reference page on the package website (#485)

- `convert_blanks_to_na()` enables the conversion of SAS blanks, i.e. `""`, into proper R `NA` values (#482).

- `derive_agegr_fda()` used to return `NA` for ages less than or equal 18. It now returns `<=18` (#441).

- A search feature has been added to the package website (#438).

- `use_ad_template()` has been refactored to no longer make use of the {usethis} package which is no longer a dependency of {admiral} (#433).

- The `drop_values_from` parameter has been removed from `derive_summary_records()` (#425).

- All derivation functions check whether the input dataset is grouped and throw an error if it is (#408).

- `call_derivation()` enables users to call the same function multiple times with some parameters being fixed across iterations and others varying (#403)

- `derive_var_atirel()` enables the derivation of the "Analysis Time Relative to Reference" (#397).

- `derive_vars_atc()` can be used to add ATC variables from FACM to ADCM (#396).

- A performance issue in `derive_vars_dt()` has been resolved (#384).

- `derive_vars_dtm_to_dt()` enables the easy conversion of datetime to date variables (#376).

- `derive_param_bmi()`, `derive_param_bsa()` and `derive_param_map()` enables the derivation of the body mass index, body surface area and mean arterial pressure parameters respectively (#368).

- `derive_param_qtc()` enables the derivation of corrected QT intervals according to the formula of Bazett, Fridericia or Sagie (#325).

- `derive_var_ontrtfl()` can now handle events with a start and end date rather than just a single assessment date (#395).

- A template for ADCM has been created (#268)

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
