# Changelog

## admiral (development version)

### Updates of Existing Functions

- Grading metadata amended to ensure the criteria is exactly the same
  across `SI` and `CV` units when the `UNIT_CHECK` is identical or
  missing for a particular `TERM`. This allows grading metadata to be
  combined across `SI` and `CV` units. Also, comment added about
  `"Creatinine Clearance"` and lab values less than 10 for `NCICTCAEv6`.
  `DAIDs` grading metadata updated to use `ADLB.ADT` instead of
  `ADLB.LBDT`
  ([\#2958](https://github.com/pharmaverse/admiral/issues/2958)).

### Documentation

`Lab Grading` vignette updated to add more clarity for
`"Creatinine Increased"` and baseline visits.
([\#2958](https://github.com/pharmaverse/admiral/issues/2958))

- Fix typo in
  [`compute_bsa()`](https:/pharmaverse.github.io/admiral/main/reference/compute_bsa.md)
  in `ADPPK` template and vignette.
  ([\#2963](https://github.com/pharmaverse/admiral/issues/2963))

## admiral 1.4.0

CRAN release: 2026-01-15

### New Features

- New vignette “Explore ADaM Templates” added under the “Get Started”
  section to allow users to peruse the
  [admiral](https://pharmaverse.github.io/admiral/) ADaM templates
  directly from the documentation website.
  ([\#2935](https://github.com/pharmaverse/admiral/issues/2935))

- Updated PK programming vignette to use new experimental functions for
  deriving `NFRLT`
  ([\#2927](https://github.com/pharmaverse/admiral/issues/2927)).

- Added experimental function
  [`derive_var_nfrlt()`](https:/pharmaverse.github.io/admiral/main/reference/derive_var_nfrlt.md)
  to derive `NFRLT` from timepoint and visit (e.g. `VISITDY` and
  `PCTPT`).
  ([\#2929](https://github.com/pharmaverse/admiral/issues/2929))

- Updated PK programming vignette to use new experimental functions for
  deriving `NFRLT`.
  ([\#2927](https://github.com/pharmaverse/admiral/issues/2927))

- Added experimental function
  [`convert_xxtpt_to_hours()`](https:/pharmaverse.github.io/admiral/main/reference/convert_xxtpt_to_hours.md)
  to parse timepoint values (e.g. `PCTPT`, `EGTPT`, `VSTPT`) into hours.
  ([\#2916](https://github.com/pharmaverse/admiral/issues/2916))

- New experimental `ADAB` template script available `ad_adab.R` which
  creates Anti-drug Antibody Analysis Dataset.
  ([\#2805](https://github.com/pharmaverse/admiral/issues/2805))

- Lab grading metadata for NCI-CTCAE version 6.0 is now available for
  both SI and US (Conventional) units via `atoxgr_criteria_ctcv6` and
  `atoxgr_criteria_ctcv6_uscv`. This includes grading criteria for
  Creatinine Clearance decreased.
  ([\#1858](https://github.com/pharmaverse/admiral/issues/1858))

- Lab grading metadata is now maintained in JSON format for all
  supported criteria (NCI-CTCAEv4, NCI-CTCAEv5, NCI-CTCAEv6, and DAIDS)
  in both SI and US conventional unit systems, improving maintainability
  and readability.
  ([\#1858](https://github.com/pharmaverse/admiral/issues/1858))

### Updates of Existing Functions

- The functions
  [`derive_param_bmi()`](https:/pharmaverse.github.io/admiral/main/reference/derive_param_bmi.md),
  [`derive_param_bsa()`](https:/pharmaverse.github.io/admiral/main/reference/derive_param_bsa.md),
  [`derive_param_framingham()`](https:/pharmaverse.github.io/admiral/main/reference/derive_param_framingham.md),
  [`derive_param_map()`](https:/pharmaverse.github.io/admiral/main/reference/derive_param_map.md),
  and
  [`derive_param_qtc()`](https:/pharmaverse.github.io/admiral/main/reference/derive_param_qtc.md)
  were updated such that they no longer fail if
  [admiral](https://pharmaverse.github.io/admiral/) is not loaded.
  ([\#2667](https://github.com/pharmaverse/admiral/issues/2667))

- [`slice_derivation()`](https:/pharmaverse.github.io/admiral/main/reference/slice_derivation.md)
  was updated such that it no longer fails if `args = NULL` is
  specified.
  ([\#2875](https://github.com/pharmaverse/admiral/issues/2875))

- [`derive_summary_records()`](https:/pharmaverse.github.io/admiral/main/reference/derive_summary_records.md)
  now issues an error if a non summary function is used in
  `set_values_to` which results in multiple records per by group.
  ([\#2872](https://github.com/pharmaverse/admiral/issues/2872))

- [`derive_summary_records()`](https:/pharmaverse.github.io/admiral/main/reference/derive_summary_records.md)
  was updated to allow constants to be provided under the
  `constant_values` argument, which will be present in both summary and
  missing rows.
  ([\#2668](https://github.com/pharmaverse/admiral/issues/2668))

### Breaking Changes

- In
  [`derive_var_atoxgr_dir()`](https:/pharmaverse.github.io/admiral/main/reference/derive_var_atoxgr_dir.md),
  the `abnormal_indicator` argument is deprecated and replaced with
  `low_indicator` and `high_indicator` arguments. This change enables
  independent control of baseline range indicators for low and high
  abnormalities, which is required for NCI-CTCAE version 6.0 criteria.
  The deprecated argument will continue to work with a deprecation
  message until the beginning of 2027.
  ([\#1858](https://github.com/pharmaverse/admiral/issues/1858))

- The default value of `ignore_seconds_flag` is set to `TRUE`.
  ([\#2798](https://github.com/pharmaverse/admiral/issues/2798))

- [`derive_var_merged_summary()`](https:/pharmaverse.github.io/admiral/main/reference/derive_var_merged_summary.md)
  is deprecated and will be replaced by
  [`derive_vars_merged_summary()`](https:/pharmaverse.github.io/admiral/main/reference/derive_vars_merged_summary.md).
  This is just a rename of the function to align with our programming
  conventions, i.e. functions that can derive multiple variables have
  `_vars_` in the name rather than `_var_`.
  ([\#2874](https://github.com/pharmaverse/admiral/issues/2874))

- The following function arguments are entering the next phase of the
  [deprecation
  process](https://pharmaverse.github.io/admiraldev/articles/programming_strategy.html#deprecation):

  **Phase 1 (message)**

  - `derive_var_atoxgr_dir(abnormal_indicator = )` is deprecated and
    replaced by the new `low_indicator` and `high_indicator` arguments
    for enhanced flexibility in lab grading.
    ([\#1858](https://github.com/pharmaverse/admiral/issues/1858))
  - [`derive_var_merged_summary()`](https:/pharmaverse.github.io/admiral/main/reference/derive_var_merged_summary.md)
    is deprecated and will be replaced by
    [`derive_vars_merged_summary()`](https:/pharmaverse.github.io/admiral/main/reference/derive_vars_merged_summary.md).
    ([\#2874](https://github.com/pharmaverse/admiral/issues/2874))

  **Phase 2 (warning)**

  - [`call_user_fun()`](https:/pharmaverse.github.io/admiral/main/reference/call_user_fun.md)
    is deprecated and will have no replacement.
    ([\#2678](https://github.com/pharmaverse/admiral/issues/2678))
  - [`derive_param_extreme_record()`](https:/pharmaverse.github.io/admiral/main/reference/derive_param_extreme_record.md)
    is deprecated and replaced by
    [`derive_extreme_event()`](https:/pharmaverse.github.io/admiral/main/reference/derive_extreme_event.md).
  - [`derive_var_dthcaus()`](https:/pharmaverse.github.io/admiral/main/reference/derive_var_dthcaus.md)
    is deprecated and replaced by
    [`derive_vars_extreme_event()`](https:/pharmaverse.github.io/admiral/main/reference/derive_vars_extreme_event.md).
  - [`date_source()`](https:/pharmaverse.github.io/admiral/main/reference/date_source.md)
    is deprecated and replaced by
    [`event()`](https:/pharmaverse.github.io/admiral/main/reference/event.md).
  - [`dthcaus_source()`](https:/pharmaverse.github.io/admiral/main/reference/dthcaus_source.md)
    is deprecated and replaced by
    [`event()`](https:/pharmaverse.github.io/admiral/main/reference/event.md).
  - [`derive_var_extreme_dt()`](https:/pharmaverse.github.io/admiral/main/reference/derive_var_extreme_dt.md)
    and
    [`derive_var_extreme_dtm()`](https:/pharmaverse.github.io/admiral/main/reference/derive_var_extreme_dtm.md)
    are deprecated and replaced by  
    [`derive_vars_extreme_event()`](https:/pharmaverse.github.io/admiral/main/reference/derive_vars_extreme_event.md).
  - [`get_summary_records()`](https:/pharmaverse.github.io/admiral/main/reference/get_summary_records.md)
    is deprecated. Please use
    [`derive_summary_records()`](https:/pharmaverse.github.io/admiral/main/reference/derive_summary_records.md)
    with the `dataset_add` argument and without the `dataset` argument.

  **Phase 3 (error)**

  No functions or arguments in this Phase

  **Phase 4 (removed)**

  No functions or arguments in this Phase

### Documentation

- New experimental vignette “Creating an ADAB ADaM”, which describes how
  to create an anti-drug antibody analysis dataset.
  ([\#2873](https://github.com/pharmaverse/admiral/issues/2873))

- The lab grading vignette has been updated with examples for NCI-CTCAE
  version 6.0 criteria, including usage of the new `low_indicator` and
  `high_indicator` arguments.
  ([\#1858](https://github.com/pharmaverse/admiral/issues/1858))

- The “Ask AI” widget was added to the bottom right of each page. It
  enables users to ask questions about
  [admiral](https://pharmaverse.github.io/admiral/) and
  [admiraldev](https://pharmaverse.github.io/admiraldev/) and receive
  answers from an LLM. It is trained on the documentation of both
  packages and provided by
  [kapa.ai](https://docs.kapa.ai/kapa-for-open-source).
  ([\#2887](https://github.com/pharmaverse/admiral/issues/2887))

- The BDS Findings vignette was updated to move derivation of `ASEQ`
  after any new rows.
  ([\#2780](https://github.com/pharmaverse/admiral/issues/2780))

- `ADLBHY` template was updated to keep `PARAM` in final dataset.
  ([\#2804](https://github.com/pharmaverse/admiral/issues/2804))

- A link to the [{admiral}
  ecosystem](https://pharmaverse.org/e2eclinical/adam/) page was added
  to the README sidebar and main text.
  ([\#2881](https://github.com/pharmaverse/admiral/issues/2881))

- The ADSL template and vignette were updated to add derivation of
  analysis age (`AAGE`/`AAGEU`) using
  [`derive_vars_aage()`](https:/pharmaverse.github.io/admiral/main/reference/derive_vars_aage.md).
  This includes deriving birth date (`BRTHDT`) from birth date character
  variable (`BRTHDTC`) using
  [`derive_vars_dt()`](https:/pharmaverse.github.io/admiral/main/reference/derive_vars_dt.md).
  ([\#2584](https://github.com/pharmaverse/admiral/issues/2584))

- Standardized variable notation across documentation to use `--` for
  SDTM variables (e.g., `--DTC`) and `*` for ADaM variables (e.g.,
  `*DTM`, `*DT`).
  ([\#2757](https://github.com/pharmaverse/admiral/issues/2757))

- For
  [`derive_vars_dt()`](https:/pharmaverse.github.io/admiral/main/reference/derive_vars_dt.md),
  clarify `highest_imputation` definitions and expand date examples.
  ([\#2841](https://github.com/pharmaverse/admiral/issues/2841))

- The documentation was enhanced:
  ([\#2585](https://github.com/pharmaverse/admiral/issues/2585))

  - For
    [`derive_locf_records()`](https:/pharmaverse.github.io/admiral/main/reference/derive_locf_records.md),
    each example now has a title (which is also shown in the TOC) and a
    description, improving readability.
    ([\#2889](https://github.com/pharmaverse/admiral/issues/2889))
  - For
    [`derive_vars_cat()`](https:/pharmaverse.github.io/admiral/main/reference/derive_vars_cat.md),
    each example now has a title (which is also shown in the TOC) and a
    description, improving readability.
    ([\#2701](https://github.com/pharmaverse/admiral/issues/2701))

### Various

Developer Notes

- Updated the [lintr](https://lintr.r-lib.org) preferences to use the
  shared [admiraldev](https://pharmaverse.github.io/admiraldev/)
  configurations.
  ([\#2863](https://github.com/pharmaverse/admiral/issues/2863)) and
  ([\#2913](https://github.com/pharmaverse/admiral/issues/2913))
- To reduce the size of the package:
  ([\#2944](https://github.com/pharmaverse/admiral/issues/2944))
  - The compression method in `data-raw/adlb_grading/atoxgr_sources.R`
    was updated to `"bzip2"` everywhere.
  - The `admiral_adlb` dataset was updated so that it only contains the
    `"AST"`, `"ALT"` and `"BILI"` parameters.

## admiral 1.3.1

CRAN release: 2025-07-29

### Documentation

- The ADSL template and vignette were updated to make example derivation
  of `SAFFL` CDISC-compliant.
  ([\#2782](https://github.com/pharmaverse/admiral/issues/2782))

- In the function documentation references to vignettes were updated to
  meet CRAN requirements.
  ([\#2788](https://github.com/pharmaverse/admiral/issues/2788))

- The minimum [dplyr](https://dplyr.tidyverse.org) version was updated
  to 1.1.1.
  ([\#2788](https://github.com/pharmaverse/admiral/issues/2788))

## admiral 1.3.0

CRAN release: 2025-06-25

### New Features

- The documentation was enhanced:
  ([\#2585](https://github.com/pharmaverse/admiral/issues/2585))

  - The default value of an argument is now displayed in the argument
    description.
  - For some complex functions each example has now a title, which is
    also shown in the TOC, and a description. This enabled adding more
    examples without losing readability. The output of the structured
    examples used for complex functions is displayed in the help pages
    in RStudio. The following existing functions (as well as any new
    functions added in this release) received this enhancement:
    - [`derive_extreme_event()`](https:/pharmaverse.github.io/admiral/main/reference/derive_extreme_event.md)
      ([\#2735](https://github.com/pharmaverse/admiral/issues/2735))
    - [`derive_extreme_records()`](https:/pharmaverse.github.io/admiral/main/reference/derive_extreme_records.md)
      ([\#2585](https://github.com/pharmaverse/admiral/issues/2585))
    - [`derive_param_computed()`](https:/pharmaverse.github.io/admiral/main/reference/derive_param_computed.md)
      ([\#2701](https://github.com/pharmaverse/admiral/issues/2701))
    - [`derive_param_tte()`](https:/pharmaverse.github.io/admiral/main/reference/derive_param_tte.md)
      ([\#2704](https://github.com/pharmaverse/admiral/issues/2704))
    - [`derive_summary_records()`](https:/pharmaverse.github.io/admiral/main/reference/derive_summary_records.md)
      ([\#2707](https://github.com/pharmaverse/admiral/issues/2707))
    - [`derive_var_extreme_flag()`](https:/pharmaverse.github.io/admiral/main/reference/derive_var_extreme_flag.md)
      ([\#2752](https://github.com/pharmaverse/admiral/issues/2752))
    - [`derive_var_joined_exist_flag()`](https:/pharmaverse.github.io/admiral/main/reference/derive_var_joined_exist_flag.md)
      ([\#2729](https://github.com/pharmaverse/admiral/issues/2729))
    - [`derive_var_trtemfl()`](https:/pharmaverse.github.io/admiral/main/reference/derive_var_trtemfl.md)
      ([\#2746](https://github.com/pharmaverse/admiral/issues/2746))
    - [`derive_vars_crit_flag()`](https:/pharmaverse.github.io/admiral/main/reference/derive_vars_crit_flag.md)
      ([\#2744](https://github.com/pharmaverse/admiral/issues/2744))
    - [`derive_vars_dt()`](https:/pharmaverse.github.io/admiral/main/reference/derive_vars_dt.md)
      ([\#2715](https://github.com/pharmaverse/admiral/issues/2715))
    - [`derive_vars_dtm()`](https:/pharmaverse.github.io/admiral/main/reference/derive_vars_dtm.md)
      ([\#2715](https://github.com/pharmaverse/admiral/issues/2715))
    - [`derive_vars_joined()`](https:/pharmaverse.github.io/admiral/main/reference/derive_vars_joined.md)
      ([\#2727](https://github.com/pharmaverse/admiral/issues/2727))
    - [`derive_vars_merged()`](https:/pharmaverse.github.io/admiral/main/reference/derive_vars_merged.md)
      ([\#2727](https://github.com/pharmaverse/admiral/issues/2727))
    - [`filter_joined()`](https:/pharmaverse.github.io/admiral/main/reference/filter_joined.md)
      ([\#2729](https://github.com/pharmaverse/admiral/issues/2729))

- New function
  [`derive_vars_joined_summary()`](https:/pharmaverse.github.io/admiral/main/reference/derive_vars_joined_summary.md)
  to derive summary variables from selected records of an additional
  dataset where the selection depends on variables from both the input
  dataset and the additional dataset. For example, the cumulative dose
  up to each adverse event in `ADAE` can be derived with the new
  function.
  ([\#2652](https://github.com/pharmaverse/admiral/issues/2652))

- New lab grading metadata for US (Conventional) units for the three
  grading criteria `admiral` already produces for SI units
  ([\#2557](https://github.com/pharmaverse/admiral/issues/2557)).

  - `atoxgr_criteria_ctcv4_uscv` (NCI-CTCAEv4 criteria)
  - `atoxgr_criteria_ctcv5_uscv` (NCI-CTCAEv5 criteria)
  - `atoxgr_criteria_daids_uscv` (DAIDs criteria)

### Updates of Existing Functions

- In
  [`derive_vars_joined()`](https:/pharmaverse.github.io/admiral/main/reference/derive_vars_joined.md),
  the `filter_add` argument is now correctly applied for all join types,
  fixing an issue where it was ignored when `join_type != "all"`.
  ([\#2682](https://github.com/pharmaverse/admiral/issues/2682))

- The function
  [`extract_duplicate_records()`](https:/pharmaverse.github.io/admiral/main/reference/extract_duplicate_records.md)
  was updated to consider all variables in the input dataset for the by
  group if the `by_vars` argument is omitted entirely.
  ([\#2644](https://github.com/pharmaverse/admiral/issues/2644))

- In
  [`slice_derivation()`](https:/pharmaverse.github.io/admiral/main/reference/slice_derivation.md),
  previously the derivation is not called for empty subsets, however
  this can lead to issues when the input dataset is empty. Now the
  derivation is called for all subsets.
  ([\#2645](https://github.com/pharmaverse/admiral/issues/2645))

- The examples section for the function
  [`derive_var_trtemfl()`](https:/pharmaverse.github.io/admiral/main/reference/derive_var_trtemfl.md)
  was enhanced to include a showcasing of all scenarios discussed in the
  following [PHUSE White Paper on Treatment-Emergent
  AEs](https://phuse.s3.eu-central-1.amazonaws.com/Deliverables/Safety+Analytics/WP-087+Recommended+Definition+of++Treatment-Emergent+Adverse+Events+in+Clinical+Trials+.pdf).
  ([\#2455](https://github.com/pharmaverse/admiral/issues/2455))

- [`derive_var_atoxgr_dir()`](https:/pharmaverse.github.io/admiral/main/reference/derive_var_atoxgr_dir.md)
  updated to handle more than one unit in grading metadata. Related to
  providing US (Conventional) units for grading
  ([\#2557](https://github.com/pharmaverse/admiral/issues/2557)).

- NCICTCAEv4 and NCICTCAEv5 grading criteria (`atoxgr_criteria_ctcv4`,
  `atoxgr_criteria_ctcv4_uscv`, `atoxgr_criteria_ctcv5`,
  `atoxgr_criteria_ctcv5_uscv`), updated to add terms `"Acidosis"` and
  `"Alkalosis"`
  ([\#2768](https://github.com/pharmaverse/admiral/issues/2768)).

- The functions
  [`derive_summary_records()`](https:/pharmaverse.github.io/admiral/main/reference/derive_summary_records.md)
  and
  [`derive_param_exist_flag()`](https:/pharmaverse.github.io/admiral/main/reference/derive_param_exist_flag.md)
  were updated to fix an issue where if a variable was in both
  `dataset_add` and `dataset_ref`, it was added to the new records even
  if it was not in `by_vars`.
  ([\#2664](https://github.com/pharmaverse/admiral/issues/2664))

- [`derive_var_atoxgr_dir()`](https:/pharmaverse.github.io/admiral/main/reference/derive_var_atoxgr_dir.md)
  updated to handle more than one unit in grading metadata. Related to
  providing US (Conventional) units for grading.
  ([\#2557](https://github.com/pharmaverse/admiral/issues/2557))

- The background checks in
  [`derive_summary_records()`](https:/pharmaverse.github.io/admiral/main/reference/derive_summary_records.md)
  were too restrictive: `by_vars` were expected in `dataset` although
  the code did not require it. This requirement has therefore been
  dropped.
  ([\#2686](https://github.com/pharmaverse/admiral/issues/2686))

- The functions
  [`derive_vars_joined()`](https:/pharmaverse.github.io/admiral/main/reference/derive_vars_joined.md),
  [`derive_var_joined_exist_flag()`](https:/pharmaverse.github.io/admiral/main/reference/derive_var_joined_exist_flag.md),
  and
  [`filter_joined()`](https:/pharmaverse.github.io/admiral/main/reference/filter_joined.md)
  produce correct results now when they are used with
  `join_type = "before"` or `join_type = "after"` and `dataset` and
  `dataset_add` differ or the `filter_add` argument is used.
  ([\#2863](https://github.com/pharmaverse/admiral/issues/2863))

- The function
  [`derive_locf_records()`](https:/pharmaverse.github.io/admiral/main/reference/derive_locf_records.md)
  was updated to include two new arguments: `id_vars_ref` and
  `imputation`. The `id_vars_ref` argument allows users to select the
  variables to group by in the reference dataset (`dataset_ref`) when
  determining which observations to add to the input dataset. The
  `imputation` argument lets users decide whether to update
  `analysis_var` when its value is `NA` (“update” and “update_add”), or
  to add a new observation instead (“add”).
  ([\#2694](https://github.com/pharmaverse/admiral/issues/2694))
  ([\#2680](https://github.com/pharmaverse/admiral/issues/2680))
  ([\#2717](https://github.com/pharmaverse/admiral/issues/2717))

- [`derive_vars_dt()`](https:/pharmaverse.github.io/admiral/main/reference/derive_vars_dt.md),
  [`derive_vars_dtm()`](https:/pharmaverse.github.io/admiral/main/reference/derive_vars_dtm.md),
  [`impute_dtc_dt()`](https:/pharmaverse.github.io/admiral/main/reference/impute_dtc_dt.md),
  [`impute_dtc_dtm()`](https:/pharmaverse.github.io/admiral/main/reference/impute_dtc_dtm.md),
  [`convert_dtc_to_dt()`](https:/pharmaverse.github.io/admiral/main/reference/convert_dtc_to_dt.md),
  &
  [`convert_dtc_to_dtm()`](https:/pharmaverse.github.io/admiral/main/reference/convert_dtc_to_dtm.md)
  and related functions will now throw an error instead of a warning
  when `highest_imputation = "Y"` but neither `min_date` (when
  `date_imputation = "first"`) nor `max_dates` (when
  `date_imputation = "last"`) are specified.
  ([\#2654](https://github.com/pharmaverse/admiral/issues/2654))

- [`create_query_data()`](https:/pharmaverse.github.io/admiral/main/reference/create_query_data.md)
  no longer issues warnings when
  [`basket_select()`](https:/pharmaverse.github.io/admiral/main/reference/basket_select.md)
  objects with custom arguments of length greater than one are used.
  ([\#2751](https://github.com/pharmaverse/admiral/issues/2751))

- The `order` argument in
  [`derive_var_joined_exist_flag()`](https:/pharmaverse.github.io/admiral/main/reference/derive_var_joined_exist_flag.md)
  and
  [`filter_joined()`](https:/pharmaverse.github.io/admiral/main/reference/filter_joined.md)
  is now optional unless `join_type = "after"`, `join_type = "before"`,
  `first_cond_lower`, `first_cond_upper`, or `tmp_obs_nr_var` are
  specified.
  ([\#2729](https://github.com/pharmaverse/admiral/issues/2729))

### Breaking Changes

- [`derive_vars_dtm()`](https:/pharmaverse.github.io/admiral/main/reference/derive_vars_dtm.md)
  issues a message alerting users to a coming change in `admiral 1.4.0`
  where the default behavior of `ignore_seconds_flag` will be changed
  from `FALSE` to `TRUE`.
  ([\#2661](https://github.com/pharmaverse/admiral/issues/2661))

- Lab grading metadata
  [`atoxgr_criteria_ctcv4()`](https:/pharmaverse.github.io/admiral/main/reference/atoxgr_criteria_ctcv4.md),
  [`atoxgr_criteria_ctcv5()`](https:/pharmaverse.github.io/admiral/main/reference/atoxgr_criteria_ctcv5.md)
  and
  [`atoxgr_criteria_daids()`](https:/pharmaverse.github.io/admiral/main/reference/atoxgr_criteria_daids.md)
  variable `SI_UNIT_CHECK` renamed to `UNIT_CHECK`.
  ([\#2557](https://github.com/pharmaverse/admiral/issues/2557))

- The values of the variable specified for `tmp_obs_nr_var` in
  [`derive_vars_joined()`](https:/pharmaverse.github.io/admiral/main/reference/derive_vars_joined.md),
  [`derive_var_joined_exist_flag()`](https:/pharmaverse.github.io/admiral/main/reference/derive_var_joined_exist_flag.md),
  [`filter_joined()`](https:/pharmaverse.github.io/admiral/main/reference/filter_joined.md)
  are now populated differently if there are multiple records in
  `dataset` or `dataset_add` for the same values of `by_vars` and
  `order`. Before each of these records was assigned a different value,
  i.e., the variable (together with `by_vars`) was a unique identifier.
  Now the value is the same for all these records.
  ([\#2683](https://github.com/pharmaverse/admiral/issues/2683))

- The following function arguments are entering the next phase of the
  [deprecation
  process](https://pharmaverse.github.io/admiraldev/articles/programming_strategy.html#deprecation):
  ([\#2487](https://github.com/pharmaverse/admiral/issues/2487))
  ([\#2595](https://github.com/pharmaverse/admiral/issues/2595))

  **Phase 1 (message)**

  - [`call_user_fun()`](https:/pharmaverse.github.io/admiral/main/reference/call_user_fun.md)
    is deprecated and will have no replacement.
    ([\#2678](https://github.com/pharmaverse/admiral/issues/2678))
  - [`derive_param_extreme_record()`](https:/pharmaverse.github.io/admiral/main/reference/derive_param_extreme_record.md)
    is deprecated and replaced by
    [`derive_extreme_event()`](https:/pharmaverse.github.io/admiral/main/reference/derive_extreme_event.md)
  - [`derive_var_dthcaus()`](https:/pharmaverse.github.io/admiral/main/reference/derive_var_dthcaus.md)
    is deprecated and replaced by
    [`derive_vars_extreme_event()`](https:/pharmaverse.github.io/admiral/main/reference/derive_vars_extreme_event.md)
  - [`date_source()`](https:/pharmaverse.github.io/admiral/main/reference/date_source.md)
    is deprecated and replaced by
    [`event()`](https:/pharmaverse.github.io/admiral/main/reference/event.md)
  - [`dthcaus_source()`](https:/pharmaverse.github.io/admiral/main/reference/dthcaus_source.md)
    is deprecated and replaced by
    [`event()`](https:/pharmaverse.github.io/admiral/main/reference/event.md)
  - [`derive_var_extreme_dt()`](https:/pharmaverse.github.io/admiral/main/reference/derive_var_extreme_dt.md)
    and
    [`derive_var_extreme_dtm()`](https:/pharmaverse.github.io/admiral/main/reference/derive_var_extreme_dtm.md)
    are deprecated and replaced by  
    [`derive_vars_extreme_event()`](https:/pharmaverse.github.io/admiral/main/reference/derive_vars_extreme_event.md)
  - [`get_summary_records()`](https:/pharmaverse.github.io/admiral/main/reference/get_summary_records.md)
    is deprecated. Please use
    [`derive_summary_records()`](https:/pharmaverse.github.io/admiral/main/reference/derive_summary_records.md)
    with the `dataset_add` argument and without the `dataset` argument.

  **Phase 2 (warning)**

  No functions or arguments in this Phase

  **Phase 3 (error)**

  No functions or arguments in this Phase

  **Phase 4 (removed)**

  No functions or arguments in this Phase

### Documentation

- Improved documentation, error messages, and argument assertions of
  [`derive_vars_dt()`](https:/pharmaverse.github.io/admiral/main/reference/derive_vars_dt.md),
  [`derive_vars_dtm()`](https:/pharmaverse.github.io/admiral/main/reference/derive_vars_dtm.md),
  [`impute_dtc_dt()`](https:/pharmaverse.github.io/admiral/main/reference/impute_dtc_dt.md),
  [`impute_dtc_dtm()`](https:/pharmaverse.github.io/admiral/main/reference/impute_dtc_dtm.md),
  [`convert_dtc_to_dt()`](https:/pharmaverse.github.io/admiral/main/reference/convert_dtc_to_dt.md),
  &
  [`convert_dtc_to_dtm()`](https:/pharmaverse.github.io/admiral/main/reference/convert_dtc_to_dtm.md).
  ([\#2654](https://github.com/pharmaverse/admiral/issues/2654))

- Added an example to the
  [`derive_vars_transposed()`](https:/pharmaverse.github.io/admiral/main/reference/derive_vars_transposed.md)
  reference page to showcase how duplicates-related errors can arise
  when records in `dataset_merge` are not uniquely identified.
  ([\#2609](https://github.com/pharmaverse/admiral/issues/2609))

- Default value of `type` in
  [`derive_vars_aage()`](https:/pharmaverse.github.io/admiral/main/reference/derive_vars_aage.md)
  is now shown as `interval` to match the function behavior.
  ([\#2685](https://github.com/pharmaverse/admiral/issues/2685))

- The “Lab Grading” vignette was updated to correct some typos and make
  text easier to read
  ([\#2623](https://github.com/pharmaverse/admiral/issues/2623)) and
  updated to include new metadata for grading using US (Conventional)
  units. ([\#2557](https://github.com/pharmaverse/admiral/issues/2557))

- The “BDS Time-to-Event” vignette was updated to include `SRCSEQ`
  consistently.
  ([\#2658](https://github.com/pharmaverse/admiral/issues/2658))

- The template for ADAE and the OCCDS vignette were updated to include
  an example of the `DOSEON` and `DOSEU` variables
  ([\#2737](https://github.com/pharmaverse/admiral/issues/2737)).

- The template for ADPC and vignette were updated to include an example
  of using `DTYPE` for imputed records
  ([\#2657](https://github.com/pharmaverse/admiral/issues/2657)).

- The “Higher Order Functions” vignette was updated to showcase an
  example of two higher order functions used in combination. The
  documentation for each of the higher order functions was also
  corrected by removing the stated requirement that the `derivation`
  takes a `dataset` argument
  ([\#2656](https://github.com/pharmaverse/admiral/issues/2656)).

- The ‘Details’ section of the
  [`derive_var_analysis_ratio()`](https:/pharmaverse.github.io/admiral/main/reference/derive_var_analysis_ratio.md)
  function and the ‘Derive Analysis Ratio’ section of the “Creating a
  BDS Finding ADaM” vignette were updated to include references to
  `R2AyHI` and `R2AyLO`.
  ([\#2548](https://github.com/pharmaverse/admiral/issues/2548))

- The
  [`derive_basetype_records()`](https:/pharmaverse.github.io/admiral/main/reference/derive_basetype_records.md)
  documentation was updated to clarify `BASETYPE` derivations.
  ([\#2545](https://github.com/pharmaverse/admiral/issues/2545))

- The
  [`derive_extreme_records()`](https:/pharmaverse.github.io/admiral/main/reference/derive_extreme_records.md)
  documentation was updated to clarify which variables are populated
  from `dataset_ref` for the new observations.
  ([\#2664](https://github.com/pharmaverse/admiral/issues/2664))

- The ‘Assign `PARAMCD`, `PARAM`, `PARAMN`, `PARCAT1`’ section of the
  “Creating a BDS Finding ADaM” vignette was updated to clarify `PARAM`
  to `PARCAT1` mapping.
  ([\#2547](https://github.com/pharmaverse/admiral/issues/2547))

- The package documentation for (1) all CRAN-released versions starting
  from [admiral](https://pharmaverse.github.io/admiral/) 1.0.0 up until
  the current version and (2) the latest development version (listed
  under “main”) are all now accessible using the “Versions” selector in
  the toolbar.
  ([\#2766](https://github.com/pharmaverse/admiral/issues/2766))

### Various

Developer Notes

- Removed CODEOWNERS file from repo
  ([\#2674](https://github.com/pharmaverse/admiral/issues/2674))

## admiral 1.2.0

CRAN release: 2025-01-15

### New Features

- New function
  [`derive_vars_cat()`](https:/pharmaverse.github.io/admiral/main/reference/derive_vars_cat.md)
  for deriving pairs of variables or more, e.g.  `AVALCATy` &
  `AVALCAyN`.
  ([\#2480](https://github.com/pharmaverse/admiral/issues/2480))
- New function
  [`derive_vars_crit_flag()`](https:/pharmaverse.github.io/admiral/main/reference/derive_vars_crit_flag.md)
  for deriving criterion flag variables (`CRITy`, `CRITyFL`, `CRITyFN`).
  ([\#2468](https://github.com/pharmaverse/admiral/issues/2468))
- New function
  [`transform_range()`](https:/pharmaverse.github.io/admiral/main/reference/transform_range.md)
  to transform values from a source range to a target range.
  ([\#2571](https://github.com/pharmaverse/admiral/issues/2571))

### Updates of Existing Functions

- Added `"message"` as option for `check_type` argument in
  [`derive_var_obs_number()`](https:/pharmaverse.github.io/admiral/main/reference/derive_var_obs_number.md)
  function.
  ([\#2481](https://github.com/pharmaverse/admiral/issues/2481))

- Added `"message"` as option for `check_type` argument in
  [`filter_extreme()`](https:/pharmaverse.github.io/admiral/main/reference/filter_extreme.md)
  function.
  ([\#2481](https://github.com/pharmaverse/admiral/issues/2481))

- Users can now specify how duplicate records are handled in
  [`derive_param_tte()`](https:/pharmaverse.github.io/admiral/main/reference/derive_param_tte.md)
  using the `check_type` argument, with options including `"error"`,
  `"warning"`, `"message"`, or `"none"`, allowing for greater
  flexibility in managing duplicate data scenarios.
  ([\#2481](https://github.com/pharmaverse/admiral/issues/2481))

- The `order` argument has been added to
  [`event_source()`](https:/pharmaverse.github.io/admiral/main/reference/event_source.md)
  and
  [`censor_source()`](https:/pharmaverse.github.io/admiral/main/reference/censor_source.md)
  and  
  defaulted to `NULL` to allow specifying variables in addition to the
  date variable. This can be used to ensure the uniqueness of the select
  records if there is more than one record per date.
  ([\#2481](https://github.com/pharmaverse/admiral/issues/2481))

- NCICTCAEv5 grading criteria (`atoxgr_criteria_ctcv5`):

  - fix for `TERM = "INR increased"`, criteria was wrongly using
    `x ULN`, for first part of criteria for grades 1 to 3. For example,
    `">2.5 x ULN"` changed to `">2.5"` for grade 3.
    ([\#2534](https://github.com/pharmaverse/admiral/issues/2534)).
  - when looking at abnormal baseline we now use `BNRIND` instead of
    comparing `BASE` with `ANRHI`, as `ANRHI` may differ within a
    subject and lab test due to data from different lab vendors. This
    effects 5 terms, namely, `Alanine aminotransferase increased`,
    `Alkaline phosphatase increased`,
    `Aspartate aminotransferase increased`, `Blood bilirubin increased`
    and `GGT Increased`.
    ([\#249](https://github.com/pharmaverse/admiral/issues/249))
  - [`derive_var_atoxgr_dir()`](https:/pharmaverse.github.io/admiral/main/reference/derive_var_atoxgr_dir.md):
    new argument `abnormal_indicator` to pass in value of `BNRIND` to
    indicate lab test is abnormal. This is only used for the 5 lab tests
    described above.
    ([\#249](https://github.com/pharmaverse/admiral/issues/249))

- The `keep_nas` argument of
  [`derive_param_computed()`](https:/pharmaverse.github.io/admiral/main/reference/derive_param_computed.md)
  was enhanced such that it is now possible to specify a list of
  variables for which `NA`s are acceptable. I.e., records are added even
  if some of the specified variables are `NA`.
  ([\#2510](https://github.com/pharmaverse/admiral/issues/2510))

- [`derive_param_tte()`](https:/pharmaverse.github.io/admiral/main/reference/derive_param_tte.md)
  now provides a useful error message if in `event_conditions` or
  `censor_conditions` a dataset is referenced which is not specified in
  `source_datasets`.
  ([\#2519](https://github.com/pharmaverse/admiral/issues/2519))

- The
  [`derive_param_qtc()`](https:/pharmaverse.github.io/admiral/main/reference/derive_param_qtc.md)
  function accepts now both `"ms"` and `"msec"` as unit of the input
  parameters.
  ([\#2513](https://github.com/pharmaverse/admiral/issues/2513))

- In
  [`derive_vars_query()`](https:/pharmaverse.github.io/admiral/main/reference/derive_vars_query.md)
  the error message was improved for the cases that some of the
  requested query variables are already present in the input dataset or
  that the queries dataset contains duplicates.
  ([\#2543](https://github.com/pharmaverse/admiral/issues/2543))

- [`derive_vars_atc()`](https:/pharmaverse.github.io/admiral/main/reference/derive_vars_atc.md)
  and
  [`create_single_dose_dataset()`](https:/pharmaverse.github.io/admiral/main/reference/create_single_dose_dataset.md)
  `by_vars` argument updated to use `get_admiral_option("subject_keys")`
  instead of `USUBJID` or `STUDYID` in `bds_exposure.Rmd`.
  ([\#2501](https://github.com/pharmaverse/admiral/issues/2501))

- The test scripts, R, and markdown files for
  [`create_single_dose_dataset()`](https:/pharmaverse.github.io/admiral/main/reference/create_single_dose_dataset.md)
  and `occds.Rmd` updated to include a `STUDYID` column because of
  `get_admiral_option("subject_keys")` update above.
  ([\#2501](https://github.com/pharmaverse/admiral/issues/2501))

- Update
  [`derive_vars_period()`](https:/pharmaverse.github.io/admiral/main/reference/derive_vars_period.md)
  to make it work when there is only one new variable.
  ([\#2582](https://github.com/pharmaverse/admiral/issues/2582))

- A check was added to
  [`derive_vars_transposed()`](https:/pharmaverse.github.io/admiral/main/reference/derive_vars_transposed.md)
  and
  [`derive_vars_atc()`](https:/pharmaverse.github.io/admiral/main/reference/derive_vars_atc.md)
  which stops execution if the records in `dataset_merge` or
  `dataset_facm` respectively are not unique.
  ([\#2563](https://github.com/pharmaverse/admiral/issues/2563))

- The functions
  [`derive_vars_joined()`](https:/pharmaverse.github.io/admiral/main/reference/derive_vars_joined.md),
  [`derive_var_joined_exist_flag()`](https:/pharmaverse.github.io/admiral/main/reference/derive_var_joined_exist_flag.md),
  [`derive_extreme_event()`](https:/pharmaverse.github.io/admiral/main/reference/derive_extreme_event.md),
  and
  [`filter_joined()`](https:/pharmaverse.github.io/admiral/main/reference/filter_joined.md)
  were updated to reduce their memory consumption. As the new code
  increases the run-time, it is not used by default. To enable it the
  new admiral option `save_memory` has to be set to `TRUE`.
  ([\#2590](https://github.com/pharmaverse/admiral/issues/2590))

- The function
  [`compute_egfr()`](https:/pharmaverse.github.io/admiral/main/reference/compute_egfr.md)
  updated to allow missing values for sex which result in missing values
  for output.
  ([\#2612](https://github.com/pharmaverse/admiral/issues/2612))

### Breaking Changes

- The following function arguments are entering the next phase of the
  [deprecation
  process](https://pharmaverse.github.io/admiraldev/articles/programming_strategy.html#deprecation):
  ([\#2487](https://github.com/pharmaverse/admiral/issues/2487))
  ([\#2595](https://github.com/pharmaverse/admiral/issues/2595))

  **Phase 1 (message)**

  - [`derive_param_extreme_record()`](https:/pharmaverse.github.io/admiral/main/reference/derive_param_extreme_record.md)
    is deprecated and replaced by
    [`derive_extreme_event()`](https:/pharmaverse.github.io/admiral/main/reference/derive_extreme_event.md)
  - [`derive_var_dthcaus()`](https:/pharmaverse.github.io/admiral/main/reference/derive_var_dthcaus.md)
    is deprecated and replaced by
    [`derive_vars_extreme_event()`](https:/pharmaverse.github.io/admiral/main/reference/derive_vars_extreme_event.md)
  - [`date_source()`](https:/pharmaverse.github.io/admiral/main/reference/date_source.md)
    is deprecated and replaced by
    [`event()`](https:/pharmaverse.github.io/admiral/main/reference/event.md)
  - [`dthcaus_source()`](https:/pharmaverse.github.io/admiral/main/reference/dthcaus_source.md)
    is deprecated and replaced by
    [`event()`](https:/pharmaverse.github.io/admiral/main/reference/event.md)
  - [`derive_var_extreme_dt()`](https:/pharmaverse.github.io/admiral/main/reference/derive_var_extreme_dt.md)
    and
    [`derive_var_extreme_dtm()`](https:/pharmaverse.github.io/admiral/main/reference/derive_var_extreme_dtm.md)
    are deprecated and replaced by
    [`derive_vars_extreme_event()`](https:/pharmaverse.github.io/admiral/main/reference/derive_vars_extreme_event.md)
  - [`get_summary_records()`](https:/pharmaverse.github.io/admiral/main/reference/get_summary_records.md)
    is deprecated. Please use
    [`derive_summary_records()`](https:/pharmaverse.github.io/admiral/main/reference/derive_summary_records.md)
    with the `dataset_add` argument and without the `dataset` argument.

  **Phase 2 (warning)**

  No functions or arguments in this Phase

  **Phase 3 (error)**

  No functions or arguments in this Phase

  **Phase 4 (removed)**

  - `consolidate_metadata(check_keys)`
  - Removed at v1.1.1 `compute_egfr(wt)`
  - Removed at v1.1.1 `derive_expected_records(dataset_expected_obs)`
  - Removed at v1.1.1 `derive_locf_records(dataset_expected_obs)`
  - `derive_extreme_event(ignore_event_order)`
  - `derive_vars_merged(match_flag)`
  - `derive_var_merged_summary(new_var, analysis_var, summary_fun)`
  - Removed at v1.1.1
    `derive_param_computed(analysis_value, analysis_var)`
  - `derive_param_exposure(filter, analysis_var, summary_fun)`
  - `derive_summary_records(filter)`
  - Removed at v1.1.1 `derive_extreme_records(filter)`
  - `derive_var_joined_exist_flag(first_cond, filter)`
  - `event_joined(first_cond)`
  - `filter_joined(first_cond, filter)`
  - In
    [`get_summary_records()`](https:/pharmaverse.github.io/admiral/main/reference/get_summary_records.md),
    previously deprecated formal arguments `analysis_var` and  
    `summary_fun` now removed from function, documentation, tests etc.
    ([\#2521](https://github.com/pharmaverse/admiral/issues/2521))

### Documentation

- [`derive_locf_records()`](https:/pharmaverse.github.io/admiral/main/reference/derive_locf_records.md)
  documentation example was fixed to display LOCF records.
  ([\#2461](https://github.com/pharmaverse/admiral/issues/2461))

- The “Find my function” and “Presentation Archive”” links were made
  more prominent in the website navigation bar.
  ([\#2536](https://github.com/pharmaverse/admiral/issues/2536))

- [`derive_var_joined_exist_flag()`](https:/pharmaverse.github.io/admiral/main/reference/derive_var_joined_exist_flag.md)
  documentation updated with extra examples.
  ([\#2523](https://github.com/pharmaverse/admiral/issues/2523))

- Updated the Cheat Sheet to be in line with the 1.2 release of
  [admiral](https://pharmaverse.github.io/admiral/).
  ([\#2458](https://github.com/pharmaverse/admiral/issues/2458))

- In the
  [`derive_param_tte()`](https:/pharmaverse.github.io/admiral/main/reference/derive_param_tte.md)
  documentation is was clarified which event/censoring is selected if
  there is more than one at the same date (for events the first one
  specified in `event_conditions` and for censoring the last one in
  `censor_conditions`).
  ([\#2639](https://github.com/pharmaverse/admiral/issues/2639))

### Various

- Replace use of `data("sdtm")` with `sdtm <- pharmaverse::sdtm` in
  templates and vignettes.
  ([\#2498](https://github.com/pharmaverse/admiral/issues/2498))
- Remove
  [`dthcaus_source()`](https:/pharmaverse.github.io/admiral/main/reference/dthcaus_source.md)
  calls in `ADSL` template because they are deprecated.
  ([\#2517](https://github.com/pharmaverse/admiral/issues/2517))
- Update `ADEG` template to flag `ABLFL` and `ANL01FL` based on
  `DTYPE == "AVERAGE"` records.
  ([\#2561](https://github.com/pharmaverse/admiral/issues/2561))

Developer Notes

- Created unit tests for developer internal function
  `restricted_imputed_dtc_dt()`
  ([\#2495](https://github.com/pharmaverse/admiral/issues/2495))
- Adopted `data-raw/data` R Package Convention
  ([\#2427](https://github.com/pharmaverse/admiral/issues/2427),
  [\#2584](https://github.com/pharmaverse/admiral/issues/2584))
- [`compute_bsa()`](https:/pharmaverse.github.io/admiral/main/reference/compute_bsa.md)
  now uses the more common (but equivalent) version of the DuBois-DuBois
  formula for BSA. The results have not changed.
  ([\#2532](https://github.com/pharmaverse/admiral/issues/2532))  
- Removed `.devcontainer` file (codespace)
  ([\#2524](https://github.com/pharmaverse/admiral/issues/2524))
- Restructured `derive_adeg_parms.R` and `derive_advs_parms.R` and
  related test files for easier reference
  ([\#2551](https://github.com/pharmaverse/admiral/issues/2551))

## admiral 1.1.1

CRAN release: 2024-06-17

- [`derive_extreme_event()`](https:/pharmaverse.github.io/admiral/main/reference/derive_extreme_event.md)
  was fixed such that `check_type = "none"` is accepted again.
  ([\#2462](https://github.com/pharmaverse/admiral/issues/2462))

## admiral 1.1.0

CRAN release: 2024-06-07

### New Features

- Added helper functions to
  [`derive_vars_merged()`](https:/pharmaverse.github.io/admiral/main/reference/derive_vars_merged.md)
  ([`get_flagged_records()`](https:/pharmaverse.github.io/admiral/main/reference/get_flagged_records.md))
  and
  [`derive_vars_query()`](https:/pharmaverse.github.io/admiral/main/reference/derive_vars_query.md)
  ([`get_vars_query()`](https:/pharmaverse.github.io/admiral/main/reference/get_vars_query.md))
  so that those can be called independently as per user’s request.
  [`derive_vars_merged()`](https:/pharmaverse.github.io/admiral/main/reference/derive_vars_merged.md)
  function call results are not impacted by the change
  ([\#2441](https://github.com/pharmaverse/admiral/issues/2441)).

- Error Messaging has been made more “user-friendly”.
  ([\#2372](https://github.com/pharmaverse/admiral/issues/2372))

- New
  [`country_code_lookup()`](https:/pharmaverse.github.io/admiral/main/reference/country_code_lookup.md)
  metadata added to decode countries based on [ISO 3166
  codes](https://www.iso.org/iso-3166-country-codes.html).
  ([\#2388](https://github.com/pharmaverse/admiral/issues/2388))

### Updates of Existing Functions

- `group_var` (optional) parameter is added to
  [`derive_var_trtemfl()`](https:/pharmaverse.github.io/admiral/main/reference/derive_var_trtemfl.md)
  to derive `TRTEMFL` for AE data if the data are collected as one
  episode of AE with multiple lines.
  ([\#2302](https://github.com/pharmaverse/admiral/issues/2302))

- Templates for ADPC, ADPPK and ADPP are updated to handle urine
  records.
  ([\#2392](https://github.com/pharmaverse/admiral/issues/2392))

- [`create_single_dose_dataset()`](https:/pharmaverse.github.io/admiral/main/reference/create_single_dose_dataset.md)
  has been updated to error if the `lookup_table` contains duplicates.
  ([\#2247](https://github.com/pharmaverse/admiral/issues/2247))

- [`derive_vars_merged()`](https:/pharmaverse.github.io/admiral/main/reference/derive_vars_merged.md)
  and
  [`derive_vars_transposed()`](https:/pharmaverse.github.io/admiral/main/reference/derive_vars_transposed.md)
  have a `relationship` argument added (the same as found in
  `dplyr::*_join()` functions) for users to specify what type of join
  (one-to-one, one-to-many, etc.) should take place.
  ([\#2247](https://github.com/pharmaverse/admiral/issues/2247))

- [`basket_select()`](https:/pharmaverse.github.io/admiral/main/reference/basket_select.md)
  function updated to add `...` argument to allow other qualifiers to be
  passed to user-defined function specified in `get_terms_fun()`
  argument for function
  [`create_query_data()`](https:/pharmaverse.github.io/admiral/main/reference/create_query_data.md).
  ([\#2265](https://github.com/pharmaverse/admiral/issues/2265))

- Messaging updated for
  [`derive_extreme_event()`](https:/pharmaverse.github.io/admiral/main/reference/derive_extreme_event.md)
  to improve clarity around duplicates.
  [\#2405](https://github.com/pharmaverse/admiral/issues/2405)

- The `id_vars` argument was added to
  [`derive_vars_transposed()`](https:/pharmaverse.github.io/admiral/main/reference/derive_vars_transposed.md)
  and
  [`derive_vars_atc()`](https:/pharmaverse.github.io/admiral/main/reference/derive_vars_atc.md)
  to allow additional variables, beyond those in `by_vars`, to uniquely
  identify records in the `dataset_merge` argument.
  ([\#2325](https://github.com/pharmaverse/admiral/issues/2325))

- Update PK Programming vignette and templates for ADPC and ADPPK for
  the nominal time formula `NFRLT` to reduce duplicate records in dose
  expansion with
  [`create_single_dose_dataset()`](https:/pharmaverse.github.io/admiral/main/reference/create_single_dose_dataset.md).
  ([\#2426](https://github.com/pharmaverse/admiral/issues/2426))

- Template for ADSL updated so that `EOSSTT` is assigned as `"ONGOING"`
  when no study completion rows exist yet in DS.
  ([\#2436](https://github.com/pharmaverse/admiral/issues/2436))

- The
  [`slice_derivation()`](https:/pharmaverse.github.io/admiral/main/reference/slice_derivation.md)
  function was updated such that it works now when called in a function
  where objects from the function environment are used.
  ([\#2244](https://github.com/pharmaverse/admiral/issues/2244))

### Breaking Changes

- The following function arguments are entering the next phase of the
  deprecation process:
  ([\#2299](https://github.com/pharmaverse/admiral/issues/2299))

  - `compute_egfr(wt)`
  - `consolidate_metadata(check_keys)`
  - `derive_expected_records(dataset_expected_obs)`
  - `derive_locf_records(dataset_expected_obs)`
  - `derive_extreme_event(ignore_event_order)`
  - `derive_vars_merged(match_flag, new_var, analysis_var, summary_fun)`
  - `derive_param_computed(analysis_value, analysis_var)`
  - `derive_param_exposure(filter, analysis_var, summary_fun)`
  - `derive_summary_records(filter)`
  - `derive_extreme_records(filter)`
  - `derive_var_joined_exist_flag(first_cond, filter)`
  - `event_joined(first_cond)`
  - `filter_joined(first_cond, filter)`

- The following function arguments have reached the end of the
  deprecation process and been removed:
  ([\#2299](https://github.com/pharmaverse/admiral/issues/2299))

  - `dthcaus_source(traceability_vars)`
  - `date_source(traceability_vars)`
  - `derive_var_ontrtfl(span_period)`
  - `derive_var_shift(na_val)`
  - `derive_vars_aage(unit)`

### Documentation

- Documentation for
  [`derive_extreme_event()`](https:/pharmaverse.github.io/admiral/main/reference/derive_extreme_event.md)
  has been updated to include a description for the value of
  [`derive_extreme_event()`](https:/pharmaverse.github.io/admiral/main/reference/derive_extreme_event.md)
  when `keep_source_vars = NULL`.
  ([\#2398](https://github.com/pharmaverse/admiral/issues/2398))

- The “Visit and Period Variables” vignette was updated and refactored
  to include example code to create a period reference dataset.
  ([\#2321](https://github.com/pharmaverse/admiral/issues/2321))

- The documentation of
  [`derive_vars_merged()`](https:/pharmaverse.github.io/admiral/main/reference/derive_vars_merged.md)
  function is updated to describe that the `check_type` argument is
  ignored (an error is issued) if `order` is not specified.
  ([\#2326](https://github.com/pharmaverse/admiral/issues/2326))

- The “User Guides” section has been reorganized. A new “Programming
  Concepts and Conventions” vignette was also added to provide more
  context and information around common
  [admiral](https://pharmaverse.github.io/admiral/) behaviors and ways
  of working.
  ([\#2395](https://github.com/pharmaverse/admiral/issues/2395))

- The “Get Started” section has been revamped, placing greater focus on
  material that may help users familiarize themselves with
  [admiral](https://pharmaverse.github.io/admiral/). There are now new
  sections showcasing the various types of
  [admiral](https://pharmaverse.github.io/admiral/) functions and some
  of the more advanced topics have been moved to the new “Programming
  Concepts and Conventions” vignette.
  ([\#2395](https://github.com/pharmaverse/admiral/issues/2395))

- The Examples section of
  [`derive_param_computed()`](https:/pharmaverse.github.io/admiral/main/reference/derive_param_computed.md)
  now contains a new item showcasing how to create a derived parameter
  in the case that a variable contributing to the derived parameter has
  some/all of its values missing.
  ([\#2338](https://github.com/pharmaverse/admiral/issues/2338))

### Various

- Templates and vignettes do not add or populate `AVALC` for
  BDS-findings datasets where the information contained in `AVALC` would
  be redundant with `AVAL`.
  ([\#2442](https://github.com/pharmaverse/admiral/issues/2442))

- The function
  [`dplyr::transmute()`](https://dplyr.tidyverse.org/reference/transmute.html)
  is superseded in favor of `dplyr::mutate(.keep = "none")`.
  Consequently, all the admiral functions that utilized the former have
  been updated accordingly.
  ([\#2274](https://github.com/pharmaverse/admiral/issues/2274))

- The templates for ADPP and ADPC are updated for missing variables
  ([\#2308](https://github.com/pharmaverse/admiral/issues/2308)) and to
  make `ATPT` and `ATPTN` consistent.
  ([\#2328](https://github.com/pharmaverse/admiral/issues/2328))

- ADLB template updated to make `PARAM` consistent for `PARAMCD` values
  `"BASO"` and `"LYMPH"`.
  ([\#2327](https://github.com/pharmaverse/admiral/issues/2327))

Developer Notes

- In the previous version, `renv` was the default framework used to
  manage package dependencies. Now, we use `devtools` as our main
  package manager (some changes also occurred for [admiralci
  workflows](https://github.com/pharmaverse/admiralci)). There is a
  possibility to get package dependency versions used for the workflows
  to ensure local reproducibility. For this, you need to go under the
  latest action summary in your current PR. You can see a deps artifact.
  For each version of R used for `R CMD CHECKS` jobs, there is an
  associated renv.lock file (under the deps artifact).
- Splitting out `R` and `test` files for date/time functions for
  cyclomatic complexity refactor
  ([\#2340](https://github.com/pharmaverse/admiral/issues/2340))([\#2339](https://github.com/pharmaverse/admiral/issues/2339))
- Created three unit tests for
  [`get_summary_records()`](https:/pharmaverse.github.io/admiral/main/reference/get_summary_records.md).
  ([\#2304](https://github.com/pharmaverse/admiral/issues/2304))
- Created unit tests for developer internal function
  [`get_imputation_target_date()`](https:/pharmaverse.github.io/admiral/main/reference/get_imputation_target_date.md)
  ([\#2378](https://github.com/pharmaverse/admiral/issues/2378))
- Modified date/time unit tests to use unified example
  ([\#2424](https://github.com/pharmaverse/admiral/issues/2424))

## admiral 1.0.2

CRAN release: 2024-03-05

- Fix bug in
  [`derive_param_tte()`](https:/pharmaverse.github.io/admiral/main/reference/derive_param_tte.md)
  where argument `dataset` populated and `PARAMCD` in `set_values_to`
  argument is an expression. Previously, there was a check early in
  function to see if `PARAMCD` defined in `set_values_to` argument,
  already existed in the dataset passed into `dataset` argument. If
  `PARAMCD` was not an expression i.e. `PARAMCD = "XYZ"` then check
  worked. However, if `PARAMCD` to be created was an expression, and
  wasn’t resolved yet, this caused an ERROR. The check has been moved to
  near the end of the function, where `PARAMCD` is resolved in the
  dataset holding the new parameters.
  ([\#2336](https://github.com/pharmaverse/admiral/issues/2336))

## admiral 1.0.1

CRAN release: 2024-02-01

- Fix bug in
  [`derive_vars_query()`](https:/pharmaverse.github.io/admiral/main/reference/derive_vars_query.md)
  where if AE terms were in mixed case no terms are flagged.
  ([\#2311](https://github.com/pharmaverse/admiral/issues/2311))

## admiral 1.0.0

CRAN release: 2023-12-15

### New Features

- The new function
  [`derive_vars_extreme_event()`](https:/pharmaverse.github.io/admiral/main/reference/derive_vars_extreme_event.md),
  which works as
  [`derive_extreme_event()`](https:/pharmaverse.github.io/admiral/main/reference/derive_extreme_event.md)
  but adds variables instead of a parameter.
  ([\#2138](https://github.com/pharmaverse/admiral/issues/2138))

- The new function
  [`derive_var_merged_ef_msrc()`](https:/pharmaverse.github.io/admiral/main/reference/derive_var_merged_ef_msrc.md)
  is provided to add a flag indicating if one of the conditions in one
  of multiple source datasets is fulfilled.
  ([\#1728](https://github.com/pharmaverse/admiral/issues/1728))

- New global option created `signif_digits` and added to
  [`set_admiral_options()`](https:/pharmaverse.github.io/admiral/main/reference/set_admiral_options.md)
  to handle floating point issue, the value is set to `15`, and is used
  with the `base R` function
  [`signif()`](https://rdrr.io/r/base/Round.html) when comparing 2
  numeric values. This is implemented in `admiral` functions
  [`derive_var_atoxgr_dir()`](https:/pharmaverse.github.io/admiral/main/reference/derive_var_atoxgr_dir.md)
  and
  [`derive_var_anrind()`](https:/pharmaverse.github.io/admiral/main/reference/derive_var_anrind.md).
  ([\#2134](https://github.com/pharmaverse/admiral/issues/2134))

  For more information, please see blog: [How admiral handles floating
  points](https://pharmaverse.github.io/blog/posts/2023-10-30_floating_point/floating_point.html)

- The new function
  [`derive_vars_computed()`](https:/pharmaverse.github.io/admiral/main/reference/derive_vars_computed.md)
  is provided which has the same functionality as
  [`derive_param_computed()`](https:/pharmaverse.github.io/admiral/main/reference/derive_param_computed.md)
  but instead of adding the computed values as a new parameter, adds it
  as a new variable.
  ([\#2178](https://github.com/pharmaverse/admiral/issues/2178))

### Updates of Existing Functions

- Fixed a bug in
  [`compute_tmf()`](https:/pharmaverse.github.io/admiral/main/reference/compute_tmf.md)
  where the time imputation flag was being incorrectly populated when
  any of the existing time components (hour, minute and/or second) of
  the date character vector (`'--DTC'`), was imputed.
  ([\#2146](https://github.com/pharmaverse/admiral/issues/2146))

- [`derive_extreme_records()`](https:/pharmaverse.github.io/admiral/main/reference/derive_extreme_records.md),
  [`derive_var_extreme_flag()`](https:/pharmaverse.github.io/admiral/main/reference/derive_var_extreme_flag.md),[`derive_vars_joined()`](https:/pharmaverse.github.io/admiral/main/reference/derive_vars_joined.md)
  and
  [`derive_vars_merged()`](https:/pharmaverse.github.io/admiral/main/reference/derive_vars_merged.md)
  were enhanced with the arguments `true_value` and `false_value` to
  align with preexisting functions that had similar functionality.
  ([\#2125](https://github.com/pharmaverse/admiral/issues/2125))

- [`restrict_derivation()`](https:/pharmaverse.github.io/admiral/main/reference/restrict_derivation.md)
  now allows [dplyr](https://dplyr.tidyverse.org) functions like
  [`mutate()`](https://dplyr.tidyverse.org/reference/mutate.html) in the
  `derivation` argument.
  ([\#2143](https://github.com/pharmaverse/admiral/issues/2143))

- [`derive_summary_records()`](https:/pharmaverse.github.io/admiral/main/reference/derive_summary_records.md),
  [`derive_var_merged_summary()`](https:/pharmaverse.github.io/admiral/main/reference/derive_var_merged_summary.md),
  and
  [`get_summary_records()`](https:/pharmaverse.github.io/admiral/main/reference/get_summary_records.md)
  were enhanced such that more than one summary variable can be derived,
  e.g., `AVAL` as the sum and `ADT` as the maximum of the contributing
  records.
  ([\#1792](https://github.com/pharmaverse/admiral/issues/1792))

- [`derive_summary_records()`](https:/pharmaverse.github.io/admiral/main/reference/derive_summary_records.md)
  was enhanced with the following arguments: `dataset_add` (required),
  `dataset_ref` (optional), `missing_values` (optional). These arguments
  respectively, generate summary variables from additional datasets,
  retain/add specific records from a reference dataset, and impute
  user-defined missing values. Note that `dataset_add` can be set to the
  same value as `dataset` if a different additional dataset is not
  required.
  [`derive_param_exposure()`](https:/pharmaverse.github.io/admiral/main/reference/derive_param_exposure.md)
  was enhanced with `dataset_add` as well.
  ([\#2142](https://github.com/pharmaverse/admiral/issues/2142))

- The `missing_values` argument was added to
  [`derive_var_merged_summary()`](https:/pharmaverse.github.io/admiral/main/reference/derive_var_merged_summary.md).
  It allows to define values for by groups, e.g., subjects which are not
  in the additional dataset.
  ([\#2230](https://github.com/pharmaverse/admiral/issues/2230))

- The argument `dataset` is now optional for
  [`derive_summary_records()`](https:/pharmaverse.github.io/admiral/main/reference/derive_summary_records.md)
  and
  [`derive_param_exposure()`](https:/pharmaverse.github.io/admiral/main/reference/derive_param_exposure.md).
  ([\#2142](https://github.com/pharmaverse/admiral/issues/2142))

- The “joined” functions
  ([`derive_vars_joined()`](https:/pharmaverse.github.io/admiral/main/reference/derive_vars_joined.md),
  [`derive_var_joined_exist_flag()`](https:/pharmaverse.github.io/admiral/main/reference/derive_var_joined_exist_flag.md),
  [`filter_joined()`](https:/pharmaverse.github.io/admiral/main/reference/filter_joined.md),
  and
  [`event_joined()`](https:/pharmaverse.github.io/admiral/main/reference/event_joined.md))
  were unified:
  ([\#2126](https://github.com/pharmaverse/admiral/issues/2126))

  - The `dataset_add` and `filter_add` arguments were added to
    [`derive_var_joined_exist_flag()`](https:/pharmaverse.github.io/admiral/main/reference/derive_var_joined_exist_flag.md)
    and
    [`filter_joined()`](https:/pharmaverse.github.io/admiral/main/reference/filter_joined.md).
  - The `filter` argument was renamed to `filter_join` in
    [`derive_var_joined_exist_flag()`](https:/pharmaverse.github.io/admiral/main/reference/derive_var_joined_exist_flag.md)
    and
    [`filter_joined()`](https:/pharmaverse.github.io/admiral/main/reference/filter_joined.md).
  - The `tmp_obs_nr_var`, the `join_type`, the `first_cond_lower`, and
    the `first_cond_upper` arguments were added to
    [`derive_vars_joined()`](https:/pharmaverse.github.io/admiral/main/reference/derive_vars_joined.md).
  - In
    [`derive_var_joined_exist_flag()`](https:/pharmaverse.github.io/admiral/main/reference/derive_var_joined_exist_flag.md),
    [`filter_joined()`](https:/pharmaverse.github.io/admiral/main/reference/filter_joined.md),
    and
    [`event_joined()`](https:/pharmaverse.github.io/admiral/main/reference/event_joined.md)
    the `first_cond` argument was renamed to `first_cond_upper` and the
    `first_cond_lower` argument was added.
  - In all “joined” functions the `filter_add` argument is applied to
    the additional dataset grouped by `by_vars` and the `filter_join`
    argument is applied to the joined dataset grouped by the
    observations from the input dataset. I.e., summary functions like
    [`all()`](https://rdrr.io/r/base/all.html) or
    [`any()`](https://rdrr.io/r/base/any.html) can be used.

- The `tmp_event_nr_var` argument was added to
  [`derive_extreme_records()`](https:/pharmaverse.github.io/admiral/main/reference/derive_extreme_records.md)
  to allow more control of the selection of records. It creates a
  temporary variable for the event number, which can be used in `order`.
  ([\#2140](https://github.com/pharmaverse/admiral/issues/2140))

- `signif_dig` argument added to both
  [`derive_var_atoxgr_dir()`](https:/pharmaverse.github.io/admiral/main/reference/derive_var_atoxgr_dir.md)
  and
  [`derive_var_anrind()`](https:/pharmaverse.github.io/admiral/main/reference/derive_var_anrind.md)
  functions with default value set to general option `signif_digits`.
  The new argument to these functions handles any floating point issues.
  ([\#2134](https://github.com/pharmaverse/admiral/issues/2134))

- Fixed a bug in
  [`derive_vars_period()`](https:/pharmaverse.github.io/admiral/main/reference/derive_vars_period.md)
  where the function was throwing an error whenever `dataset_ref`
  contained variables that were neither key variables, nor `APERIOD`,
  `ASPER`, `APHASEN`, nor mentioned in the `new_vars` argument.
  ([\#2231](https://github.com/pharmaverse/admiral/issues/2231))

- [`compute_duration()`](https:/pharmaverse.github.io/admiral/main/reference/compute_duration.md),
  [`derive_vars_duration()`](https:/pharmaverse.github.io/admiral/main/reference/derive_vars_duration.md),
  [`derive_vars_aage()`](https:/pharmaverse.github.io/admiral/main/reference/derive_vars_aage.md)
  now accepts more terms for the `in_unit`, `out_unit`, and `age_unit`
  arguments
  ([\#2255](https://github.com/pharmaverse/admiral/issues/2255))

- Updated the unit test for
  [`derive_var_obs_number()`](https:/pharmaverse.github.io/admiral/main/reference/derive_var_obs_number.md).
  The new test checked the derivation of the default and customized
  `new_var`, sorting with the the missing value and expected conditions.
  ([\#2260](https://github.com/pharmaverse/admiral/issues/2260))

- The check for existence of `TERMNUM`/`TERMCHAR` in queries dataset is
  now less strict depending on values of `SRCVAR` for
  [`derive_vars_query()`](https:/pharmaverse.github.io/admiral/main/reference/derive_vars_query.md)
  ([\#2264](https://github.com/pharmaverse/admiral/issues/2264))

- DAIDS grading criteria fixed for `Grade = 0` for
  `TERM = "Absolute Lymphocyte Count, Low"`, criteria was
  `AVAL <= 0.65`, now corrected to `AVAL >= 0.65`
  ([\#2284](https://github.com/pharmaverse/admiral/issues/2284)).

- A bug in
  [`derive_extreme_event()`](https:/pharmaverse.github.io/admiral/main/reference/derive_extreme_event.md)
  was fixed. The `condition` field is no longer ignored if `mode` is
  specified for
  [`event()`](https:/pharmaverse.github.io/admiral/main/reference/event.md)
  ([\#2291](https://github.com/pharmaverse/admiral/issues/2291)).

- A bug in
  [`derive_vars_joined()`](https:/pharmaverse.github.io/admiral/main/reference/derive_vars_joined.md)
  was fixed. The function no longer fails if renaming is used in
  `by_vars` and `new_vars` is not specified
  ([\#2289](https://github.com/pharmaverse/admiral/issues/2289)).

### Breaking Changes

- [admiral](https://pharmaverse.github.io/admiral/) now only supports R
  \>= 4.0.0

- In
  [`derive_extreme_records()`](https:/pharmaverse.github.io/admiral/main/reference/derive_extreme_records.md)
  the `dataset_add` argument is now mandatory.
  ([\#2139](https://github.com/pharmaverse/admiral/issues/2139))

- In
  [`derive_summary_records()`](https:/pharmaverse.github.io/admiral/main/reference/derive_summary_records.md)
  and
  [`get_summary_records()`](https:/pharmaverse.github.io/admiral/main/reference/get_summary_records.md)
  the arguments `analysis_var` and `summary_fun` were deprecated in
  favor of `set_values_to`.
  ([\#1792](https://github.com/pharmaverse/admiral/issues/1792))

- In
  [`derive_summary_records()`](https:/pharmaverse.github.io/admiral/main/reference/derive_summary_records.md)
  and
  [`derive_param_exposure()`](https:/pharmaverse.github.io/admiral/main/reference/derive_param_exposure.md)
  the argument `filter` was renamed to `filter_add`
  ([\#2142](https://github.com/pharmaverse/admiral/issues/2142))

- In
  [`derive_var_merged_summary()`](https:/pharmaverse.github.io/admiral/main/reference/derive_var_merged_summary.md)
  the arguments `new_var`, `analysis_var`, and `summary_fun` were
  deprecated in favor of `new_vars`.
  ([\#1792](https://github.com/pharmaverse/admiral/issues/1792))

- In
  [`derive_vars_merged()`](https:/pharmaverse.github.io/admiral/main/reference/derive_vars_merged.md),
  the argument `match_flag` was renamed to `exist_flag`
  ([\#2125](https://github.com/pharmaverse/admiral/issues/2125))

- The default value for the `false_value` argument in
  [`derive_extreme_records()`](https:/pharmaverse.github.io/admiral/main/reference/derive_extreme_records.md)
  was changed to `NA_character_`
  ([\#2125](https://github.com/pharmaverse/admiral/issues/2125))

- In
  [`consolidate_metadata()`](https:/pharmaverse.github.io/admiral/main/reference/consolidate_metadata.md),
  the argument `check_keys` was renamed to `check_type` to align with
  other functions
  ([\#2184](https://github.com/pharmaverse/admiral/issues/2184))

- In
  [`filter_joined()`](https:/pharmaverse.github.io/admiral/main/reference/filter_joined.md)
  and
  [`derive_var_joined_exist_flag()`](https:/pharmaverse.github.io/admiral/main/reference/derive_var_joined_exist_flag.md)
  ([\#2126](https://github.com/pharmaverse/admiral/issues/2126))

  - the `first_cond` argument was deprecated in favor of
    `first_cond_upper` and
  - the `filter` argument was deprecated in favor of `filter_join`.

- In
  [`event_joined()`](https:/pharmaverse.github.io/admiral/main/reference/event_joined.md)
  the `first_cond` argument was deprecated in favor of
  `first_cond_upper`.
  ([\#2126](https://github.com/pharmaverse/admiral/issues/2126))

- In
  [`derive_extreme_event()`](https:/pharmaverse.github.io/admiral/main/reference/derive_extreme_event.md),
  the `ignore_event_order` argument was deprecated and the selection of
  the records was changed to allow more control. Before, the records
  were selected first by event and then by `order`. Now they are
  selected by `order` only, but the event number can be added to it.

  To achieve the old behavior update

      order = exprs(my_order_var),
      ignore_event_order = FALSE,

  to

      tmp_event_nr_var = event_nr,
      order = exprs(event_nr, my_order_var),

  and

      order = exprs(my_order_var),
      ignore_event_order = TRUE,

  to

      order = exprs(my_order_var),

- [`create_query_data()`](https:/pharmaverse.github.io/admiral/main/reference/create_query_data.md)
  and
  [`derive_vars_query()`](https:/pharmaverse.github.io/admiral/main/reference/derive_vars_query.md)
  were updated to rename variables in query data set as follows:
  ([\#2186](https://github.com/pharmaverse/admiral/issues/2186))

  - `TERMNAME` to `TERMCHAR`
  - `TERMID` to `TERMNUM`

  Users need to adjust their `get_terms()` function accordingly.

- The following functions, which were deprecated in previous
  [admiral](https://pharmaverse.github.io/admiral/) versions, have been
  removed:
  ([\#2098](https://github.com/pharmaverse/admiral/issues/2098))

  - `derive_param_extreme_event()`
  - `derive_vars_last_dose()`
  - `derive_var_last_dose_amt()`
  - `derive_var_last_dose_date()`
  - `derive_var_last_dose_grp()`
  - `derive_var_basetype()`
  - `derive_var_merged_cat()`
  - `derive_var_merged_character()`
  - `derive_var_confirmation_flag()`

- The following function arguments are entering the next phase of the
  deprecation process:
  ([\#2098](https://github.com/pharmaverse/admiral/issues/2098))

  - `compute_egfr(wt)`
  - `derive_extreme_records(filter)`
  - `derive_param_computed(analysis_value, analysis_var)`
  - `derive_var_shift(na_val)`
  - `derive_expected_records(dataset_expected_obs)`
  - `derive_var_ontrtfl(span_period)`

- The
  [`derive_param_extreme_record()`](https:/pharmaverse.github.io/admiral/main/reference/derive_param_extreme_record.md)
  function has been superseded in favor of
  [`derive_extreme_event()`](https:/pharmaverse.github.io/admiral/main/reference/derive_extreme_event.md).
  ([\#2141](https://github.com/pharmaverse/admiral/issues/2141))

- The functions
  [`derive_var_dthcaus()`](https:/pharmaverse.github.io/admiral/main/reference/derive_var_dthcaus.md),
  [`derive_var_extreme_dt()`](https:/pharmaverse.github.io/admiral/main/reference/derive_var_extreme_dt.md),
  and
  [`derive_var_extreme_dtm()`](https:/pharmaverse.github.io/admiral/main/reference/derive_var_extreme_dtm.md)
  have been superseded in favor of
  [`derive_vars_extreme_event()`](https:/pharmaverse.github.io/admiral/main/reference/derive_vars_extreme_event.md).
  ([\#2138](https://github.com/pharmaverse/admiral/issues/2138))

### Documentation

- The documentation of the `by_vars` and `constant_by_vars` argument was
  improved and unified across all functions where it is used.
  ([\#2137](https://github.com/pharmaverse/admiral/issues/2137))

- The functions
  [`assert_db_requirements()`](https:/pharmaverse.github.io/admiral/main/reference/assert_db_requirements.md),
  [`assert_terms()`](https:/pharmaverse.github.io/admiral/main/reference/assert_terms.md),
  `assert_valid_queries()`, `extend_source_datasets()`,
  [`filter_date_sources()`](https:/pharmaverse.github.io/admiral/main/reference/filter_date_sources.md),
  [`validate_basket_select()`](https:/pharmaverse.github.io/admiral/main/reference/validate_basket_select.md),
  [`validate_query()`](https:/pharmaverse.github.io/admiral/main/reference/validate_query.md)
  are no longer exported and have had documentation removed.
  ([\#2220](https://github.com/pharmaverse/admiral/issues/2220))

- The function
  [`extract_duplicate_records()`](https:/pharmaverse.github.io/admiral/main/reference/extract_duplicate_records.md)
  has been re-classified as an `internal` function, which means that the
  function still appears in our help pages but not on our website.
  ([\#2220](https://github.com/pharmaverse/admiral/issues/2220))

- The “Generic Functions” vignette (now “Generic Derivations”) was
  rewritten. Now it provides a more complete overview of the generic
  derivations, describe the common concepts, and makes it easier to find
  the appropriate function.
  ([\#2230](https://github.com/pharmaverse/admiral/issues/2230))

- A way to standardize roxygen labels and descriptions for function
  arguments was implemented and tested.
  ([\#2034](https://github.com/pharmaverse/admiral/issues/2034))

- A link to published CDISC Population PK (ADPPK) implementation guide
  was added.
  ([\#2161](https://github.com/pharmaverse/admiral/issues/2161))

- Removed Deprecation section in Reference tab. Added new Superseded
  section in Reference tab.
  ([\#2174](https://github.com/pharmaverse/admiral/issues/2174))

- Added a link to the previous versions of the website to the navigation
  bar. ([\#2205](https://github.com/pharmaverse/admiral/issues/2205))

- The meaning of `date_imputation = "mid"` was clarified in the
  documentation of the imputation functions,
  e.g. [`derive_vars_dtm()`](https:/pharmaverse.github.io/admiral/main/reference/derive_vars_dtm.md).
  ([\#2222](https://github.com/pharmaverse/admiral/issues/2222))

- Added an example derivation of `DTHCGR1` to the ADSL vignette.
  ([\#2218](https://github.com/pharmaverse/admiral/issues/2218))

- Moved Development Process from
  [admiraldev](https://pharmaverse.github.io/admiraldev/) to
  Contribution Model in the `admiral` website, updated GitHub strategy.
  ([\#2196](https://github.com/pharmaverse/admiral/issues/2196))

- Added new drop downs in Get Started navigation bar- Getting Started,
  Admiral Discovery, Cheat Sheet. Community removed from top.
  ([\#2217](https://github.com/pharmaverse/admiral/issues/2217))

- All “Example Script(s)” sections in the User Guide vignettes were
  updated to point the user towards using `use_ad_template("ADaM")`
  rather than linking to the template in the code repository.
  ([\#2239](https://github.com/pharmaverse/admiral/issues/2239))

- Handling of `NA` values was added to the documentation of the `order`
  argument for all functions.
  ([\#2230](https://github.com/pharmaverse/admiral/issues/2230),
  [\#2257](https://github.com/pharmaverse/admiral/issues/2257))

### Various

- Website now has button/links to Slack channel and GitHub Issues.
  ([\#2127](https://github.com/pharmaverse/admiral/issues/2127))

- Added example derivations of `DTHCAUS` and `DTHCGR1` to the ADSL
  template.
  ([\#2218](https://github.com/pharmaverse/admiral/issues/2218))

- Cheat Sheet now added to website front page
  ([\#2130](https://github.com/pharmaverse/admiral/issues/2130))

## admiral 0.12.3

CRAN release: 2023-10-18

- Fixed a bug in
  [`derive_var_dthcaus()`](https:/pharmaverse.github.io/admiral/main/reference/derive_var_dthcaus.md)
  where if a subject has observations in more than one of the sources,
  the one from the last source was selected regardless of the date. Now
  the function works as described in its documentation.
  ([\#2154](https://github.com/pharmaverse/admiral/issues/2154))

## admiral 0.12.2

CRAN release: 2023-10-06

- A unit test for
  [`derive_param_computed()`](https:/pharmaverse.github.io/admiral/main/reference/derive_param_computed.md)
  was modified in anticipation of major user-facing changes to R version
  4.4 ([\#2147](https://github.com/pharmaverse/admiral/issues/2147))

## admiral 0.12.1

CRAN release: 2023-09-25

- [`derive_extreme_records()`](https:/pharmaverse.github.io/admiral/main/reference/derive_extreme_records.md)
  no longer fails if `dataset_add` is specified and a variable specified
  for `order` is not in `dataset`.
  ([\#2113](https://github.com/pharmaverse/admiral/issues/2113))

- The `type` argument in
  [`compute_duration()`](https:/pharmaverse.github.io/admiral/main/reference/compute_duration.md)
  changed the underlying default behavior in
  [`derive_vars_duration()`](https:/pharmaverse.github.io/admiral/main/reference/derive_vars_duration.md)
  without allowing the user to toggle between `"duration"` and
  `"interval"` as originally intended. This was fixed by adding the
  `type` argument for
  [`derive_vars_duration()`](https:/pharmaverse.github.io/admiral/main/reference/derive_vars_duration.md)
  and a wrapper function
  [`derive_vars_aage()`](https:/pharmaverse.github.io/admiral/main/reference/derive_vars_aage.md)
  such that it gets passed through
  [`compute_duration()`](https:/pharmaverse.github.io/admiral/main/reference/compute_duration.md)
  appropriately
  ([\#2112](https://github.com/pharmaverse/admiral/issues/2112))

- Template `ad_adpp.R` updated to replace
  [`left_join()`](https://dplyr.tidyverse.org/reference/mutate-joins.html)
  with
  [`derive_vars_merged()`](https:/pharmaverse.github.io/admiral/main/reference/derive_vars_merged.md)
  ([\#2109](https://github.com/pharmaverse/admiral/issues/2109)).

## admiral 0.12.0

CRAN release: 2023-09-12

### New Features

- [`event_joined()`](https:/pharmaverse.github.io/admiral/main/reference/event_joined.md)
  events were added. They can be specified for the `events` argument in
  [`derive_extreme_event()`](https:/pharmaverse.github.io/admiral/main/reference/derive_extreme_event.md).
  This allows to define events based on more than one observation, e.g.,
  events which need to be confirmed by a second assessment.
  ([\#1960](https://github.com/pharmaverse/admiral/issues/1960))

- `atoxgr_criteria_daids.rda` added, which holds metadata for [Division
  of AIDS (DAIDS) Table for Grading the Severity of Adult and Pediatric
  Adverse
  Events](https://rsc.niaid.nih.gov/sites/default/files/daidsgradingcorrectedv21.pdf).
  You can find additional documentation here
  [`atoxgr_criteria_daids()`](https:/pharmaverse.github.io/admiral/main/reference/atoxgr_criteria_daids.md)

### Updates of Existing Functions

- The functions
  [`derive_param_bmi()`](https:/pharmaverse.github.io/admiral/main/reference/derive_param_bmi.md)
  and
  [`derive_param_bsa()`](https:/pharmaverse.github.io/admiral/main/reference/derive_param_bsa.md)
  are updated to have the option of producing more values at visits when
  only weight is collected
  ([\#1228](https://github.com/pharmaverse/admiral/issues/1228)).

- The functions
  [`derive_var_age_years()`](https:/pharmaverse.github.io/admiral/main/reference/derive_var_age_years.md)
  and
  [`compute_age_years()`](https:/pharmaverse.github.io/admiral/main/reference/compute_age_years.md)
  are updated to return an `NA` age in the case that the age unit is
  missing.
  ([\#2001](https://github.com/pharmaverse/admiral/issues/2001)) The
  argument `unit` for
  [`derive_vars_aage()`](https:/pharmaverse.github.io/admiral/main/reference/derive_vars_aage.md)
  is also changed to `age_unit` for consistency between these
  age-related functions.
  ([\#2025](https://github.com/pharmaverse/admiral/issues/2025))

- The
  [`derive_var_ontrtfl()`](https:/pharmaverse.github.io/admiral/main/reference/derive_var_ontrtfl.md)
  function has been updated to allow for the column passed in
  `ref_end_date` to contain `NA` values. Previously, if the end date was
  `NA`, the row would never be flagged. Now, an `NA` value is
  interpreted as the treatment being ongoing, for example.
  ([\#1984](https://github.com/pharmaverse/admiral/issues/1984))

- The function
  [`derive_var_extreme_flag()`](https:/pharmaverse.github.io/admiral/main/reference/derive_var_extreme_flag.md)
  has a new function argument, `flag_all` that additionally flags all
  records if the first or last record is not unique.
  ([\#1979](https://github.com/pharmaverse/admiral/issues/1979))

- The function
  [`derive_param_computed()`](https:/pharmaverse.github.io/admiral/main/reference/derive_param_computed.md)
  was enhanced:
  ([\#1968](https://github.com/pharmaverse/admiral/issues/1968))

  - The `analysis_value` and `analysis_var` arguments were deprecated in
    favor of `set_values_to`. This enables users to compute more than
    one variable.
  - The `keep_nas` argument was added. If it is set to `TRUE`,
    observations are created even if values contributing to the computed
    values are `NA`.

- The function
  [`derive_vars_dy()`](https:/pharmaverse.github.io/admiral/main/reference/derive_vars_dy.md)
  is updated to avoid potential error when the input `dataset` with
  columns ending with `temp`.
  ([\#2012](https://github.com/pharmaverse/admiral/issues/2012))

- Argument `keep_source_vars` was added to
  [`derive_extreme_records()`](https:/pharmaverse.github.io/admiral/main/reference/derive_extreme_records.md)
  which specifies which variables in the new observations should be
  kept. ([\#1697](https://github.com/pharmaverse/admiral/issues/1697))

- Templates, vignettes, and other uses of `{admiral.test}` SDTM data are
  updated to use
  [pharmaversesdtm](https://pharmaverse.github.io/pharmaversesdtm/)
  instead.
  ([\#2040](https://github.com/pharmaverse/admiral/issues/2040))

- The `traceability_vars` argument in
  [`date_source()`](https:/pharmaverse.github.io/admiral/main/reference/date_source.md)
  and `dthcaus_source` were deprecated in favor of `set_values_to`. The
  [`date_source()`](https:/pharmaverse.github.io/admiral/main/reference/date_source.md)
  function creates a date_source object as input for
  [`derive_var_extreme_dt()`](https:/pharmaverse.github.io/admiral/main/reference/derive_var_extreme_dt.md)
  and
  [`derive_var_extreme_dtm()`](https:/pharmaverse.github.io/admiral/main/reference/derive_var_extreme_dtm.md),users
  can now define the traceability variables by assigning those variables
  to the `set_values_to`argument.Similarly, the `dthcaus_source` creates
  a dthcaus_source Object.
  ([\#2068](https://github.com/pharmaverse/admiral/issues/2068))

- [`derive_extreme_event()`](https:/pharmaverse.github.io/admiral/main/reference/derive_extreme_event.md)
  was enhanced
  ([\#1960](https://github.com/pharmaverse/admiral/issues/1960)):

  - [`event_joined()`](https:/pharmaverse.github.io/admiral/main/reference/event_joined.md)
    events can be specified for the `events` argument. This allows to
    define events based on more than one observation, e.g., events which
    need to be confirmed by a second assessment.

  - The `source_datasets` argument was added to the function and the
    `dataset_name` field to
    [`event()`](https:/pharmaverse.github.io/admiral/main/reference/event.md).
    It can be used to define events based on a different dataset than
    the input dataset.

  - The `keep_source_vars` argument was added to the function and the
    `keep_source_vars` field to
    [`event()`](https:/pharmaverse.github.io/admiral/main/reference/event.md).
    It allows to select which variables should be kept for the selected
    observations.

  - The `mode` and `order` field were added to
    [`event()`](https:/pharmaverse.github.io/admiral/main/reference/event.md).
    They allow to select the first or last observation per by group if
    there are multiple observation fulfilling the event condition.

  - The `ignore_event_order` argument was added.

  - The `description` field was added to
    [`event()`](https:/pharmaverse.github.io/admiral/main/reference/event.md).
    It can be used to provide a description of the event in plain
    language.

- [`derive_var_atoxgr_dir()`](https:/pharmaverse.github.io/admiral/main/reference/derive_var_atoxgr_dir.md)
  was enhanced
  ([\#1859](https://github.com/pharmaverse/admiral/issues/1859)):

  - Can now select `atoxgr_criteria_daids` in argument `meta_criteria`
    to create `ATOXGRL` and `ATOXGRH` based on [Division of AIDS (DAIDS)
    Table for Grading the Severity of Adult and Pediatric Adverse
    Events](https://rsc.niaid.nih.gov/sites/default/files/daidsgradingcorrectedv21.pdf)

  - New argument `signif_dig` added to control the number of significant
    digits to use when comparing 2 numeric values.
    (<https://github.com/pharmaverse/admiral/pull/2060>)

### Breaking Changes

- The `compute_duration(type)` argument added the `"duration"` type
  calculation, and this is the new default (previously `"interval"`
  differences were returned). See function help file for details on the
  difference between `"duration"` and `"interval"` calculations.
  ([\#1875](https://github.com/pharmaverse/admiral/issues/1875))

- The following functions, which were deprecated in previous
  [admiral](https://pharmaverse.github.io/admiral/) versions, have been
  removed:
  ([\#1950](https://github.com/pharmaverse/admiral/issues/1950))

  - `derive_var_disposition_status()`
  - `derive_vars_disposition_reason()`
  - `format_eoxxstt_default()`
  - `format_reason_default()`
  - `derive_var_worst_flag()`

- The following functions have been deprecated from previous
  [admiral](https://pharmaverse.github.io/admiral/) versions using the
  next phase of the deprecation process:
  ([\#1950](https://github.com/pharmaverse/admiral/issues/1950))

  - `derive_param_extreme_event()`
  - `derive_vars_last_dose()`
  - `derive_var_last_dose_amt()`
  - `derive_var_last_dose_date()`
  - `derive_var_last_dose_grp()`
  - `derive_var_basetype()`
  - `derive_var_merged_cat()`
  - `derive_var_merged_character()`

- The arguments `dataset_adsl` in the function
  [`derive_param_exist_flag()`](https:/pharmaverse.github.io/admiral/main/reference/derive_param_exist_flag.md)
  and `subject_keys` have been deprecated versions using the next phase
  of the deprecation process.
  ([\#1950](https://github.com/pharmaverse/admiral/issues/1950))

- The argument `wt` in the function
  [`compute_egfr()`](https:/pharmaverse.github.io/admiral/main/reference/compute_egfr.md)
  was deprecated in favor of `weight` using the first phase of the
  deprecation process.
  ([\#2020](https://github.com/pharmaverse/admiral/issues/2020))

- The `filter` argument in
  [`derive_extreme_records()`](https:/pharmaverse.github.io/admiral/main/reference/derive_extreme_records.md)
  was deprecated in favor of the `filter_add` using the next phase of
  the deprecation process.
  ([\#1950](https://github.com/pharmaverse/admiral/issues/1950))

- The `analysis_value` and `analysis_var` arguments in
  [`derive_param_computed()`](https:/pharmaverse.github.io/admiral/main/reference/derive_param_computed.md)
  were deprecated in favor of `set_values_to`
  ([\#1968](https://github.com/pharmaverse/admiral/issues/1968)).

- The `na_val` argument in
  [`derive_var_shift()`](https:/pharmaverse.github.io/admiral/main/reference/derive_var_shift.md)
  has been deprecated in favor of `missing_value` using the first phase
  of the deprecation process.
  ([\#2014](https://github.com/pharmaverse/admiral/issues/2014))

- The `dataset_expected_obs` argument in
  [`derive_expected_records()`](https:/pharmaverse.github.io/admiral/main/reference/derive_expected_records.md)
  and
  [`derive_locf_records()`](https:/pharmaverse.github.io/admiral/main/reference/derive_locf_records.md)
  has been deprecated in favor of `dataset_ref`.
  ([\#2037](https://github.com/pharmaverse/admiral/issues/2037))

- The `span_period` argument in
  [`derive_var_ontrtfl()`](https:/pharmaverse.github.io/admiral/main/reference/derive_var_ontrtfl.md)
  has been updated to only accept `TRUE` or `FALSE`, where is previously
  accepted `"Y"` and `NULL`.
  ([\#2033](https://github.com/pharmaverse/admiral/issues/2033))

### Documentation

- Non-exported utility and print functions were previously listed on the
  admiral website reference page. They have been removed.
  ([\#2049](https://github.com/pharmaverse/admiral/issues/2049),
  [\#2050](https://github.com/pharmaverse/admiral/issues/2050))

- The description of the argument `reference_date` in the function
  [`derive_vars_dy()`](https:/pharmaverse.github.io/admiral/main/reference/derive_vars_dy.md)
  has been clarified to make it agnostic to start/end selection.
  ([\#2027](https://github.com/pharmaverse/admiral/issues/2027))

- Date and Time Imputation User Guide/Vignette has section on preserving
  partial dates updated
  ([\#2028](https://github.com/pharmaverse/admiral/issues/2028))

### Various

- The list of package authors/contributors has been reformatted so that
  those who are actively maintaining the code base are now marked as
  *authors*, whereas those who made a significant contribution in the
  past are now down as *contributors*. All other acknowledgments have
  been moved to README section
  ([\#1941](https://github.com/pharmaverse/admiral/issues/1941)).

- [`derive_vars_joined()`](https:/pharmaverse.github.io/admiral/main/reference/derive_vars_joined.md)
  had two bugs with regards to duplicates messaging and when `new_vars`
  was set to `NULL` that have now been addressed
  ([\#1966](https://github.com/pharmaverse/admiral/issues/1966)).

- [`compute_dtf()`](https:/pharmaverse.github.io/admiral/main/reference/compute_dtf.md)
  had a bug with regards to imputing days to full date-time character
  strings.
  ([\#2042](https://github.com/pharmaverse/admiral/issues/2042))

## admiral 0.11.1

CRAN release: 2023-07-06

- Fix bug in
  [`derive_param_tte()`](https:/pharmaverse.github.io/admiral/main/reference/derive_param_tte.md).
  ([\#1962](https://github.com/pharmaverse/admiral/issues/1962))
- Get Started page now points to correct article.
  ([\#1969](https://github.com/pharmaverse/admiral/issues/1969))

## admiral 0.11.0

CRAN release: 2023-06-08

### New Features

- In the function
  [`derive_var_anrind()`](https:/pharmaverse.github.io/admiral/main/reference/derive_var_anrind.md),
  added argument `use_a1hia1lo` to turn the usage of `A1HI` and `A1LO`
  off and on, with the default being off.
  ([\#1795](https://github.com/pharmaverse/admiral/issues/1795))

- Added a “Report a bug” link to
  [admiral](https://pharmaverse.github.io/admiral/) website.
  ([\#1836](https://github.com/pharmaverse/admiral/issues/1836))

- New function
  [`compute_age_years()`](https:/pharmaverse.github.io/admiral/main/reference/compute_age_years.md)
  for converting a vector of age values to years.
  ([\#1794](https://github.com/pharmaverse/admiral/issues/1794))

- New functions
  [`filter_exist()`](https:/pharmaverse.github.io/admiral/main/reference/filter_exist.md)
  and
  [`filter_not_exist()`](https:/pharmaverse.github.io/admiral/main/reference/filter_not_exist.md)
  for selecting records from a dataset dependent on the existence of the
  corresponding by groups in a filtered source dataset.
  ([\#1699](https://github.com/pharmaverse/admiral/issues/1699))

- New function
  [`derive_param_extreme_record()`](https:/pharmaverse.github.io/admiral/main/reference/derive_param_extreme_record.md)
  that adds parameter based on the first or last record from multiple
  sources.
  ([\#1822](https://github.com/pharmaverse/admiral/issues/1822))

- New ADPPK template script available `ad_adppk.R` which creates
  Population PK Analysis Dataset based on forthcoming CDISC
  Implementation guide.
  ([\#1772](https://github.com/pharmaverse/admiral/issues/1772))

- New function
  [`compute_egfr()`](https:/pharmaverse.github.io/admiral/main/reference/compute_egfr.md)
  for calculating Estimated Glomerular Filtration Rate (eGFR) and
  Creatinine Clearance for Kidney Function.
  ([\#1826](https://github.com/pharmaverse/admiral/issues/1826))

### Updates of Existing Functions

- [`derive_extreme_records()`](https:/pharmaverse.github.io/admiral/main/reference/derive_extreme_records.md)
  was enhanced such that it includes the functionality of
  `derive_param_extreme_event()`.
  ([\#1725](https://github.com/pharmaverse/admiral/issues/1725))

- For the `set_values_to` argument expressions are accepted now. For
  example, `set_values_to = exprs(PARAMCD = str_to_upper(QSTESTCD))`.
  This affects
  [`censor_source()`](https:/pharmaverse.github.io/admiral/main/reference/censor_source.md),
  [`derive_expected_records()`](https:/pharmaverse.github.io/admiral/main/reference/derive_expected_records.md),
  [`derive_extreme_event()`](https:/pharmaverse.github.io/admiral/main/reference/derive_extreme_event.md),
  [`derive_extreme_records()`](https:/pharmaverse.github.io/admiral/main/reference/derive_extreme_records.md),
  [`derive_param_bmi()`](https:/pharmaverse.github.io/admiral/main/reference/derive_param_bmi.md),
  [`derive_param_bsa()`](https:/pharmaverse.github.io/admiral/main/reference/derive_param_bsa.md),
  [`derive_param_computed()`](https:/pharmaverse.github.io/admiral/main/reference/derive_param_computed.md),
  [`derive_param_doseint()`](https:/pharmaverse.github.io/admiral/main/reference/derive_param_doseint.md),
  [`derive_param_exposure()`](https:/pharmaverse.github.io/admiral/main/reference/derive_param_exposure.md),
  [`derive_param_framingham()`](https:/pharmaverse.github.io/admiral/main/reference/derive_param_framingham.md),
  [`derive_param_map()`](https:/pharmaverse.github.io/admiral/main/reference/derive_param_map.md),
  [`derive_param_exist_flag()`](https:/pharmaverse.github.io/admiral/main/reference/derive_param_exist_flag.md),
  `derive_param_extreme_event()`,
  [`derive_param_qtc()`](https:/pharmaverse.github.io/admiral/main/reference/derive_param_qtc.md),
  [`derive_param_rr()`](https:/pharmaverse.github.io/admiral/main/reference/derive_param_rr.md),
  [`derive_param_wbc_abs()`](https:/pharmaverse.github.io/admiral/main/reference/derive_param_wbc_abs.md),
  [`derive_summary_records()`](https:/pharmaverse.github.io/admiral/main/reference/derive_summary_records.md),
  [`event_source()`](https:/pharmaverse.github.io/admiral/main/reference/event_source.md),
  [`get_summary_records()`](https:/pharmaverse.github.io/admiral/main/reference/get_summary_records.md).
  ([\#1727](https://github.com/pharmaverse/admiral/issues/1727))

- For the `order` argument expressions are accepted now.
  ([\#1727](https://github.com/pharmaverse/admiral/issues/1727))

- [`derive_vars_merged()`](https:/pharmaverse.github.io/admiral/main/reference/derive_vars_merged.md)
  updates:
  ([\#1727](https://github.com/pharmaverse/admiral/issues/1727))

  - The `missing_values` argument to assign values to the new variables
    for non-matching observations was added.
  - The `new_vars` argument accepts expressions now.

- [`derive_vars_joined()`](https:/pharmaverse.github.io/admiral/main/reference/derive_vars_joined.md)
  updates:
  ([\#1727](https://github.com/pharmaverse/admiral/issues/1727))

  - The `missing_values` argument to assign values to the new variables
    for non-matching observations was added.
  - The `new_vars` and the `join_vars` argument accept expressions now.

- The `date` field of
  [`date_source()`](https:/pharmaverse.github.io/admiral/main/reference/date_source.md)
  accepts expressions now. This affects
  [`derive_var_extreme_dt()`](https:/pharmaverse.github.io/admiral/main/reference/derive_var_extreme_dt.md)
  and
  [`derive_var_extreme_dtm()`](https:/pharmaverse.github.io/admiral/main/reference/derive_var_extreme_dtm.md).
  ([\#1727](https://github.com/pharmaverse/admiral/issues/1727))

- The `date` and `dthcaus` field of
  [`dthcaus_source()`](https:/pharmaverse.github.io/admiral/main/reference/dthcaus_source.md)
  accept expressions now. This affects
  [`derive_var_dthcaus()`](https:/pharmaverse.github.io/admiral/main/reference/derive_var_dthcaus.md).
  ([\#1727](https://github.com/pharmaverse/admiral/issues/1727))

- The `date` field of
  [`event_source()`](https:/pharmaverse.github.io/admiral/main/reference/event_source.md)
  and
  [`censor_source()`](https:/pharmaverse.github.io/admiral/main/reference/censor_source.md)
  accepts expressions now. This affects
  [`derive_param_tte()`](https:/pharmaverse.github.io/admiral/main/reference/derive_param_tte.md).
  ([\#1727](https://github.com/pharmaverse/admiral/issues/1727))

- The
  [`derive_param_computed()`](https:/pharmaverse.github.io/admiral/main/reference/derive_param_computed.md)
  function was enhanced:
  ([\#1873](https://github.com/pharmaverse/admiral/issues/1873))

  - The new `dataset_add` argument allows to consider parameters from a
    different dataset than the input dataset.
  - The new `analysis_var` argument allows to specify the variable to be
    populated, e.g., `AVALC`.
  - For `parameters` and `constant_parameters` a list of expressions can
    be specified now. This allows to create temporary parameter codes,
    e.g., if SDTM data is used as input.
  - The `analysis_value` argument was enhanced such that any variable of
    the form `<variable>.<parameter>` can be used, e.g.,
    `QSORRES.CHSF13`.

### Breaking Changes

- [`create_query_data()`](https:/pharmaverse.github.io/admiral/main/reference/create_query_data.md)
  and
  [`derive_vars_query()`](https:/pharmaverse.github.io/admiral/main/reference/derive_vars_query.md)
  updated to rename variables in query data set as follows:
  ([\#1907](https://github.com/pharmaverse/admiral/issues/1907))

  - `VAR_PREFIX` to `PREFIX`
  - `QUERY_NAME` to `GRPNAME`
  - `QUERY_ID` to `GRPID`
  - `QUERY_SCOPE` to `SCOPE`
  - `QUERY_SCOPE_NUM` to `SCOPEN`
  - `TERM_LEVEL` to `SRCVAR`
  - `TERM_NAME` to `TERMNAME`
  - `TERM_ID` to `TERMID`

  Users need to adjust their `get_terms()` function accordingly.

- The `aval_fun` argument of
  [`derive_param_exist_flag()`](https:/pharmaverse.github.io/admiral/main/reference/derive_param_exist_flag.md)
  was deprecated in favor of the `set_values_to` argument.
  ([\#1727](https://github.com/pharmaverse/admiral/issues/1727))

- `derive_var_merged_cat()` and `derive_var_merged_character()` have
  been deprecated in favor of
  [`derive_vars_merged()`](https:/pharmaverse.github.io/admiral/main/reference/derive_vars_merged.md).
  ([\#1727](https://github.com/pharmaverse/admiral/issues/1727))

- The following functions, which were deprecated in previous
  [admiral](https://pharmaverse.github.io/admiral/) versions, have been
  removed:
  ([\#1747](https://github.com/pharmaverse/admiral/issues/1747))

  - `derive_vars_merged_dt()`
  - `derive_vars_merged_dtm()`
  - `derive_var_agegr_ema()`
  - `derive_var_agegr_fda()`
  - `derive_param_first_event()`
  - `derive_derived_param()`
  - `derive_var_confirmation_flag()`
  - `filter_confirmation()`

- The following functions have been deprecated from previous
  [admiral](https://pharmaverse.github.io/admiral/) versions using the
  next phase of the deprecation process:
  ([\#1747](https://github.com/pharmaverse/admiral/issues/1747))

  - `derive_var_disposition_status()`
  - `derive_vars_disposition_reason()`
  - `format_eoxxstt_default()`
  - `format_reason_default()`
  - `derive_var_worst_flag()`

- `derive_param_extreme_event()` was deprecated in favor of
  [`derive_extreme_records()`](https:/pharmaverse.github.io/admiral/main/reference/derive_extreme_records.md).
  ([\#1725](https://github.com/pharmaverse/admiral/issues/1725))

- The `filter` argument in
  [`derive_extreme_records()`](https:/pharmaverse.github.io/admiral/main/reference/derive_extreme_records.md)
  was deprecated in favor of the `filter_add` argument.
  ([\#1725](https://github.com/pharmaverse/admiral/issues/1725))

- `derive_vars_last_dose()`, `derive_var_last_dose_amt()`,
  `derive_var_last_dose_date()`, `derive_var_last_dose_grp()`, were
  deprecated in favor of
  [`derive_vars_joined()`](https:/pharmaverse.github.io/admiral/main/reference/derive_vars_joined.md).
  ([\#1797](https://github.com/pharmaverse/admiral/issues/1797))

- `derive_var_basetype()` was deprecated in favor of
  [`derive_basetype_records()`](https:/pharmaverse.github.io/admiral/main/reference/derive_basetype_records.md).
  ([\#1796](https://github.com/pharmaverse/admiral/issues/1796))

- In the function
  [`derive_param_exist_flag()`](https:/pharmaverse.github.io/admiral/main/reference/derive_param_exist_flag.md)
  the arguments `dataset_adsl` and `subject_keys` have been renamed to
  `dataset_ref` and `by_vars` respectively.
  ([\#1793](https://github.com/pharmaverse/admiral/issues/1793))

### Documentation

- Updated example dataset to trigger deterioration flags in the vignette
  “Creating Questionnaire ADaMs”.
  ([\#1853](https://github.com/pharmaverse/admiral/issues/1853),
  [\#1854](https://github.com/pharmaverse/admiral/issues/1854))

- Updated PK Programming vignette to include new Population PK Template
  `ad_adppk.R`.
  ([\#1772](https://github.com/pharmaverse/admiral/issues/1772))

- Updated “Lab Grading” Vignette to link to grading metadata available
  in [admiral](https://pharmaverse.github.io/admiral/) and clarify how
  abnormal baseline values are assigned in NCI-CTCAEv5.
  ([\#1863](https://github.com/pharmaverse/admiral/issues/1863))

- Updated “Visit and Period Variables” Vignette to add more detail about
  Study Specific Code that is required.
  ([\#1831](https://github.com/pharmaverse/admiral/issues/1831))

- Increased documentation for those functions which are regarded as
  wrapper functions.
  ([\#1726](https://github.com/pharmaverse/admiral/issues/1726))

- Examples in function documentation no longer rely on
  [`library(admiral.test)`](https://rdrr.io/r/base/library.html).
  ([\#1752](https://github.com/pharmaverse/admiral/issues/1752))

- Conferences where [admiral](https://pharmaverse.github.io/admiral/)
  was presented were updated on the `README.md`.
  ([\#1890](https://github.com/pharmaverse/admiral/issues/1890))

### Various

- [`vars()`](https://dplyr.tidyverse.org/reference/vars.html) which was
  used in the admiral function calls that expected a list of quosures
  has been removed. The admiral option `force_admiral_vars` was removed
  as well.
  ([\#1694](https://github.com/pharmaverse/admiral/issues/1694))

- [`derive_vars_dtm()`](https:/pharmaverse.github.io/admiral/main/reference/derive_vars_dtm.md)
  and
  [`derive_vars_dt()`](https:/pharmaverse.github.io/admiral/main/reference/derive_vars_dt.md)
  had a bug pertaining to imputations associated with supplying both
  `min_dates` and `max_dates` that has now been resolved.
  ([\#1843](https://github.com/pharmaverse/admiral/issues/1843))

- Examples for
  [`derive_var_extreme_flag()`](https:/pharmaverse.github.io/admiral/main/reference/derive_var_extreme_flag.md)
  were reworked to reduce runtime that occasionally led to failing CI
  check. ([\#1780](https://github.com/pharmaverse/admiral/issues/1780))

- [`create_period_dataset()`](https:/pharmaverse.github.io/admiral/main/reference/create_period_dataset.md)
  had a bug that led to an error when both DT and DTM columns existed.
  ([\#1845](https://github.com/pharmaverse/admiral/issues/1845))

- External functions are now consistently imported via namespace.
  `package::function()` calls have been removed from `admiral`
  functions.
  ([\#1842](https://github.com/pharmaverse/admiral/issues/1842))

- [`restrict_derivation()`](https:/pharmaverse.github.io/admiral/main/reference/restrict_derivation.md)
  had a bug which led to failure if the `derivation` argument was not in
  the global environment.
  ([\#1765](https://github.com/pharmaverse/admiral/issues/1765))

## admiral 0.10.2

CRAN release: 2023-04-25

- Changing package maintainer from Thomas Neitmann to Ben Straub.
  ([\#1848](https://github.com/pharmaverse/admiral/issues/1848))

## admiral 0.10.1

CRAN release: 2023-03-14

- Fix checks on
  [`derive_vars_dtm()`](https:/pharmaverse.github.io/admiral/main/reference/derive_vars_dtm.md)
  and
  [`derive_vars_dt()`](https:/pharmaverse.github.io/admiral/main/reference/derive_vars_dt.md)
  that were too restrictive.
  ([\#1810](https://github.com/pharmaverse/admiral/issues/1810))

## admiral 0.10.0

CRAN release: 2023-03-08

### New Features

- Using [testthat](https://testthat.r-lib.org) 3rd edition for unit
  testing. This is stricter in that messages must be addressed and
  deprecated functions throw errors.
  ([\#1754](https://github.com/pharmaverse/admiral/issues/1754))

- New function
  [`consolidate_metadata()`](https:/pharmaverse.github.io/admiral/main/reference/consolidate_metadata.md)
  for consolidating multiple meta datasets into a single one.
  ([\#1479](https://github.com/pharmaverse/admiral/issues/1479))

- New function
  [`compute_scale()`](https:/pharmaverse.github.io/admiral/main/reference/compute_scale.md)
  for computing the average of a vector and transforming the result from
  a source to a target range.
  ([\#1692](https://github.com/pharmaverse/admiral/issues/1692))

- New ADPC template script available `ad_adpc.R` which creates PK
  Concentration Analysis Dataset
  ([\#849](https://github.com/pharmaverse/admiral/issues/849)). This
  script includes formatting suitable for Non-Compartmental Analysis
  (ADNCA). ([\#851](https://github.com/pharmaverse/admiral/issues/851))

- New function
  [`derive_expected_records()`](https:/pharmaverse.github.io/admiral/main/reference/derive_expected_records.md)
  for adding expected records.
  ([\#1729](https://github.com/pharmaverse/admiral/issues/1729))

- New function
  [`derive_extreme_event()`](https:/pharmaverse.github.io/admiral/main/reference/derive_extreme_event.md)
  for adding the worst or best observation for each by group as new
  records.
  ([\#1755](https://github.com/pharmaverse/admiral/issues/1755))

### Updates of Existing Functions

- Arguments `analysis_var`, `keep_vars` were added to
  [`derive_locf_records()`](https:/pharmaverse.github.io/admiral/main/reference/derive_locf_records.md),  
  `analysis_var` allows to specify analysis variable, `keep_vars` keeps
  variables that need carrying the last observation forward other than
  `analysis_var` (e.g., `PARAMN`, `VISITNUM`).
  ([\#1636](https://github.com/pharmaverse/admiral/issues/1636))

- The function
  [`create_single_dose_dataset()`](https:/pharmaverse.github.io/admiral/main/reference/create_single_dose_dataset.md)
  adds support for expanding relative nominal time (e.g. NFRLT) used in
  Pharmacokinetic (PK) analyses. The new parameter `nominal_time`
  defaults as `NULL` and does not change the normal operation of the
  function. If a `nominal_time` is specified such as NFRLT (Nominal
  Relative Time from First Dose) then the nominal time is incremented by
  the interval specified in `EXDOSFRQ` for example for “QD” records the
  NFRLT is incremented by 24 hours, e.g. 0, 24, 48…
  ([\#1640](https://github.com/pharmaverse/admiral/issues/1640))

- [`create_single_dose_dataset()`](https:/pharmaverse.github.io/admiral/main/reference/create_single_dose_dataset.md)
  is also updated for values of `EXDOSFRQ` with units in days but
  expected values less than 24 hours, such as “BID”, “TID”, and “QID”.
  Previously these values of `EXDOSFRQ` may result in duplicate records
  where the day values are incremented but the time values are not.
  ([\#1643](https://github.com/pharmaverse/admiral/issues/1643))

- The function `derive_var_confirmation_flag()` and
  `filter_confirmation()` gained the `tmp_obs_nr_var` argument. It helps
  flagging or selecting consecutive observations or the first or last
  observation in a by group.
  ([\#1724](https://github.com/pharmaverse/admiral/issues/1724))

- The functions
  [`derive_vars_merged()`](https:/pharmaverse.github.io/admiral/main/reference/derive_vars_merged.md),
  `derive_var_merged_cat()`, `derive_var_merged_character()`,
  [`derive_var_merged_exist_flag()`](https:/pharmaverse.github.io/admiral/main/reference/derive_var_merged_exist_flag.md),
  [`derive_var_merged_summary()`](https:/pharmaverse.github.io/admiral/main/reference/derive_var_merged_summary.md),
  and
  [`derive_vars_merged_lookup()`](https:/pharmaverse.github.io/admiral/main/reference/derive_vars_merged_lookup.md)
  were updated to allow renaming in the argument `by_vars`.
  ([\#1680](https://github.com/pharmaverse/admiral/issues/1680))

- The units “min” and “sec” are added as valid values of `out_unit` in
  [`compute_duration()`](https:/pharmaverse.github.io/admiral/main/reference/compute_duration.md)
  and
  [`derive_vars_duration()`](https:/pharmaverse.github.io/admiral/main/reference/derive_vars_duration.md).
  ([\#1647](https://github.com/pharmaverse/admiral/issues/1647))

- The function
  [`derive_vars_query()`](https:/pharmaverse.github.io/admiral/main/reference/derive_vars_query.md)
  now includes a consistency check for `QUERY_SCOPE` and
  `QUERY_SCOPE_NUM` values.
  ([\#652](https://github.com/pharmaverse/admiral/issues/652))

- Argument `new_var` in `derive_param_extreme_event()` is made optional.
  ([\#1630](https://github.com/pharmaverse/admiral/issues/1630))

- `derive_vars_last_dose()` no longer fails if `USUBJID` is not included
  in the input dataset.
  ([\#1787](https://github.com/pharmaverse/admiral/issues/1787))

### Breaking Changes

- All function arguments which expected a list of quosures created by
  [`vars()`](https://dplyr.tidyverse.org/reference/vars.html) are now
  expecting a list of expressions created by
  [`exprs()`](https:/pharmaverse.github.io/admiral/main/reference/reexport-exprs.md).
  For example, instead of `by_vars = vars(STUDYID, USUBJID)`
  `by_vars = exprs(STUDYID, USUBJID)` must be used now.

  To enable running old scripts using
  [`vars()`](https://dplyr.tidyverse.org/reference/vars.html) in the
  admiral function calls admiral redefines the
  [`vars()`](https://dplyr.tidyverse.org/reference/vars.html) function
  such that it returns a list of expressions. This can be disabled by
  the admiral option `force_admiral_vars` (see
  [`set_admiral_options()`](https:/pharmaverse.github.io/admiral/main/reference/set_admiral_options.md)).
  Please note that this is a temporary solution and will be removed in a
  future admiral release.
  ([\#1627](https://github.com/pharmaverse/admiral/issues/1627))

- Function
  [`derive_param_tte()`](https:/pharmaverse.github.io/admiral/main/reference/derive_param_tte.md)
  has been updated such that only observations are added for subjects
  who have both an event or censoring and an observation in
  `dataset_adsl`.
  ([\#1576](https://github.com/pharmaverse/admiral/issues/1576))

- Function `derive_var_disposition_status()` has been deprecated, please
  use `derive_var_merged_cat()` instead.
  ([\#1681](https://github.com/pharmaverse/admiral/issues/1681))

- Function `derive_var_worst_flag()` has been deprecated, in favor of
  [`slice_derivation()`](https:/pharmaverse.github.io/admiral/main/reference/slice_derivation.md)/[`derive_var_extreme_flag()`](https:/pharmaverse.github.io/admiral/main/reference/derive_var_extreme_flag.md).
  ([\#1682](https://github.com/pharmaverse/admiral/issues/1682))

- Function `derive_vars_disposition_reason()` has been deprecated, in
  favor of
  [`derive_vars_merged()`](https:/pharmaverse.github.io/admiral/main/reference/derive_vars_merged.md).
  ([\#1683](https://github.com/pharmaverse/admiral/issues/1683))

- The following functions have been deprecated from previous
  [admiral](https://pharmaverse.github.io/admiral/) versions using the
  next phase of the deprecation process:
  ([\#1712](https://github.com/pharmaverse/admiral/issues/1712))

  - `derive_derived_param()`
  - `derive_param_first_event()`
  - `derive_vars_merged_dt()`
  - `derive_vars_merged_dtm()`
  - `derive_var_agegr_ema()`
  - `derive_var_agegr_fda()`

- The following functions, which were deprecated in previous
  [admiral](https://pharmaverse.github.io/admiral/) versions, have been
  removed:
  ([\#1712](https://github.com/pharmaverse/admiral/issues/1712))

  - `derive_var_ady()`
  - `derive_var_aendy()`
  - `derive_var_astdy()`
  - `derive_var_atirel()`
  - `derive_vars_suppqual()`
  - `smq_select()`
  - `sdg_select()`

- The following parameters, which were deprecated in previous
  [admiral](https://pharmaverse.github.io/admiral/) versions, have been
  removed:
  ([\#1712](https://github.com/pharmaverse/admiral/issues/1712))

  - `meddra_version`, `whodd_version`, `get_smq_fun` and `get_sdg_fun`
    from the
    [`create_query_data()`](https:/pharmaverse.github.io/admiral/main/reference/create_query_data.md)
    function
  - `date_imputation`, `time_imputation` and `preserve` parameters from
    [`date_source()`](https:/pharmaverse.github.io/admiral/main/reference/date_source.md)
    function
  - `filter` parameter from
    [`derive_var_extreme_flag()`](https:/pharmaverse.github.io/admiral/main/reference/derive_var_extreme_flag.md)

- `ADLB` metadata data set called `atoxgr_criteria_ctcv5` updated to
  remove unit check for `HYPERURICEMIA` as grade criteria based on
  `ANRHI` only. This metadata holds criteria for lab grading based on
  [Common Terminology Criteria for Adverse Events (CTCAE)
  v5.0](https://dctd.cancer.gov/research/ctep-trials/for-sites/adverse-events#ctep-ctcae).
  ([\#1650](https://github.com/pharmaverse/admiral/issues/1650))

- Renamed `derive_var_confirmation_flag()` and `filter_confirmation()`
  to
  [`derive_var_joined_exist_flag()`](https:/pharmaverse.github.io/admiral/main/reference/derive_var_joined_exist_flag.md)
  and
  [`filter_joined()`](https:/pharmaverse.github.io/admiral/main/reference/filter_joined.md)
  respectively.
  ([\#1738](https://github.com/pharmaverse/admiral/issues/1738))

### Documentation

- New vignette “Creating a PK NCA ADaM (ADPC/ADNCA)”.
  ([\#1639](https://github.com/pharmaverse/admiral/issues/1639))

- New vignette “Hy’s Law Implementation”.
  ([\#1637](https://github.com/pharmaverse/admiral/issues/1637))

- New vignette “Creating Questionnaire ADaMs”.
  ([\#1715](https://github.com/pharmaverse/admiral/issues/1715))

- The expected value for the `derivation` argument of
  [`restrict_derivation()`](https:/pharmaverse.github.io/admiral/main/reference/restrict_derivation.md),
  [`slice_derivation()`](https:/pharmaverse.github.io/admiral/main/reference/slice_derivation.md),
  and
  [`call_derivation()`](https:/pharmaverse.github.io/admiral/main/reference/call_derivation.md)
  is described now.
  ([\#1698](https://github.com/pharmaverse/admiral/issues/1698))

- Removed authors from function documentation, as we will now only be
  tracking an overall list of authors for admiral.
  ([\#1673](https://github.com/pharmaverse/admiral/issues/1673))

- Added an imputation example for `create_single_source_dataset()` in
  function documentation.
  ([\#1408](https://github.com/pharmaverse/admiral/issues/1408),
  [\#1760](https://github.com/pharmaverse/admiral/issues/1760))

- Updates to examples for
  [`derive_var_age_years()`](https:/pharmaverse.github.io/admiral/main/reference/derive_var_age_years.md)
  and
  [`derive_vars_duration()`](https:/pharmaverse.github.io/admiral/main/reference/derive_vars_duration.md).
  ([\#1620](https://github.com/pharmaverse/admiral/issues/1620),
  [\#1634](https://github.com/pharmaverse/admiral/issues/1634))

- Increased the level of documentation for
  [`derive_var_age_years()`](https:/pharmaverse.github.io/admiral/main/reference/derive_var_age_years.md)
  to describe the data type of the newly created `new_var` column.
  ([\#970](https://github.com/pharmaverse/admiral/issues/970))

### Various

- Functions
  [`derive_vars_dtm()`](https:/pharmaverse.github.io/admiral/main/reference/derive_vars_dtm.md)
  and
  [`derive_vars_dt()`](https:/pharmaverse.github.io/admiral/main/reference/derive_vars_dt.md)
  had a bug pertaining to imputations associated with `NA` values that
  has now been fixed.
  ([\#1646](https://github.com/pharmaverse/admiral/issues/1646))

## admiral 0.9.1

CRAN release: 2022-12-23

- Implement changes to
  [`if_else()`](https://dplyr.tidyverse.org/reference/if_else.html) from
  the release of `dplyr` version 1.1.0, which affects
  [`derive_vars_dtm()`](https:/pharmaverse.github.io/admiral/main/reference/derive_vars_dtm.md)
  and `and compute_tmf()`.
  ([\#1641](https://github.com/pharmaverse/admiral/issues/1641))

## admiral 0.9.0

CRAN release: 2022-12-06

### New Features

- The new function
  [`derive_vars_joined()`](https:/pharmaverse.github.io/admiral/main/reference/derive_vars_joined.md)
  adds variables from an additional dataset. The selection of the
  observations can depend on variables from both datasets. This can be
  used for adding `AVISIT`, `AWLO`, `AWHI` based on time windows and
  `ADY` or deriving the lowest value (nadir) before the current
  observation.
  ([\#1448](https://github.com/pharmaverse/admiral/issues/1448))

- New function
  [`derive_var_trtemfl()`](https:/pharmaverse.github.io/admiral/main/reference/derive_var_trtemfl.md)
  for deriving treatment emergent flags.
  ([\#989](https://github.com/pharmaverse/admiral/issues/989))

- The new function
  [`chr2vars()`](https:/pharmaverse.github.io/admiral/main/reference/chr2vars.md)
  turns a character vector into a list of quosures.
  ([\#1448](https://github.com/pharmaverse/admiral/issues/1448))

- New function
  [`derive_var_relative_flag()`](https:/pharmaverse.github.io/admiral/main/reference/derive_var_relative_flag.md)
  for flagging observations before or after a condition is fulfilled.
  ([\#1453](https://github.com/pharmaverse/admiral/issues/1453))

- New functions
  [`get_admiral_option()`](https:/pharmaverse.github.io/admiral/main/reference/get_admiral_option.md)
  and
  [`set_admiral_options()`](https:/pharmaverse.github.io/admiral/main/reference/set_admiral_options.md)
  to allow more flexibility on common function inputs; e.g. like
  `subject_keys` to avoid several find and replace instances of
  `vars(STUDYID, USUBJID)`.
  ([\#1338](https://github.com/pharmaverse/admiral/issues/1338))

- The new function
  [`create_period_dataset()`](https:/pharmaverse.github.io/admiral/main/reference/create_period_dataset.md)
  for creating a reference dataset for subperiods, periods, or phases
  from the ADSL dataset was added. The reference dataset can be used to
  create subperiod, period, and phase variables in OCCDS and BDS
  datasets.
  ([\#1477](https://github.com/pharmaverse/admiral/issues/1477))

- The new function
  [`derive_vars_period()`](https:/pharmaverse.github.io/admiral/main/reference/derive_vars_period.md)
  adds subperiod, period, or phase variables to ADSL. The values for the
  new variables are provided by a period reference dataset.
  ([\#1477](https://github.com/pharmaverse/admiral/issues/1477))

- New function
  [`derive_var_merged_summary()`](https:/pharmaverse.github.io/admiral/main/reference/derive_var_merged_summary.md)
  adds a variable of summarized values to the input dataset.
  ([\#1564](https://github.com/pharmaverse/admiral/issues/1564))

- A [`print()`](https://rdrr.io/r/base/print.html) method was added for
  all S3 objects defined by admiral, e.g.,
  [`date_source()`](https:/pharmaverse.github.io/admiral/main/reference/date_source.md),
  [`dthcaus_source()`](https:/pharmaverse.github.io/admiral/main/reference/dthcaus_source.md),
  … ([\#858](https://github.com/pharmaverse/admiral/issues/858))

- New metadata data set called `atoxgr_criteria_ctcv5` which holds
  criteria for lab grading based on [Common Terminology Criteria for
  Adverse Events (CTCAE)
  v5.0](https://dctd.cancer.gov/research/ctep-trials/for-sites/adverse-events#ctep-ctcae).

- Removed the `{assertthat}` dependency in
  [admiral](https://pharmaverse.github.io/admiral/).
  ([\#1392](https://github.com/pharmaverse/admiral/issues/1392))

- Removed R Version 3.6 check in CI/CD workflows in favor of the three
  most recent versions: 4.0, 4.1 and 4.2.
  ([\#1556](https://github.com/pharmaverse/admiral/issues/1556))

- The new function
  [`derive_locf_records()`](https:/pharmaverse.github.io/admiral/main/reference/derive_locf_records.md)
  adds LOCF records as new observations. This can be used when the input
  dataset does not contain observations for missed visits/time points or
  when `AVAL` is `NA` for particular visits/time points.
  ([\#1316](https://github.com/pharmaverse/admiral/issues/1316))

- New function
  [`convert_na_to_blanks()`](https:/pharmaverse.github.io/admiral/main/reference/convert_na_to_blanks.md)
  to convert character `NA` to blanks.
  ([\#1624](https://github.com/pharmaverse/admiral/issues/1624))

### Updates of Existing Functions

- Function `derive_param_first_event()` has been replaced by a more
  generalized `derive_param_extreme_event()` function with new argument
  `mode` allowing for the selection of either the `"first"` or `"last"`
  event record according to the conditions provided. Also the `date_var`
  argument has been replaced with the `order` argument instead. In
  addition, three new arguments `new_var`, `true_value`, and
  `false_value` have been added to allow the user to choose what
  variable is used to indicate whether an event happened, and the values
  it is given.
  ([\#1317](https://github.com/pharmaverse/admiral/issues/1317),
  [\#1242](https://github.com/pharmaverse/admiral/issues/1242))

- Argument `ignore_time_for_ref_end_date` was added to
  [`derive_var_ontrtfl()`](https:/pharmaverse.github.io/admiral/main/reference/derive_var_ontrtfl.md),
  which controls if time is considered for the condition if `start_date`
  is after `ref_end_date` + `ref_end_window` days.
  ([\#989](https://github.com/pharmaverse/admiral/issues/989))

- [`derive_var_atoxgr_dir()`](https:/pharmaverse.github.io/admiral/main/reference/derive_var_atoxgr_dir.md)
  default value of `atoxgr_criteria_ctcv4` removed for parameter
  `meta_criteria`. Can now also choose `atoxgr_criteria_ctcv5` for
  parameter `meta_criteria`, to implement NCI-CTCAEv5 grading criteria.

- *Environment* objects were consolidated into a single
  `admiral_environment` object under `R/admiral__environment.R`.
  ([\#1572](https://github.com/pharmaverse/admiral/issues/1572))

- The default value of the `keep_source_vars` argument in
  [`create_single_dose_dataset()`](https:/pharmaverse.github.io/admiral/main/reference/create_single_dose_dataset.md)
  was updated such that it takes the values of the other arguments into
  account and the `start_datetime` and `end_datetime` arguments are
  optional now.
  ([\#1598](https://github.com/pharmaverse/admiral/issues/1598))

- Function
  [`create_query_data()`](https:/pharmaverse.github.io/admiral/main/reference/create_query_data.md)
  has been updated such that the dictionary version is stored in the
  output dataset.
  ([\#1337](https://github.com/pharmaverse/admiral/issues/1337))

### Breaking Changes

- Function `derive_param_first_event()` has been deprecated. Please use
  `derive_param_extreme_event()` with the `order` argument instead of
  the `date_var` argument.
  ([\#1317](https://github.com/pharmaverse/admiral/issues/1317))

- Functions `smq_select()` and `sdg_select()` have been deprecated and
  replaced with
  [`basket_select()`](https:/pharmaverse.github.io/admiral/main/reference/basket_select.md).
  In the
  [`create_query_data()`](https:/pharmaverse.github.io/admiral/main/reference/create_query_data.md)
  function, `meddra_version` and `whodd_version` argument has been
  replaced by `version` and `get_smq_fun` and `get_sdg_fun` argument by
  `get_terms_fun`.
  ([\#1597](https://github.com/pharmaverse/admiral/issues/1597))

### Documentation

- New vignette “Generic Functions”.
  ([\#734](https://github.com/pharmaverse/admiral/issues/734))
- New vignette “Visit and Period Variables”.
  ([\#1478](https://github.com/pharmaverse/admiral/issues/1478))

### Various

- Function
  [`derive_param_tte()`](https:/pharmaverse.github.io/admiral/main/reference/derive_param_tte.md)
  had a bug that set `ADT` to `NA` when `start_date` was missing, which
  has now been fixed.
  ([\#1540](https://github.com/pharmaverse/admiral/issues/1540))

- Function
  [`derive_vars_merged()`](https:/pharmaverse.github.io/admiral/main/reference/derive_vars_merged.md)
  had an improperly formatted error message which has been corrected.
  ([\#1473](https://github.com/pharmaverse/admiral/issues/1473))

- Templates now save datasets as `.rds` instead of `.rda`.
  ([\#1501](https://github.com/pharmaverse/admiral/issues/1501))

- Function
  [`create_single_dose_dataset()`](https:/pharmaverse.github.io/admiral/main/reference/create_single_dose_dataset.md)
  no longer fails if the input dataset contains observations with dose
  frequency `"ONCE"`.
  ([\#1375](https://github.com/pharmaverse/admiral/issues/1375))

## admiral 0.8.4

CRAN release: 2022-10-14

- Fixed a bug where a recent update to `{lifecylce}` caused several
  `admiral` tests to break
  ([\#1500](https://github.com/pharmaverse/admiral/issues/1500))

## admiral 0.8.3

CRAN release: 2022-10-07

- Second attempt to address issue where CRAN identified a failing test
  when “a strict Latin-1\* locale” is used
  ([\#1469](https://github.com/pharmaverse/admiral/issues/1469))
- Fixed a bug in
  [`derive_vars_duration()`](https:/pharmaverse.github.io/admiral/main/reference/derive_vars_duration.md)
  that surfaced after changes in R-devel
  ([\#1486](https://github.com/pharmaverse/admiral/issues/1486))

## admiral 0.8.2

CRAN release: 2022-09-29

- Fixed an issue where CRAN identified a failing test when “a strict
  Latin-1\* locale” is used
  ([\#1469](https://github.com/pharmaverse/admiral/issues/1469))

## admiral 0.8.1

CRAN release: 2022-09-20

- [`derive_var_extreme_dt()`](https:/pharmaverse.github.io/admiral/main/reference/derive_var_extreme_dt.md)
  and
  [`derive_var_extreme_dtm()`](https:/pharmaverse.github.io/admiral/main/reference/derive_var_extreme_dtm.md)
  were updated such that source observations where the date is `NA` are
  excluded
  ([\#1419](https://github.com/pharmaverse/admiral/issues/1419))

## admiral 0.8.0

CRAN release: 2022-09-05

### New Features

- [`get_summary_records()`](https:/pharmaverse.github.io/admiral/main/reference/get_summary_records.md)
  creates summary records e.g. derive analysis value (`AVAL`) from
  multiple records, only keeping the derived observations
  ([\#525](https://github.com/pharmaverse/admiral/issues/525))

- [`derive_param_framingham()`](https:/pharmaverse.github.io/admiral/main/reference/derive_param_framingham.md)
  adds a Parameter for Framingham Heart Study Cardiovascular Disease
  10-Year Risk Score
  ([\#977](https://github.com/pharmaverse/admiral/issues/977))

- [`compute_qual_imputation()`](https:/pharmaverse.github.io/admiral/main/reference/compute_qual_imputation.md)
  imputes values when qualifier exists in character result
  ([\#976](https://github.com/pharmaverse/admiral/issues/976))

- [`derive_vars_merged_lookup()`](https:/pharmaverse.github.io/admiral/main/reference/derive_vars_merged_lookup.md)
  maps lookup tables
  ([\#940](https://github.com/pharmaverse/admiral/issues/940))

- `filter_confirmation()` filters out confirmed observations
  ([\#1292](https://github.com/pharmaverse/admiral/issues/1292))
  including supporting functions
  [`count_vals()`](https:/pharmaverse.github.io/admiral/main/reference/count_vals.md),
  [`min_cond()`](https:/pharmaverse.github.io/admiral/main/reference/min_cond.md),
  and
  [`max_cond()`](https:/pharmaverse.github.io/admiral/main/reference/max_cond.md).

- `derive_var_confirmation_flag()` derives a flag which depends on other
  observations of the input dataset
  ([\#1293](https://github.com/pharmaverse/admiral/issues/1293))

- [`derive_var_atoxgr()`](https:/pharmaverse.github.io/admiral/main/reference/derive_var_atoxgr.md)
  derives lab toxicity/severity grade `ATOXGR` from `ATOXGRL` and
  `ATOXGRH`. `ATOXGRL` holds toxicity/severity grade for low lab values,
  and `ATOXGRH` holds toxicity/severity grade for high lab values.

- [`derive_var_atoxgr_dir()`](https:/pharmaverse.github.io/admiral/main/reference/derive_var_atoxgr_dir.md)
  derives lab toxicity/severity grade for low lab values (`ATOXGRL`) or
  for high lab values (`ATOXGRH`). The grading is created from metadata.

- New metadata data set called `atoxgr_criteria_ctcv4` which holds
  criteria for lab grading based on [Common Terminology Criteria for
  Adverse Events (CTCAE)
  v4.0](https://dctd.cancer.gov/research/ctep-trials/trial-development#ctcae-and-ctep-codes)

### Updates of Existing Functions

- [`list_tte_source_objects()`](https:/pharmaverse.github.io/admiral/main/reference/list_tte_source_objects.md)
  gains a `package` parameter and is now exported
  ([\#1212](https://github.com/pharmaverse/admiral/issues/1212))

- [`list_all_templates()`](https:/pharmaverse.github.io/admiral/main/reference/list_all_templates.md)
  and
  [`use_ad_template()`](https:/pharmaverse.github.io/admiral/main/reference/use_ad_template.md)
  gain a `package` parameter which can be used to indicate in which
  package to look for templates
  ([\#1205](https://github.com/pharmaverse/admiral/issues/1205))

- Randomization Date `RANDDT` variable added to ADSL template and
  vignette
  ([\#1126](https://github.com/pharmaverse/admiral/issues/1126))

- Renamed `derive_derived_param()` to
  [`derive_param_computed()`](https:/pharmaverse.github.io/admiral/main/reference/derive_param_computed.md)
  and added a deprecation notice
  ([\#1229](https://github.com/pharmaverse/admiral/issues/1229))

- [`derive_vars_duration()`](https:/pharmaverse.github.io/admiral/main/reference/derive_vars_duration.md)
  updated to not display units when there is missing duration
  ([\#1207](https://github.com/pharmaverse/admiral/issues/1207))

- `value_var` parameter added to
  [`derive_vars_atc()`](https:/pharmaverse.github.io/admiral/main/reference/derive_vars_atc.md)
  ([\#1120](https://github.com/pharmaverse/admiral/issues/1120))

- `format_eoxxstt_default()` - Updated the default value of EOSSTT for
  screen failure patients
  ([\#885](https://github.com/pharmaverse/admiral/issues/885))

- The imputation functions
  ([`derive_vars_dtm()`](https:/pharmaverse.github.io/admiral/main/reference/derive_vars_dtm.md),
  [`derive_vars_dt()`](https:/pharmaverse.github.io/admiral/main/reference/derive_vars_dt.md),
  [`convert_dtc_to_dtm()`](https:/pharmaverse.github.io/admiral/main/reference/convert_dtc_to_dtm.md),
  [`convert_dtc_to_dt()`](https:/pharmaverse.github.io/admiral/main/reference/convert_dtc_to_dt.md))
  have been enhanced to address users feedback
  ([\#1300](https://github.com/pharmaverse/admiral/issues/1300)):

  - Partial dates with missing components in the middle like
    `"2003-12-15T-:15:18"`, `"2003-12-15T13:-:19"`, `"2020-07--T00:00"`
    are handled now.

  - The control of the level of imputation has been refined by adding
    the `highest_imputation` argument. For example,
    `highest_imputation = "D"` requests imputation for day and time but
    not for year and month.

    (For the `date_imputation` and the `time_imputation` argument `NULL`
    is no longer a permitted value.)

  - It is now possible to impute completely missing dates by specifying
    `highest_imputation = "Y"` and the `min_dates` or `max_dates`
    argument.

- `order` parameter added to
  [`dthcaus_source()`](https:/pharmaverse.github.io/admiral/main/reference/dthcaus_source.md)
  which allows an additional character vector to be used for sorting the
  `dataset`, `derive_vars_dthcaus()` updated to process additional
  parameter
  ([\#1125](https://github.com/pharmaverse/admiral/issues/1125)).

- [`create_single_dose_dataset()`](https:/pharmaverse.github.io/admiral/main/reference/create_single_dose_dataset.md)
  Fixed bug where `ASTDTM` and `AENDTM` were not updated when
  `start_date = ASTDT` and `end_date = AENDT`. The function has been
  amended to now require `start_datetime` and `end_datetime` parameters
  in addition to `start_date` and `end_date`.The `keep_source_vars` has
  been added to specify the variables to be retained from the source
  dataset ([\#1224](https://github.com/pharmaverse/admiral/issues/1224))

### Breaking Changes

- Moved all developer-facing functions and vignettes to
  [admiraldev](https://pharmaverse.github.io/admiraldev/).
  [admiraldev](https://pharmaverse.github.io/admiraldev/) is now a
  dependency of [admiral](https://pharmaverse.github.io/admiral/)
  ([\#1231](https://github.com/pharmaverse/admiral/issues/1231))

- All ADaM datasets but `admiral_adsl` have been removed from the
  package ([\#1234](https://github.com/pharmaverse/admiral/issues/1234))

- `derive_var_agegr_ema()` and `derive_var_agegr_fda()` have been
  deprecated
  ([\#1333](https://github.com/pharmaverse/admiral/issues/1333))

- Imputation related arguments have been deprecated for all functions
  except the imputation functions themselves
  ([\#1299](https://github.com/pharmaverse/admiral/issues/1299)). I.e.,
  if a derivation like last known alive date is based on dates, DTC
  variables have to be converted to numeric date or datetime variables
  in a preprocessing step. For examples see the [ADSL
  vignette](https://pharmaverse.github.io/admiral/cran-release/articles/adsl.html).
  The following arguments were deprecated:

  - `date_imputation`, `time_imputation`, and `preserve` in
    [`date_source()`](https:/pharmaverse.github.io/admiral/main/reference/date_source.md)

  The following arguments no longer accept DTC variables:

  - `date` in
    [`date_source()`](https:/pharmaverse.github.io/admiral/main/reference/date_source.md),
    [`dthcaus_source()`](https:/pharmaverse.github.io/admiral/main/reference/dthcaus_source.md),
    [`censor_source()`](https:/pharmaverse.github.io/admiral/main/reference/censor_source.md),
    and
    [`event_source()`](https:/pharmaverse.github.io/admiral/main/reference/event_source.md)
  - `dose_date` and `analysis_date` in `derive_vars_last_dose()`,
    `derive_var_last_dose_amt()`, `derive_var_last_dose_date()`,
    `derive_var_last_dose_grp()`

  The following functions were deprecated:

  - `derive_vars_merged_dt()`
  - `derive_vars_merged_dtm()`

- For the `date_imputation` and the `time_imputation` argument of the
  imputation functions
  ([`derive_vars_dtm()`](https:/pharmaverse.github.io/admiral/main/reference/derive_vars_dtm.md),
  [`derive_vars_dt()`](https:/pharmaverse.github.io/admiral/main/reference/derive_vars_dt.md),
  [`convert_dtc_to_dtm()`](https:/pharmaverse.github.io/admiral/main/reference/convert_dtc_to_dtm.md),
  [`convert_dtc_to_dt()`](https:/pharmaverse.github.io/admiral/main/reference/convert_dtc_to_dt.md))
  `NULL` is no longer a permitted value. The level of imputation can be
  controlled by the `highest_imputation` argument now.

- The following functions, which were deprecated in previous {admiral}
  versions, have been removed:

  - `derive_var_disposition_dt()`
  - `derive_var_lstalvdt()`
  - `lstalvdt_source()`
  - `derive_var_trtedtm()`
  - `derive_var_trtsdtm()`

- The following functions and parameters, which were deprecated in
  previous {admiral} versions, are now defunct and will output an ERROR
  if used:

  - `derive_var_ady()`
  - `derive_var_aendy()`
  - `derive_var_astdy()`
  - `derive_var_atirel()`
  - `filter` parameter in
    [`derive_var_extreme_flag()`](https:/pharmaverse.github.io/admiral/main/reference/derive_var_extreme_flag.md)
    and `derive_var_worst_flag()`

### Documentation

- New ADMH template script can be accessed using
  `admiral::use_ad_template("admh")`
  ([\#502](https://github.com/pharmaverse/admiral/issues/502))

- New vignette “Higher Order Functions”
  ([\#1047](https://github.com/pharmaverse/admiral/issues/1047))

- New vignette “Lab Grading”
  ([\#1369](https://github.com/pharmaverse/admiral/issues/1369))

- Fixed `derive_var_disposition_status()` argument to render correctly
  ([\#1268](https://github.com/pharmaverse/admiral/issues/1268))

- Added link to [pharmaverse YouTube
  channel](https://www.youtube.com/channel/UCxQFEv8HNqM01DXzdQLCy6Q) to
  README

### Various

- Restructured Reference page and updated **all** functions to use
  `family` tag in roxygen headers for finding similar functions.
  ([\#1105](https://github.com/pharmaverse/admiral/issues/1105))

- Rename “Articles” page on website to “User Guides” and moved developer
  vignettes to [admiraldev](https://pharmaverse.github.io/admiraldev/)
  website ([\#1356](https://github.com/pharmaverse/admiral/issues/1356))

## admiral 0.7.1

CRAN release: 2022-07-18

- `derive_vars_last_dose()` no longer fails when a variable renamed in
  `new_vars` is supplied to the `dose_date` parameter
  ([\#1206](https://github.com/pharmaverse/admiral/issues/1206))

- [`derive_vars_duration()`](https:/pharmaverse.github.io/admiral/main/reference/derive_vars_duration.md)
  updated to not display units when there is missing duration
  ([\#1207](https://github.com/pharmaverse/admiral/issues/1207))

- `derive_param_first_event()` was updated
  ([\#1214](https://github.com/pharmaverse/admiral/issues/1214)) such
  that

  - `AVAL` is derived instead of `AVALN` and
  - all variables from the source dataset are kept.

- [`create_single_dose_dataset()`](https:/pharmaverse.github.io/admiral/main/reference/create_single_dose_dataset.md)
  Fixed bug where ASTDTM and AENDTM were not updated when
  `start_date=ASTDT` and `end_date=AENDT`. The function has been amended
  to now require start_datetime and end_datetime parameters in addition
  to start_date and end_date.The keep_source_vars has been added to
  specify the variables to be retained from the source dataset.

- [`slice_derivation()`](https:/pharmaverse.github.io/admiral/main/reference/slice_derivation.md)
  was updated such that it no longer fails if a slice is empty
  ([\#1309](https://github.com/pharmaverse/admiral/issues/1309))

## admiral 0.7.0

CRAN release: 2022-05-31

### New Features

- Updates to date/time imputation functions
  ([\#761](https://github.com/pharmaverse/admiral/issues/761)):

  - [`convert_date_to_dtm()`](https:/pharmaverse.github.io/admiral/main/reference/convert_date_to_dtm.md)
    and
    [`convert_dtc_to_dtm()`](https:/pharmaverse.github.io/admiral/main/reference/convert_dtc_to_dtm.md)
    now have time_imputation = “00:00:00” as default

  - [`derive_vars_dt()`](https:/pharmaverse.github.io/admiral/main/reference/derive_vars_dt.md)now
    has flag_imputation = “auto” as default

- New functions for merging variables
  ([\#607](https://github.com/pharmaverse/admiral/issues/607)):

  - [`derive_vars_merged()`](https:/pharmaverse.github.io/admiral/main/reference/derive_vars_merged.md) -
    Merge Variables from a Dataset to the Input Dataset
  - `derive_vars_merged_dt()` - Merge a (Imputed) Date Variable
  - `derive_vars_merged_dtm()` - Merge a (Imputed) Datetime Variable
  - `derive_var_merged_cat()` - Merge a Categorization Variable
  - [`derive_var_merged_exist_flag()`](https:/pharmaverse.github.io/admiral/main/reference/derive_var_merged_exist_flag.md) -
    Merge an Existence Flag
  - `derive_var_merged_character()` - Merge a Character Variable

- [`create_query_data()`](https:/pharmaverse.github.io/admiral/main/reference/create_query_data.md)
  is provided to create the [queries
  dataset](https://pharmaverse.github.io/admiral/cran-release/articles/queries_dataset.html)
  required as input for
  [`derive_vars_query()`](https:/pharmaverse.github.io/admiral/main/reference/derive_vars_query.md)
  ([\#606](https://github.com/pharmaverse/admiral/issues/606))

- [`create_single_dose_dataset()`](https:/pharmaverse.github.io/admiral/main/reference/create_single_dose_dataset.md) -
  Derives dataset of single dose from aggregate dose information
  ([\#660](https://github.com/pharmaverse/admiral/issues/660))

- New functions for deriving first or last dates from multiple source
  datasets ([\#753](https://github.com/pharmaverse/admiral/issues/753)):

  - [`derive_var_extreme_dtm()`](https:/pharmaverse.github.io/admiral/main/reference/derive_var_extreme_dtm.md) -
    Derive First or Last Datetime from Multiple Sources
  - [`derive_var_extreme_dt()`](https:/pharmaverse.github.io/admiral/main/reference/derive_var_extreme_dt.md) -
    Derive First or Last Date from Multiple Sources

- New function
  [`derive_extreme_records()`](https:/pharmaverse.github.io/admiral/main/reference/derive_extreme_records.md)
  for adding the first or last observation within each by group to the
  dataset ([\#1042](https://github.com/pharmaverse/admiral/issues/1042))

- New function `derive_param_first_event()`: Add a new parameter for the
  first event occurring in a dataset.
  ([\#1063](https://github.com/pharmaverse/admiral/issues/1063))

- New function
  [`derive_param_exist_flag()`](https:/pharmaverse.github.io/admiral/main/reference/derive_param_exist_flag.md):
  Add a new parameter indicating that a certain event exists in a
  dataset.
  ([\#1064](https://github.com/pharmaverse/admiral/issues/1064))

- New high order functions
  ([\#701](https://github.com/pharmaverse/admiral/issues/701)):

  - [`restrict_derivation()`](https:/pharmaverse.github.io/admiral/main/reference/restrict_derivation.md) -
    Execute a derivation on a subset of the input dataset
  - [`slice_derivation()`](https:/pharmaverse.github.io/admiral/main/reference/slice_derivation.md) -
    The input dataset is split into slices (subsets) and for each slice
    a derivation is called separately. Some or all arguments of the
    derivation may vary depending on the slice.

- [`filter_relative()`](https:/pharmaverse.github.io/admiral/main/reference/filter_relative.md) -
  Selects observations before or after the observation where a specified
  condition is fulfilled. For example, all observations up to first
  disease progression.
  ([\#1023](https://github.com/pharmaverse/admiral/issues/1023))

#### ADLB

- New ADLB template script available `ad_adlb.R`, specific ADLB
  functions developed and [BDS Finding
  vignette](https://pharmaverse.github.io/admiral/cran-release/articles/bds_finding.html)
  has examples enhanced with ADLB functions.
  ([\#1122](https://github.com/pharmaverse/admiral/issues/1122))

- [`derive_var_shift()`](https:/pharmaverse.github.io/admiral/main/reference/derive_var_shift.md) -
  Derives a character shift variable containing concatenated shift in
  values based on user-defined pairing
  ([\#944](https://github.com/pharmaverse/admiral/issues/944))

- [`derive_var_analysis_ratio()`](https:/pharmaverse.github.io/admiral/main/reference/derive_var_analysis_ratio.md) -
  Derives a ratio variable based on user-supplied variables from a BDS
  dataset, e.g. ADLB.
  ([\#943](https://github.com/pharmaverse/admiral/issues/943))

- [`derive_param_wbc_abs()`](https:/pharmaverse.github.io/admiral/main/reference/derive_param_wbc_abs.md) -
  Adds a parameter for lab differentials converted to absolute values.
  ([\#941](https://github.com/pharmaverse/admiral/issues/941))

#### ADPP

- New ADPP template script available `ad_adpp.R` which creates
  Pharmacokinetics Parameters Analysis Dataset
  ([\#850](https://github.com/pharmaverse/admiral/issues/850))

### Updates of Existing Functions

- Datasets internal to the package have been renamed with prefix
  `admiral_`, e.g. `adsl` has been renamed to `admiral_adsl`.
  Corresponding SDTM datasets in `{admiral.test}` have also been
  renamed, e.g.`dm` to `admiral_dm`. These changes will impact examples,
  vignettes, unit tests and templates
  ([\#1108](https://github.com/pharmaverse/admiral/issues/1108) and
  [\#1088](https://github.com/pharmaverse/admiral/issues/1088))

- When
  [`derive_vars_dtm_to_tm()`](https:/pharmaverse.github.io/admiral/main/reference/derive_vars_dtm_to_tm.md)
  was called for variables created by
  [`derive_vars_dtm()`](https:/pharmaverse.github.io/admiral/main/reference/derive_vars_dtm.md)
  the function failed. This bug was fixed
  ([\#1097](https://github.com/pharmaverse/admiral/issues/1097)).

- `impute_dtc()` - Fixed imputation bug. A user setting
  `date_imputation = MID` and `preserve = FALSE` would expect the date
  `2019---07` to be imputed to `2019-06-30`, but the function was
  returning `2019-06-15`. Now returns it correctly. This bug fix also
  addresses the issue in the downstream functions
  [`derive_vars_dt()`](https:/pharmaverse.github.io/admiral/main/reference/derive_vars_dt.md)
  and
  [`derive_vars_dtm()`](https:/pharmaverse.github.io/admiral/main/reference/derive_vars_dtm.md).
  ([\#1081](https://github.com/pharmaverse/admiral/issues/1081))

- `format_eoxxstt_default()` - Updated to have a more meaningful
  parameter name i.e. the parameter that was x is now status
  ([\#911](https://github.com/pharmaverse/admiral/issues/911))

### Breaking Changes

- `derive_var_lstalvdt()` has been deprecated in favor of
  [`derive_var_extreme_dt()`](https:/pharmaverse.github.io/admiral/main/reference/derive_var_extreme_dt.md)
  ([\#753](https://github.com/pharmaverse/admiral/issues/753)).

- `derive_vars_disposition_reason()` now is updated such that the
  default is populating `DCSREASP` only when `DSDECOD` is equal to
  `'OTHER'`, which is consistent with ADaMIG_v1.3
  ([\#886](https://github.com/pharmaverse/admiral/issues/886)).

- `derive_vars_suppqual()` has been removed from {admiral} as adding
  supplementary qualifiers is now done in another package called
  [{metatools}](https://github.com/pharmaverse/metatools) in a function
  called `combine_supp()` and is available on CRAN
  ([\#950](https://github.com/pharmaverse/admiral/issues/950))

- The `filter` parameter in
  [`derive_var_extreme_flag()`](https:/pharmaverse.github.io/admiral/main/reference/derive_var_extreme_flag.md)
  and `derive_var_worst_flag()` has been deprecated in favor of
  [`restrict_derivation()`](https:/pharmaverse.github.io/admiral/main/reference/restrict_derivation.md)
  ([\#701](https://github.com/pharmaverse/admiral/issues/701)).

- The following functions and parameters, which were deprecated in
  previous {admiral} versions, have been removed
  ([\#1056](https://github.com/pharmaverse/admiral/issues/1056)):

  - `derive_agegr_ema()`
  - `derive_agegr_fda()`
  - `derive_disposition_dt()`
  - `derive_disposition_status()`
  - `derive_extreme_flag()`
  - `derive_worst_flag()`
  - `derive_obs_number()`
  - `derive_disposition_reason()`
  - `derive_var_basec()`
  - `derive_baseline()`
  - `derive_params_exposure()`
  - `derive_last_dose()`
  - `dataset` parameter in `lstalvdt_source` and `dthcaus_source`

- The following functions were deprecated in favor of
  [`derive_vars_dy()`](https:/pharmaverse.github.io/admiral/main/reference/derive_vars_dy.md)
  ([\#1076](https://github.com/pharmaverse/admiral/issues/1076)):

  - `derive_var_ady()` - Derive Analysis Study Day
  - `derive_var_aendy()` - Derive Analysis End Relative Day
  - `derive_var_astdy()` - Derive Analysis Start Relative Day

- The following functions were deprecated in favor of
  `derive_vars_merged_dtm()`
  ([\#1076](https://github.com/pharmaverse/admiral/issues/1076)):

  - `derive_var_trtedtm()` - Derive Datetime of Last Exposure to
    Treatment
  - `derive_var_trtsdtm()` - Derive Datetime of First Exposure to
    Treatment

- The `derive_var_disposition_dt()` function was deprecated in favor of
  `derive_vars_merged_dt()`
  ([\#1076](https://github.com/pharmaverse/admiral/issues/1076))

- The `derive_var_atirel()` function was deprecated, as it is deemed as
  too specific for admiral. Derivations like this can be implemented
  calling
  [`mutate()`](https://dplyr.tidyverse.org/reference/mutate.html) and
  [`case_when()`](https://dplyr.tidyverse.org/reference/case_when.html).

### Documentation

- Additional explanation added to `derive_param_*` and
  `derive_derived_param` functions regarding which variables are
  populated in the additional rows
  ([\#939](https://github.com/pharmaverse/admiral/issues/939))

- Updated `derive_var_worst_flag()` and
  [`derive_var_extreme_flag()`](https:/pharmaverse.github.io/admiral/main/reference/derive_var_extreme_flag.md)
  vignettes to clarify their purpose
  ([\#691](https://github.com/pharmaverse/admiral/issues/691))

- Added example of ASEQ derivation in ADCM to [OCCDS
  vignette](https://pharmaverse.github.io/admiral/cran-release/articles/occds.html#aseq)
  ([\#720](https://github.com/pharmaverse/admiral/issues/720))

- Examples have been added for `format_reason_default()`,
  `format_eoxxstt_default()`, `extend_source_datasets()` and
  [`filter_date_sources()`](https:/pharmaverse.github.io/admiral/main/reference/filter_date_sources.md)
  ([\#745](https://github.com/pharmaverse/admiral/issues/745))

### Various

- Naming convention of admiral.xxx packages change to admiralxxx from
  this point onwards
  ([\#968](https://github.com/pharmaverse/admiral/issues/968))

## admiral 0.6.3

CRAN release: 2022-02-17

Address [CRAN
comments](https://github.com/pharmaverse/admiral/issues/946) raised
after submitting v0.6.2
([\#946](https://github.com/pharmaverse/admiral/issues/946))

## admiral 0.6.2

Address [CRAN
comments](https://github.com/pharmaverse/admiral/issues/925) raised
after submitting v0.6.1
([\#925](https://github.com/pharmaverse/admiral/issues/925))

## admiral 0.6.1

Address [CRAN
comments](https://github.com/pharmaverse/admiral/issues/918) raised
after submitting v0.6.0
([\#918](https://github.com/pharmaverse/admiral/issues/918))

## admiral 0.6.0

### New Features

- [`derive_vars_dy()`](https:/pharmaverse.github.io/admiral/main/reference/derive_vars_dy.md)
  derives the analysis day from one or more `--DT(M)` variables
  ([\#700](https://github.com/pharmaverse/admiral/issues/700))

### Updates of Existing Functions

- The `derive_last_dose()` function has been split into a general
  function `derive_vars_last_dose()` and three wrapper functions
  `derive_var_last_dose_amt()`, `derive_var_last_dose_date()`, and
  `derive_var_last_dose_grp()`
  ([\#385](https://github.com/pharmaverse/admiral/issues/385))

- [`derive_var_ontrtfl()`](https:/pharmaverse.github.io/admiral/main/reference/derive_var_ontrtfl.md)
  now has a `new_var` parameter to support the derivation of `ONTRxxFL`
  and `ONTRTwFL` variables
  ([\#721](https://github.com/pharmaverse/admiral/issues/721))

- [`derive_vars_dtm()`](https:/pharmaverse.github.io/admiral/main/reference/derive_vars_dtm.md),
  `derive_var_disposition` and `derive_var_lstalvdt` now have `preserve`
  argument. A user can preserve partial dates when doing date
  imputation, e.g. `2019---07` would become `2019-06-07` by setting
  `preserve` to `TRUE` when doing date_imputation
  ([\#592](https://github.com/pharmaverse/admiral/issues/592))

- [`derive_vars_dtm()`](https:/pharmaverse.github.io/admiral/main/reference/derive_vars_dtm.md)
  now has `ignore_seconds_flag` argument so users can suppress `"S"`
  flag if seconds are not recorded in the data
  ([\#589](https://github.com/pharmaverse/admiral/issues/589))

### Breaking Changes

- `derive_agegr_ema()`, `derive_agegr_fda()`, `derive_disposition_dt()`,
  `derive_disposition_status()`,`derive_extreme_flag()`,
  `derive_worst_flag()`, `derive_obs_number()`,
  `derive_disposition_reason()` have been deprecated and renamed in
  favor of `derive_var_agegr_ema()`, `derive_var_agegr_fda()`,
  `derive_var_disposition_dt()`, `derive_var_disposition_status()`,
  [`derive_var_extreme_flag()`](https:/pharmaverse.github.io/admiral/main/reference/derive_var_extreme_flag.md),
  `derive_var_worst_flag()`, `derive_var_last_dose()`,
  [`derive_var_obs_number()`](https:/pharmaverse.github.io/admiral/main/reference/derive_var_obs_number.md),
  and `derive_vars_disposition_reason()` respectively
  ([\#738](https://github.com/pharmaverse/admiral/issues/738))

- `derive_var_basec()` and `derive_baseline()` have been deprecated in
  favor of the extended
  [`derive_var_base()`](https:/pharmaverse.github.io/admiral/main/reference/derive_var_base.md)
  function ([\#695](https://github.com/pharmaverse/admiral/issues/695))

- `derive_params_exposure()` has been deprecated and renamed as
  [`derive_param_exposure()`](https:/pharmaverse.github.io/admiral/main/reference/derive_param_exposure.md)
  ([\#722](https://github.com/pharmaverse/admiral/issues/722))

- The `derive_last_dose()` function has been deprecated in favor of
  `derive_var_last_dose_date()`
  ([\#385](https://github.com/pharmaverse/admiral/issues/385))

- The behavior of all functions providing the `date_imputation`
  parameter, e.g.,
  [`derive_vars_dtm()`](https:/pharmaverse.github.io/admiral/main/reference/derive_vars_dtm.md)
  and
  [`derive_vars_dt()`](https:/pharmaverse.github.io/admiral/main/reference/derive_vars_dt.md),
  has changed for `date_imputation = "mid"`. Before the date was imputed
  as June 15th if both month and day were missing. Now it is imputed as
  June 30th. For the old behavior please specify
  `date_imputation = "06-15"`. Please note the behavior has not changed
  if only the day is missing. In this case the day is imputed as `15`
  ([\#592](https://github.com/pharmaverse/admiral/issues/592))

- [`derive_var_ontrtfl()`](https:/pharmaverse.github.io/admiral/main/reference/derive_var_ontrtfl.md)
  now has a `new_var` parameter to support the derivation of `ONTRxxFL`
  and `ONTRTwFL` variables
  ([\#721](https://github.com/pharmaverse/admiral/issues/721))

- The following functions and parameters, which were deprecated in
  previous {admiral} versions, were removed
  ([\#513](https://github.com/pharmaverse/admiral/issues/513)):

  - `derive_aage()`, `derive_duration()`, `derive_query_vars()`, and
    `derive_suppqual_vars()` function
  - `fns` and `filter_rows` parameters in
    [`derive_summary_records()`](https:/pharmaverse.github.io/admiral/main/reference/derive_summary_records.md)
  - `date_var` and `traceabilty_vars` parameters in
    [`dthcaus_source()`](https:/pharmaverse.github.io/admiral/main/reference/dthcaus_source.md)
  - `flag_filter` parameter in `derive_extreme_flag()`
  - `flag_filter` parameter in
    [`derive_var_extreme_flag()`](https:/pharmaverse.github.io/admiral/main/reference/derive_var_extreme_flag.md)
  - `date_var` parameter in `lstalvdt_source()`
  - `date` parameter in
    [`derive_var_ontrtfl()`](https:/pharmaverse.github.io/admiral/main/reference/derive_var_ontrtfl.md)

- `derive_var_agegr_fda()` has been updated to use ranges \<18, 18-64,
  \>=65 ([\#829](https://github.com/pharmaverse/admiral/issues/829))

### Documentation

- README and site homepage has been updated with important new section
  around expectations of {admiral}, as well as other useful references
  such as links to conference talks
  ([\#868](https://github.com/pharmaverse/admiral/issues/868) &
  [\#802](https://github.com/pharmaverse/admiral/issues/802))

- New vignette [Development
  Process](https://pharmaverse.github.io/admiral/cran-release/CONTRIBUTING.html)
  and improvements made to contribution vignettes
  ([\#765](https://github.com/pharmaverse/admiral/issues/765) &
  [\#758](https://github.com/pharmaverse/admiral/issues/758))

- Updated [Pull Request Review
  Guidance](https://pharmaverse.github.io/admiraldev/articles/pr_review_guidance.html)
  on using `task-list-completed` workflow
  ([\#817](https://github.com/pharmaverse/admiral/issues/817))

### Various

- GitHub repo moved to pharmaverse org and associated broken site links
  fixed ([\#803](https://github.com/pharmaverse/admiral/issues/803) &
  [\#820](https://github.com/pharmaverse/admiral/issues/820))

- Examples have been added for `format_reason_default`,
  `format_eoxxstt_default`, `extend_source_datasets` and
  `filter_date_sources`
  ([\#745](https://github.com/pharmaverse/admiral/issues/745))

## admiral 0.5.0

- The first truly open source release licensed under Apache 2.0
  ([\#680](https://github.com/pharmaverse/admiral/issues/680))

- New vignette [Contributing to
  admiral](https://pharmaverse.github.io/admiral/cran-release/CONTRIBUTING.html)
  ([\#679](https://github.com/pharmaverse/admiral/issues/679))

- New vignette [Unit Test
  Guidance](https://pharmaverse.github.io/admiraldev/articles/unit_test_guidance.html)
  ([\#679](https://github.com/pharmaverse/admiral/issues/679))

- Broken links in README have been fixed
  ([\#564](https://github.com/pharmaverse/admiral/issues/564))

## admiral 0.4.0

### New Features

#### General

- [`derive_vars_dtm_to_tm()`](https:/pharmaverse.github.io/admiral/main/reference/derive_vars_dtm_to_tm.md)
  enables the easy conversion of datetime to time variables
  ([\#551](https://github.com/pharmaverse/admiral/issues/551))

- [`derive_var_age_years()`](https:/pharmaverse.github.io/admiral/main/reference/derive_var_age_years.md)
  derives age in years from a variable providing the age in different
  units ([\#569](https://github.com/pharmaverse/admiral/issues/569))

#### BDS

- [`derive_param_tte()`](https:/pharmaverse.github.io/admiral/main/reference/derive_param_tte.md)
  derives time-to-event-parameters
  ([\#546](https://github.com/pharmaverse/admiral/issues/546))

- For common time-to-event endpoints [event and censoring source
  objects](https://pharmaverse.github.io/admiral/cran-release/reference/index.html#section-pre-defined-time-to-event-sources)
  are provided
  ([\#612](https://github.com/pharmaverse/admiral/issues/612))

#### Developer

- [`assert_list_element()`](https://pharmaverse.github.io/admiraldev/reference/assert_list_element.html)
  checks if an element of a list of lists/classes fulfills a condition

- [`assert_one_to_one()`](https://pharmaverse.github.io/admiraldev/reference/assert_one_to_one.html)
  checks if there is a one to one mapping between two lists of variables

- [`negate_vars()`](https:/pharmaverse.github.io/admiral/main/reference/negate_vars.md)
  negates a list of variables to remove them from a dataset with
  [`select()`](https://dplyr.tidyverse.org/reference/select.html)

### Updates of Existing Functions

- Unit checks in `derive_param_*()` functions are no longer case
  sensitive ([\#631](https://github.com/pharmaverse/admiral/issues/631))

- `derive_agegr_ema()` and `derive_agegr_fda()` gain a `age_unit`
  parameter used to specify the unit of the input age
  ([\#569](https://github.com/pharmaverse/admiral/issues/569))

### Breaking Changes

- All SDTM datasets have been moved to the {admiral.test} package
  ([\#639](https://github.com/pharmaverse/admiral/issues/639))

- The `min_dates` and `max_dates` parameters of
  [`derive_vars_dt()`](https:/pharmaverse.github.io/admiral/main/reference/derive_vars_dt.md)
  and
  [`derive_vars_dtm()`](https:/pharmaverse.github.io/admiral/main/reference/derive_vars_dtm.md)
  no longer expect a [`list()`](https://rdrr.io/r/base/list.html) but
  [`vars()`](https://dplyr.tidyverse.org/reference/vars.html) as input
  ([\#405](https://github.com/pharmaverse/admiral/issues/405))

### Bug Fixes

- [`derive_vars_dtm()`](https:/pharmaverse.github.io/admiral/main/reference/derive_vars_dtm.md)
  no longer shifts the time of the input `--DTC` variable
  ([\#436](https://github.com/pharmaverse/admiral/issues/436))

- [`derive_vars_dtm()`](https:/pharmaverse.github.io/admiral/main/reference/derive_vars_dtm.md)
  Change the min_dates with max_dates in the `lapply` statement when
  computing max_dates
  ([\#687](https://github.com/pharmaverse/admiral/issues/687))

### Documentation

- New vignette [Creating a BDS Time-to-Event
  ADaM](https://pharmaverse.github.io/admiral/cran-release/articles/bds_tte.html)
  ([\#549](https://github.com/pharmaverse/admiral/issues/549))

- New vignette [Queries Dataset
  Documentation](https://pharmaverse.github.io/admiral/cran-release/articles/queries_dataset.html)
  ([\#561](https://github.com/pharmaverse/admiral/issues/561))

- New vignette [Writing
  Vignettes](https://pharmaverse.github.io/admiraldev/articles/writing_vignettes.html)
  ([\#334](https://github.com/pharmaverse/admiral/issues/334))

- New vignette [Pull Request Review
  Guidance](https://pharmaverse.github.io/admiraldev/articles/pr_review_guidance.html)
  ([\#554](https://github.com/pharmaverse/admiral/issues/554))

- A section on handling missing values when working with {admiral} has
  been added to the “Get Started” vignette
  ([\#577](https://github.com/pharmaverse/admiral/issues/577))

- Package installation instructions have been added to the README
  ([\#558](https://github.com/pharmaverse/admiral/issues/558))

- The documentation of
  [`derive_vars_dtm()`](https:/pharmaverse.github.io/admiral/main/reference/derive_vars_dtm.md)
  falsely stated that the `flag_imputation` parameter should be either
  `TRUE` or `FALSE`. It now correctly states that the possible values
  are `"time"`, `"date"` or `"auto"`
  ([\#539](https://github.com/pharmaverse/admiral/issues/539))

## admiral 0.3.0

### New Features

#### General

- [`convert_blanks_to_na()`](https:/pharmaverse.github.io/admiral/main/reference/convert_blanks_to_na.md)
  can be used to convert SAS blanks, i.e. `""`, into proper R `NA`
  values ([\#482](https://github.com/pharmaverse/admiral/issues/482))

- [`call_derivation()`](https:/pharmaverse.github.io/admiral/main/reference/call_derivation.md)
  enables users to call the same function multiple times with some
  parameters being fixed across iterations and others varying
  ([\#403](https://github.com/pharmaverse/admiral/issues/403))

- [`derive_vars_dtm_to_dt()`](https:/pharmaverse.github.io/admiral/main/reference/derive_vars_dtm_to_dt.md)
  enables the easy conversion of datetime to date variables
  ([\#376](https://github.com/pharmaverse/admiral/issues/376))

- [`derive_var_ontrtfl()`](https:/pharmaverse.github.io/admiral/main/reference/derive_var_ontrtfl.md)
  can now handle events with a start and end date rather than just a
  single assessment date
  ([\#395](https://github.com/pharmaverse/admiral/issues/395))

- `derive_worst_flag()` enables the flagging of worst records
  ([\#300](https://github.com/pharmaverse/admiral/issues/300))

#### BDS

- `derive_derived_param()` can be used to derive a new parameter from
  existing parameters in a BDS dataset
  ([\#325](https://github.com/pharmaverse/admiral/issues/325))

- [`derive_param_bmi()`](https:/pharmaverse.github.io/admiral/main/reference/derive_param_bmi.md),
  [`derive_param_bsa()`](https:/pharmaverse.github.io/admiral/main/reference/derive_param_bsa.md)
  and
  [`derive_param_map()`](https:/pharmaverse.github.io/admiral/main/reference/derive_param_map.md)
  enables the derivation of the body mass index, body surface area and
  mean arterial pressure parameters respectively
  ([\#368](https://github.com/pharmaverse/admiral/issues/368))

- [`derive_param_qtc()`](https:/pharmaverse.github.io/admiral/main/reference/derive_param_qtc.md)
  enables the derivation of corrected QT intervals according to the
  formula of Bazett, Fridericia or Sagie
  ([\#325](https://github.com/pharmaverse/admiral/issues/325))

- [`derive_param_rr()`](https:/pharmaverse.github.io/admiral/main/reference/derive_param_rr.md)
  enables the derivation of the RR interval
  ([\#325](https://github.com/pharmaverse/admiral/issues/325))

- `derive_params_exposure()` enables the derivation of summary exposure
  parameters
  ([\#400](https://github.com/pharmaverse/admiral/issues/400))

- [`derive_param_doseint()`](https:/pharmaverse.github.io/admiral/main/reference/derive_param_doseint.md)
  enables the derivation of dose intensity
  ([\#179](https://github.com/pharmaverse/admiral/issues/179))

#### OCCDS

- `derive_var_atirel()` enables the derivation of the “Analysis Time
  Relative to Reference”
  ([\#397](https://github.com/pharmaverse/admiral/issues/397))

- [`derive_vars_atc()`](https:/pharmaverse.github.io/admiral/main/reference/derive_vars_atc.md)
  can be used to add ATC variables from FACM to ADCM
  ([\#396](https://github.com/pharmaverse/admiral/issues/396))

### Updates of Existing Functions

- [`derive_var_anrind()`](https:/pharmaverse.github.io/admiral/main/reference/derive_var_anrind.md)
  now checks whether the `AVAL` variable is present in the input dataset
  ([\#486](https://github.com/pharmaverse/admiral/issues/486))

- All derivation functions check whether the input dataset is grouped
  and throw an error if it is
  ([\#408](https://github.com/pharmaverse/admiral/issues/408))

- [`use_ad_template()`](https:/pharmaverse.github.io/admiral/main/reference/use_ad_template.md)
  has been refactored to no longer make use of the {usethis} package
  which is no longer a dependency of {admiral}
  ([\#433](https://github.com/pharmaverse/admiral/issues/433))

- A performance issue in
  [`derive_vars_dt()`](https:/pharmaverse.github.io/admiral/main/reference/derive_vars_dt.md)
  has been resolved
  ([\#384](https://github.com/pharmaverse/admiral/issues/384))

### Breaking Changes

- The `drop_values_from` parameter has been removed from
  [`derive_summary_records()`](https:/pharmaverse.github.io/admiral/main/reference/derive_summary_records.md)
  ([\#425](https://github.com/pharmaverse/admiral/issues/425))

- The format of the `date_imputation` parameter of
  [`derive_vars_dt()`](https:/pharmaverse.github.io/admiral/main/reference/derive_vars_dt.md)
  and
  [`derive_vars_dtm()`](https:/pharmaverse.github.io/admiral/main/reference/derive_vars_dtm.md)
  has been changed from “dd-mm” to “mm-dd”. Thus, “01-12” now refers to
  January 12th rather than December 1st
  ([\#492](https://github.com/pharmaverse/admiral/issues/492))

- Several functions have been renamed. The old names are now deprecated.
  They can still be used but a warning will be issued
  ([\#507](https://github.com/pharmaverse/admiral/issues/507))

  - `derive_aage()` -\>
    [`derive_vars_aage()`](https:/pharmaverse.github.io/admiral/main/reference/derive_vars_aage.md)
  - `derive_duration()` -\>
    [`derive_vars_duration()`](https:/pharmaverse.github.io/admiral/main/reference/derive_vars_duration.md)
  - `derive_query_vars()` -\>
    [`derive_vars_query()`](https:/pharmaverse.github.io/admiral/main/reference/derive_vars_query.md)
  - `derive_suppqual_vars()` -\> `derive_vars_suppqual()`

- The `date_var` parameter of `lstalvdt_source()` has been renamed to
  `date`

- The `filter_rows` parameter of
  [`derive_summary_records()`](https:/pharmaverse.github.io/admiral/main/reference/derive_summary_records.md)
  has been renamed to `filter`. The `fns` parameter has been deprecated
  in favor of `analysis_var` and `summary_fun`
  ([\#491](https://github.com/pharmaverse/admiral/issues/491))

- The `date_var` and `traceabilty_vars` parameters of
  [`dthcaus_source()`](https:/pharmaverse.github.io/admiral/main/reference/dthcaus_source.md)
  have been renamed to `date` and `traceability_vars`, respectively
  ([\#493](https://github.com/pharmaverse/admiral/issues/493))

- The `flag_filter` parameter of `derive_extreme_flag()` has been
  renamed to `filter`
  ([\#487](https://github.com/pharmaverse/admiral/issues/487))

### Bug Fixes

- `derive_agegr_fda()` used to return `NA` for ages less than or
  equal 18. It now returns `<=18`
  ([\#441](https://github.com/pharmaverse/admiral/issues/441))

### Documentation

- New vignette on “Date and Time Imputation” has been created
  ([\#198](https://github.com/pharmaverse/admiral/issues/198))

- A “Guidance for git Usage” has been created
  ([\#266](https://github.com/pharmaverse/admiral/issues/266))

- “OCCDS” has been added as a new section on the reference page on the
  package website
  ([\#485](https://github.com/pharmaverse/admiral/issues/485))

- The Programming Strategy has been updated
  ([\#495](https://github.com/pharmaverse/admiral/issues/495))

- A search feature has been added to the package website
  ([\#438](https://github.com/pharmaverse/admiral/issues/438))

- New template scripts for ADEX
  ([\#181](https://github.com/pharmaverse/admiral/issues/181)), ADCM
  ([\#268](https://github.com/pharmaverse/admiral/issues/268)) and ADEG
  ([\#258](https://github.com/pharmaverse/admiral/issues/258)) have been
  created

- New vignette for programming ADEX has been created
  ([\#372](https://github.com/pharmaverse/admiral/issues/372))

- A section on how to create query variables (e.g. SMQs in ADAE) has
  been added to the Occurrence datasets vignette
  ([\#370](https://github.com/pharmaverse/admiral/issues/370))

- The BDS vignette has been updated to incorporate examples of ADVS and
  ADEG specific functions
  ([\#371](https://github.com/pharmaverse/admiral/issues/371))

## admiral 0.2.1

- Fixed a critical bug in
  [`use_ad_template()`](https:/pharmaverse.github.io/admiral/main/reference/use_ad_template.md)
  that prevented the function from being usable at all
  ([\#326](https://github.com/pharmaverse/admiral/issues/326))

## admiral 0.2.0

### New Features

#### General

- Function argument checks have been completely re-written to provide
  clearer error messages to users
  ([\#263](https://github.com/pharmaverse/admiral/issues/263),
  [\#288](https://github.com/pharmaverse/admiral/issues/288))

- SDTM `SUPP--` datasets can be merged onto their parent domain using
  `derive_suppqual_vars()`
  ([\#145](https://github.com/pharmaverse/admiral/issues/145))

- In case a derivation detects duplicate records after applying a
  `filter`, the dataset of duplicate records is made available to users
  via
  [`get_duplicates_dataset()`](https:/pharmaverse.github.io/admiral/main/reference/get_duplicates_dataset.md)
  ([\#202](https://github.com/pharmaverse/admiral/issues/202))

- [`derive_vars_dt()`](https:/pharmaverse.github.io/admiral/main/reference/derive_vars_dt.md)
  and
  [`derive_vars_dtm()`](https:/pharmaverse.github.io/admiral/main/reference/derive_vars_dtm.md)
  gain a `min_dates` and `max_dates` parameter which can be used to
  ensure that the imputed date(time) is not before the `min_dates` nor
  after the `max_dates`, e.g. avoid that `AENDT` is after the data cut
  date or `ASTDT` is before the first treatment date
  ([\#158](https://github.com/pharmaverse/admiral/issues/158))

- [`use_ad_template()`](https:/pharmaverse.github.io/admiral/main/reference/use_ad_template.md)
  can be used to open a template script for an ADaM dataset; all
  available templates can be displayed using
  [`list_all_templates()`](https:/pharmaverse.github.io/admiral/main/reference/list_all_templates.md)
  ([\#110](https://github.com/pharmaverse/admiral/issues/110))

#### ADSL

- EMA and FDA defined age groupings can be derived using
  `derive_agegr_ema()` and `derive_agegr_fda()`, respectively

- Disposition Status can be derived using `derive_disposition_status()`
  ([\#92](https://github.com/pharmaverse/admiral/issues/92))

- Disposition Reason can be derived using `derive_disposition_reason()`

- Disposition Dates can be derived using `derive_disposition_dt()`
  ([\#91](https://github.com/pharmaverse/admiral/issues/91))

- Date Last Known Alive can be derived using `derive_var_lstalvdt()`
  ([\#94](https://github.com/pharmaverse/admiral/issues/94))

- Cause of Death can be derived using
  [`derive_var_dthcaus()`](https:/pharmaverse.github.io/admiral/main/reference/derive_var_dthcaus.md)
  ([\#93](https://github.com/pharmaverse/admiral/issues/93))

#### BDS

- Summary records for BDS datasets, e.g. with `DTYPE == "AVERAGE"`, can
  be derived using
  [`derive_summary_records()`](https:/pharmaverse.github.io/admiral/main/reference/derive_summary_records.md)
  ([\#177](https://github.com/pharmaverse/admiral/issues/177))

#### OCCDS

- Last Dose Date(time) can be derived using `derive_last_dose()`

### Breaking Changes

- `derive_merged_vars()` has been removed from {admiral} in favor of
  smaller special purpose functions, e.g. `derive_disposition_status()`
  ([\#167](https://github.com/pharmaverse/admiral/issues/167))

- Function arguments no longer accept expressions created with
  [`expr()`](https://rlang.r-lib.org/reference/expr.html) or
  [`exprs()`](https:/pharmaverse.github.io/admiral/main/reference/reexport-exprs.md)
  as inputs; instead filter expressions can be passed “as is” and
  multiple variables have to be wrapped inside
  [`vars()`](https://dplyr.tidyverse.org/reference/vars.html)
  ([\#187](https://github.com/pharmaverse/admiral/issues/187))

  **Old:**

  ``` r
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

  ``` r
  derive_extreme_flag(
    vs,
    new_var = LASTFL,
    by_vars = vars(USUBJID, VSTESTCD, VISIT),
    order = vars(VSTPTNUM),
    mode = "last",
    flag_filter = VISIT != "BASELINE"
  )
  ```

- `read_dap_m3()` and
  [`initialize()`](https://rdrr.io/r/methods/new.html) have been
  migrated to {admiral.roche}
  ([\#272](https://github.com/pharmaverse/admiral/issues/272))

- The `start_date` and `end_date` parameters of `derive_var_ady()`,
  `derive_var_aendy()` and `derive_var_astdy()` have been renamed to
  `reference_date` and `date`, respectively
  ([\#121](https://github.com/pharmaverse/admiral/issues/121))

### Bug Fixes

- `derive_var_basetype()` no longer drops records which do not match any
  condition defined in the `basetype` argument
  ([\#226](https://github.com/pharmaverse/admiral/issues/226))

- Join warnings like “Column `USUBJID` has different attributes on LHS
  and RHS of join when using left_join()” are no longer displayed
  ([\#271](https://github.com/pharmaverse/admiral/issues/271))

- For datetimes with time imputed to “00:00:00” the time part is now
  displayed ([\#206](https://github.com/pharmaverse/admiral/issues/206))

### Documentation

- [Frequently Asked
  Questions](https://pharmaverse.github.io/admiral/cran-release/articles/faq.html)

- [Creating
  ADSL](https://pharmaverse.github.io/admiral/cran-release/articles/adsl.html)

- [Creating a BDS Finding
  ADaM](https://pharmaverse.github.io/admiral/cran-release/articles/bds_finding.html)

- [Creating an OCCDS
  ADaM](https://pharmaverse.github.io/admiral/cran-release/articles/occds.html)
