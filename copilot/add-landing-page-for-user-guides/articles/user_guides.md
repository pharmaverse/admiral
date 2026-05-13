# User Guides

The [admiral](https://pharmaverse.github.io/admiral/) User Guides are
grouped by topic so you can quickly find conceptual guidance, end-to-end
ADaM examples, and focused methodology references.

## General Guides

### [Programming Concepts and Conventions](https:/pharmaverse.github.io/admiral/copilot/add-landing-page-for-user-guides/articles/concepts_conventions.md)

Learn the core [admiral](https://pharmaverse.github.io/admiral/)
programming principles, including expression handling, missing-value
conventions, and common workflow expectations.

### [Generic Derivations](https:/pharmaverse.github.io/admiral/copilot/add-landing-page-for-user-guides/articles/generic.md)

Get an overview of reusable derivation patterns and functions for adding
variables, parameters, and records across many ADaM use cases.

### [Higher Order Functions](https:/pharmaverse.github.io/admiral/copilot/add-landing-page-for-user-guides/articles/higher_order.md)

See how to use higher-order helpers like
[`call_derivation()`](https:/pharmaverse.github.io/admiral/copilot/add-landing-page-for-user-guides/reference/call_derivation.md),
[`restrict_derivation()`](https:/pharmaverse.github.io/admiral/copilot/add-landing-page-for-user-guides/reference/restrict_derivation.md),
and
[`slice_derivation()`](https:/pharmaverse.github.io/admiral/copilot/add-landing-page-for-user-guides/reference/slice_derivation.md)
to reduce repetitive derivation code.

## ADaM Script Guides

### [Creating a Basic ADSL](https:/pharmaverse.github.io/admiral/copilot/add-landing-page-for-user-guides/articles/adsl.md)

Follow a practical step-by-step workflow to build ADSL from SDTM
sources, from treatment and disposition variables through final
labeling.

### [Creating an OCCDS ADaM](https:/pharmaverse.github.io/admiral/copilot/add-landing-page-for-user-guides/articles/occds.md)

Walk through a full OCCDS example (ADAE-focused) covering timing,
treatment-emergent logic, flags, sequence derivation, and final dataset
assembly.

### [Creating a BDS Findings ADaM](https:/pharmaverse.github.io/admiral/copilot/add-landing-page-for-user-guides/articles/bds_finding.md)

Build a BDS findings dataset (ADVS-style) with examples for parameters,
baseline/change logic, analysis flags, categorization, and record
creation.

### [Creating a BDS Exposure ADaM](https:/pharmaverse.github.io/admiral/copilot/add-landing-page-for-user-guides/articles/bds_exposure.md)

Explore example approaches for structuring exposure analysis data in BDS
format, including one-to-one records, summaries, and parameter mappings.

### [Creating a BDS Time-to-Event ADaM](https:/pharmaverse.github.io/admiral/copilot/add-landing-page-for-user-guides/articles/bds_tte.md)

Learn how to define events and censoring sources and derive core
time-to-event analysis variables and parameters.

### [Creating Questionnaire ADaMs](https:/pharmaverse.github.io/admiral/copilot/add-landing-page-for-user-guides/articles/questionnaires.md)

Review common questionnaire derivation patterns for original items,
transformed values, and score/scale parameters rather than a single
fixed workflow.

### [Creating a PK NCA or Population PK ADaM](https:/pharmaverse.github.io/admiral/copilot/add-landing-page-for-user-guides/articles/pk_adnca.md)

Work through examples for building PK analysis datasets (ADNCA/ADPC and
ADPPK), with emphasis on relative timing and dose-reference derivations.

### [Creating an Antidrug Antibody ADaM](https:/pharmaverse.github.io/admiral/copilot/add-landing-page-for-user-guides/articles/adab_web.md)

See an experimental ADAB workflow covering data staging, parameter
computation, visit/summary assembly, and final analysis-variable
derivations.

### [Hy’s Law Implementation](https:/pharmaverse.github.io/admiral/copilot/add-landing-page-for-user-guides/articles/hys_law.md)

Use ADLB-based examples to derive criteria flags and joined records that
support identifying potential Hy’s Law drug-induced liver injury events.

## Special Topics Guides

### [Date and Time Imputation](https:/pharmaverse.github.io/admiral/copilot/add-landing-page-for-user-guides/articles/imputation.md)

Understand [admiral](https://pharmaverse.github.io/admiral/) date/time
imputation rules and apply both vector-level and dataset-level
imputation workflows for ADaM timing variables.

### [Visit and Period Variables](https:/pharmaverse.github.io/admiral/copilot/add-landing-page-for-user-guides/articles/visits_periods.md)

Learn practical strategies and helper functions for deriving visit
windows plus period, subperiod, and phase variables in study-specific
designs.

### [Queries Dataset Documentation](https:/pharmaverse.github.io/admiral/copilot/add-landing-page-for-user-guides/articles/queries_dataset.md)

Reference the required structure and content of query metadata datasets
used by
[`derive_vars_query()`](https:/pharmaverse.github.io/admiral/copilot/add-landing-page-for-user-guides/reference/derive_vars_query.md)
for grouped-event derivations.

### [Lab Grading](https:/pharmaverse.github.io/admiral/copilot/add-landing-page-for-user-guides/articles/lab_grading.md)

Understand how [admiral](https://pharmaverse.github.io/admiral/)
implements lab toxicity grading from metadata (e.g., CTCAE and DAIDS)
and how to apply grading in ADLB workflows.

### [Estimands](https:/pharmaverse.github.io/admiral/copilot/add-landing-page-for-user-guides/articles/estimands_web.md)

Review example patterns for representing intercurrent events and
implementing estimand-oriented analysis logic in
[admiral](https://pharmaverse.github.io/admiral/) datasets.

## FAQ

### [FAQ](https:/pharmaverse.github.io/admiral/copilot/add-landing-page-for-user-guides/articles/faq.md)

Find quick answers to common questions about
[admiral](https://pharmaverse.github.io/admiral/) goals, scope,
validation expectations, and practical usage.
