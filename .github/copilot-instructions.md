# GitHub Copilot Instructions - admiraldev Guidelines

This file provides context to GitHub Copilot about admiraldev programming standards
and best practices. Copilot will automatically reference these guidelines when
providing code suggestions in this repository.

**Auto-generated:** 2025-12-09 21:56:30  
**admiraldev Version:** GitHub main branch  
**Source:** https://github.com/pharmaverse/admiraldev

<e2><9a><a0><ef><b8><8f> **DO NOT EDIT MANUALLY** - Run `source('.github/scripts/sync_admiraldev_docs.R')` to update

---

## Purpose

These guidelines ensure that code in this repository follows admiral ecosystem standards:

- <e2><9c><85> Consistent programming patterns across admiral packages
- <e2><9c><85> Proper function design and documentation
- <e2><9c><85> High-quality vignettes that help users
- <e2><9c><85> Comprehensive unit tests with good coverage
- <e2><9c><85> Code that is maintainable and readable

GitHub Copilot will use these guidelines to provide better suggestions that align
with admiral best practices.

---

## Table of Contents

1. [Programming Strategy](#-rogramming-trategy)
2. [Writing Vignettes](#-riting-ignettes)
3. [Unit Test Guidance](#-nit-est-uidance)

---

# Programming Strategy

**Description:** Core programming principles and strategies for admiral packages  
**Source File:** `programming_strategy.Rmd`  
**Source:** GitHub (pharmaverse/admiraldev)  
**URL:** https://raw.githubusercontent.com/pharmaverse/admiraldev/main/vignettes/programming_strategy.Rmd

---

---
title: "Programming Strategy"
output: 
  rmarkdown::html_vignette:
    toc: true
vignette: >
  %\VignetteIndexEntry{Programming Strategy}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r setup, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
```

# Introduction

As `{admiral}` is intended to be contributed by the user community, this 
article is meant for developers that want to either expand `{admiral}` functionalities or build on top of `{admiral}`. 
In order to keep the framework robust across the whole community, 
we have defined a programming strategy that should be followed in such cases.
These contributions could include, for example, company specific derivations of ADaM datasets. 


# Functional Workflow

* Overall programming will follow a functional approach.
* We mandate the use of tidyverse (e.g. dplyr) over similar functionality existing in base R
* Each ADaM dataset is built with a set of functions and not with free flow code.
* Each ADaM dataset has a specific programming workflow.
* Each function has a specific purpose that supports the ADaM Dataset programming workflow. It could be an `{admiral}` function or a company specific function.
* Admiral functions can be re-used for company specific functions.
* Each function belongs to one category defined in keywords/family.
* Each function that is used to derive one or multiple variable(s) is required to be unit tested.
* Functions have a standard naming convention.
* Double coding is not used as a QC method (only if absolutely necessary).
* ADaMs are created with readable, submission-ready code.

# Functions in R

## Function Design

Firstly, it is important to explain how we decide on the need for new derivation functions. 

If a derivation rule or algorithm is common and highly similar across different variables/parameters
(e.g. study day or duration) then we would provide a generic function that can be used to satisfy all
the times this may be needed across different ADaMs. Similarly, if we feel that a certain derivation
could be useful beyond a single purpose we also would provide a generic function (e.g. instead of a
last known alive date function, we have an extreme date function where a user could find the last date
from a selection, or for example the first).

Otherwise, if we feel that a derivation rule is a unique need or sufficiently complex to justify then we
opt for a dedicated function for that specific variable/parameter (e.g. treatment-emergent flag for AEs).

If certain variables are closely connected (e.g. an imputed date and the corresponding imputation flag)
then a single function would provide both variables.

If something needed for ADaM could be achieved simply via an existing tidyverse function, then we do not
wrap this into an admiral function, as that would add an unnecessary extra layer for users.

The following principles are key when designing a new function:

* _**Modularity**_ - All code follows a modular approach, i.e. the steps must be clearly separated and
have a dedicated purpose. This applies to scripts creating a dataset where each module should create a
single variable or parameter. But also to complex derivations with several steps. Commenting on these
steps is key for readability.

* _**Avoid Copy and Paste**_ - If the same or very similar code is used multiple times, it should be put
into a separate function. This improves readability and maintainability and makes unit testing easier.
This should not be done for every simple programming step where tidyverse can be used. But rather for
computational functions or data checks. However, also consider not to nest too many functions.

* _**Checks**_ - Whenever a function fails, a meaningful error message must be provided with a clear
reference to the input which caused the failure. A users should not have to dig into detailed
code if they only want to apply a function.  A meaningful error message supports usability.

* _**Flexibility**_ - Functions should be as flexible as possible as long as it does not reduce the usability.
For example:

  * The source variables or newly created variables and conditions for selecting observations should not be hard-coded.

  * It is useful if an argument triggers optional steps, e.g. if the `filter` argument is specified, the input dataset
is restricted, otherwise this step is skipped.

  * However, arguments should not trigger completely different algorithms. For example `BNRIND` could be derived based
on `BASE` or based on `ANRIND`. It should not be implemented within one function as the algorithms are completely different.
If `BASE` is used, the values are categorized while if `ANRIND` is used, the values are merged from the baseline observation.

## Input, Output, and Side-effects

* The behavior of the function is only determined by its input, not by any global object,  
i.e. all input like datasets, variable names, options, <U+2026> must be provided to the function by arguments.
* It is expected that the input datasets are not grouped. If any are grouped, the function must issue an error.
* If a function requires grouping, the function must provide the `by_vars` argument.
* The output dataset must be ungrouped.
* The functions should not sort (arrange) the output dataset at the end.
* If the function needs to create temporary variables in an input dataset, names
for these variables must be generated by `get_new_tmp_var()` to avoid that
variables of the input dataset are accidentally overwritten. The temporary
variables must be removed from the output dataset by calling
`remove_tmp_vars()`.
* If developers find the need to use or create *environment* objects to achieve flexibility, use the `admiral_environment` environment object created in `admiral_environment.R`. All objects which are stored in this environment must be documented in `admiral_environment.R`. An equivalent environment object and `.R` file exist for admiraldev as well. For more details how environments work, see relevant sections on environments in [R Packages](https://r-pkgs.org) and [Advanced R](https://adv-r.hadley.nz) textbooks.
* In general, the function must not have any side-effects like creating or modifying global objects, printing, writing files, ...

## Admiral Options

* An exception is made for admiral options, see `get_admiral_option()` and
`set_admiral_options()`, where we have certain pre-defined defaults with added
flexibility to allow for user-defined defaults on *commonly used* function
arguments e.g. `subject_keys` currently pre-defined as `exprs(STUDYID,
USUBJID)`, but can be modified using `set_admiral_options(subject_keys =
exprs(...))` at the top of a script. The reasoning behind this was to relieve
the user of repeatedly changing aforementioned *commonly used* function
arguments multiple times in a script, which may be called across many admiral
functions.
* If this additional flexibility needs to be added for another *commonly used*
function argument e.g. `future_input` to be set as `exprs(...)` it can be added
as an admiral option. In the function formals define `future_input =
get_admiral_option("future_input")` then proceed to modify the body and roxygen
documentation of `set_admiral_options()`.

## Function Names

* Function names should start with a verb and use snake case, e.g. `derive_var_base()`. 

| Function name prefix                         | Description                                                                                         |
|----------------------------------------------|-----------------------------------------------------------------------------------------------------|
| `assert_` / `warn_` / `is_`                  | Functions that check other functions' inputs                                                        |
| `derive_`                                    | Functions that take a dataset as input and return a new dataset with additional rows and/or columns |
| `derive_var_` (e.g. `derive_var_trtdurd`)    | Functions which add a single variable                                                               |
| `derive_vars_` (e.g. `derive_vars_dt`)       | Functions which add multiple variables                                                              |
| `derive_param_` (e.g. `derive_param_os`)     | Functions which add a single parameter                                                              |
| `compute_` /  `calculate_` / ...             | Functions that take vectors as input and return a vector                                            |
| `create_`  /  `consolidate_`                 | Functions that create datasets without keeping the original observations                            |
| `get_`                                       | Usually utility functions that return very specific objects that get passed through other functions |
| `filter_`                                    | Functions that filter observations based on conditions associated with common clinical trial syntax |

| Function Name Suffix                         | Description                                                                                         |
|----------------------------------------------|-----------------------------------------------------------------------------------------------------|
| `_derivation` (suffix)                       | High order functions that call a user specified derivation                                          |
| `_date` / `_time` / `_dt` / `_dtc` / `_dtm`  | Functions associated with dates, times, datetimes, and their character equivalents.                 |
| `_source`                                    | Functions that create source datasets that usually will be passed through other `derive_` functions.|

| Other Common Function Name Terms             | Description                                                                                         |
|----------------------------------------------|-----------------------------------------------------------------------------------------------------|
| `_merged_` / `_joined_` / `_extreme_`        | Functions that follow the [generic function user-guide](https://pharmaverse.github.io/admiral/articles/generic.html).                                              |



Please note that the appropriate *var*/*vars* prefix should be used for all cases in which the function creates any variable(s), regardless of the presence of a `new_var` argument in the function call. 

Oftentimes when creating a new `derive_var` or `derive_param` function there may be some sort of non-trivial calculation involved that you may want to write a customized function for. This is when creating a `compute_` function becomes appropriate, such that the calculation portion is contained in one step as part of the overall `derive_` function, reducing clutter in the main function body and assisting in debugging. In addition, a `compute_` function should be implemented if the calculation could be used for more than one derivation. For example `compute_bmi()` could be used to derive a baseline BMI variable in ADSL (based on baseline weight and baseline height variables) and could also be used to derive a BMI parameter in ADVS (based on weight and height parameters). Please see `compute_age_years()` and `derive_var_age_years()` as another example. 


## Function Arguments

The default value of optional arguments should be `NULL`.

There is a recommended argument order that all contributors are asked to adhere to 
(in order to keep consistency across functions):

1. `dataset` (and any additional datasets denoted by `dataset_*`)
1. `by_vars`
1. `order`
1. `new_var` (and any related `new_var_*` arguments)
1. `filter` (and any additional filters denoted by `filter_*`)
1. all additional arguments:
    * Make sure to always mention `start_date` before `end_date` (or related).

Names of variables inside a dataset should be passed as symbols rather than
strings, i.e. `AVAL` rather than `"AVAL"`. If an argument accepts one or more
variables or expressions as input then the variables and expressions should be
wrapped inside `exprs()`.

For example:

* `new_var = TEMPBL`
* `by_vars = exprs(PARAMCD, AVISIT)`
* `filter = PARAMCD == "TEMP"`
* `order = exprs(AVISIT, desc(AESEV))`
* `new_vars = exprs(LDOSE = EXDOSE, LDOSEDT = convert_dtc_to_dt(EXSTDTC))`

Each function argument needs to be tested with `assert_` type of function. 

Each expression needs to be tested for the following 
(there are many utility functions in `{admiral}` available to the contributor):

* whether it is an expression (or a list of expressions, depending on the function)
* whether it is a valid expression (i.e. whether it evaluates without error)

## Common Function Arguments Naming Convention

The first argument of `derive_` functions should be the input dataset and it
should be named `dataset`. If more than one input dataset is required, the other
input dataset should start with `dataset_`, e.g., `dataset_ex.`

Arguments for specifying items to add should start with `new_`. If a variable is
added, the second part of the argument name should be var, if a parameter is
added, it should be `param.` For example: `new_var`, `new_var_unit`,
`new_param`.

Arguments which expect a boolean or boolean vector must start with a verb, e.g.,
`is_imputed` or `impute_date`.

Arguments which only expect one value or variable name must be a singular version of the word(s), e.g., `missing_value` or `new_var`. Arguments which expect several values or variable names (as a list, expressions, etc.) must be a plural version of the word(s), e.g., `missing_values` or `new_vars`.

## List of Common Arguments
| Argument Name   | Description                                                                                                         |
|------------------|--------------------------------------------------------------------------------------------------------------------|
| `dataset`        | The input dataset. Expects a data.frame or a tibble.                                                               |
| `dataset_ref`    | The reference dataset, e.g. ADSL. Typically includes just one observation per subject.                             |
| `dataset_add`    | An additional dataset. Used in some `derive_xx` and `filter_xx` functions to access variables from an additional dataset.            |
| `by_vars`        | Variables to group by.                                                                                             |
| `order`          | List of expressions for sorting a dataset, e.g., `exprs(PARAMCD, AVISITN, desc(AVAL))`.                            |
| `new_var`        | Name of a single variable to be added to the dataset.                                                              |
| `new_vars`       | List of variables to be added to the dataset.                                                                      |
| `new_var_unit`   | Name of the unit variable to be added. It should be the unit of the variable specified for the `new_var` argument. |
| `filter`         | Expression to filter a dataset, e.g., `PARAMCD == "TEMP"`.                                                         |
| `start_date`     | The start date of an event/interval. Expects a date object.                                                        |
| `end_date`       | The end date of an event/interval. Expects a date object.                                                          |
| `start_dtc`      | (Partial) start date/datetime in ISO 8601 format.                                                                  |
| `dtc`            | (Partial) date/datetime in ISO 8601 format.                                                                        |
| `date`           | Date of an event / interval. Expects a date object.                                                                |
| `subject_keys`   | Variables to uniquely identify a subject, defaults to `exprs(STUDYID, USUBJID)`. In function formals, use `subject_keys = get_admiral_option("subject_keys")`                                     |
| `set_values_to`  | List of variable name-value pairs. Use `process_set_values_to()` for processing the value and providing user friendly error messages. |
| `keep_source_vars` | Specifies which variables from the selected observations should be kept. The default of the argument should be `exprs(everything())`. The primary difference between `set_values_to` and `keep_source_vars` is that `keep_source_vars` only selects and retains the variables from a source dataset, so e.g. `keep_source_vars = exprs(DOMAIN)` would join + keep the `DOMAIN` variable, whereas `set_values_to` can make renaming and inline function changes such as `set_values_to = exprs(LALVDOM = DOMAIN)`                       |
| `missing_value`    | A singular value to be entered if the data is missing.                                                           |
| `missing_values`   | A named list of expressions where the names are variables in the dataset and the values are a value to be entered if the data is missing, e.g., `exprs(BASEC = "MISSING", BASE = -1)`.      |

## Source Code Formatting

All source code should be formatted according to the [tidyverse](https://style.tidyverse.org/) style guide. The
[lintr](https://github.com/jimhester/lintr) and [styler](https://github.com/r-lib/styler) packages are used to
check and enforce this. 

With regards to [lintr](https://github.com/jimhester/lintr), `{admiral}` and all
related packages should maintain consistent linting standards by ensuring that their 
`.lintr.R` configuration files use the `admiral_linters()` function. This contains 
agreed-upon preferences and conventions such as avoiding the use of `stop()` and
`warning()` in favor of `cli::abort()` and `cli::warn()`.
The `admiral_linters()` function is stored under `inst/lintr/linters.R` in `{admiraldev}` 
(so as not to expose it to users) and can be loaded within the `.lintr.R` configuration 
file with `source(system.file("lintr/linters.R", package = "admiraldev")`. An example
`.lintr.R` configuration file is shown below:

```{R, eval = F}
library(lintr)

source(system.file("lintr/linters.R", package = "admiraldev"))

linters <- admiral_linters()

exclusions <- list(
  "R/data.R" = Inf,
  "inst" = list(undesirable_function_linter = Inf),
  "vignettes" = list(undesirable_function_linter = Inf)
)
```

If there is a good case to be made for altering any of the default configurations,
this can be done by passing arguments to `admiral_linters()`, for instance:

```{R, eval = F}
admiral_linters(
  line_length = line_length_linter(80), # (dropped down to 80 from the default 100)
  object_name_linter = object_name_linter( # (activated the object name linter)
    styles = c("snake_case", "symbols", "SNAKE_CASE"),
    regexes = c("^pl__.*", "^_.*")
  )
)
```

## Comments

Comments should be added to help other readers than the author to understand the
code. There are two main cases:

- If the intention of a chunk of code is not clear, a comment should be added.
The comment should not rephrase the code but provide additional information.

    *Bad*
    
    ```
      # If AVAL equals zero, set it to 0.0001. Otherwise, do not change it
      mutate(dataset, AVAL = if_else(AVAL == 0, 0.0001, AVAL))
    ```

    *Good*
    
    ```
      # AVAL is to be displayed on a logarithmic scale.
      # Thus replace zeros by a small value to avoid gaps.
      mutate(dataset, AVAL = if_else(AVAL == 0, 0.0001, AVAL))
    ```

- For long functions (>100 lines) comments can be added to structure the code
and simplify navigation. In this case the comment should end with `----` to add
an entry to the document outline in RStudio. For example:
  ```
    # Check arguments ----
  ```

The formatting of the comments must follow the
[tidyverse](https://style.tidyverse.org/syntax.html#comments) style guide. I.e.,
the comment should start with a single `#` and a space. No decoration (except
for outline entries) must be added.

*Bad*
```
# This is a comment #

###########################
# This is another comment #
###########################

#+++++++++++++++++++++++++++++++
# This is a section comment ----
#+++++++++++++++++++++++++++++++
```

*Good*
```
# This is a comment

# This is another comment

# This is a section comment ----
```

## Input Checking

In line with the [fail-fast](https://en.wikipedia.org/wiki/Fail-fast) design principle, 
function inputs should be checked for validity 
and, if there's an invalid input, the function should stop immediately with an error. 
An exception is the case where a variable to be added by a function already exists in the input dataset: 
here only a warning should be displayed and the function should continue executing.

Inputs should be checked using custom assertion functions defined in [`R/assertions.R`](https://github.com/pharmaverse/admiraldev/blob/main/R/assertions.R).
These custom assertion functions should either return an error in case of an invalid input or return nothing.

For the most common types of input arguments like a single variable, a list of
variables, a dataset, ... functions for checking are available (see
[assertions](../reference/index.html#section-assertions)).

Arguments which expect keywords should handle them in a case-insensitive manner,
e.g., both `date_imputation = "FIRST"` and `date_imputation = "first"` should be
accepted. The `assert_character_scalar()` function helps with handling arguments
in a case-insensitive manner.

A argument should not be checked in an outer function if the argument name is the same as in the inner function. 
This rule is applicable only if both functions are part of `{admiral}`.


## Function Header (Documentation)

Every function that is exported from the package must have an accompanying
header that should be formatted according to the
[roxygen2](https://roxygen2.r-lib.org/) convention. We have also implemented a
custom roclet to enhance our documentation and examples for more complex
functions - see `rdx_roclet()`, the example function `demo_fun()`, and further
down for more details.

In addition to the standard roxygen2 tags, the `@family` and `@keywords` tags are also used. 

The family/keywords are used to categorize the function, which is used both on our website and the internal package help pages. Please see section [Categorization of functions](https://pharmaverse.github.io/admiraldev/articles/programming_strategy.html#categorization-of-functions).

An example is given below:

```{r, eval=FALSE}
#' Derive Relative Day Variables
#'
#' Adds relative day variables (`--DY`) to the dataset, e.g., `ASTDY` and
#' `AENDY`.
#'
#' @param dataset Input dataset
#'
#'   The columns specified by the `reference_date` and the `source_vars`
#'   argument are expected.
#'
#' @permitted A dataset
#'
#' @param reference_date The start date column, e.g., date of first treatment
#'
#'   A date or date-time object column is expected.
#'
#'   Refer to `derive_var_dt()` to impute and derive a date from a date
#'   character vector to a date object.
#'
#' @permitted An unquoted variable name of the input dataset
#'
#' @param source_vars A list of datetime or date variables created using
#'   `exprs()` from which dates are to be extracted. This can either be a list of
#'   date(time) variables or named `--DY` variables and corresponding --DT(M)
#'   variables e.g. `exprs(TRTSDTM, ASTDTM, AENDT)` or `exprs(TRTSDT, ASTDTM,
#'   AENDT, DEATHDY = DTHDT)`. If the source variable does not end in --DT(M), a
#'   name for the resulting `--DY` variable must be provided.
#'
#' @permitted [var_list]
#'
#' @details The relative day is derived as number of days from the reference
#'   date to the end date. If it is nonnegative, one is added. I.e., the
#'   relative day of the reference date is 1. Unless a name is explicitly
#'   specified, the name of the resulting relative day variable is generated
#'   from the source variable name by replacing DT (or DTM as appropriate) with
#'   DY.
#'
#' @returns The input dataset with `--DY` corresponding to the `--DTM` or `--DT`
#'   source variable(s) added
#'
#' @keywords der_date_time
#' @family der_date_time
#'
#' @export
#'
#' @examples
#' library(lubridate)
#' library(dplyr, warn.conflicts = FALSE)
#'
#' datain <- tribble(
#'   ~TRTSDTM,              ~ASTDTM,               ~AENDT,
#'   "2014-01-17T23:59:59", "2014-01-18T13:09:O9", "2014-01-20"
#' ) %>%
#'   mutate(
#'     TRTSDTM = as_datetime(TRTSDTM),
#'     ASTDTM = as_datetime(ASTDTM),
#'     AENDT = ymd(AENDT)
#'   )
#'
#' derive_vars_dy(
#'   datain,
#'   reference_date = TRTSDTM,
#'   source_vars = exprs(TRTSDTM, ASTDTM, AENDT)
#' )
```


The following fields are mandatory:

* `@param`: One entry per function argument. 
The following attributes should be described: expected data type (e.g.
`data.frame`, `logical`, `numeric` etc.), permitted values (if applicable),
optionality (i.e. is this a required argument). If the expected input is a
dataset then the required variables should be clearly stated. Describing the
default value becomes difficult to maintain and subject to manual error when it
is already declared in the function arguments. For the description of the
permitted values the (custom) `@permitted` tag should be used (see
`rdx_roclet()` for more details).
* `@details`: A natural-language description of the derivation used inside the function.
* `@keyword`: One applicable tag to the function - identical to family.
* `@family`: One applicable tag to the function - identical to keyword.
* `@returns`: A description of the return value of the function. 
Any newly added variable(-s) should be mentioned here.
* `@examples` or `@caption`, `@info`, `@code`: Fully self-contained examples of
how to use the function. Self-contained means that, if this code is executed in
a new R session, it will run without errors. That means any packages need to be
loaded with `library()` and any datasets needed either to be created directly
inside the example code or loaded using `pkg_name::dataset_name`, e.g., `adsl <-
admiral::admiral_adsl`. If a dataset is created in the example, it should be
done so using the function `tribble()` (specify `library(dplyr)` before calling
this function). Make sure to align columns as this ensures quick code
readability. If other functions are called in the example, please specify
`library(pkg_name)` then refer to the respective function `fun()` as opposed to
the preferred `pkg_name::fun()` notation as specified in [Unit Test
Guidance](https://pharmaverse.github.io/admiraldev/articles/unit_test_guidance.html#set-up-the-test-script).

    The `@examples` tag should be used for simple functions which require only a
    few examples and no explanation. For more complex functions, the (custom) 
    `@caption`, `@info`, and `@code` tags should be used. Please see the separate 
    vignette on 
    [Writing Custom Examples](https://pharmaverse.github.io/admiraldev/articles/writing_custom_examples.html)
    for detailed guidance on how these are constructed, and [`derive_extreme_records.R` in             admiral](https://github.com/pharmaverse/admiral/blob/main/R/derive_extreme_records.R) 
    for an example of this in action.

Copying descriptions should be avoided as it makes the documentation hard to
maintain. For example if the same argument with the same description is used by
more than one function, the argument should be described for one function and
the other functions should use `@inheritParams <function name where the
argument is described>`.

Please note that if `@inheritParams func_first` is used in the header of the
`func_second()` function, those argument descriptions of `func_first()` are
included in the documentation of `func_second()` for which

- the argument is offered by `func_second()` and
- no `@param` tag for the argument is included in the header of
`func_second()`.

The order of the `@param` tags should be the same as in the function definition.
The `@inheritParams` tags should be after the `@param`. This does not affect the
order of the argument description in the rendered documentation but makes it
easier to maintain the headers.

Variable names, expressions, functions, and any other code must be enclosed
which backticks. This will render it as code.

For functions which derive a specific CDISC variable, the title must state the 
label of the variable without the variable name. The variable should be stated 
in the description.

To avoid confusion the term "parameter" should be used for CDISC parameters
only. For function arguments the term "argument" should be used.

## Categorization of Functions

The functions are categorized by keywords and families within the roxygen header. Categorization is important 
as the `admiral` user-facing functions base totals above 125 and is growing!  However, to ease the burden for developers, we have decided that
the keywords and families should be identical in the roxygen header, which are specified via the `@keywords` and `@family` fields. 
To reiterate, each function must use the **same keyword and family**. Also, please note that the keywords and families are case-sensitive.

### `@keywords`

The keywords allows for the reference page to be easily organized when using certain 
`pgkdown` functions. For example, using the function `has_keyword(der_bds_gen)` in the `_pkgdown.yml` file while building
the website will collect all the BDS General Derivation functions and display them in alphabetical order on the Reference Page in a section called
BDS-Specific.

### `@family`

The families allow for similar functions to be displayed in the **See Also** section of a function's documentation. For example, a user looking at
`derive_vars_dy()` function documentation might be interested in other Date/Time functions.  Using the `@family` tag `der_date_time` will display
all the Date/Time functions available in admiral to the user in the **See Also** section of `derive_vars_dy()` function documentation. Please take a look at the
function documentation for `derive_vars_dy()` to see the family tag in action.

Below are the list of available keyword/family tags to be used in `admiral` functions. If you think an additional keyword/family tag should be added, then please
add an issue in GitHub for discussion.


| Keyword/Family                                                                    | Description                                                                                                              |
|-----------------------------------------------------------------------------------|--------------------------------------------------------------------------------------------------------------------------|
| `com_date_time`                                                                   | Date/Time Computation Functions that returns a vector                                                                    | 	
| `com_bds_findings`                                                                | BDS-Findings Functions that returns a vector                                                                             |
| `create_aux`                                                                      | Functions for Creating Auxiliary Datasets |
| `datasets`                                                                        | Example datasets used within admiral                                                                                     |
| `der_gen`                                                                         | General Derivation Functions that can be used for any ADaM.                                                              |
| `der_date_time`                                                                   | Date/Time Derivation Function                                                                                            |
| `der_bds_gen`                                                                     | Basic Data Structure (BDS) Functions that can be used across different BDS ADaM (adex, advs, adlb, etc)                  |
| `der_bds_findings`                                                                | Basic Data Structure (BDS) Functions specific to the BDS-Findings ADaMs                                                  |
| `der_prm_bds_findings`                                                            | BDS-Findings Functions for adding Parameters                                                                             |
| `der_adsl`                                                                        | Functions that can only be used for creating ADSL.                                                                       |
| `der_tte`                                                                         | Function used only for creating a Time to Event (TTE) Dataset                                                            |
| `der_occds`                                                                       | OCCDS specific derivation of helper Functions                                                                            |
| `der_prm_tte`                                                                     | TTE Functions for adding Parameters to TTE Dataset                                                                       |
| `deprecated`                                                                      | Function which will be removed from admiral after next release.  See  [Deprecation Guidance](#deprecation).                                                                                                               |
| `metadata`                                                                        | Auxiliary datasets providing definitions as input for derivations, e.g. grading criteria or dose frequencies         |
| `utils_ds_chk`                                                                    | Utilities for Dataset Checking                                                                                           |	
| `utils_fil`                                                                       | Utilities for Filtering Observations                                                                                     |
| `utils_fmt`                                                                       | Utilities for Formatting Observations                                                                                    |	
| `utils_print`                                                                       | Utilities for Printing Objects in the Console                                                                                    |	
| `utils_help` 	                                                                    | Utilities used within Derivation functions                                                                               |	
| `utils_examples`                                                                  | Utilities used for examples and template scripts                                                                         |
| `source_specifications` 	                                                        | Source Objects                                                                                                             |	
| `other_advanced` 	                                                                | Other Advanced Functions                                                                                                   |	
| `high_order_function`                                                             |	Higher Order Functions                                                                                                   |	                                                                     |	
| `internal`                                                                        | Internal functions only available to admiral developers                                                                  |
|                                                                                   |                                                                                                                         |
| `assertion`*                                                                       | Asserts a certain type and gives warning, error to user        |
| `warning`                                                                         | Provides custom warnings to user       |
| `what`                                                                            | A function that ...      |
| `is`                                                                              | A function that ...      |
| `get`                                                                             | A function that ...      |

**NOTE:** It is strongly encouraged that each `@keyword` and `@family` are to be identical.  This eases the burden of development and maintenance for admiral functions. If you need to use multiple keywords or families, please reach out to the core development team for discussion.  


# Missing values

Missing values (`NA`s) need to be explicitly shown.

Regarding character vectors converted from SAS files: SAS treats missing character values as blank. 
Those are imported into R as empty strings (`""`) although in nature they are missing values (`NA`). 
All empty strings that originate like this need to be converted to proper R missing values `NA`.

# File Structuring

Organizing functions into files is more of an art than a science.
Thus, there are no hard rules but just recommendations.
First and foremost, there are two extremes that should be avoided:
putting each function into its own file and putting all functions into a single file.
Apart from that the following recommendations should be taken into consideration when deciding upon file structuring:

- If a function is very long (together with its documentation), store it in a separate file
- If some functions are documented together, put them into one file
- If some functions have some sort of commonality or relevance with one another (like `dplyr::bind_rows()` and `dplyr::bind_cols()`), put them into one file
- Store functions together with their helpers and methods
- Have no more than 1000 lines in a single file, unless necessary (exceptions are, for example, classes with methods)

It is the responsibility of both the author of a new function and reviewer to ensure that these recommendations are put into practice.


# R Package Dependencies

Package dependencies have to be documented in the `DESCRIPTION` file. 
If a package is used only in examples and/or unit tests then it should be listed in `Suggests`, otherwise in `Imports`.

Functions from other packages have to be explicitly imported by using the `@importFrom` tag in the `R/admiral-package.R` file. 
To import the `if_else()` and `mutate()` function from `dplyr` the following line would have to be included in that file:
`#' @importFrom dplyr if_else mutate`.
By using the `@importFrom` tag, it is easier to track all of our dependencies in one place and improves code readability. 

Some of these functions become critically important while using admiral and
should be included as an export. This applies to functions which are frequently
called within `{admiral }`function calls like `rlang::exprs()`, `dplyr::desc()`
or the pipe operator `dplyr::%>%`. To export these functions, the following R
code should be included in the `R/reexports.R` file using the format:

```
#' @export
pkg_name::fun
```

# Metadata

Functions should only perform the derivation logic and not add any kind of metadata, e.g. labels.


# Unit Testing

A function requires a set of unit tests to verify it produces the expected result.
See [Writing Unit Tests in {admiral}](https://pharmaverse.github.io/admiraldev/articles/unit_test_guidance.html#writing-unit-tests-in-admiral) for details.

# Deprecation

The below deprecation strategy provides stability to users while allowing admiral developers
the ability to remove and update the code base in the coming days.

- **Phase 1:** In the release where the identified function or argument is to
be deprecated, there will be a message issued when using the function or argument
using `deprecate_inform()`. This message will appear to the user for at least
_one year_. Templates, vignettes and any internal calls should be updated to use the new recommended function/argument.   

- **Phase 2:** After at least _one year_ and in the closet next release, a warning will be 
issued when using the function or argument using `deprecate_warn()`. This warning 
message will appear for at least _one year_.   

- **Phase 3:** After at least _one year_ and in the closest next release, an error will be thrown
when using the function or argument using `deprecate_stop()` and follow similar process
for Phase 1 and Phase 2.

- **Phase 4:** Finally after three years from the time of being identified for deprecation, the
function or argument will be completely removed from `{admiral}`.

**NB:** Major/Minor release make the most sense for deprecation updates. However, if
a release cycle becomes multiple years, then patch releases should be considered to
help keep `{admiral}` neat and tidy!

**NB:** Take care with the `NEWS.md` entries around deprecation as the person continuing this
process might not be you!

## Documentation

If a function or argument is removed, the documentation must be updated to
indicate the function or the argument is now deprecated and which new
function/argument should be used instead.

The documentation will be updated at Phase 1:

+ the description level for a function will have a lifecycle badge added, 
+ the `@keywords` and `@family` roxygen tags will be replaced with `deprecated`

  ```{r, eval=FALSE}
#' Title of the function
#'
#' @description
#' `r lifecycle::badge("deprecated")`
#'
#' This function is *deprecated*, please use `new_fun()` instead.
#' .
#' @family deprecated
#' @keywords deprecated
  ``` 


Example for documentation at the argument level

+ the `@param` level for a argument.

  ```{r, eval = FALSE}
#' @param old_param `r lifecycle::badge("deprecated")` Please use `new_param` instead.
  ```

The documentation will be further updated at Phase 3:

+ the `@examples` section should be removed.

## Handling of Messages, Warnings and Errors

When a function or argument is deprecated, the function must be updated to issue
a message, warning or error using `deprecate_inform()`, `deprecate_warn()` or 
`deprecate_stop()`, respectively, as described above.

There should be a test case added in the test file of the function that checks
whether this message/warning/error is issued as appropriate when using the deprecated
function or argument.

### Function

**Phase 1:** At the start of this phase the call to `deprecate_inform()` will 
appear as:

```r
fun_xxx <- function(dataset, some_param, other_param) {
  deprecate_inform(
    when = "x.y.z",
    what = "fun_xxx()",
    with = "new_fun_xxx()",
    details = c(
      x = "This message will turn into a warning {at the beginning of 20XX}.",
      i = "See admiral's deprecation guidance:
              https://pharmaverse.github.io/admiraldev/dev/articles/programming_strategy.html#deprecation"
    )
  )
  
  new_fun_xxx(
    dataset = dataset,
    some_param = some_param,
    other_param = other_param
  )
}
```
__NB:__ Please adjust phrase `{at the beginning of 20XX}` to relevant timeline.

The code of the deprecated function should be replaced with a call to the new
function which should be used instead.

**Phase 2:** At the start of this phase the call to `deprecate_warn()` will 
appear as:

```r
fun_xxx <- function(dataset, some_param, other_param) {
  deprecate_warn(
    when = "x.y.z",
    what = "fun_xxx()",
    with = "new_fun_xxx()",
    details = c(
      x = "This message will turn into a error {at the beginning of 20XX}.",
      i = "See admiral's deprecation guidance:
              https://pharmaverse.github.io/admiraldev/dev/articles/programming_strategy.html#deprecation"
    )
  )

  new_fun_xxx(
    dataset = dataset,
    some_param = some_param,
    other_param = other_param
  )
}
```
__NB:__ Please adjust phrase `{at the beginning of 20XX}` to relevant timeline.

**Phase 3:** At the start of this phase the call to `deprecate_stop()` will 
appear as:

```r
fun_xxx <- function(dataset, some_param, other_param) {
  deprecate_stop(
    when = "x.y.z",
    what = "fun_xxx()",
    with = "new_fun_xxx()"
  )

  new_fun_xxx(
    dataset = dataset,
    some_param = some_param,
    other_param = other_param
  )
}
```

**Phase 4:** Function should be removed from the package.

### Argument

**Phase 1:** If the argument is renamed or replaced, a **message** must be issued and the
new argument takes the value of the old argument until the next phase. Note:
arguments which are not passed as `exprs()` argument (e.g. `new_var = VAR1` or
`filter = AVAL > 10`) will need to be quoted.

``` 
  if (!missing(old_param)) {
    deprecate_inform("x.y.z", "fun_xxx(old_param = )", "fun_xxx(new_param = )")
    # old_param is given using exprs()
    new_param <- old_param
    # old_param is NOT given using exprs()
    new_param <- enexpr(old_param)
  }
```

**Phase 2:** If the argument is renamed or replaced, a **warning** must be issued and the
new argument takes the value of the old argument until the next phase Note:
arguments which are not passed as `exprs()` argument (e.g. `new_var = VAR1` or
`filter = AVAL > 10`) will need to be quoted.

``` 
  if (!missing(old_param)) {
    deprecate_warn("x.y.z", "fun_xxx(old_param = )", "fun_xxx(new_param = )")
    # old_param is given using exprs()
    new_param <- old_param
    # old_param is NOT given using exprs()
    new_param <- enexpr(old_param)
  }
```

**Phase 3:**  If an argument is removed and is not replaced, an **error** must be generated:

```
  if (!missing(old_param)) {
    deprecate_stop("x.y.z", "fun_xxx(old_param = )", "fun_xxx(new_param = )")
  }
```

**Phase 4:**  All mentions of the argument are completely removed from admiral.  

## Unit Testing

Unit tests for deprecated functions and arguments must be added to the test 
file [^1] of the function to ensure that a message, warning, or error is issued.

[^1]: For example, if `derive_var_example()` is going to be deprecated and
it is defined in `examples.R`, the unit tests are in
`tests/testthat/test-examples.R`.

The unit-test should follow the corresponding format, per the [unit test
guidance](https://pharmaverse.github.io/admiraldev/articles/unit_test_guidance.html#writing-unit-tests-in-admiral).

### For Deprecated Functions that Issue a Message (Phase 1)

* Please put tests for deprecation at the top of test file to make finding the specific test
easier for next phase of deprecation.
* Tests that call multiple functions with deprecation messages can be wrapped using 
parentheses and curly brackets, e.g. `expect_snapshot({})`.
* You can use `withr::local_options(list(lifecycle_verbosity = "quiet"))` to suppress
the deprecation messages in already created tests.  

    ```
    ## Test 1: deprecation message if function is called ----
    test_that("derive_var_example() Test #: deprecation message if function is called", {
      expect_snapshot({
        ae <- date_source(...)
        ...
        derive_var_example(...)
      })
    })

    ## Test 2: Test of function argument 1 ----
    test_that("derive_var_example() Test 2: Test of function argument 1", {
      withr::local_options(list(lifecycle_verbosity = "quiet"))
      ...
    })
```

### For Deprecated Functions that Issue a Warning (Phase 2)

The snapshot of the deprecation message test must be updated because instead of
a message a warning is issued now.

### For Deprecated Functions that Issue an Error (Phase 3)

A unit test like the following must be added.
```
## Test #: error if function is called ----
test_that("derive_var_example() Test #: deprecation error if function is called", {
  expect_error(
    derive_var_example(),
    class = "lifecycle_error_deprecated"
  )
})
```
When writing the unit test, check that the error has the right class, i.e.,
`"lifecycle_error_deprecated"`.

Other unit tests of the deprecated function must be removed.

# Best Practices and Hints

Please take the following list as recommendation and try to adhere to its rules if possible.

* Arguments in function calls should be named except for the first parameter 
(e.g. `assert_data_frame(dataset, required_vars = exprs(var1, var2), optional = TRUE)`).
* `dplyr::if_else()` should be used when there are only two conditions. 
Try to always set the `missing` argument whenever appropriate.

## How Quoting is used

* Some admiral arguments require selecting one particular option like `mode`, e.g. `mode = "last"`. Use quotation marks to capture these. The expected assertion function corresponding to these arguments is `assert_character_scalar()/assert_character_vector()`.
* Many admiral arguments require capturing an expression, typically encased in a `exprs()` statement, which are to be evaluated _later_ inside the function body, see arguments like `new_vars`, e.g. `new_vars = exprs(TRTSDTM = EXSTDTM)`. Oftentimes, the assertion function corresponding to these are `assert_expr()/assert_expr_list()`. These arguments are unquoted by using `!!!`.
* Some admiral arguments like `new_var` or `filter` which expect a _single_ variable or expression are not quoted in the call. In the function body, it has to be quoted by using `enexpr().` Usually this is combined with the assertion, e.g., `new_var <- assert_symbol(enexpr(new_var))`. These arguments are unquoted by using `!!`.
* Keep in mind `!!` is a one-to-one replacement and `!!!` is a one-to-many replacement. Please see [this chapter](https://adv-r.hadley.nz/quasiquotation.html) in the Advanced R textbook for more details.

## Standardizing Text Used to Label and Describe Arguments

In the following [PR](https://github.com/pharmaverse/admiral/pull/2065/files), you will find an example of how the function argument `dataset` was able to be standardized such that the Label and Description of said function argument was aligned across the codebase. Please see the changes to the file `derive_adeg_params.R` for further details. 

The benefits of having a programmatic way to write documentation is that if any changes need to be made, _making the modification on the corresponding function, in this case, `roxygen_param_dataset()`, scales across the codebase, can be tested, and is less prone to user-error such as typos or grammar mistakes_.

These functions are implemented in `roxygen2.R` and the naming convention for each argument will be as follows `roxygen_param_xxx()`, where "xxx" is the be replaced with the argument name.

# R and Package Versions for Development

* The choice of R Version and Package versions are not set in stone.  However, a common development environment is important to establish when working across multiple companies and multiple developers. We currently recommend developers work with the latest R version and latest available packages. However, this will deviate over time as developers come and go from `{admiral}`. We actually see this as a positive, i.e. the deviations between developers, as this introduces a bit of random stress testing to our code base. 
* GitHub allows us through the Actions/Workflows to test `{admiral}` under several versions of R as well as several versions of dependent R packages needed for `{admiral}`. Currently we test `{admiral}` against the two latest R Versions and the closest snapshots of packages to those R versions.  You can view this workflow and others on our [admiralci GitHub Repository](https://github.com/pharmaverse/admiralci).

---

# Writing Vignettes

**Description:** Guidelines for creating comprehensive and useful vignettes  
**Source File:** `writing_vignettes.Rmd`  
**Source:** GitHub (pharmaverse/admiraldev)  
**URL:** https://raw.githubusercontent.com/pharmaverse/admiraldev/main/vignettes/writing_vignettes.Rmd

---

---
title: "Writing Vignettes"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Writing Vignettes}
  %\VignetteEngine{knitr::rmarkdown}
---

```{r setup, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
```

# Introduction

This guidance ensures consistency across the vignettes in the `{admiral}` package in terms of content, structure and code execution. 
As a general rule, the vignette workflow defined in [r-pkgs.org](https://r-pkgs.org/vignettes.html) is followed.


# Metadata

Each vignette in `{admiral}` should start with the following metadata.

```
---
title: "<Title>"
output:  rmarkdown::html_vignette:
vignette: >
  %\VignetteIndexEntry{<Title>}
  %\VignetteEngine{knitr::rmarkdown}
---
```

The `<Title>` of the vignette should be meaningful. 

# Markdown

## Default Options

If any chunks are used within the vignette, the following options should be set after the metadata 
to ensure that the chunks are rendered consistently across all vignettes.

    `r ''````{r setup, include=FALSE}
    knitr::opts_chunk$set(
      collapse = TRUE,
      comment = "#>"
    )
    ```


## Format Sections

### Table of Contents

Headings must be title case and start from Heading 1:

```
# Heading 1
This is the description of my first section.

## Heading 1.1
This is the description of my first sub-section within my first section.

## Heading 1.2
This is the description of my second sub-section within my first section.

```

The first section gives a brief introduction of the vignette. The last sub-section 
of the introduction should describe the packages required to run the `{admiral}` functions. 
The `{admiral}` package should be described first.

```
# Introduction

This is the introduction of my vignette.

## Required Packages
```

    `r ''````{r, warning=FALSE, message=FALSE}
    library(admiral)
    # <all packages used in the vignette>
    ```

The `warning=FALSE` and `message=FALSE` options are there to prevent the usual messages:

<span style="color: red;font-family: textmate, monospace; font-size = 10pt;" >Attaching package: 'xxxx'</br>The following objects are masked from 'package:yyyyy'</br>  fun1, fun2</span>

### Conventions

#### General Conventions

+ Any new vignette must be added to the `_pkgdown.yml` file in the relevant section.

+ Any variable name, dataset name, function, argument name must be quoted with backticks: e.g.


```
  The `date` parameter of the `derive_my_var()` function expects a date variable, e.g., `ADT`.
```

+ Functions must also end with `()`.

+ Variables and datasets name are expected to be uppercase.

+ All the codes must be described, executed and the output result displayed once the code is executed.
Use:
 

      `r ''````{r}
      #<code>
      ```

+ Any output created must clearly show what the function has derived. 
It should at least show the variables/records used as input to the function and the derived 
variables/records. If a dataset must be displayed, it will be formatted using the `dataset_vignette()`
function so that it is displayed consistently across the vignettes.E.g.


  + Description and execution of the code used to derive the variable/record
  


```{r, include=FALSE}
library(lubridate)
library(dplyr)
library(DT)

data(vs)
vs <- tribble(
  ~USUBJID, ~VSTESTCD, ~VISIT, ~VSDTC,
  "01-701-1015", "WEIGHT", "BASELINE", "2013-12-26",
  "01-701-1015", "WEIGHT", "WEEK 2", "2014-01-14",
  "01-701-1016", "WEIGHT", "BASELINE", "2014-01-06",
  "01-701-1016", "WEIGHT", "WEEK 2", "2014-01-22",
)
```

```{r, eval=FALSE, echo=TRUE}
vs1 <- vs %>%
  derive_vars_dt(
    new_vars_prefix = "A",
    dtc = VSDTC,
    date_imputation = "first"
  )
```


  + Output dataset formatted using `dataset_vignette()`...

```{r, eval=FALSE, echo=FALSE}
vs1 <- vs %>%
  derive_vars_dt(
    new_vars_prefix = "A",
    dtc = VSDTC,
    date_imputation = "first"
  )

dataset_vignette(
  vs1,
  display_vars = exprs(USUBJID, VSTESTCD, VISIT, VSDTC, ADT),
  filter = VSTESTCD == "WEIGHT"
)
```


Note: The call to get the formatted dataset would be: 
 

      ```{r, eval=FALSE, echo=TRUE}
dataset_vignette(
  vs1,
  display_vars = exprs(USUBJID, VSTESTCD, VISIT, VSDTC, ADT),
  filter = VSTESTCD == "WEIGHT"
)
      ```
      
Displaying many big datasets on a vignette, may require a long time to load
the page and may cause messages from the browser that the page is not
responsive. In this case the number of displayed observations should be
restricted either by restricting the source datasets at the beginning of the
vignette or in the call of `dataset_vignette()` if only some calls are
affected.

#### Using Footnotes

Footnotes can be useful to add context without adding clutter to the primary subject matter of the vignette being written.

To use footnotes, add a caret and an identifier inside brackets e.g. `([^1])`. The identifiers can be numbers or words, but they can't contain spaces or tabs.

The following markdown text will render as follows:

```
Here is an example [^1]

[^1]: Cool note about the example.
```

Here is an example [^2]

[^2]: Cool note about the example.


#### Conventions for ADaM Workflow

For vignettes describing an ADaM workflow, 

+ the second section will summarize the programming workflow. The first sub-section within this 
workflow will always describe the data read to demonstrate the use of the `{admiral}` functions,

+ Each sub-section within the programming workflow should be tagged (e.g. [Step1] (#step)), so that 
the user can go on the relevant section from the programming workflow (in addition to the Table of
contents). Don't use a tag with a number but use a meaningful name (e.g. do not use `(#link1)`,
but use `(#this_action)`)

+ the last section should link to a template script. 


```
# Programming Workflow

* [Read in Data](#readdata)
* [Derive/Impute End and Start Analysis Date/time and Relative Day](#datetime)
* ...
* [Assign `ASEQ`](#aseq)

## Read in Data {#readdata}
## Derive/Impute End and Start Analysis Date/time and Relative Day {#datetime}
## ...
## Assign `ASEQ` {#aseq}

# Another Section

# Example Script
ADaM | Sample Code
---- | --------------
ADxx | [ad_adxx.R](https://github.com/pharmaverse/admiral/blob/main/inst/templates/ad_adxx.R){target="_blank"}

```

+ `ADSL` variables

  All `ADSL` variables which are required for any derivation are merged to the SDTM dataset before the
  first derivation. 
  These `ADSL` variables have to be added to the by-variables for derivations which add observations. 
  This ensures that the `ADSL` variables are populated `for` the new observations. 
  A `adsl_vars` variable should be created at the beginning of the script and added to the `by_vars` 
  parameter for derivations which add observations.

  `ADSL` variables which should be in the final dataset but are not required for any derivation are
  merged to the dataset after the last derivation.



---

# Unit Test Guidance

**Description:** Best practices for writing unit tests in admiral packages  
**Source File:** `unit_test_guidance.Rmd`  
**Source:** GitHub (pharmaverse/admiraldev)  
**URL:** https://raw.githubusercontent.com/pharmaverse/admiraldev/main/vignettes/unit_test_guidance.Rmd

---

---
title: "Unit Test Guidance"
output: 
  rmarkdown::html_vignette:
    toc: true
vignette: >
  %\VignetteIndexEntry{Unit Test Guidance}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r setup, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
```


# Why Write Unit Tests?

## Unit Tests Become a Safety Net for Developers

A comprehensive suite of unit tests can act as a safety net for developers. By 
frequently running the tests, they can assure their recent modifications to the 
code haven't broken anything. In other words, unit tests help prevent regressions.

## Unit Tests Can Contribute to Higher Code Quality

Since unit tests act as a safety net, developers become more confident when 
changing the code. They can refactor the code without fear of breaking things, 
driving the general quality of the code base up.

## Unit Tests Can Contribute to Better Application Architecture

If you can add unit tests easily to a code base, that's usually a good sign 
regarding the quality of the app's architecture. So, the drive to write testable
code can be an incentive for better architecture.

## Detects Code Smells in your Codebase

If ease of adding unit tests to a code base is a good sign, the opposite is also
true. Having a hard time creating unit tests for a given piece of code might be 
a sign of code smells in the code<U+2014>e.g. functions that are too complex.

# Writing Good Unit Tests

## Tests Should Be Fast

If they're slow, developers won't run them as often as they should. That defeats
the whole purpose of having a suite of unit tests in the first place, which is 
to boost the developers' confidence to make changes to the code. The tests can't
work as the safety net if they're not run often.

## Tests Should Be Simple

There are several techniques we can apply to have a high degree of confidence in
the correctness of our tests. One of those is to keep your tests with low 
cyclomatic complexity. Cyclomatic complexity is a code metric that indicates the
number of possible execution paths a given method can follow. A piece of code 
with lower complexity is easier to understand and maintain, which means developers 
are less likely to introduce bugs when working on it. We can measure the 
cyclomatic complexity of your tests (using, for instance, a linter tool) and do 
your best to keep it low.

## Test Shouldn't Duplicate Implementation Logic

If the same person wrote both the test and the implementation, it's possible they
made the same errors in both places. Since the tests mirror the implementation,
they might still pass, and the implementation could be wrong, but the tests might
fool you into thinking otherwise. Resist the urge to make your tests fancy, keep
them simple, and your testing suite will be better for it.

## Tests Should Be Readable

This best practice overlaps a little bit with the one about keeping your tests 
simple. If tests are hard to read, developers are more likely to misunderstand 
them and introduce bugs. Test cases could be used as a form of documentation, so
they obviously need to be readable.

## Running Unit Tests Part of the Build Process

Automate the whole process of running the unit tests and taking some action when
they fail. Your build process should execute your unit tests and mark the build 
as broken when the tests fail.

# Writing Unit Tests in {admiral}

## Plan your Unit Tests

Start by considering the derivation rule you are testing and the possible 
arguments/flexibilities of your function code. Then plan which scenarios you will
test. These can either involve generating different input test cases or feeding 
them into different calls of your function.

## Test coverage

Unit tests should cover the functionality of the function.
If another function `g()` is called within a function `f()`, the unit tests of `f()` should not test the functionality of `g()`. 
This should be tested by the unit tests of `g()`, i.e. unit tests should be added at the lowest level.

## Tests Should be Robust to Cover Realistic Data Scenarios

For generating input test cases, it can be helpful to consider regular cases 
(expected common data scenarios), boundary cases (where data points are close or
equal), and special cases (uncommon but valid data scenarios, e.g. missing or 
special characters). Although you will never cover every single eventuality of 
possible input data (no reliability testing method ever gives 100% certainty), 
you do need to give confidence that the code is robust enough to work across most
data scenarios.

## Testing Should Cover Possible Arguments

For the different calls of your function, consider how the user might apply your
function and test a variety of possible calls, whilst still remembering the tips
above that tests should be fast and simple. 
This is only needed in cases where the complexity and level of flexibility of 
your function justifies it, e.g. see the test script: https://github.com/pharmaverse/admiral/blob/main/tests/testthat/test-derive_var_extreme_flag.R.

## Exported Functions

Don't forget to add a unit test for each exported function.

## Snapshot Testing

Standard unit tests are not always convenient to record the expected behavior with code.
Some challenges include:

  - Output that is large, making it painful to define the reference output, and bloating the size of the test file and making it hard to navigate.
  - Text output that includes many characters like quotes and newlines that require special handling in a string.  
  - Binary formats like plots or images, which are very difficult to describe in code: i.e. the plot looks right, the error message is useful to a human, the print method uses color effectively.

For these situations, testthat provides an alternative mechanism: snapshot tests.
Snapshot tests record results in a separate human readable file and records the results, including output, messages, warnings, and errors.
Review the [{testthat} snapshot vignette](https://testthat.r-lib.org/articles/snapshotting.html) for details.

## Set up the Test Script

Within the `tests/testthat` folder of the project, add a script with the naming
convention `test-<script_containing_function>.R`., the unit test script can be
created from the console also, as follows:

```
usethis::use_test("<script_containing_function>")
```
the testing framework used is testthat and has the following format :

```
## Test 1: <Explanation of the test> ----
test_that("<function_name> Test 1: <Explanation of the test>", {
  
  input <- dplyr::tribble(
    ~inputvar1, ~inputvar2, ...
    <Add Test Data Scenarios>
    ...
  )

  expected_output <- mutate(input, outputvar = c(<Add Expected Outputs>))

  expect_dfs_equal(<function name>(input), expected_output)
  
})
```

For example, if you are testing a function called `my_new_func` that is contained 
in script `all_funcs.R` then from console use:

```
usethis::use_test("all_funcs")
```

Open the newly created file `test-all_funcs.R`  and use the following format:

```
# my_new_func ----
## Test 1: <Explanation of the test> ----
test_that("my_new_func Test 1: <Explanation of the test>", {
  
  input <- dplyr::tribble(
    ~inputvar1, ~inputvar2, ...
    <Add Test Data Scenarios>
    ...
  )

  expected_output <- mutate(input, outputvar = c(<Add Expected Outputs>))

  expect_dfs_equal(<function name>(input), expected_output)
})
```
**Note**: When comparing datasets in `{admiral}` we use function `expect_dfs_equal()`. 

The input and expected output for the unit tests must follow the following rules:

* Input and output should be as simple as possible.
* Values should be hard-coded whenever possible.
* If values need to be derived, only unit tested functions can be used.

In contrast to the [Programming Strategy](https://pharmaverse.github.io/admiraldev/articles/programming_strategy.html#function-header-documentation) documentation for function examples, test files should not include `library(pkg_name)` calls. 
If a dataset needs to be created for testing purposes, it should be done so using the function `tribble()` from the `tibble` package with the following command `dplyr::tribble(<data here>)`. 
Furthermore, if other functions need to be called, it should also be done using `pkg_name::fun()`notation. 
Make sure to align columns as well. This ensures quick code readability.

Ensure you give a meaningful explanation of the test in the testthat call, as 
these will be compiled in the package validation report. Having the name of the
function and test ID included in title will also help with traceability.

The comments ending with `----` create entries in the TOC in RStudio.

```{r echo=FALSE, out.width='120%'}
knitr::include_graphics("./unit_test_toc.png")
```


## Addin `pharmaverse4devs::format_test_that_file()`

To ease the burden on developers for writing and adding tests we have developed an Addin for formatting test_that test files according to admiral programming standards. The Addin will add and update comments as well as number or re-numbers the tests. To access the Addin, be sure to install the {pharmaverse4devs} from Github. 

To install the latest development version of the package directly from
GitHub use the following code:
```r
    if (!requireNamespace("remotes", quietly = TRUE)) {
      install.packages("remotes")
    }

    remotes::install_github("pharmaverse/pharmaverse4devs")
```

Then use the Addin button and select the "Format
test_that test file" as seen in the image.  Be sure to have the test-file open and selected when calling the Addin.

```{r echo=FALSE, out.width='120%'}
knitr::include_graphics("./unit_test_format_tests.png")
```

The Addin will perform the following:

- Updates or adds the number of the tests in the comments and in the
`test_that()` call
- Updates the comments based on the description provided in the `test_that()`
call
- Updates the function name in the `test_that()` call. The function name is
extracted from the last `# <function name> ----` comment before the
`test_that()` call. If a test file tests more than one function, such comments
should be added before the first test of each function. If a test files tests a
single function only, the comments can be omitted. In this case the addin
determines the function name from the file name by stripping of the "test-"
prefix and the ".R" suffix.

When writing new unit tests, just provide a description in the `test_that()`
call and if necessary the function name in a `# <function name> ----` comment:
```
# derive_vars_merged ----
test_that( "works if it merges all variables", {
  actual <- derive_vars_merged(advs,
    dataset_add = adsl,
    by_vars = exprs(STUDYID, USUBJID)
  )

# convert_dtm_to_dtc ----
test_that("works if dtm is in correct format", {
  expect_equal(
    convert_dtm_to_dtc(as.POSIXct("2022-04-05 15:34:07 UTC")),
    "2022-04-05T15:34:07"
  )
})

test_that("Error is thrown if dtm is not in correct format", {
  expect_error(
    convert_dtm_to_dtc("2022-04-05T15:26:14"),
    "lubridate::is.instant(dtm) is not TRUE",
    fixed = TRUE
  )
})
```

Call the addin and get:
```
# derive_vars_merged ----
## Test 1: derive_vars_merged ----
test_that( "derive_vars_merged Test 1: it merges all variables", {
  actual <- derive_vars_merged(advs,
    dataset_add = adsl,
    by_vars = exprs(STUDYID, USUBJID)
  )

# convert_dtm_to_dtc ----
## Test 2: works if dtm is in correct format ----
test_that("convert_dtm_to_dtc Test 2: works if dtm is in correct format", {
  expect_equal(
    convert_dtm_to_dtc(as.POSIXct("2022-04-05 15:34:07 UTC")),
    "2022-04-05T15:34:07"
  )
})

## Test 3: Error is thrown if dtm is not in correct format ----
test_that("convert_dtm_to_dtc Test 3: Error is thrown if dtm is not in correct format", {
  expect_error(
    convert_dtm_to_dtc("2022-04-05T15:26:14"),
    "lubridate::is.instant(dtm) is not TRUE",
    fixed = TRUE
  )
})
```

Once you have tested your unit test program, you can run all unit tests from
the console, as follows.

```
devtools::test()
```

For running just the tests of the current file call

```
devtools::test_file()
```

## Automation of Unit Tests

When a user actions a pull request in {admiral} GitHub repo, the unit tests are 
automatically run and pull request will be denied if any unit tests fail.

## Flow Chart

```{r echo=FALSE, out.width='120%'}
knitr::include_graphics("./unit_test_guidance.png")
```

---

---

## Keeping This Document Updated

**Generated:** 2025-12-09 21:56:31  
**admiraldev Version:** GitHub main branch  
**Vignettes Included:** 3 of 3  
**Output File:** `/home/jeffreyd/admiral/.github/copilot-instructions.md`

### Manual Update
```r
source('.github/scripts/sync_admiraldev_docs.R')
```

### Automatic Updates

Consider setting up a monthly cron job or GitHub Action to keep these
guidelines synchronized with the latest admiraldev version.

See `.github/scripts/sync_admiraldev_docs.R` for setup instructions.

---

## How GitHub Copilot Uses This File

GitHub Copilot automatically reads files in `.github/` directories and uses them
as context when providing code suggestions. By keeping admiraldev guidelines here:

1. **Copilot suggests code** that follows admiral conventions
2. **Function names and patterns** match admiral ecosystem standards  
3. **Documentation style** aligns with admiral expectations
4. **Test structures** follow unit test guidance
5. **Vignette examples** use recommended patterns

You don't need to do anything special - just having this file in `.github/`
is enough for Copilot to use it.

---

*This is an automated sync of admiraldev documentation for GitHub Copilot context.*
