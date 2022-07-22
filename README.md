<!-- Please do not edit the README.md file as it is auto-generated. Only edit the README.Rmd file -->

# admiraldev - ADaM in R Asset Library Development Utilities

<img src="man/figures/admiraldev.png" align="right" width="200" style="margin-left:50px;"/>

<!-- badges: start -->

[<img src="http://pharmaverse.org/shields/admiral.svg">](https://pharmaverse.org)
[![CRAN
status](https://www.r-pkg.org/badges/version/admiral)](https://CRAN.R-project.org/package=admiral)
[![R-CMD-check](https://github.com/pharmaverse/admiraltemplate/actions/workflows/R-CMD-check.yml/badge.svg)](https://github.com/pharmaverse/admiraltemplate/actions/workflows/R-CMD-check.yml)
[![Test
Coverage](https://raw.githubusercontent.com/pharmaverse/admiraltemplate/badges/main/test-coverage.svg)](https://github.com/pharmaverse/admiraltemplate/actions/workflows/code-coverage.yml)
<!-- badges: end -->

## Purpose

Tools for developing functions and maintaining a healthy codebase within
the family of admiral R packages.

## Installation

To install the latest development version of the package directly from
GitHub use the following code:

    if (!requireNamespace("remotes", quietly = TRUE)) {
      install.packages("remotes")
    }

    remotes::install_github("pharmaverse/admiral.test", ref = "devel") # This is a required dependency of {admiral}
    remotes::install_github("pharmaverse/admiral", ref = "devel")
    remotes::install_github("pharmaverse/admiraldev", ref = "devel")

## Structure of package

    ## R
    ## ├── addin_format_testthat.R
    ## ├── admiraldev-package.R
    ## ├── assertions.R
    ## ├── compat_friendly_type.R
    ## ├── dataset_vignette.R
    ## ├── dev_utilities.R
    ## ├── expect_dfs_equal.R
    ## ├── get.R
    ## ├── global.R
    ## ├── is.R
    ## ├── quo.R
    ## ├── quote.R
    ## ├── reexports.R
    ## ├── warnings.R
    ## └── what.R

## Keyword/Family

<table>
<colgroup>
<col style="width: 40%" />
<col style="width: 59%" />
</colgroup>
<thead>
<tr class="header">
<th>Keyword</th>
<th>Description</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td><code>assertion</code></td>
<td>Asserts a certain type and gives warning, error to user</td>
</tr>
<tr class="even">
<td><code>warning</code></td>
<td>Provides custom warnings to user</td>
</tr>
<tr class="odd">
<td><code>what</code></td>
<td>A function that …</td>
</tr>
<tr class="even">
<td><code>is</code></td>
<td>A function that …</td>
</tr>
<tr class="odd">
<td><code>get</code></td>
<td>A function that …</td>
</tr>
</tbody>
</table>

## Release Schedule

Releases are done at the end of every of quarter on the last Monday.
Pull Requests will be frozen the week before a release. The `admiral`
family has several downstream and upstream dependencies and so this
release shall be done in three Phases.

-   Phase 1 is for `admiraldev` and `admiral.test`, which feed into all
    `admiral` packages
-   Phase 2 is only for core `admiral`
-   Phase 3 is extension packages, e.g. `admiralonco`

<table>
<colgroup>
<col style="width: 24%" />
<col style="width: 25%" />
<col style="width: 24%" />
<col style="width: 25%" />
</colgroup>
<thead>
<tr class="header">
<th>Release Schedule</th>
<th>Phase 1- Date and Packages</th>
<th>Phase 2- Date and Packages</th>
<th>Phase 3- Date and Packages</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td>Q3-2022</td>
<td>August 29th</td>
<td>September 5th</td>
<td>September 12th</td>
</tr>
<tr class="even">
<td></td>
<td><code>admiraldev</code></td>
<td><code>admiral</code></td>
<td><code>admiralonco</code></td>
</tr>
<tr class="odd">
<td></td>
<td><code>admiral.test</code></td>
<td></td>
<td><code>admiralroche</code></td>
</tr>
<tr class="even">
<td>Q4-2022</td>
<td>November 28th</td>
<td>December 5th</td>
<td>December 12th</td>
</tr>
<tr class="odd">
<td></td>
<td><code>admiraldev</code></td>
<td><code>admiral</code></td>
<td><code>admiralonco</code></td>
</tr>
<tr class="even">
<td></td>
<td><code>admiral.test</code></td>
<td></td>
<td><code>admiralroche</code></td>
</tr>
<tr class="odd">
<td>Q1-2023</td>
<td>February 27th</td>
<td>March 6th</td>
<td>March 13th</td>
</tr>
<tr class="even">
<td></td>
<td><code>admiraldev</code></td>
<td><code>admiral</code></td>
<td><code>admiralonco</code></td>
</tr>
<tr class="odd">
<td></td>
<td><code>admiral.test</code></td>
<td></td>
<td><code>admiralroche</code></td>
</tr>
<tr class="even">
<td>Q2-2023</td>
<td>May 29th</td>
<td>June 5th</td>
<td>June 12th</td>
</tr>
<tr class="odd">
<td></td>
<td><code>admiraldev</code></td>
<td><code>admiral</code></td>
<td><code>admiralonco</code></td>
</tr>
<tr class="even">
<td></td>
<td><code>admiral.test</code></td>
<td></td>
<td><code>admiralroche</code></td>
</tr>
</tbody>
</table>
