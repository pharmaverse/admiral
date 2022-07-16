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

QStandardPaths: XDG\_RUNTIME\_DIR not set, defaulting to
‘/tmp/runtime-r590548’ TypeError: Attempting to change the setter of an
unconfigurable property. TypeError: Attempting to change the setter of
an unconfigurable property.
![](README_files/figure-markdown_strict/unnamed-chunk-2-1.png)
