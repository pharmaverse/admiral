# Output Boilerplate for the "Add Labels and Attributes" Vignette Section

Outputs the standard markdown text used in ADaM-specific vignettes for
the "Add Labels and Attributes" section. This function is intended to be
called inside an R Markdown code chunk with `results='asis'` and
`echo=FALSE`.

## Usage

``` r
admiral_add_labels_attrs_section(header_lvl = "##")
```

## Arguments

- header_lvl:

  The markdown header level for the section heading. Must be a character
  string consisting only of hash marks, e.g. `"#"`, `"##"`, or `"###"`.
  The default is `"##"`.

  Default value

  :   `"##"`

## Value

No return value. The function outputs text directly to the console or
output stream via [`cat()`](https://rdrr.io/r/base/cat.html), intended
for use in R Markdown documents with `results='asis'`.

## Details

The outputted section describes how to add variable labels and other
metadata to ADaM datasets as a final step in the derivation process,
using the `{metacore}`, `{metatools}`, and `{xportr}` packages.

## See also

Utilities used for examples and template scripts:
[`list_all_templates()`](https:/pharmaverse.github.io/admiral/copilot/feature-boiler-plate-language-function/reference/list_all_templates.md),
[`use_ad_template()`](https:/pharmaverse.github.io/admiral/copilot/feature-boiler-plate-language-function/reference/use_ad_template.md)

## Examples

``` r
admiral_add_labels_attrs_section()
#> ## Add Labels and Attributes {#attributes}
#> 
#> Note that attributes may not be preserved in some cases after processing
#> with `{admiral}`. The recommended approach is to apply variable labels
#> and other metadata as a final step in your data derivation process using
#> packages like:
#> 
#> -   [metacore](https://atorus-research.github.io/metacore/): establish a
#>     common foundation for the use of metadata within an R session.
#> 
#> -   [metatools](https://pharmaverse.github.io/metatools/): enable the
#>     use of metacore objects. Metatools can be used to build datasets or
#>     enhance columns in existing datasets as well as checking datasets
#>     against the metadata.
#> 
#> -   [xportr](https://atorus-research.github.io/xportr/): functionality
#>     to associate all metadata information to a local R data frame,
#>     perform data set level validation checks and convert into a
#>     [transport v5
#>     file(xpt)](https://documentation.sas.com/doc/en/pgmsascdc/9.4_3.5/movefile/n1xbwdre0giahfn11c99yjkpi2yb.htm).
#> 
#> NOTE: Together with `{admiral}` these packages comprise an End to End
#> pipeline under the umbrella of the
#> [pharmaverse](https://github.com/pharmaverse). An example of applying
#> metadata and performing associated checks can be found at the [pharmaverse
#> E2E example](https://pharmaverse.github.io/examples/adam/adsl).

admiral_add_labels_attrs_section(header_lvl = "#")
#> # Add Labels and Attributes {#attributes}
#> 
#> Note that attributes may not be preserved in some cases after processing
#> with `{admiral}`. The recommended approach is to apply variable labels
#> and other metadata as a final step in your data derivation process using
#> packages like:
#> 
#> -   [metacore](https://atorus-research.github.io/metacore/): establish a
#>     common foundation for the use of metadata within an R session.
#> 
#> -   [metatools](https://pharmaverse.github.io/metatools/): enable the
#>     use of metacore objects. Metatools can be used to build datasets or
#>     enhance columns in existing datasets as well as checking datasets
#>     against the metadata.
#> 
#> -   [xportr](https://atorus-research.github.io/xportr/): functionality
#>     to associate all metadata information to a local R data frame,
#>     perform data set level validation checks and convert into a
#>     [transport v5
#>     file(xpt)](https://documentation.sas.com/doc/en/pgmsascdc/9.4_3.5/movefile/n1xbwdre0giahfn11c99yjkpi2yb.htm).
#> 
#> NOTE: Together with `{admiral}` these packages comprise an End to End
#> pipeline under the umbrella of the
#> [pharmaverse](https://github.com/pharmaverse). An example of applying
#> metadata and performing associated checks can be found at the [pharmaverse
#> E2E example](https://pharmaverse.github.io/examples/adam/adsl).
```
