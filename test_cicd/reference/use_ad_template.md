# Open an ADaM Template Script

Open an ADaM Template Script

## Usage

``` r
use_ad_template(
  adam_name = "adsl",
  save_path = paste0("./", adam_name, ".R"),
  package = "admiral",
  overwrite = FALSE,
  open = interactive()
)
```

## Arguments

- adam_name:

  An ADaM dataset name. You can use any of the available dataset names
  `"ADAB"`, `"ADAE"`, `"ADCM"`, `"ADEG"`, `"ADEX"`, `"ADLB"`,
  `"ADLBHY"`, `"ADMH"`, `"ADPC"`, `"ADPP"`, `"ADPPK"`, `"ADSL"`,
  `"ADVS"`. The dataset name is case-insensitive. The default dataset
  name is `"ADSL"`.

  Default value

  :   `"adsl"`

- save_path:

  Path to save the script.

  Default value

  :   `paste0("./", adam_name, ".R")`

- package:

  The R package in which to look for templates. By default `"admiral"`.

  Default value

  :   `"admiral"`

- overwrite:

  Whether to overwrite an existing file named `save_path`.

  Default value

  :   `FALSE`

- open:

  Whether to open the script right away.

  Default value

  :   [`interactive()`](https://rdrr.io/r/base/interactive.html)

## Value

No return values, called for side effects

## Details

Running without any arguments such as `use_ad_template()` auto-generates
`adsl.R` in the current path. Use
[`list_all_templates()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/list_all_templates.md)
to discover which templates are available.

## See also

Utilities used for examples and template scripts:
[`list_all_templates()`](https:/pharmaverse.github.io/admiral/test_cicd/reference/list_all_templates.md)

## Examples

``` r
if (interactive()) {
  use_ad_template("adsl")
}
```
