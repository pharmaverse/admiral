# List all `tte_source` Objects Available in a Package

List all `tte_source` Objects Available in a Package

## Usage

``` r
list_tte_source_objects(package = "admiral")
```

## Arguments

- package:

  The name of the package in which to search for `tte_source` objects

  Default value

  :   `"admiral"`

## Value

A `data.frame` where each row corresponds to one `tte_source` object or
`NULL` if `package` does not contain any `tte_source` objects

## See also

Other Advanced Functions:
[`params()`](https:/pharmaverse.github.io/admiral/v1.4.1/reference/params.md)

## Examples

``` r
list_tte_source_objects()
#>              object dataset_name                                        filter
#> 1      ae_gr3_event         adae                TRTEMFL == "Y" & ATOXGR == "3"
#> 2       ae_wd_event         adae    TRTEMFL == "Y" & AEACN == "DRUG WITHDRAWN"
#> 3     ae_gr35_event         adae TRTEMFL == "Y" & ATOXGR %in% c("3", "4", "5")
#> 4  lastalive_censor         adsl                                          NULL
#> 5      ae_gr1_event         adae                TRTEMFL == "Y" & ATOXGR == "1"
#> 6      ae_ser_event         adae                 TRTEMFL == "Y" & AESER == "Y"
#> 7      ae_gr2_event         adae                TRTEMFL == "Y" & ATOXGR == "2"
#> 8          ae_event         adae                                TRTEMFL == "Y"
#> 9      ae_gr4_event         adae                TRTEMFL == "Y" & ATOXGR == "4"
#> 10     ae_gr5_event         adae                TRTEMFL == "Y" & ATOXGR == "5"
#> 11     ae_sev_event         adae            TRTEMFL == "Y" & AESEV == "SEVERE"
#> 12      death_event         adsl                                  DTHFL == "Y"
#>        date censor
#> 1     ASTDT      0
#> 2     ASTDT      0
#> 3     ASTDT      0
#> 4  LSTALVDT      1
#> 5     ASTDT      0
#> 6     ASTDT      0
#> 7     ASTDT      0
#> 8     ASTDT      0
#> 9     ASTDT      0
#> 10    ASTDT      0
#> 11    ASTDT      0
#> 12    DTHDT      0
#>                                                                                                 set_values_to
#> 1                     EVNTDESC: "GRADE 3 ADVERSE EVENT"<br>SRCDOM: "ADAE"<br>SRCVAR: "ASTDT"<br>SRCSEQ: AESEQ
#> 2  EVNTDESC: "ADVERSE EVENT LEADING TO DRUG WITHDRAWAL"<br>SRCDOM: "ADAE"<br>SRCVAR: "ASTDT"<br>SRCSEQ: AESEQ
#> 3                   EVNTDESC: "GRADE 3-5 ADVERSE EVENT"<br>SRCDOM: "ADAE"<br>SRCVAR: "ASTDT"<br>SRCSEQ: AESEQ
#> 4                                                   EVNTDESC: "ALIVE"<br>SRCDOM: "ADSL"<br>SRCVAR: "LSTALVDT"
#> 5                     EVNTDESC: "GRADE 1 ADVERSE EVENT"<br>SRCDOM: "ADAE"<br>SRCVAR: "ASTDT"<br>SRCSEQ: AESEQ
#> 6                     EVNTDESC: "SERIOUS ADVERSE EVENT"<br>SRCDOM: "ADAE"<br>SRCVAR: "ASTDT"<br>SRCSEQ: AESEQ
#> 7                     EVNTDESC: "GRADE 2 ADVERSE EVENT"<br>SRCDOM: "ADAE"<br>SRCVAR: "ASTDT"<br>SRCSEQ: AESEQ
#> 8                             EVNTDESC: "ADVERSE EVENT"<br>SRCDOM: "ADAE"<br>SRCVAR: "ASTDT"<br>SRCSEQ: AESEQ
#> 9                     EVNTDESC: "GRADE 4 ADVERSE EVENT"<br>SRCDOM: "ADAE"<br>SRCVAR: "ASTDT"<br>SRCSEQ: AESEQ
#> 10                    EVNTDESC: "GRADE 5 ADVERSE EVENT"<br>SRCDOM: "ADAE"<br>SRCVAR: "ASTDT"<br>SRCSEQ: AESEQ
#> 11                     EVNTDESC: "SEVERE ADVERSE EVENT"<br>SRCDOM: "ADAE"<br>SRCVAR: "ASTDT"<br>SRCSEQ: AESEQ
#> 12                                                     EVNTDESC: "DEATH"<br>SRCDOM: "ADSL"<br>SRCVAR: "DTHDT"
```
