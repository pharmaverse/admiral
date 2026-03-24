# Pre-Defined Time-to-Event Source Objects

These pre-defined `tte_source` objects can be used as input to
[`derive_param_tte()`](https:/pharmaverse.github.io/admiral/main/reference/derive_param_tte.md).

## Usage

``` r
death_event

lastalive_censor

ae_event

ae_ser_event

ae_gr1_event

ae_gr2_event

ae_gr3_event

ae_gr4_event

ae_gr5_event

ae_gr35_event

ae_sev_event

ae_wd_event
```

## Details

To see the definition of the various objects simply print the object in
the R console, e.g. `print(death_event)`. For details of how to use
these objects please refer to
[`derive_param_tte()`](https:/pharmaverse.github.io/admiral/main/reference/derive_param_tte.md).

## See also

[`derive_param_tte()`](https:/pharmaverse.github.io/admiral/main/reference/derive_param_tte.md),
[`tte_source()`](https:/pharmaverse.github.io/admiral/main/reference/tte_source.md),
[`event_source()`](https:/pharmaverse.github.io/admiral/main/reference/event_source.md),
[`censor_source()`](https:/pharmaverse.github.io/admiral/main/reference/censor_source.md)

Source Objects:
[`basket_select()`](https:/pharmaverse.github.io/admiral/main/reference/basket_select.md),
[`censor_source()`](https:/pharmaverse.github.io/admiral/main/reference/censor_source.md),
[`event()`](https:/pharmaverse.github.io/admiral/main/reference/event.md),
[`event_joined()`](https:/pharmaverse.github.io/admiral/main/reference/event_joined.md),
[`event_source()`](https:/pharmaverse.github.io/admiral/main/reference/event_source.md),
[`flag_event()`](https:/pharmaverse.github.io/admiral/main/reference/flag_event.md),
[`query()`](https:/pharmaverse.github.io/admiral/main/reference/query.md),
[`records_source()`](https:/pharmaverse.github.io/admiral/main/reference/records_source.md),
[`tte_source()`](https:/pharmaverse.github.io/admiral/main/reference/tte_source.md)

## Examples

``` r
# This shows the definition of all pre-defined `tte_source` objects that ship
# with {admiral}
for (obj in list_tte_source_objects()$object) {
  cat(obj, "\n")
  print(get(obj))
  cat("\n")
}
#> ae_gr3_event 
#> <event_source> object
#> dataset_name: "adae"
#> filter: TRTEMFL == "Y" & ATOXGR == "3"
#> date: ASTDT
#> censor: 0
#> set_values_to:
#>   EVNTDESC: "GRADE 3 ADVERSE EVENT"
#>   SRCDOM: "ADAE"
#>   SRCVAR: "ASTDT"
#>   SRCSEQ: AESEQ
#> order: NULL
#> 
#> ae_wd_event 
#> <event_source> object
#> dataset_name: "adae"
#> filter: TRTEMFL == "Y" & AEACN == "DRUG WITHDRAWN"
#> date: ASTDT
#> censor: 0
#> set_values_to:
#>   EVNTDESC: "ADVERSE EVENT LEADING TO DRUG WITHDRAWAL"
#>   SRCDOM: "ADAE"
#>   SRCVAR: "ASTDT"
#>   SRCSEQ: AESEQ
#> order: NULL
#> 
#> ae_gr35_event 
#> <event_source> object
#> dataset_name: "adae"
#> filter: TRTEMFL == "Y" & ATOXGR %in% c("3", "4", "5")
#> date: ASTDT
#> censor: 0
#> set_values_to:
#>   EVNTDESC: "GRADE 3-5 ADVERSE EVENT"
#>   SRCDOM: "ADAE"
#>   SRCVAR: "ASTDT"
#>   SRCSEQ: AESEQ
#> order: NULL
#> 
#> lastalive_censor 
#> <censor_source> object
#> dataset_name: "adsl"
#> filter: NULL
#> date: LSTALVDT
#> censor: 1
#> set_values_to:
#>   EVNTDESC: "ALIVE"
#>   SRCDOM: "ADSL"
#>   SRCVAR: "LSTALVDT"
#> order: NULL
#> 
#> ae_gr1_event 
#> <event_source> object
#> dataset_name: "adae"
#> filter: TRTEMFL == "Y" & ATOXGR == "1"
#> date: ASTDT
#> censor: 0
#> set_values_to:
#>   EVNTDESC: "GRADE 1 ADVERSE EVENT"
#>   SRCDOM: "ADAE"
#>   SRCVAR: "ASTDT"
#>   SRCSEQ: AESEQ
#> order: NULL
#> 
#> ae_ser_event 
#> <event_source> object
#> dataset_name: "adae"
#> filter: TRTEMFL == "Y" & AESER == "Y"
#> date: ASTDT
#> censor: 0
#> set_values_to:
#>   EVNTDESC: "SERIOUS ADVERSE EVENT"
#>   SRCDOM: "ADAE"
#>   SRCVAR: "ASTDT"
#>   SRCSEQ: AESEQ
#> order: NULL
#> 
#> ae_gr2_event 
#> <event_source> object
#> dataset_name: "adae"
#> filter: TRTEMFL == "Y" & ATOXGR == "2"
#> date: ASTDT
#> censor: 0
#> set_values_to:
#>   EVNTDESC: "GRADE 2 ADVERSE EVENT"
#>   SRCDOM: "ADAE"
#>   SRCVAR: "ASTDT"
#>   SRCSEQ: AESEQ
#> order: NULL
#> 
#> ae_event 
#> <event_source> object
#> dataset_name: "adae"
#> filter: TRTEMFL == "Y"
#> date: ASTDT
#> censor: 0
#> set_values_to:
#>   EVNTDESC: "ADVERSE EVENT"
#>   SRCDOM: "ADAE"
#>   SRCVAR: "ASTDT"
#>   SRCSEQ: AESEQ
#> order: NULL
#> 
#> ae_gr4_event 
#> <event_source> object
#> dataset_name: "adae"
#> filter: TRTEMFL == "Y" & ATOXGR == "4"
#> date: ASTDT
#> censor: 0
#> set_values_to:
#>   EVNTDESC: "GRADE 4 ADVERSE EVENT"
#>   SRCDOM: "ADAE"
#>   SRCVAR: "ASTDT"
#>   SRCSEQ: AESEQ
#> order: NULL
#> 
#> ae_gr5_event 
#> <event_source> object
#> dataset_name: "adae"
#> filter: TRTEMFL == "Y" & ATOXGR == "5"
#> date: ASTDT
#> censor: 0
#> set_values_to:
#>   EVNTDESC: "GRADE 5 ADVERSE EVENT"
#>   SRCDOM: "ADAE"
#>   SRCVAR: "ASTDT"
#>   SRCSEQ: AESEQ
#> order: NULL
#> 
#> ae_sev_event 
#> <event_source> object
#> dataset_name: "adae"
#> filter: TRTEMFL == "Y" & AESEV == "SEVERE"
#> date: ASTDT
#> censor: 0
#> set_values_to:
#>   EVNTDESC: "SEVERE ADVERSE EVENT"
#>   SRCDOM: "ADAE"
#>   SRCVAR: "ASTDT"
#>   SRCSEQ: AESEQ
#> order: NULL
#> 
#> death_event 
#> <event_source> object
#> dataset_name: "adsl"
#> filter: DTHFL == "Y"
#> date: DTHDT
#> censor: 0
#> set_values_to:
#>   EVNTDESC: "DEATH"
#>   SRCDOM: "ADSL"
#>   SRCVAR: "DTHDT"
#> order: NULL
#> 
```
