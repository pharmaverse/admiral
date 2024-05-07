# pharmaversesdtm 0.2.0

## New Features

- Following Vaccine Specific SDTM datasets have been added. (#4)

    - `ce_vaccine`
    - `dm_vaccine`
    - `ex_vaccine`
    - `face_vaccine`
    - `is_vaccine`
    - `vs_vaccine`
    - `suppce_vaccine`
    - `suppdm_vaccine`
    - `suppex_vaccine`
    - `suppface_vaccine`
 
- Oncology response data for iRECIST criteria (`rs_onco_irecist`) was added. (#32)

- `get_terms()` now expects `TERMCHAR` instead of `TERMNAME` in alignment with [this](https://github.com/pharmaverse/admiral/issues/2186) `{admiral}` issue. (#76)

# pharmaversesdtm 0.1.1

## Documentation

 - Fixed redirected links on website for CRAN release. 

# pharmaversesdtm 0.1.0

## New Features

 - Ophthalmology variants of `ex` and `qs` SDTM datasets added. (#15)
 - Migrate data and function `get_terms()` from `admiral.test`. (#1, #49)
 - Oncology datasets `tu_onco_recist`, `tr_onco_recist`, and `rs_onco_recist`
 using RECIST 1.1 response criteria. The datasets contain just a few patients.
 They are intended for vignettes and examples of ADaM datasets creation. (#33)

