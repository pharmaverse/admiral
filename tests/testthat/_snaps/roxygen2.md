# roxygen_param_by_vars Test 1: Text variations

    Code
      roxygen_param_by_vars()
    Output
      [1] "*Permitted Values*: list of variables created by `exprs()` \ne.g. `exprs(USUBJID, VISIT)`"

---

    Code
      roxygen_param_by_vars(rename = TRUE)
    Output
      [1] "Variables can be renamed by naming the element, i.e. \n`by_vars = exprs(<name in input dataset> = <name in additional dataset>)`, similar to the `dplyr` joins.\n \n*Permitted Values*: list of variables created by `exprs()` \ne.g. `exprs(USUBJID, VISIT)`"

