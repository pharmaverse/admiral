# roxygen_param_by_vars Test 1: Text variations

    Code
      roxygen_param_by_vars()
    Output
      [1] ""

---

    Code
      roxygen_param_by_vars(rename = TRUE)
    Output
      [1] "Variables can be renamed by naming the element, i.e. \n`by_vars = exprs(<name in input dataset> = <name in additional dataset>)`, similar to the `dplyr` joins."

