# roxygen_param_by_vars Test 1: Text variations

    Code
      roxygen_param_by_vars()
    Output
      [1] "*Permitted Values*: list of variables created by exprs() \n \ne.g. exprs(USUBJID, VISIT)"

---

    Code
      roxygen_param_by_vars(additional_dataset = "additional_dataset", rename = TRUE)
    Output
      [1] "Variables from `additional_dataset` can be renamed by naming the element, I.e. \nby_vars = exprs(<name in input dataset> = <name in additional dataset>),similar to the dplyr joins.\n \n*Permitted Values*: list of variables created by exprs() \n \ne.g. exprs(USUBJID, VISIT)"

---

    Code
      roxygen_param_by_vars(additional_dataset = "additional_dataset", unique = TRUE)
    Output
      [1] "Variables must be a unique key of the selected observations in `additional_dataset`. \n \n*Permitted Values*: list of variables created by exprs() \n \ne.g. exprs(USUBJID, VISIT)"

---

    Code
      roxygen_param_by_vars(additional_dataset = "additional_dataset", unique = TRUE,
        rename = TRUE)
    Output
      [1] "Variables must be a unique key of the selected observations in `additional_dataset`. \n \nVariables from `additional_dataset` can be renamed by naming the element, I.e. \nby_vars = exprs(<name in input dataset> = <name in additional dataset>),similar to the dplyr joins.\n \n*Permitted Values*: list of variables created by exprs() \n \ne.g. exprs(USUBJID, VISIT)"

