# derive_vars_computed Test 2: no new variables added if filtered dataset is empty

    Code
      result <- derive_vars_computed(dataset = adsl, dataset_add = advs, by_vars = exprs(
        USUBJID), parameters = c("WEIGHT"), constant_by_vars = exprs(USUBJID),
      constant_parameters = c("HEIGHT"), new_vars = exprs(BMIBL = compute_bmi(height = AVAL.HEIGHT,
        weight = AVAL.WEIGHT)), filter_add = ABLFL == "")
    Condition
      Warning:
      The input dataset does not contain any observations fulfilling the filter condition (`ABLFL == ""`) for the parameter codes (`PARAMCD`): WEIGHT.
      i No new observations were added.

# derive_vars_computed Test 3: no new variables are added if a parameter is missing

    Code
      result <- derive_vars_computed(dataset = adsl, dataset_add = advs, by_vars = exprs(
        STUDYID, USUBJID), parameters = c("WEIGHT"), constant_by_vars = exprs(STUDYID,
        USUBJID), constant_parameters = c("HEIGHT"), new_vars = exprs(BMIBL = compute_bmi(
        height = AVAL.HEIGHT, weight = AVAL.WEIGHT)), filter_add = ABLFL == "Y")
    Condition
      Warning:
      The input dataset does not contain any observations fulfilling the filter condition (`ABLFL == "Y"`) for the parameter codes (`PARAMCD`): HEIGHT.
      i No new observations were added.

