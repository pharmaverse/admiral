list(
  rdx_permitted_values = list(
    char_scalar = "a character scalar, i.e., a character vector of length one",
    condition = "an unquoted condition, e.g., `AVISIT == \"BASELINE\"`",
    dataset = "a dataset, i.e., a `data.frame` or tibble",
    dataset_list = "named list of datasets, e.g., `list(adsl = adsl, ae = ae)`",
    date = "a date variable",
    expr_list_summary = "list of named expressions created by `exprs()`, e.g., `exprs(CUMDOSA = sum(AVAL, na.rm = TRUE), AVALU = \"ml\")`",
    mode = "`\"first\"`, `\"last\"`",
    msg_type = "`\"none\"`, `\"message\"`, `\"warning\"`, `\"error\"`",
    source_list = "a list of source objects, e.g., `list(pd, death)`",
    var_list = "list of variables created by `exprs()`, e.g., `exprs(USUBJID, VISIT)`",
    expr_list_formula = "list of named expressions created by a formula using `exprs()`, e.g., `exprs(AVALC = VSSTRESC, AVAL = yn_to_numeric(AVALC))`" #nolint
  )
)
