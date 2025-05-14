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
    order_optional = "list of expressions created by `exprs()`, e.g., `exprs(ADT, desc(AVAL))` or `NULL`",
    source_list = "a list of source objects, e.g., `list(pd, death)`",
    symbol = "an unquoted symbol, e.g., `AVAL`",
    var_list = "list of variables created by `exprs()`, e.g., `exprs(USUBJID, VISIT)`",
    var_list_rename = "list of (optionally named) variables created by `exprs()`, e.g., `exprs(USUBJID, ADY = ASTDY)`",
    var_expr_list = "list of variables or named expressions created by `exprs()`, e.g., `exprs(EXSTDY, EXSTDTM = convert_dtc_to_dtm(EXSTDTC))`"
  )
)
