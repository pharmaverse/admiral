list(
  rdx_permitted_values = list(
    char_scalar = "a character scalar, i.e., a character vector of length one",
    condition = "an unquoted condition, e.g., `AVISIT == \"BASELINE\"`",
    dataset = "a dataset, i.e., a `data.frame` or tibble",
    dataset_list = "named list of datasets, e.g., `list(adsl = adsl, ae = ae)`",
    date = "a date variable",
    mode = "`\"first\"`, `\"last\"`",
    msg_type = "`\"none\"`, `\"message\"`, `\"warning\"`, `\"error\"`",
    source_list = "a list of source objects, e.g., `list(pd, death)`",
    symbol = "an unquoted symbol, e.g., `AVAL`",
    var_list = "list of variables created by `exprs()`, e.g., `exprs(USUBJID, VISIT)`"
  )
)
