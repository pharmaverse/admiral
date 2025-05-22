list(
  rdx_permitted_values = list(
    boolean = "`\"TRUE\"`, `\"FALSE\"`",
    char_scalar = "a character scalar, i.e., a character vector of length one",
    condition = "an unquoted condition, e.g., `AVISIT == \"BASELINE\"`",
    dataset = "a dataset, i.e., a `data.frame` or tibble",
    dataset_list = "named list of datasets, e.g., `list(adsl = adsl, ae = ae)`",
    date = "a date variable",
    event = "an `event()` or `event_joined()` object",
    expr_list = "list of expressions created by `exprs()`, e.g., `exprs(BASEC = \"MISSING\", BASE = -1)`",
    expr_list_formula = "list of named expressions created by a formula using `exprs()`, e.g., `exprs(AVALC = VSSTRESC, AVAL = yn_to_numeric(AVALC))`", # nolint
    expr_list_summary = "list of named expressions created by `exprs()`, e.g., `exprs(CUMDOSA = sum(AVAL, na.rm = TRUE), AVALU = \"ml\")`",
    flag_event = "a list of `flag_event()` objects",
    join = "`\"before\"`, `\"after\"`, `\"all\"`",
    merge_rel = "`\"one-to-one\"`, `\"many-to-one\"`",
    mode = "`\"first\"`, `\"last\"`",
    msg = "a console message to be printed, e.g. `\"Attention\"` or for longer messages use `paste(\"Line 1\", \"Line 2\")`",
    msg_type = "`\"none\"`, `\"message\"`, `\"warning\"`, `\"error\"`",
    source_list = "a list of source objects, e.g., `list(pd, death)`",
    var = "a variable",
    var_list = "list of variables created by `exprs()`, e.g., `exprs(USUBJID, VISIT)`"
  )
)
