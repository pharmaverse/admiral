list(
  rdx_permitted_values = list(
    boolean = "`TRUE`, `FALSE`",
    char_scalar = "a character scalar, i.e., a character vector of length one",
    condition = "an unquoted condition, e.g., `AVISIT == \"BASELINE\"`",
    dataset = "a dataset, i.e., a `data.frame` or tibble",
    dataset_list = "named list of datasets, e.g., `list(adsl = adsl, ae = ae)`",
    date = "a date or datetime variable",
    date_chr = "a character date variable",
    date_chr_vector = "a character date vector",
    date_imp = "`\"first\"`, `\"mid\"`, `\"last\"`, or user-defined",
    date_high_imp = "`\"Y\"` (year, highest level), `\"M\"` (month), `\"D\"`
    (day), `\"n\"` (none, lowest level)",
    date_list = "a list of dates, e.g. `list(ymd_hms(\"2021-07-01T04:03:01\"), ymd_hms(\"2022-05-12T13:57:23\"))`",
    date_flag_imp = "`\"auto\"`, `\"date\"` or `\"none\"`",
    date_time_flag_imp = "`\"auto\"`, `\"date\"`,`\"time\"`, `\"both\"` or `\"none\"`",
    date_time_high_imp = "`\"Y\"` (year, highest level), `\"M\"` (month), `\"D\"`
    (day), `\"h\"` (hour), `\"m\"` (minute), `\"s\"` (second), `\"n\"` (none, lowest
    level)",
    event = "an `event()` or `event_joined()` object",
    expr_list_formula = "list of named expressions created by a formula using `exprs()`, e.g., `exprs(AVALC = VSSTRESC, AVAL = yn_to_numeric(AVALC))`",
    expr_list_summary = "list of named expressions created by `exprs()`, e.g., `exprs(CUMDOSA = sum(AVAL, na.rm = TRUE), AVALU = \"ml\")`",
    flag_event = "a list of `flag_event()` objects",
    join_type = "`\"before\"`, `\"after\"`, `\"all\"`",
    merge_rel = "`\"one-to-one\"`, `\"many-to-one\"`",
    mode = "`\"first\"`, `\"last\"`",
    msg = "a console message to be printed, e.g. `\"Attention\"` or for longer messages use `paste(\"Line 1\", \"Line 2\")`",
    msg_type = "`\"none\"`, `\"message\"`, `\"warning\"`, `\"error\"`",
    order_optional = "list of expressions created by `exprs()`, e.g., `exprs(ADT, desc(AVAL))` or `NULL`",
    pos_int = "a positive integer, e.g. `2` or `5`",
    source_list = "a list of source objects, e.g., `list(pd, death)`",
    time_imp = "`\"first\"`, `\"last\"`, or user-defined",
    var = "an unquoted symbol, e.g., `AVAL`",
    var_list = "list of variables created by `exprs()`, e.g., `exprs(USUBJID, VISIT)`",
    var_list_rename = "list of (optionally named) variables created by `exprs()`, e.g., `exprs(USUBJID, ADY = ASTDY)`",
    var_list_tidyselect = "list of variables or tidyselect expressions created by `exprs()`, e.g., `exprs(DTHDT, starts_with(\"AST\"))` or `exprs(everything)`",
    var_expr_list = "list of variables or named expressions created by `exprs()`, e.g., `exprs(EXSTDY, EXSTDTM = convert_dtc_to_dtm(EXSTDTC))`"
  )
)
