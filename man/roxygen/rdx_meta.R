list(
  rdx_permitted_values = list(
    char_scalar = "a character scalar, i.e., a character vector of length one",
    condition = "an unquoted condition, e.g., `AVISIT == \"BASELINE\"`",
    dataset = "a dataset, i.e., a `data.frame` or tibble",
    mode = "`\"first\"`, `\"last\"`",
    msg_type = "`\"none\"`, `\"message\"`, `\"warning\"`, `\"error\"`",
    var_list = "list of variables created by `exprs()`, e.g., `exprs(USUBJID, VISIT)`"
  )
)
