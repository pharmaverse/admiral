#' @keywords internal
#' @family internal
#' @import admiraldev
#' @importFrom cli cli_abort ansi_collapse cli_div cli_inform cli_text cli_warn
#' @importFrom dplyr across arrange bind_cols bind_rows case_when coalesce
#'             desc distinct ends_with everything filter first full_join
#'             group_by group_by_at if_else mutate n n_distinct na_if pull
#'             rename rename_with row_number select semi_join slice starts_with
#'             summarise summarise_all tibble tribble ungroup union lag
#' @importFrom hms as_hms
#' @importFrom lifecycle deprecate_warn deprecate_stop deprecated
#' @importFrom lubridate %--% as_datetime ceiling_date date days duration
#'             floor_date hour hours is.Date is.instant is.POSIXct minute
#'             minutes rollback second time_length weeks ymd ymd_hms years
#' @importFrom magrittr %>%
#' @importFrom purrr compose discard every flatten keep map map_chr map_dbl
#'             map_if map_lgl map2 modify_at modify_if pmap reduce transpose
#'             walk
#' @importFrom rlang := abort arg_match as_data_mask as_function as_label
#'   as_name as_string call2 call_name caller_env current_env .data enexpr
#'   enexprs eval_bare eval_tidy expr expr_interp exec expr_label exprs f_lhs
#'   f_rhs  inform is_call is_expression is_missing is_named list2
#'   new_environment new_formula parse_expr parse_exprs set_names sym syms
#'   type_of warn
#' @importFrom stats setNames
#' @importFrom stringr str_c str_count str_detect str_extract str_glue
#'             str_length str_locate str_locate_all str_match str_remove
#'             str_remove_all str_replace str_replace_all str_split str_starts
#'             str_sub str_subset str_to_lower str_to_title str_to_upper
#'             str_trim
#' @importFrom tidyr crossing drop_na fill nest pivot_longer pivot_wider unnest
#' @importFrom tidyselect all_of any_of contains matches vars_select
#' @importFrom utils capture.output file.edit str
#'
"_PACKAGE"
