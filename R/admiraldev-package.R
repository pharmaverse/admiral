#' @keywords internal
#' @importFrom dplyr across arrange bind_rows case_when desc ends_with filter
#'   full_join group_by if_else mutate n pull rename row_number select slice
#'   starts_with transmute ungroup n_distinct union distinct summarise coalesce
#'   bind_cols na_if tibble
#' @importFrom magrittr %>%
#' @importFrom rlang := abort arg_match as_function as_label as_name as_string
#'   call2 caller_env call_name current_env .data enexpr enquo eval_bare
#'   eval_tidy expr expr_interp expr_label exprs f_lhs f_rhs inform missing_arg
#'   is_bare_formula is_call is_character is_expression is_formula is_integerish
#'   is_logical is_missing is_quosure is_symbol is_symbolic new_formula
#'   parse_expr parse_exprs quo quo_is_missing quo_is_null set_names sym syms
#'   type_of warn
#' @importFrom utils capture.output str
#' @importFrom purrr map map2 map_chr map_lgl reduce walk keep map_if transpose
#'             flatten every modify_at modify_if reduce compose
#' @importFrom stringr str_c str_detect str_extract str_glue str_match
#'             str_remove str_remove_all str_replace str_trim str_to_lower str_subset
#'             str_to_title str_to_upper
#' @importFrom lubridate as_datetime ceiling_date date days duration floor_date is.Date is.instant
#'             time_length %--% ymd ymd_hms weeks years hours minutes
#' @importFrom tidyr drop_na nest pivot_longer pivot_wider unnest
#' @importFrom tidyselect all_of contains vars_select
#' @importFrom hms as_hms
#' @importFrom lifecycle deprecate_warn deprecated deprecate_stop
"_PACKAGE"
