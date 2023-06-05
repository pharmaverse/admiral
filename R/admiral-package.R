#' @keywords internal
#' @family internal
#' @import admiraldev
#' @importFrom dplyr across arrange bind_rows case_when desc ends_with
#'   everything filter full_join group_by if_else mutate n pull rename
#'   rename_with row_number select slice semi_join starts_with transmute ungroup
#'   n_distinct union distinct summarise coalesce bind_cols na_if tibble tribble
#'   summarise_all group_by_at first
#' @importFrom magrittr %>%
#' @importFrom rlang := abort arg_match as_function as_label as_name as_string
#'   call2 caller_env call_name current_env .data enexpr eval_bare eval_tidy
#'   expr expr_interp expr_label exprs f_lhs f_rhs inform is_call is_expression
#'   is_missing new_formula parse_expr parse_exprs set_names sym syms type_of
#'   warn as_data_mask list2
#' @importFrom utils capture.output str file.edit
#' @importFrom purrr map map2 map_chr map_lgl reduce walk keep map_if transpose
#'             flatten every modify_at modify_if reduce compose pmap
#' @importFrom stringr str_c str_count str_detect str_extract str_glue
#'   str_length str_locate str_locate_all str_match str_remove str_remove_all
#'   str_replace str_replace_all str_split str_starts str_sub str_subset
#'   str_trim str_to_lower str_to_title str_to_upper
#' @importFrom lubridate as_datetime ceiling_date date days duration floor_date is.Date is.instant
#'             rollback time_length %--% ymd ymd_hms weeks years hours minutes is.POSIXct
#' @importFrom tidyr crossing drop_na fill nest pivot_longer pivot_wider unnest
#' @importFrom tidyselect all_of any_of contains matches vars_select
#' @importFrom hms as_hms
#' @importFrom lifecycle deprecate_warn deprecated deprecate_stop
"_PACKAGE"
