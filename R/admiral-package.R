#' @keywords internal
#' @family internal
#' @importFrom dplyr arrange bind_rows case_when desc ends_with filter full_join group_by
#'             if_else mutate mutate_at mutate_if n pull rename rename_at row_number select slice
#'             semi_join starts_with transmute ungroup vars n_distinct union distinct
#'             summarise_at summarise coalesce bind_cols na_if tibble
#' @importFrom magrittr %>%
#' @importFrom rlang := abort arg_match as_function as_label as_string call2
#'   caller_env call_name current_env .data enexpr enquo eval_bare eval_tidy
#'   expr expr_interp expr_label f_lhs f_rhs inform is_bare_formula is_call
#'   is_character is_formula is_integerish is_logical is_quosure is_quosures
#'   is_symbol new_formula parse_expr parse_exprs quo quo_get_expr quo_is_call
#'   quo_is_missing quo_is_null quo_is_symbol quos quo_squash quo_text set_names
#'   sym syms type_of warn quo_set_env quo_get_env
#' @importFrom utils capture.output str
#' @importFrom purrr map map2 map_chr map_lgl reduce walk keep map_if transpose
#'             flatten every modify_at modify_if reduce compose
#' @importFrom stringr str_c str_detect str_extract str_glue str_match
#'   str_remove str_remove_all str_replace str_sub str_subset str_trim
#'   str_to_lower str_to_title str_to_upper str_length str_locate
#' @importFrom assertthat assert_that is.number on_failure<-
#' @importFrom lubridate as_datetime ceiling_date date days duration floor_date is.Date is.instant
#'             rollback time_length %--% ymd ymd_hms weeks years hours minutes
#' @importFrom tidyr drop_na nest pivot_longer pivot_wider unnest
#' @importFrom tidyselect all_of contains vars_select
#' @importFrom hms as_hms
#' @importFrom lifecycle deprecate_warn deprecated deprecate_stop
#' @importFrom readxl read_excel
#' @importFrom admiraldev left_join assert_character_scalar assert_character_vector assert_data_frame
#'             assert_expr assert_filter_cond assert_function assert_function_param
#'             assert_has_variables assert_integer_scalar assert_list_element assert_list_of
#'             assert_logical_scalar assert_named_exprs assert_numeric_vector assert_one_to_one assert_order_vars
#'             assert_param_does_not_exist assert_s3_class assert_symbol assert_unit assert_vars assert_varval_list
#'             get_constant_vars get_duplicates get_source_vars is_auto is_date is_named is_order_vars is_timeunit
#'             is_valid_date_entry is_valid_day is_valid_dtc is_valid_hour is_valid_sec_min is_valid_month is_valid_time_entry what_is_it
#'             suppress_warning warn_if_incomplete_dtc warn_if_inconsistent_list
#'             warn_if_invalid_dtc warn_if_vars_exist expect_dfs_equal backquote enumerate quo_c squote dquote
#'             dataset_vignette arg_name convert_dtm_to_dtc extract_vars `%notin%` replace_values_by_names
"_PACKAGE"
