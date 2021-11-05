#' @keywords internal
#' @import cdiscpilot
#' @importFrom dplyr arrange bind_rows case_when desc ends_with filter full_join group_by
#'             if_else mutate mutate_at mutate_if n pull rename rename_at row_number select slice
#'             starts_with transmute ungroup vars n_distinct union distinct
#'             summarise_at summarise coalesce bind_cols na_if tibble
#' @importFrom magrittr %>%
#' @importFrom rlang := abort arg_match as_function as_string call2 caller_env
#'             call_name current_env .data enexpr enquo eval_bare eval_tidy expr
#'             expr_interp expr_label exprs f_lhs f_rhs friendly_type inform
#'             is_bare_formula is_call is_character is_formula is_integerish
#'             is_quosure is_quosures is_symbol new_formula parse_exprs quo
#'             quo_get_expr quo_is_call quo_is_missing quo_is_null quo_is_symbol
#'             quos quo_squash quo_text set_names sym syms type_of warn
#'             quo_set_env quo_get_env
#' @importFrom utils capture.output
#' @importFrom purrr map map2 map_chr map_lgl reduce walk keep map_if transpose
#'             flatten every modify_at modify_if reduce compose
#' @importFrom stringr str_c str_detect str_extract str_remove str_remove_all
#'             str_replace str_trim str_to_upper str_glue
#' @importFrom assertthat assert_that is.number on_failure<-
#' @importFrom lubridate as_datetime ceiling_date date days duration floor_date is.Date is.instant
#'             time_length %--% ymd ymd_hms
#' @importFrom tidyr spread gather drop_na
#' @importFrom tidyselect contains vars_select
"_PACKAGE"
