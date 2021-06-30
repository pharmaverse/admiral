#' @keywords internal
#' @importFrom dplyr anti_join arrange bind_rows case_when desc ends_with full_join filter group_by
#'             if_else mutate mutate_at n pull rename row_number select slice
#'             starts_with transmute ungroup vars n_distinct union bind_rows distinct
#'             summarise_at summarise coalesce bind_cols na_if tibble
#' @importFrom magrittr %>%
#' @importFrom rlang := abort arg_match as_string enquo expr exprs inform
#'             quo_is_null sym syms warn is_symbol is_quosures is_integerish
#'             is_call eval_bare caller_env friendly_type new_formula f_lhs
#'             f_rhs expr_interp as_function quo call_name is_formula quo_squash
#'             is_bare_formula is_quosure call2 set_names enexpr quo_get_expr
#'             quo_is_symbol is_character eval_tidy quo_is_call quo_text quo_is_null
#'             .data
#' @importFrom utils capture.output
#' @importFrom purrr map map2 map_chr map_lgl reduce walk keep map_if transpose
#'             flatten every modify_at modify_if reduce compose
#' @importFrom stringr str_c str_detect str_extract str_remove str_remove_all str_trim
#'             str_to_upper str_glue
#' @importFrom assertthat assert_that is.number on_failure<-
#' @importFrom lubridate ceiling_date date days duration floor_date is.Date is.instant time_length
#'             ymd ymd_hms %--%
#' @importFrom tidyr spread
"_PACKAGE"
