#' @keywords internal
#' @importFrom dplyr anti_join arrange bind_rows case_when desc ends_with full_join filter group_by
#'             if_else mutate n left_join pull rename row_number select slice starts_with
#'             transmute ungroup vars n_distinct union bind_rows distinct
#'             summarise_at summarise
#' @importFrom magrittr %>%
#' @importFrom rlang := abort arg_match as_string enquo expr exprs inform quo_is_null sym syms warn
#'             is_symbol is_quosures is_integerish is_call eval_bare caller_env
#'             friendly_type new_formula f_lhs f_rhs expr_interp as_function
#'             quo call_name is_formula quo_squash is_bare_formula is_quosure
#'             call2 set_names enexpr
#' @importFrom utils capture.output
#' @importFrom purrr map map2 map_chr map_lgl reduce keep map_if transpose
#' @importFrom stringr str_c str_detect str_remove str_trim str_glue
#' @importFrom assertthat assert_that is.number on_failure<-
#' @importFrom lubridate ceiling_date days duration floor_date time_length %--%
#'             ymd ymd_hms is.instant
"_PACKAGE"
