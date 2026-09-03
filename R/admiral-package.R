#' @keywords internal
#' @family internal
#' @import admiraldev
#' @importFrom cli
#'   ansi_collapse cli_abort cli_div cli_inform cli_text cli_warn qty
#'
#' @importFrom dplyr
#'   across anti_join arrange between bind_cols bind_rows case_when coalesce
#'   cross_join desc distinct ends_with everything filter filter_out first
#'   full_join group_by group_by_at group_split if_else inner_join lag
#'   left_join mutate n n_distinct na_if pull recode_values reframe rename
#'   rename_with replace_when row_number select semi_join slice starts_with
#'   summarise summarise_all tibble tribble ungroup union
#'
#' @importFrom hms as_hms
#'
#' @importFrom lifecycle deprecate_warn deprecate_stop deprecated
#'
#' @importFrom lubridate
#'   %--% as_datetime ceiling_date date days duration floor_date hour hours
#'   is.Date is.instant is.POSIXct minute minutes rollback second time_length
#'   weeks ymd ymd_hms years
#'
#' @importFrom magrittr %>%
#'
#' @importFrom purrr
#'   discard every flatten keep map map_chr map_dbl map_if map_lgl map2 pmap
#'   reduce walk
#'
#' @importFrom rlang
#'   := arg_match as_data_mask as_label as_name as_string call2 caller_env
#'   cnd_muffle cnd_signal current_env .data enexpr enexprs eval_tidy expr
#'   exec expr_label exprs inform is_call is_expression is_missing is_named
#'   list2 new_environment parse_expr parse_exprs set_names sym syms zap
#'
#' @importFrom stats setNames
#'
#' @importFrom stringr
#'   fixed str_c str_count str_detect str_extract str_glue str_length
#'   str_locate str_locate_all str_match str_remove str_remove_all str_replace
#'   str_replace_all str_split str_starts str_sub str_subset str_to_lower
#'   str_to_title str_to_upper str_trim regex
#'
#' @importFrom tidyr crossing drop_na fill nest pivot_longer pivot_wider unnest
#'
#' @importFrom tidyselect all_of any_of contains matches
#'
#' @importFrom utils capture.output object.size str
#'
"_PACKAGE"
