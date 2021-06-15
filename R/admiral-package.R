#' @keywords internal
#' @importFrom dplyr anti_join arrange bind_rows case_when desc distinct ends_with full_join filter
#'             group_by if_else inner_join mutate mutate_at n left_join pull rename row_number
#'             select slice starts_with summarise transmute ungroup vars
#' @importFrom magrittr %>%
#' @importFrom rlang := .data abort arg_match as_string enquo exprs inform is_call
#'             is_symbol is_quosure is_quosures quo_get_expr quo_text quo_is_call
#'             quo_is_null quo_is_symbol sym syms warn
#' @importFrom utils capture.output
#' @importFrom purrr map map2 map_chr map_lgl reduce walk
#' @importFrom stringr str_c str_detect str_remove str_to_upper str_trim
#' @importFrom assertthat assert_that is.number on_failure<-
#' @importFrom lubridate ceiling_date date days duration floor_date is.Date is.instant time_length %--%
#'             ymd ymd_hms
"_PACKAGE"
