#' @keywords internal
#' @importFrom dplyr anti_join arrange bind_rows case_when desc ends_with full_join filter group_by
#'             if_else mutate n left_join pull rename row_number select slice starts_with
#'             transmute ungroup vars
#' @importFrom magrittr %>%
#' @importFrom rlang := abort arg_match as_string enquo exprs inform is_call
#'             is_symbol quo_get_expr quo_is_null quo_is_symbol sym syms warn
#' @importFrom utils capture.output
#' @importFrom purrr map map2 map_chr map_lgl reduce walk
#' @importFrom stringr str_c str_detect str_remove str_to_upper str_trim
#' @importFrom assertthat assert_that is.number on_failure<-
#' @importFrom lubridate ceiling_date days duration floor_date time_length %--%
#'             ymd ymd_hms is.instant
"_PACKAGE"
