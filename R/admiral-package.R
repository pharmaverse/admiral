#' @keywords internal
#' @importFrom dplyr anti_join arrange bind_rows case_when desc ends_with filter group_by if_else
#'             mutate n left_join pull rename row_number select slice starts_with
#'             transmute ungroup
#' @importFrom magrittr %>%
#' @importFrom rlang := abort arg_match as_string enquo expr exprs inform quo_is_null sym syms warn
#' @importFrom utils capture.output
#' @importFrom purrr map2 map_chr map_lgl
#' @importFrom stringr str_c str_detect
#' @importFrom assertthat assert_that is.number on_failure<-
#' @importFrom lubridate ceiling_date days duration floor_date time_length %--%
#'             ymd ymd_hms is.instant
"_PACKAGE"
