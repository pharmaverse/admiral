#' Derive ATC Class Variables
#'
#' @description Add Anatomical Therapeutic Chemical class variables from `FACM` to `ADCM`.
#'
#' **Note:** This is a wrapper function for the more generic `derive_vars_transposed()`.
#'
#' @param dataset
#'   `r roxygen_param_dataset(expected_vars = c("by_vars"))`
#'
#' @param dataset_facm FACM dataset
#'
#'   The variables specified by the `by_vars`, `id_vars`, and `value_var`
#'   arguments and `FATESTCD` are required. The variables `by_vars`, `id_vars`,
#'   and `FATESTCD` must be a unique key.
#'
#' @param by_vars Grouping variables
#'
#'  Keys used to merge `dataset_facm` with `dataset`.
#'
#' @param id_vars ID variables
#'
#'  Variables (excluding by_vars) that uniquely identify each observation in `dataset_merge`.
#'
#' `r roxygen_param_by_vars()`
#'
#' @param value_var The variable of `dataset_facm` containing the values of the
#'   transposed variables
#'
#' @return The input dataset with ATC variables added
#'
#' @seealso [derive_vars_transposed()]
#'
#' @family der_occds
#' @keywords der_occds
#'
#' @export
#'
#' @examples
#' library(tibble)
#'
#' cm <- tribble(
#'   ~STUDYID,  ~USUBJID,       ~CMGRPID, ~CMREFID,  ~CMDECOD,
#'   "STUDY01", "BP40257-1001", "14",     "1192056", "PARACETAMOL",
#'   "STUDY01", "BP40257-1001", "18",     "2007001", "SOLUMEDROL",
#'   "STUDY01", "BP40257-1002", "19",     "2791596", "SPIRONOLACTONE"
#' )
#' facm <- tribble(
#'   ~STUDYID,  ~USUBJID,       ~FAGRPID, ~FAREFID,  ~FATESTCD,  ~FASTRESC,
#'   "STUDY01", "BP40257-1001", "1",      "1192056", "CMATC1CD", "N",
#'   "STUDY01", "BP40257-1001", "1",      "1192056", "CMATC2CD", "N02",
#'   "STUDY01", "BP40257-1001", "1",      "1192056", "CMATC3CD", "N02B",
#'   "STUDY01", "BP40257-1001", "1",      "1192056", "CMATC4CD", "N02BE",
#'   "STUDY01", "BP40257-1001", "1",      "2007001", "CMATC1CD", "D",
#'   "STUDY01", "BP40257-1001", "1",      "2007001", "CMATC2CD", "D10",
#'   "STUDY01", "BP40257-1001", "1",      "2007001", "CMATC3CD", "D10A",
#'   "STUDY01", "BP40257-1001", "1",      "2007001", "CMATC4CD", "D10AA",
#'   "STUDY01", "BP40257-1001", "2",      "2007001", "CMATC1CD", "D",
#'   "STUDY01", "BP40257-1001", "2",      "2007001", "CMATC2CD", "D07",
#'   "STUDY01", "BP40257-1001", "2",      "2007001", "CMATC3CD", "D07A",
#'   "STUDY01", "BP40257-1001", "2",      "2007001", "CMATC4CD", "D07AA",
#'   "STUDY01", "BP40257-1001", "3",      "2007001", "CMATC1CD", "H",
#'   "STUDY01", "BP40257-1001", "3",      "2007001", "CMATC2CD", "H02",
#'   "STUDY01", "BP40257-1001", "3",      "2007001", "CMATC3CD", "H02A",
#'   "STUDY01", "BP40257-1001", "3",      "2007001", "CMATC4CD", "H02AB",
#'   "STUDY01", "BP40257-1002", "1",      "2791596", "CMATC1CD", "C",
#'   "STUDY01", "BP40257-1002", "1",      "2791596", "CMATC2CD", "C03",
#'   "STUDY01", "BP40257-1002", "1",      "2791596", "CMATC3CD", "C03D",
#'   "STUDY01", "BP40257-1002", "1",      "2791596", "CMATC4CD", "C03DA"
#' )
#'
#' derive_vars_atc(cm, facm, id_vars = exprs(FAGRPID))
derive_vars_atc <- function(dataset,
                            dataset_facm,
                            by_vars = exprs(
                              !!!get_admiral_option("subject_keys"),
                              CMREFID = FAREFID
                            ),
                            id_vars = NULL,
                            value_var = FASTRESC) {
  value_var <- assert_symbol(enexpr(value_var))
  assert_vars(by_vars)
  assert_vars(id_vars, optional = TRUE)
  assert_data_frame(dataset, required_vars = replace_values_by_names(by_vars))
  assert_data_frame(
    dataset_facm,
    required_vars = exprs(!!!by_vars, !!value_var, !!!id_vars, FATESTCD)
  )

  tryCatch(
    data_transposed <- derive_vars_transposed(
      dataset,
      select(dataset_facm, !!!unname(by_vars), !!value_var, !!!id_vars, FATESTCD),
      by_vars = by_vars,
      id_vars = id_vars,
      key_var = FATESTCD,
      value_var = !!value_var,
      filter = str_detect(FATESTCD, "^CMATC[1-4](CD)?$")
    ),
    merge_duplicates = function(cnd) {
      cnd$message <- str_replace(cnd$message, "dataset_merge", "dataset_facm")
      cnd$body[[1]] <- "Please check data and `by_vars` and `id_vars` arguments."
      cnd_signal(cnd)
    }
  )
  data_transposed %>%
    select(-starts_with("FA")) %>%
    rename_with(.fn = ~ str_remove(.x, "^CM"), .cols = starts_with("CMATC"))
}
