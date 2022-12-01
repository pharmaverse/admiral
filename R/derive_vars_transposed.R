#' Derive Variables by Transposing and Merging a Second Dataset
#'
#' Adds variables from a vertical dataset after transposing it into a wide one.
#'
#' @param dataset Input dataset
#'
#'   The variables specified by the `by_vars` parameter are required
#'
#' @param dataset_merge Dataset to transpose and merge
#'
#'   The variables specified by the `by_vars`, `key_var` and `value_var` parameters
#'   are expected
#'
#' @param by_vars Keys used to merge `dataset_merge` with `dataset`
#'
#' @param key_var The variable of `dataset_merge` containing the names of the
#'   transposed variables
#'
#' @param value_var The variable of `dataset_merge` containing the values of the
#'   transposed variables
#'
#' @param filter Expression used to restrict the records of `dataset_merge` prior to transposing
#'
#' @details
#' After filtering `dataset_merge` based upon the condition provided in `filter`, this
#' dataset is transposed and subsequently merged onto `dataset` using `by_vars` as
#' keys.
#'
#' @author Thomas Neitmann
#'
#' @return The input dataset with transposed variables from `dataset_merge` added
#'
#' @family der_gen
#' @keywords der_gen
#'
#' @export
#'
#' @examples
#' library(tibble)
#' library(dplyr, warn.conflicts = FALSE)
#'
#' cm <- tribble(
#'   ~USUBJID, ~CMGRPID, ~CMREFID, ~CMDECOD,
#'   "BP40257-1001", "14", "1192056", "PARACETAMOL",
#'   "BP40257-1001", "18", "2007001", "SOLUMEDROL",
#'   "BP40257-1002", "19", "2791596", "SPIRONOLACTONE"
#' )
#' facm <- tribble(
#'   ~USUBJID, ~FAGRPID, ~FAREFID, ~FATESTCD, ~FASTRESC,
#'   "BP40257-1001", "1", "1192056", "CMATC1CD", "N",
#'   "BP40257-1001", "1", "1192056", "CMATC2CD", "N02",
#'   "BP40257-1001", "1", "1192056", "CMATC3CD", "N02B",
#'   "BP40257-1001", "1", "1192056", "CMATC4CD", "N02BE",
#'   "BP40257-1001", "1", "2007001", "CMATC1CD", "D",
#'   "BP40257-1001", "1", "2007001", "CMATC2CD", "D10",
#'   "BP40257-1001", "1", "2007001", "CMATC3CD", "D10A",
#'   "BP40257-1001", "1", "2007001", "CMATC4CD", "D10AA",
#'   "BP40257-1001", "2", "2007001", "CMATC1CD", "D",
#'   "BP40257-1001", "2", "2007001", "CMATC2CD", "D07",
#'   "BP40257-1001", "2", "2007001", "CMATC3CD", "D07A",
#'   "BP40257-1001", "2", "2007001", "CMATC4CD", "D07AA",
#'   "BP40257-1001", "3", "2007001", "CMATC1CD", "H",
#'   "BP40257-1001", "3", "2007001", "CMATC2CD", "H02",
#'   "BP40257-1001", "3", "2007001", "CMATC3CD", "H02A",
#'   "BP40257-1001", "3", "2007001", "CMATC4CD", "H02AB",
#'   "BP40257-1002", "1", "2791596", "CMATC1CD", "C",
#'   "BP40257-1002", "1", "2791596", "CMATC2CD", "C03",
#'   "BP40257-1002", "1", "2791596", "CMATC3CD", "C03D",
#'   "BP40257-1002", "1", "2791596", "CMATC4CD", "C03DA"
#' )
#'
#' cm %>%
#'   derive_vars_transposed(
#'     facm,
#'     by_vars = vars(USUBJID, CMREFID = FAREFID),
#'     key_var = FATESTCD,
#'     value_var = FASTRESC
#'   ) %>%
#'   select(USUBJID, CMDECOD, starts_with("CMATC"))
derive_vars_transposed <- function(dataset,
                                   dataset_merge,
                                   by_vars,
                                   key_var,
                                   value_var,
                                   filter = NULL) {
  key_var <- assert_symbol(enquo(key_var))
  value_var <- assert_symbol(enquo(value_var))
  filter <- assert_filter_cond(enquo(filter), optional = TRUE)
  assert_vars(by_vars)
  assert_data_frame(dataset, required_vars = replace_values_by_names(by_vars))
  assert_data_frame(dataset_merge, required_vars = quo_c(by_vars, key_var, value_var))

  dataset_transposed <- dataset_merge %>%
    filter_if(filter) %>%
    pivot_wider(names_from = !!key_var, values_from = !!value_var)

  left_join(dataset, dataset_transposed, by = vars2chr(by_vars))
}

#' Derive ATC Class Variables
#'
#' Add Anatomical Therapeutic Chemical class variables from `FACM` to `ADCM`
#'
#' @param dataset Input dataset
#'
#'   The variables specified by the `by_vars` parameter are required
#'
#' @param dataset_facm FACM dataset
#'
#'   The variables specified by the `by_vars` and `value_var` parameters,
#'   `FAGRPID` and `FATESTCD` are required
#'
#' @param by_vars Keys used to merge `dataset_facm` with `dataset`
#'
#'   *Permitted Values:* list of variables
#'
#' @param value_var The variable of `dataset_facm` containing the values of the
#'   transposed variables
#'
#'   Default: `FASTRESC`
#'
#' @author Thomas Neitmann
#'
#' @return The input dataset with ATC variables added
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
#'   ~USUBJID, ~CMGRPID, ~CMREFID, ~CMDECOD,
#'   "BP40257-1001", "14", "1192056", "PARACETAMOL",
#'   "BP40257-1001", "18", "2007001", "SOLUMEDROL",
#'   "BP40257-1002", "19", "2791596", "SPIRONOLACTONE"
#' )
#' facm <- tribble(
#'   ~USUBJID, ~FAGRPID, ~FAREFID, ~FATESTCD, ~FASTRESC,
#'   "BP40257-1001", "1", "1192056", "CMATC1CD", "N",
#'   "BP40257-1001", "1", "1192056", "CMATC2CD", "N02",
#'   "BP40257-1001", "1", "1192056", "CMATC3CD", "N02B",
#'   "BP40257-1001", "1", "1192056", "CMATC4CD", "N02BE",
#'   "BP40257-1001", "1", "2007001", "CMATC1CD", "D",
#'   "BP40257-1001", "1", "2007001", "CMATC2CD", "D10",
#'   "BP40257-1001", "1", "2007001", "CMATC3CD", "D10A",
#'   "BP40257-1001", "1", "2007001", "CMATC4CD", "D10AA",
#'   "BP40257-1001", "2", "2007001", "CMATC1CD", "D",
#'   "BP40257-1001", "2", "2007001", "CMATC2CD", "D07",
#'   "BP40257-1001", "2", "2007001", "CMATC3CD", "D07A",
#'   "BP40257-1001", "2", "2007001", "CMATC4CD", "D07AA",
#'   "BP40257-1001", "3", "2007001", "CMATC1CD", "H",
#'   "BP40257-1001", "3", "2007001", "CMATC2CD", "H02",
#'   "BP40257-1001", "3", "2007001", "CMATC3CD", "H02A",
#'   "BP40257-1001", "3", "2007001", "CMATC4CD", "H02AB",
#'   "BP40257-1002", "1", "2791596", "CMATC1CD", "C",
#'   "BP40257-1002", "1", "2791596", "CMATC2CD", "C03",
#'   "BP40257-1002", "1", "2791596", "CMATC3CD", "C03D",
#'   "BP40257-1002", "1", "2791596", "CMATC4CD", "C03DA"
#' )
#'
#' derive_vars_atc(cm, facm)
derive_vars_atc <- function(dataset,
                            dataset_facm,
                            by_vars = vars(USUBJID, CMREFID = FAREFID),
                            value_var = FASTRESC) {
  value_var <- assert_symbol(enquo(value_var))
  assert_vars(by_vars)
  assert_data_frame(dataset, required_vars = replace_values_by_names(by_vars))
  assert_data_frame(dataset_facm, required_vars = vars(!!!by_vars, !!value_var, FAGRPID, FATESTCD))

  dataset %>%
    derive_vars_transposed(
      select(dataset_facm, !!!unname(by_vars), !!value_var, FAGRPID, FATESTCD),
      by_vars = by_vars,
      key_var = FATESTCD,
      value_var = !!value_var,
      filter = str_detect(FATESTCD, "^CMATC[1-4](CD)?$")
    ) %>%
    select(-starts_with("FA")) %>%
    rename_at(vars(starts_with("CMATC")), ~ str_remove(.x, "^CM"))
}
