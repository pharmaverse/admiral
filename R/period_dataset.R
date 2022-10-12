#' @author Stefan Bundfuss
#'
#' @export
#'
#' @examples
#' library(tibble)
create_period_dataset <- function(dataset,
                                  new_vars,
                                  subject_keys = vars(STUDYID, USUBJID)) {
  assert_vars(new_vars, expect_names = TRUE)
  assert_vars(subject_keys)
  assert_data_frame(dataset, required_vars = subject_keys)

  new_vars_names <- names(new_vars)
  new_vars_chr <- vars2chr(new_vars)
  cols <- str_replace(new_vars_chr, "xx", "\\\\d\\\\d") %>%
    str_replace("w", "\\\\d")
  names_pattern <- str_replace(new_vars_chr, "xx", "\\\\d\\\\d") %>%
    str_replace("w", "\\\\d") %>%
    str_replace("(\\w+)\\\\d", "(\\1)\\\\d") %>%
    str_replace_all("((\\\\d)+)", "(\\1)")
  mode <- case_when(
    str_detect(new_vars_chr, "\\w+xx\\w+w\\w*") ~ "subperiod",
    str_detect(new_vars_chr, "\\w+xx\\w*") ~ "period",
    str_detect(new_vars_chr, "\\w+w\\w*") ~ "phase",
    TRUE ~ "none"
  ) %>% unique()
  if (any(mode == "none")) {
    abort(
      paste(
        "The right hand side values of `new_vars` have to be CDISC style subperiod, period, or phase variables.",
        "I.e., they must contain the xx or w fragment, e.g., APxxSDT, PxxSwSDT, or PHwSDT.",
        sep = "\n")
    )
  }
  if (length(mode) > 1) {
    abort(
      paste0(
        "More than one type of subperiod, period, or phase variables is specified for `new_vars`:\n",
        if_else("subperiod" %in% mode, paste0("subperiod: ", enumerate(new_vars_chr[mode == "subperiod"]), "\n"), ""),
        if_else("period" %in% mode, paste0("period: ", enumerate(new_vars_chr[mode == "period"]), "\n"), ""),
        if_else("phase" %in% mode, paste0("phase: ", enumerate(new_vars_chr[mode == "phase"]), "\n"), "")
      )
    )
  }
  prefix <- syms(str_match(cols, "(\\w+)\\\\")[,2])
  num_var_chr <- c(
    subperiod = "ASPER",
    period = "APERIOD",
    phase = "APHASEN")
  num_var <- syms(num_var_chr)
  period_ref <- vector("list", length(new_vars))
  for (i in seq_along(new_vars)) {
    if (!any(str_detect(colnames(dataset), names_pattern[[i]]))) {
      abort(paste(
        "No variables of the form",
        new_vars_chr[[i]],
        "were found in the input dataset."
        )
      )
    }
    if (mode == "subperiod") {
      period_ref[[i]] <- pivot_longer(
        select(dataset,!!!subject_keys, matches(cols[[i]])),
        matches(cols[[i]]),
        names_to = c(".value", "APERIOD", num_var_chr[[mode]]),
        names_pattern = names_pattern[[i]]
      ) %>%
        rename(!!sym(new_vars_names[[i]]) := !!prefix[[i]]) %>%
        mutate(
          APERIOD = as.integer(APERIOD),
          !!num_var[[mode]] := as.integer(!!num_var[[mode]])) %>%
        filter(!is.na(!!sym(new_vars_names[[i]])))
      by_vars <- vars(APERIOD, !!sym(num_var[[mode]]))
    } else {
      period_ref[[i]] <- pivot_longer(
        select(dataset, !!!subject_keys, matches(cols[[i]])),
        matches(cols[[i]]),
        names_to = c(".value", num_var_chr[[mode]]),
        names_pattern = names_pattern[[i]]
      ) %>%
        rename(!!sym(new_vars_names[[i]]) := !!prefix[[i]]) %>%
        mutate(!!num_var[[mode]] := as.integer(!!num_var[[mode]])) %>%
        filter(!is.na(!!sym(new_vars_names[[i]])))
      by_vars <- vars(!!sym(num_var[[mode]]))
    }
    if (i == 1) {
      period_ref_final <-  period_ref[[1]]
    } else {
      period_ref_final <- derive_vars_merged(
        period_ref_final,
        dataset_add = period_ref[[i]],
        by_vars = quo_c(subject_keys, by_vars)
      )
    }
  }
  period_ref_final
}
