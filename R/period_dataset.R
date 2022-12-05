#' Create a Reference Dataset for Subperiods, Periods, or Phases
#'
#' The function creates a reference dataset for subperiods, periods, or phases
#' from the `ADSL` dataset. The reference dataset can be used to derive
#' subperiod, period, or phase variables like `ASPER`, `ASPRSDT`, `ASPREDT`,
#' `APERIOD`, `APERSDT`, `APEREDT`, `TRTA`, `APHASEN`, `PHSDTM`, `PHEDTM`, ...
#' in OCCDS and BDS datasets.
#'
#' @param dataset ADSL dataset
#'
#'   The variables specified by `new_vars` and `subject_keys` are expected. For
#'   each element of `new_vars` at least one variable of the form of the right
#'   hand side value must be available in the dataset.
#'
#' @param new_vars New variables
#'
#'   A named list of variables like `vars(PHSDT = PHwSDT, PHEDT = PHwEDT, APHASE
#'   = APHASEw)` is expected. The left hand side of the elements defines a
#'   variable of the output dataset, the right hand side defines the source
#'   variables from the ADSL dataset in CDISC notation.
#'
#'   If the lower case letter "w"  is used it refers to a phase variable, if the
#'   lower case letters "xx" are used it refers to a period variable, and if
#'   both "xx" and "w" are used it refers to a subperiod variable.
#'
#'   Only one type must be used, e.g., all right hand side values must refer to
#'   period variables. It is not allowed to mix for example period and subperiod
#'   variables. If period *and* subperiod variables are required, separate
#'   reference datasets must be created.
#'
#' @param subject_keys Variables to uniquely identify a subject
#'
#'   A list of quosures where the expressions are symbols as returned by
#'   `vars()` is expected.
#'
#' @author Stefan Bundfuss
#'
#' @details For each subject and each subperiod/period/phase where at least one
#'   of the source variable is not `NA` an observation is added to the output
#'   dataset.
#'
#'   Depending on the type of the source variable (subperiod, period, or phase)
#'   the variable `ASPER`, `APERIOD`, or `APHASEN` is added and set to the
#'   number of the subperiod, period, or phase.
#'
#'   The variables specified for `new_vars` (left hand side) are added to the
#'   output dataset and set to the value of the source variable (right hand
#'   side).
#'
#' @return A period reference dataset (see "Details" section)
#'
#' @seealso [derive_vars_period()]
#'
#' @keywords create_aux
#' @family create_aux
#'
#' @export
#'
#' @examples
#' library(tibble)
#' library(dplyr, warn.conflicts = FALSE)
#' library(lubridate)
#'
#' # Create reference dataset for periods
#' adsl <- tribble(
#'   ~USUBJID, ~AP01SDT,     ~AP01EDT,     ~AP02SDT,     ~AP02EDT,     ~TRT01A, ~TRT02A,
#'   "1",      "2021-01-04", "2021-02-06", "2021-02-07", "2021-03-07", "A",     "B",
#'   "2",      "2021-02-02", "2021-03-02", "2021-03-03", "2021-04-01", "B",     "A",
#' ) %>%
#'   mutate(
#'     across(matches("AP\\d\\d[ES]DT"), ymd)
#'   ) %>%
#'   mutate(
#'     STUDYID = "xyz"
#'   )
#'
#' create_period_dataset(
#'   adsl,
#'   new_vars = vars(APERSDT = APxxSDT, APEREDT = APxxEDT, TRTA = TRTxxA)
#' )
#'
#' # Create reference dataset for phases
#' adsl <- tribble(
#'   ~USUBJID, ~PH1SDT,      ~PH1EDT,      ~PH2SDT,      ~PH2EDT,      ~APHASE1,    ~APHASE2,
#'   "1",      "2021-01-04", "2021-02-06", "2021-02-07", "2021-03-07", "TREATMENT", "FUP",
#'   "2",      "2021-02-02", "2021-03-02", NA,           NA,           "TREATMENT", NA
#' ) %>%
#'   mutate(
#'     across(matches("PH\\d[ES]DT"), ymd)
#'   ) %>%
#'   mutate(
#'     STUDYID = "xyz"
#'   )
#'
#' create_period_dataset(
#'   adsl,
#'   new_vars = vars(PHSDT = PHwSDT, PHEDT = PHwEDT, APHASE = APHASEw)
#' )
#'
#' # Create reference datasets for subperiods
#' adsl <- tribble(
#'   ~USUBJID, ~P01S1SDT,    ~P01S1EDT,    ~P01S2SDT,    ~P01S2EDT,    ~P02S1SDT,    ~P02S1EDT,
#'   "1",      "2021-01-04", "2021-01-19", "2021-01-20", "2021-02-06", "2021-02-07", "2021-03-07",
#'   "2",      "2021-02-02", "2021-03-02", NA,           NA,           "2021-03-03", "2021-04-01"
#' ) %>%
#'   mutate(
#'     across(matches("P\\d\\dS\\d[ES]DT"), ymd)
#'   ) %>%
#'   mutate(
#'     STUDYID = "xyz"
#'   )
#'
#' create_period_dataset(
#'   adsl,
#'   new_vars = vars(ASPRSDT = PxxSwSDT, ASPREDT = PxxSwEDT)
#' )
create_period_dataset <- function(dataset,
                                  new_vars,
                                  subject_keys = get_admiral_option("subject_keys")) {
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
        paste0(
          "The right hand side values of `new_vars` have to be CDISC style ",
          "subperiod, period, or phase variables."
        ),
        "I.e., they must contain the xx or w fragment, e.g., APxxSDT, PxxSwSDT, or PHwSDT.",
        sep = "\n"
      )
    )
  }
  if (length(mode) > 1) {
    abort(
      paste0(
        "More than one type of subperiod, period, or phase variables ",
        "is specified for `new_vars`:\n",
        if_else(
          "subperiod" %in% mode,
          paste0("subperiod: ", enumerate(new_vars_chr[mode == "subperiod"]), "\n"),
          ""
        ),
        if_else(
          "period" %in% mode,
          paste0("period: ", enumerate(new_vars_chr[mode == "period"]), "\n"),
          ""
        ),
        if_else(
          "phase" %in% mode,
          paste0("phase: ", enumerate(new_vars_chr[mode == "phase"]), "\n"),
          ""
        )
      )
    )
  }
  prefix <- syms(str_match(cols, "(\\w+)\\\\")[, 2])
  num_var_chr <- c(
    subperiod = "ASPER",
    period = "APERIOD",
    phase = "APHASEN"
  )
  num_var <- syms(num_var_chr)
  period_ref <- vector("list", length(new_vars))
  for (i in seq_along(new_vars)) {
    if (!any(str_detect(colnames(dataset), names_pattern[[i]]))) {
      abort(paste(
        "No variables of the form",
        new_vars_chr[[i]],
        "were found in the input dataset."
      ))
    }
    if (mode == "subperiod") {
      period_ref[[i]] <- pivot_longer(
        select(dataset, !!!subject_keys, matches(cols[[i]])),
        matches(cols[[i]]),
        names_to = c(".value", "APERIOD", num_var_chr[[mode]]),
        names_pattern = names_pattern[[i]]
      ) %>%
        rename(!!sym(new_vars_names[[i]]) := !!prefix[[i]]) %>%
        mutate(
          APERIOD = as.integer(APERIOD),
          !!num_var[[mode]] := as.integer(!!num_var[[mode]])
        ) %>%
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
      period_ref_final <- period_ref[[1]]
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

#' Add Subperiod, Period, or Phase Variables to ADSL
#'
#' The function adds subperiod, period, or phase variables like `P01S1SDT`,
#' `P01S2SDT`, `AP01SDTM`, `AP02SDTM`, `TRT01A`, `TRT02A`, `PH1SDT`, `PH2SDT`,
#' ... to the input dataset. The values of the variables are defined by a period
#' reference dataset which has one observations per patient and subperiod,
#' period, or phase.
#'
#' @param dataset ADSL dataset
#'
#'   The variables specified by `subject_keys` are expected.
#'
#' @param dataset_ref Period reference dataset
#'
#'   The variables specified by `new_vars` and `subject_keys` are expected.
#'
#'   If subperiod variables are requested, `APERIOD` and `ASPER` are expected.
#'   If period variables are requested. `APERIOD` is expected. If phase
#'   variables are requested, `APHASEN` is expected.
#'
#' @param new_vars New variables
#'
#'   A named list of variables like `vars(PHwSDT = PHSDT, PHwEDT = PHEDT,
#'   APHASEw = APHASE)` is expected. The left hand side of the elements defines
#'   a set of variables (in CDISC notation) to be added to the output dataset.
#'   The right hand side defines the source variable from the period reference
#'   dataset.
#'
#'   If the lower case letter "w"  is used it refers to a phase variable, if the
#'   lower case letters "xx" are used it refers to a period variable, and if
#'   both "xx" and "w" are used it refers to a subperiod variable.
#'
#'   Only one type must be used, e.g., all left hand side values must refer to
#'   period variables. It is not allowed to mix for example period and subperiod
#'   variables. If period *and* subperiod variables are required, separate calls
#'   must be used.
#'
#' @param subject_keys Variables to uniquely identify a subject
#'
#'   A list of quosures where the expressions are symbols as returned by
#'   `vars()` is expected.
#'
#' @author Stefan Bundfuss
#'
#' @details For each subperiod/period/phase in the period reference dataset and
#'   each element in `new_vars` a variable (LHS value of `new_vars`) is added to
#'   the output dataset and set to the value of the source variable (RHS value
#'   of `new_vars`.
#'
#' @return The input dataset with subperiod/period/phase variables added (see
#'   "Details" section)
#'
#' @seealso [create_period_dataset()]
#'
#' @keywords der_adsl
#' @family der_adsl
#'
#' @export
#'
#' @examples
#' library(tibble)
#' library(dplyr, warn.conflicts = FALSE)
#' library(lubridate)
#'
#' adsl <- tibble(STUDYID = "xyz", USUBJID = c("1", "2"))
#'
#' # Add period variables to ADSL
#' period_ref <- tribble(
#'   ~USUBJID, ~APERIOD, ~APERSDT,     ~APEREDT,
#'   "1",             1, "2021-01-04", "2021-02-06",
#'   "1",             2, "2021-02-07", "2021-03-07",
#'   "2",             1, "2021-02-02", "2021-03-02",
#'   "2",             2, "2021-03-03", "2021-04-01"
#' ) %>%
#'   mutate(
#'     STUDYID = "xyz",
#'     APERIOD = as.integer(APERIOD),
#'     across(matches("APER[ES]DT"), ymd)
#'   )
#'
#' derive_vars_period(
#'   adsl,
#'   dataset_ref = period_ref,
#'   new_vars = vars(APxxSDT = APERSDT, APxxEDT = APEREDT)
#' ) %>%
#'   select(STUDYID, USUBJID, AP01SDT, AP01EDT, AP02SDT, AP02EDT)
#'
#' # Add phase variables to ADSL
#' phase_ref <- tribble(
#'   ~USUBJID, ~APHASEN, ~PHSDT,       ~PHEDT,       ~APHASE,
#'   "1",             1, "2021-01-04", "2021-02-06", "TREATMENT",
#'   "1",             2, "2021-02-07", "2021-03-07", "FUP",
#'   "2",             1, "2021-02-02", "2021-03-02", "TREATMENT"
#' ) %>%
#'   mutate(
#'     STUDYID = "xyz",
#'     APHASEN = as.integer(APHASEN),
#'     across(matches("PH[ES]DT"), ymd)
#'   )
#'
#' derive_vars_period(
#'   adsl,
#'   dataset_ref = phase_ref,
#'   new_vars = vars(PHwSDT = PHSDT, PHwEDT = PHEDT, APHASEw = APHASE)
#' ) %>%
#'   select(STUDYID, USUBJID, PH1SDT, PH1EDT, PH2SDT, PH2EDT, APHASE1, APHASE2)
#'
#' # Add subperiod variables to ADSL
#' subperiod_ref <- tribble(
#'   ~USUBJID, ~APERIOD, ~ASPER, ~ASPRSDT,     ~ASPREDT,
#'   "1",             1,      1, "2021-01-04", "2021-01-19",
#'   "1",             1,      2, "2021-01-20", "2021-02-06",
#'   "1",             2,      1, "2021-02-07", "2021-03-07",
#'   "2",             1,      1, "2021-02-02", "2021-03-02",
#'   "2",             2,      1, "2021-03-03", "2021-04-01"
#' ) %>%
#'   mutate(
#'     STUDYID = "xyz",
#'     APERIOD = as.integer(APERIOD),
#'     ASPER = as.integer(ASPER),
#'     across(matches("ASPR[ES]DT"), ymd)
#'   )
#'
#' derive_vars_period(
#'   adsl,
#'   dataset_ref = subperiod_ref,
#'   new_vars = vars(PxxSwSDT = ASPRSDT, PxxSwEDT = ASPREDT)
#' ) %>%
#'   select(STUDYID, USUBJID, P01S1SDT, P01S1EDT, P01S2SDT, P01S2EDT, P02S1SDT, P02S1EDT)
derive_vars_period <- function(dataset,
                               dataset_ref,
                               new_vars,
                               subject_keys = get_admiral_option("subject_keys")) {
  assert_vars(new_vars, expect_names = TRUE)
  assert_vars(subject_keys)
  assert_data_frame(dataset, required_vars = subject_keys)
  assert_data_frame(dataset_ref, required_vars = subject_keys)

  new_vars_names <- names(new_vars)
  new_vars_chr <- vars2chr(new_vars)
  mode <- case_when(
    str_detect(new_vars_names, "\\w+xx\\w+w\\w*") ~ "subperiod",
    str_detect(new_vars_names, "\\w+xx\\w*") ~ "period",
    str_detect(new_vars_names, "\\w+w\\w*") ~ "phase",
    TRUE ~ "none"
  ) %>% unique()
  if (any(mode == "none")) {
    abort(
      paste(
        paste0(
          "The left hand side values of `new_vars` have to be CDISC style ",
          "subperiod, period, or phase variables."
        ),
        "I.e., they must contain the xx or w fragment, e.g., APxxSDT, PxxSwSDT, or PHwSDT.",
        sep = "\n"
      )
    )
  }
  if (length(mode) > 1) {
    abort(
      paste0(
        "More than one type of subperiod, period, or phase variables ",
        "is specified for `new_vars`:\n",
        if_else(
          "subperiod" %in% mode,
          paste0("subperiod: ", enumerate(new_vars_names[mode == "subperiod"]), "\n"),
          ""
        ),
        if_else(
          "period" %in% mode,
          paste0("period: ", enumerate(new_vars_names[mode == "period"]), "\n"),
          ""
        ),
        if_else(
          "phase" %in% mode,
          paste0("phase: ", enumerate(new_vars_names[mode == "phase"]), "\n"),
          ""
        )
      )
    )
  }
  if (mode == "subperiod") {
    id_vars <- vars(APERIOD, ASPER)
  } else if (mode == "period") {
    id_vars <- vars(APERIOD)
  } else {
    id_vars <- vars(APHASEN)
  }
  assert_data_frame(dataset_ref, required_vars = quo_c(subject_keys, new_vars, id_vars))

  ref_wide <- pivot_wider(
    dataset_ref,
    names_from = vars2chr(id_vars),
    values_from = unname(new_vars_chr)
  )

  # pivot_wider creates columns like APERSDT_1, APERSDT_2, ...
  # these need to be renamed to variables like AP01SDT, AP02SDT, ...
  rename_arg <- colnames(select(ref_wide, !!!negate_vars(subject_keys)))
  split_names <- str_match(rename_arg, "(\\w+?)_(\\d{1,2})_?(\\d)?")
  source_vars <- names(new_vars_chr)
  names(source_vars) <- new_vars_chr
  index <- split_names[, 3]
  if (mode == "phase") {
    names_rename_arg <- str_replace(source_vars[split_names[, 2]], "w", index)
  } else {
    # add leading zero for "xx" fragment
    index <- if_else(str_length(index) == 1, paste0("0", index), index)
    names_rename_arg <- str_replace(source_vars[split_names[, 2]], "xx", index)
    if (mode == "subperiod") {
      index2 <- split_names[, 4]
      names_rename_arg <- str_replace(names_rename_arg, "w", index2)
    }
  }
  names(rename_arg) <- names_rename_arg

  derive_vars_merged(
    dataset,
    dataset_add = ref_wide,
    by_vars = subject_keys
  ) %>% rename(all_of(rename_arg))
}
