#' Join Supplementary Qualifier Variables into the Parent SDTM Domain
#'
#' The SDTM does not allow any new variables beside ones assigned to each SDTM
#' domain. So, Supplemental Qualifier is introduced to supplement each SDTM
#' domain to contain non standard variables. `dataset_suppqual` can be either
#' a single SUPPQUAL dataset or separate supplementary data sets (SUPP) such
#' as SUPPDM, SUPPAE, and SUPPEX. When a `dataset_suppqual` is a single SUPPQUAL
#' dataset, specify two character`domain` value.
#'
#' `derive_vars_suppqual()` expects `USUBJID`, `RDOMAIN`, `IDVAR`, `IDVARVAL`,
#'  `QNAM`, `QLABEL`, and `QVAL` variables to exist in `dataset_suppqual`.
#'
#' @param dataset A SDTM domain data set.
#'
#' @param dataset_suppqual A Supplemental Qualifier (SUPPQUAL) data set.
#'
#'
#' @param domain Two letter domain value. Used when supplemental data set is
#'   common across multiple SDTM domain.
#'
#' @return A data frame with SUPPQUAL variables appended to parent data set.
#'
#' @author Vignesh Thanikachalam
#'
#' @export
#'
#' @examples
#' ## The following example includes selected variables from AE and SUPPAE
#' ## datasets for a rash whose locations are the face, neck, and chest.
#' ae <- tibble::tribble(
#'   ~STUDYID,   ~DOMAIN, ~USUBJID,   ~AESEQ, ~AETERM,  ~AELOC,
#'   "1234-005", "AE",    "XYZ-1001",      1, "RASH",  "MULTIPLE",
#'   "1234-005", "AE",    "XYZ-1002",      1, "NAUSEA", "",
#' )
#' suppae <- tibble::tribble(
#'   ~STUDYID,   ~RDOMAIN, ~USUBJID,    ~IDVAR,  ~IDVARVAL, ~QNAM,     ~QLABEL,     ~QVAL,
#'   "1234-005", "AE",     "XYZ-1001", "AESEQ", "1",        "AELOC1", "Location 1", "FACE",
#'   "1234-005", "AE",     "XYZ-1001", "AESEQ", "1",        "AELOC2", "Location 2", "NECK",
#'   "1234-005", "AE",     "XYZ-1001", "AESEQ", "1",        "AELOC3", "Location 3", "CHEST",
#' )
#' derive_vars_suppqual(ae, suppae)
#'
#' ## The following example included subjects with multiple/other specific race.
#' dm <- tibble::tribble(
#'   ~STUDYID, ~DOMAIN, ~USUBJID, ~RACE,
#'   "ABC",    "DM",    "001",    "OTHER",
#'   "ABC",    "DM",    "002",    "MULTIPLE",
#'   "ABC",    "DM",    "003",    NA,
#'   "ABC",    "DM",    "004",    "ASIAN"
#' )
#' suppdm <- tibble::tribble(
#'   ~STUDYID, ~RDOMAIN, ~USUBJID, ~IDVAR, ~IDVARVAL, ~QNAM,     ~QLABEL,       ~QVAL,
#'   "ABC",   "DM",      "001",     "",     "",       "RACEOTH", "Race, Other", "BRAZILIAN",
#'   "ABC",   "DM",      "002",     "",     "",       "RACE1"  , "Race 1",      "AMERICAN",
#'   "ABC",   "DM",      "002",     "",     "",       "RACE2"  , "Race 2",      "OTHER",
#'   "ABC",   "DM",      "002",     "",     "",       "RACEOTH", "Race, Other", "ABORIGINE"
#' )
#' derive_vars_suppqual(dm, suppdm)
derive_vars_suppqual <- function(dataset, dataset_suppqual, domain = NULL) {
  assert_data_frame(dataset)
  assert_data_frame(
    dataset_suppqual,
    required_vars = vars(USUBJID, RDOMAIN, IDVAR, IDVARVAL, QNAM, QLABEL, QVAL)
  )
  assert_character_scalar(domain, optional = TRUE)

  if (!is.null(domain)) {
    dataset_suppqual <- filter(dataset_suppqual, .data$RDOMAIN == domain)
  }

  assert_is_supp_domain(dataset, dataset_suppqual, .domain = domain)

  ## Get unique value across QNAM & IDVAR
  supp_unique <- dataset_suppqual %>%
    distinct(.data$QNAM, .data$IDVAR, .data$QLABEL)

  assert_supp_idvar(supp_unique)

  supp_unique_list <- transpose(supp_unique)

  ## Loop IDVAR/QNAM
  pivoted <- lapply(supp_unique_list, function(.x) {
    if (!.x$IDVAR %in% c(NA, "")) {
      supp <- dataset_suppqual %>%
        filter(.data$IDVAR == .x$IDVAR, .data$QNAM == .x$QNAM) %>%
        select(USUBJID, IDVARVAL, QNAM, QVAL) %>%
        spread(key = QNAM, value = QVAL) %>%
        rename(!! .x$IDVAR := "IDVARVAL") %>%
        # Convert IDVAR to match parent domain type
        modify_at(.x$IDVAR, function(x) {
          if (is.numeric(dataset[[.x$IDVAR]])) as.numeric(x) else x
        })
    } else {
      supp <- dataset_suppqual %>%
        filter(.data$QNAM == .x$QNAM) %>%
        select(USUBJID, QNAM, QVAL) %>%
        spread(key = QNAM, value = QVAL)
    }

    # Set label
    attr(supp[[.x$QNAM]], "label") <- .x$QLABEL
    supp
  })

  for (i in seq_along(pivoted)) {
    parent_nms <- names(dataset)
    supp_nms <- names(pivoted[[i]])

    new_var <- supp_unique_list[[i]]$QNAM
    by_vars <- intersect(parent_nms, setdiff(supp_nms, new_var))

    ## The following step helps parent dataset undistrubed:
    ##
    ## - Join with parent domain and retain only supp variables. This helps to
    ##   match observation and finally bind with original parent dataset
    join_var <- dataset %>%
      left_join(pivoted[[i]], by = by_vars) %>%
      select(dplyr::matches(new_var))

    ## Sometimes QNAM might have more than one IDVAR. For example, in AE domain
    ## QNAM = 'XXXX' and IDVAR = c('AEID', 'AESEQ'). Join all values into single
    ## values using [dplyr::coalesce]
    if (ncol(join_var) > 1L) {
      values <- join_var %>%
        map(~ na_if(., "")) %>%
        reduce(coalesce)

      join_var <- tibble(!! new_var := values)
    }

    rep_var <- intersect(parent_nms, names(join_var))

    if (length(rep_var) > 0L) {
      dataset <- select(dataset, -!! sym(rep_var))
    }

    dataset <- bind_cols(dataset, join_var)
  }

  dataset
}

#' Join Supplementary Qualifier Variables into the Parent SDTM Domain
#'
#' `derive_suppqual_vars()` was renamed to `derive_vars_suppqual()` to create a
#' more consistent API.
#'
#' @keywords internal
#'
#' @export
derive_suppqual_vars <- function(dataset, dataset_suppqual, domain = NULL) {
  deprecate_warn("0.3.0", "derive_suppqual_vars()", "derive_vars_suppqual()")
  derive_vars_suppqual(dataset, dataset_suppqual, domain)
}
