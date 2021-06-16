#' Combine Supplemental SDTM domain with Parent SDTM domain
#'
#' Transposes Supplemental domain and joins the
#' supplemental SDTM domain to the parent domain by the
#' IDVAR and IDVARVAL.
#'
#' @param dataset Parent SDTM dataset
#'
#' @param dataset_supp Supplemental SDTM Dataset to join to Parent
#'
#' @author Alice Ehmann
#'
#' @return The input dataset with variables from the supplemental dataset added
#'
#' @keywords derivation adam
#'
#' @export
#'
#' @examples
#' library(dplyr, warn.conflict = FALSE)
#' data("vs")
#' data("suppvs")
#'
#' # join supplemental domain with parent domain
#' derive_supplemental_vars(
#'   vs,
#'   suppvs
#' )
#'
derive_supplemental_vars <- function(dataset, dataset_supp) {


  #' # =========================================================================
  # Transpose SUPP domain onto Parent
  # Current limitations:
  # Assumes IDVAR is the same for all observations, if this is not true,
  # workaround would be to split the SUPP domain into multiple SUPPS prior to transposing
  # =========================================================================


  domain <- dataset
  suppdomain <- dataset_supp

  if (toupper(deparse(substitute(ds))) == "DM") {
    suppdomain <- mutate(suppdomain, IDVAR = "DMSEQ",
                         IDVARVAL = 1)

    domain <- mutate(domain, DMSEQ = 1)
  }


  # Get IDVAR information
  idvar <- as.symbol(as.character(distinct(suppdomain, IDVAR)))
  idvar_position <- which(names(domain) == as.character(idvar))
  idvar_label <- attr(domain[[idvar_position]], 'label')

  # Transpose SUPP domain
  supp_t1 <- pivot_wider(suppdomain, id_cols=c(USUBJID, IDVARVAL), names_from=QNAM, values_from=QVAL)

  # Make IDVARVAL appropriate type for join
  # Add IDVAR label
  if (typeof(domain[[idvar_position]]) == "character") {
    supp_t2 <- mutate(supp_t1, !!idvar := structure(IDVARVAL, label=idvar_label))
  }
  else {
    supp_t2 <- mutate(supp_t1, !!idvar := structure(as.numeric(IDVARVAL), label=idvar_label))
  }

  supp_t <- select(supp_t2, -IDVARVAL)

  # Add QNAM labels
  qnam <- distinct(suppdomain, QNAM, QLABEL)
  for (i in 1:nrow(qnam)) {
    labs <- as.character(qnam[i,2])
    names(labs) <- qnam[i,1]
    supp_t <- supp_t %>%
      add_labels(labs)
  }

  # Join SUPP to Parent
  ds_suppds <- domain %>%
    left_join(supp_t, by=c("USUBJID", as.character(idvar)))

  if (toupper(deparse(substitute(ds))) == "DM") {
    ds_suppds <- select(ds_suppds, -DMSEQ)
  }

  ds_suppds

}



#' Join supplementary qualifier variables into the parent domain of SDTM
#'
#' The SDTM does not allow any new variables beside ones assigned to each SDTM
#' domain. So, Supplemental Qualifier is introduced to supplement each SDTM
#' domain to contain non standard variables. The supplemental data set is either
#' a single SUPPQUAL dataset or separate supplementary data sets (SUPP) such
#' as SUPPDM, SUPPAE, and SUPPEX.
#'
#' `derive_suppqual_vars` expects USUBJID, RDOMAIN, IDVAR, IDVARVAL, QNAM,
#'  QLABEL, and QVAL variables to exist in SUPPQUAL data set.
#'
#' @param parent_dataset A SDTM domain data set.
#' @param suppqual_dataset A Supplemental Qualifier (SUPPQUAL) data set.
#' @param .domain Two letter domain value. Used when supplemental data set is
#'   common across multiple SDTM domain.
#'
#' @return A data frame with SUPPQUAL variables appended to parent data set.
#' @export
#'
#' @examples
#' ## Sample AE data set ----
#' ae <- tibble::tribble(
#'   ~STUDYID  , ~DOMAIN, ~USUBJID  , ~AESEQ, ~AETERM,
#'   '1234-005', 'AE'   , 'XYZ-1001',      1, 'Nausea'
#' )
#'
#' ## Sample Supplemental Qualifiers for AE (SUPPAE) ----
#' suppae <- tibble::tribble(
#'   ~STUDYID  ,  ~RDOMAIN, ~USUBJID  , ~IDVAR , ~IDVARVAL, ~QNAM    , ~QLABEL            , ~QVAL,
#'   '1234-005', 'AE'     , 'XYZ-1001', 'AESEQ', '1'      , 'AECLINT', 'Clinical Interest', 'Y',
#' )
#'
#' \dontrun{
#'   derive_suppqual_vars(ae, suppae)
#' }
derive_suppqual_vars <- function(parent_dataset, suppqual_dataset, .domain = NULL) {

  ## Convert variable names to upper case
  parent <- as_tibble(set_names(parent_dataset, toupper), .name_repair = "minimal")
  supp <- as_tibble(set_names(suppqual_dataset, toupper), .name_repair = "minimal")

  assert_has_variables(
    supp,
    c("USUBJID", "RDOMAIN", "IDVAR", "IDVARVAL", "QNAM", "QLABEL", "QVAL")
  )

  if (!is.null(.domain)) {
    supp <- filter(supp, .data$RDOMAIN == .domain)
  }

  check_supp_domain(parent, supp, .domain = .domain)

  ## Get unique value across QNAM & IDVAR
  supp_unique <- distinct(supp, .data$QNAM, .data$IDVAR, .data$QLABEL)
  check_supp_idvar(supp_unique)

  supp_unique_list <- transpose(supp_unique)

  ## Loop IDVAR/QNAM
  pivoted <- lapply(supp_unique_list, function(.x) {
    if (!.x$IDVAR %in% c(NA, "")) {
      supp <- supp %>%
        filter(.data$IDVAR == .x$IDVAR) %>%
        tidyr::pivot_wider(
          id_cols = c("USUBJID", "IDVARVAL"),
          names_from = "QNAM",
          values_from = "QVAL"
        ) %>%
        rename(!! .x$IDVAR := "IDVARVAL") %>%
        # Convert IDVAR to match parent domain type
        modify_at(.x$IDVAR, function(x) {
          if (is.numeric(parent[[.x$IDVAR]])) as.numeric(x)
        })
    } else {
      supp <- supp %>%
        filter(.data$QNAM == .x$QNAM) %>%
        tidyr::pivot_wider(
          id_cols = "USUBJID",
          names_from = "QNAM",
          values_from = "QVAL"
        )
    }

    # Set label
    attr(supp[[.x$QNAM]], "label") <- .x$QLABEL
    supp
  })

  for (i in seq_along(pivoted)) {
    parent_nms <- names(parent)
    supp_nms <- names(pivoted[[i]])

    by_vars <- intersect(parent_nms, supp_nms)
    new_var <- setdiff(supp_nms, parent_nms)

    ## The following step helps parent dataset undistrubed:
    ##
    ## - Join with parent domain and retain only supp variables. This helps to
    ##   match observation and finally bind with original parent dataset
    join_var <- parent %>%
      left_join(pivoted[[i]], by = by_vars) %>%
      select(!! sym(new_var))

    ## Sometimes QNAM might have more than one IDVAR. For example, in AE domain
    ## QNAM = 'XXXX' and IDVAR = c('AEID', 'AESEQ'). Join all values into single
    ## values using [dplyr::coalesce]
    dup <- match(new_var, names(parent_dataset), nomatch = 0L)

    if (dup > 0L) {
      values <- map(
        list(parent[[new_var]], unlist(join_var, use.names = FALSE)),
        ~ na_if(., "")
      ) %>%
        reduce(coalesce) %>%
        modify_if(is.character, ~ tidyr::replace_na(., ""))

      join_var <- tibble(!! new_var := values)
    }

    parent_dataset <- bind_cols(parent_dataset, join_var)
  }

  parent_dataset
}

check_supp_idvar <- function(x) {
  x <- unclass(x)
  dup <- duplicated(x$QNAM)
  if (any(dup)) {
    message(
      paste0(
        str_glue("More than one IDVAR = '{x$IDVAR[dup]}' for a QNAM = \\
                 '{x$QNAM[dup]}'."),
        collapse = "\n") )
  }
}

check_supp_domain <- function(parent, supp, .domain = NULL) {
  parent <- unique(parent$DOMAIN)
  supp <- unique(supp$RDOMAIN)

  if (!is.null(.domain)) {
    if (is.na(match(.domain, supp))) {
      abort(str_glue("Can't find the domain `{.domain}` in `supp_dataset`."))
    }
  }

  if (!supp %in% parent) {
    abort("DOMAIN of `parent_dataset` and RDOMAIN of `supp_dataset` is not related.")
  }
}
