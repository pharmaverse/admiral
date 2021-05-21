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
