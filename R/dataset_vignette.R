#' Output a dataset in a vignette in the admiral format
#'
#' Output a dataset in a vignette with the pre-specified format
#' By default, the variables specified in `display_vars`
#'
#' @param dataset Dataset to output in the vignette
#'
#' @param display_vars Variables selected to demonstrate the outcome of the derivation.
#'
#' Permitted Values: list of variables
#'
#' Default is NULL
#'
#' If `display_vars` is not NULL, only the selected variables are visible in the vignette while the
#' other variables are hidden. They can be made visible by clicking the`Choose the columns to
#'  display` button.
#'
#' @param filter Filter condition
#'
#'   The specified condition is applied to the dataset before it is displayed.
#'
#'   *Permitted Values:* a condition
#'
#' @keywords vignette
#'
#' @examples
#' library(DT)
#' library(admiral)
#'
#' data("dm")
#' dataset_vignette(dm)
#' dataset_vignette(dm, display_vars=vars(USUBJID, RFSTDTC, DTHDTC), filter= ARMCD == "Pbo")


dataset_vignette<-function(dataset, display_vars=NULL, filter=NULL) {

  display_vars<-assert_vars(display_vars, optional = TRUE)
  assert_data_frame(dataset,required_vars = display_vars)
  filter <- assert_filter_cond(enquo(filter), optional = TRUE)

  out<-dataset %>%
    filter_if(filter)%>%
    mutate_if(is.character,as.factor)

  if (!is.null(display_vars)){
    hide_columns<-which(!(colnames(out) %in% vars2chr(display_vars)))
    cols_to_hide<-list(list(targets = hide_columns-1, visible = FALSE))
  }
  else{
    cols_to_hide<-list()
  }

  datatable(
    out,
    rownames = FALSE,
    filter = 'top',
    extensions = c('Buttons', 'ColReorder', 'Scroller'),
    options = list(
      columnDefs =cols_to_hide,
      searchHighlight = TRUE,
      searching = TRUE,
      pageLength = 5,
      lengthMenu = c(5, 10, 15, 20, 50, 100),
      dom = 'Bfrtip',
      buttons = list(list(extend = "colvis",
                          text = "Choose the columns to display",
                          scroller = T,
                          collectionLayout= 'fixed two-column' )),
      colReorder = TRUE)
  )
}




