#' Output a Dataset in a Vignette in the admiral Format
#'
#' Output a dataset in a vignette with the pre-specified admiral format.
#'
#' @param dataset Dataset to output in the vignette
#'
#' @param display_vars Variables selected to demonstrate the outcome of the derivation
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
#' The specified condition is applied to the dataset before it is displayed.
#'
#' Permitted Values: a condition
#'
#' @return A HTML table
#'
#' @keywords dev_utility
#'
#' @export
#'
dataset_vignette <- function(dataset, display_vars = NULL, filter = NULL) {
  display_vars <- assert_vars(display_vars, optional = TRUE)
  assert_data_frame(dataset, required_vars = display_vars)
  filter <- assert_filter_cond(enexpr(filter), optional = TRUE)

  out <- dataset %>%
    filter_if(filter) %>%
    mutate(across(where(is.character), as.factor))

  # Create a short markdown table when this function is called outside {pkgdown}
  if (!identical(Sys.getenv("IN_PKGDOWN"), "true")) {
    if (is.null(display_vars)) {
      return(knitr::kable(utils::head(out, 10)))
    } else {
      return(knitr::kable(utils::head(select(out, !!!display_vars), 10)))
    }
  }

  if (!is.null(display_vars)) {
    hide_columns <- which(!(colnames(out) %in% vars2chr(display_vars)))
    cols_to_hide <- list(list(targets = hide_columns - 1, visible = FALSE))
  } else {
    cols_to_hide <- list()
  }
  htmltools::tagList(
    htmltools::htmlDependency(
      name = "dt-scroll",
      version = "1.0.0",
      src = "www",
      stylesheet = "style.css",
      package = "admiraldev"
    ),
    DT::datatable(
      out,
      rownames = FALSE,
      filter = "top",
      height = "auto",
      width = "auto",
      extensions = c("Buttons", "ColReorder", "Scroller"),
      options = list(
        columnDefs = cols_to_hide,
        searchHighlight = TRUE,
        searching = TRUE,
        pageLength = 5,
        lengthMenu = c(5, 10, 15, 20, 50, 100),
        dom = "<Bfr<\"dt-scroll\"t>ipl>",
        buttons = list(list(
          extend = "colvis",
          text = "Choose the columns to display",
          scroller = T,
          collectionLayout = "fixed two-column"
        )),
        colReorder = TRUE
      )
    )
  )
}
