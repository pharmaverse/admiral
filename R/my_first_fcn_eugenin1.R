#' Dummy Function for Onboarding
#'
#' @description
#' Returns "Welcome to the admiral family!"
#'
#' @details
#' This is a function that returns "Welcome to the admiral family!"
#' independent of its inputs.
#'
#' @param anything Element in any format the user wants.
#'
#' @return A "Welcome to the admiral family!" string.
#'
#' @examples
#' welcome_fun(anything = "a")
#' welcome_fun(anything = 2)
#' welcome_fun(anything = NULL)
#' welcome_fun(data.frame(1:10, 2:11))

welcome_fun <- function(anything = NULL){

  if(is.null(anything)){

      cat("Welcome to the admiral family!")

      } else {cat("Welcome to the admiral family!")}

}
