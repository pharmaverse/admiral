#' Environment Objects
#'
#' @details
#' Once in a while, we may encounter "locked binding for 'xxx'." errors
#' during the development process while building out functions. This may arise because
#' we want to create dynamic data/objects based on user-inputs that need modification
#' at points in time after the package has been loaded. To manage such data or objects,
#' R has a data structure known as an 'environment'. These environment objects are created
#' at build time, but can be populated with values after the package has been loaded and
#' update those values over the course of an R session. For more details how environments work,
#' see relevant sections on environments in R Packages and Advanced R textbooks for more details.
#' @noRd
admiraldev_environment <- new.env(parent = emptyenv())
# See respective ...R page for usage

# assertions.R ----
## assert_one_to_one
admiraldev_environment$many_to_one <- NULL
admiraldev_environment$one_to_many <- NULL

# datasets.R ----
## get_dataset
# Function above is used to retrieve many_to_one and one_to_many
