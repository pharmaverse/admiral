#' Derive Analysis Age
#'
#' Derives analysis age (`AAGE`) and analysis age unit (`AAGEU`)
#'
#' @param dataset Input dataset
#'
#'   The columns specified by the `start_date` and the `end_date` parameter are
#'   expected.
#'
#' @param start_date The start date
#'
#'   A date or date-time object is expected.
#'
#'   Default: `BRTHDT`
#'
#' @param end_date The end date
#'
#'   A date or date-time object is expected.
#'
#'   Default: `RANDDT`
#'
#' @param unit Unit
#'
#'   The age is derived in the specified unit
#'
#'   Default: 'years'
#'
#'   Permitted Values: 'years', 'months', 'days', 'hours', 'minutes', 'seconds'
#'
#' @details The age is derived as the integer part of the duration from start to
#'   end date in the specified unit.
#'
#' @author Stefan Bundfuss
#'
#' @return The input dataset with ``AAGE`` and ``AAGEU`` added
#'
#' @export
#'
#' @seealso [derive_duration()]
#'
#' @examples
#' data <- tibble::tribble(
#'   ~BRTHDT, ~RANDDT,
#'   lubridate::ymd("1984-09-06"), lubridate::ymd("2020-02-24")
#' )
#'
#' derive_aage(data)
derive_aage <- function(dataset,
                        start_date = BRTHDT,
                        end_date = RANDDT,
                        unit = "years") {
  derive_duration(
    dataset,
    new_var = AAGE,
    new_var_unit = AAGEU,
    start_date = !!enquo(start_date),
    end_date = !!enquo(end_date),
    out_unit = unit,
    add_one = FALSE,
    trunc_out = TRUE
  )
}


#' Calculate age groups
#'
#' Calculates age groups by given breaks and returns either factor or numeric.
#'
#' @param x Numeric vector providing the age.
#' @param breaks Numeric vector providing the breaks for the age groups.
#' @param labels Character vector providing the labels of the age groups.
#' The length of `labels` needs to be greater by one than the length of `breaks`.
#' @param leq Should intervals be defined as "lower or equal" or "greater or equal"?
#' Defaults to the first.
#'
#' @return Factor (for `compute_agegr`) or numeric (for `compute_agegrn`) with defined age groups.
#'
#' @author Ondrej Slama
#'
#' @export
#'
#' @examples
#' compute_agegr(1:10, 5)
#' # provide arbitrary labels
#' compute_agegr(1:10, 5, labels = c("Less than 5 years old",
#'                                   "More or equal to 5 years old"))
#' # change the interval borders
#' compute_agegr(1:10, 5, leq = TRUE)
#' # add more breaks
#' compute_agegr(0:100, c(30, 60))
#' # compare factor with numeric output
#' compute_agegr(1:10, 5)
#' compute_agegrn(1:10, 5)
compute_agegr <- function(x, breaks, labels = NULL, leq = TRUE) {
  assert_that(
    is.numeric(x),
    length(x) > 0,
    is.numeric(breaks),
    length(breaks) > 0,
    is.null(labels) | is.character(labels),
    is.logical(leq)
  )

  # determine functions for comparing numbers
  lower <- `if`(leq, `<=`, `<`)
  lbl_lower <- `if`(leq, "<=", "<")

  greater <- `if`(leq, `>`, `>=`)
  lbl_greater <- `if`(leq, ">", ">=")

  # use ifelse for binary outcome
  if (length(breaks) == 1) {
    # define labels if not given
    if (is.null(labels)) {
      labels <- paste0(c(lbl_lower, lbl_greater), breaks)
    } else {
      assert_that(length(labels) == 2)
    }
    out <- factor(ifelse(lower(x, breaks), labels[1], labels[2]),
                  levels = labels)
  } else {
    # otherwise use cut function
    # define labels if not given
    if (is.null(labels)) {
      lbr <- `if`(leq, "(", "[")
      rbr <- `if`(leq, "]", ")")
      labels <- c(
        paste0(lbl_lower, breaks[1]),
        paste0(lbr, breaks[1:(length(breaks) - 1)], ", ", breaks[2:length(breaks)], rbr),
        paste0(lbl_greater, breaks[length(breaks)])
      )
    } else {
      assert_that(length(labels) == length(breaks) + 1)
    }
    out <- cut(x, breaks = c(-Inf, breaks, Inf), labels = labels, include.lowest = T, right = leq)
  }

  return(out)
}

#' @describeIn compute_agegr Calculate age groups and return the number of the group.
#' @export
compute_agegrn <- function(x, breaks, labels = NULL, leq = FALSE) {
  as.integer(do.call(compute_agegr, as.list(environment(NULL))))
}
