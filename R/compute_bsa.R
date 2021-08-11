#' Derive BSA (Body Surface Area)
#'
#' Derives BSA from HEIGHT and WEIGHT making use of the specified derivation formula
#'
#' @param height_code HEIGHT parameter code
#'
#'   The observations where `PARAMCD` equals the specified value are considered
#'   as the HEIGHT assessments. It is expected that HEIGHT is measured in cm.
#'
#'   Permitted Values: character value
#'
#' @param weight_code WEIGHT parameter code
#'
#'   The observations where `PARAMCD` equals the specified value are considered
#'   as the WEIGHT assessments. It is expected that WEIGHT is measured in kg.
#'
#'   Permitted Values: character value
#'
#' @param method Derivation method to use
#'
#'   The derivation method, e.g. Mosteller will use sqrt(height(cm) * weight(kg)) / 3600
#'
#'   Permitted Values: character value
#'
#' @details The BSA (Body Surface Area) is derived with the specified formula:
#'          "Mosteller", "DuBois-DuBois", "Haycock", "Gehan-George", "Boyd",
#'          ""Fujimoto" or "Takahira".
#'
#' @author Eric Simms
#'
#' @return The BSA (Body Surface Area) in m^2.
#'
#' @keywords computation adam BSA
#'
#' @export
#'
#' @examples
#' # derive BSA by the Mosteller method
#' compute_bsa(
#'   height = 170,
#'   weight = 75,
#'   method = "Mosteller")
#'
#' # derive BSA by the DuBois & DuBois method
#' compute_bsa(
#'   height = 170,
#'   weight = 75,
#'   method = "DuBois-DuBois")
#'
# # derive BSA by the Boyd method
# compute_bsa(
#   height = c(170, 185),
#   weight = c(75, 90),
#   method = "Boyd")

compute_bsa <- function(height = height,
                        weight = weight,
                        method = "Mosteller"
                        ) {
  # Checks
  #assert_is_numeric(height, weight)
  assert_character_scalar(method, values = c("Mosteller", "DuBois-DuBois", "Haycock",
                                             "Gehan-George", "Boyd", "Fujimoto",
                                             "Takahira"))

  # Derivation
  if (method == "Mosteller") {
      bsa <- sqrt(height * weight / 3600)
  } else if (method == "DuBois-DuBois") {
      # Note: the DuBois & DuBois formula expects the value of height in meters; we need to convert from cm.
      bsa <- 0.20247 * (height/100) ^ 0.725 * weight ^ 0.425
  } else if (method == "Haycock") {
      bsa <- 0.024265 * height ^ 0.3964 * weight ^ 0.5378
  } else if (method == "Gehan-George") {
      bsa <- 0.0235 * height ^ 0.42246 * weight ^ 0.51456
  } else if (method == "Boyd") {
      # Note: the Boyd formula expects the value of weight in grams; we need to convert from kg.
      bsa <- 0.0003207 * (height ^ 0.3) *
                  (1000 * weight) ^ (0.7285 - (0.0188 * log10(1000 * weight)))
  } else if (method == "Fujimoto") {
      bsa <- 0.008883 * height ^ 0.663 * weight ^ 0.444
  } else if (method == "Takahira") {
      bsa <- 0.007241 * height ^ 0.725 * weight ^ 0.425
  }

  bsa <- round(bsa, 2)
  bsa
}

