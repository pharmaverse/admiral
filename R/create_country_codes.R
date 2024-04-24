#' Country Code Lookup
#'
#' @description
#' These pre-defined country codes are sourced from
#' [ISO Standards](https://www.iso.org/iso-3166-country-codes.html).
#'
#' @details
#'
#' `country_code` is the 3-letter county code commonly found in the ADSL COUNTRY variable.
#' `country_name` is the country long name corresponding to to the 3-letter code.
#' `country_number` is the numeric code corresponding to an alphabetic sorting of
#' the 3-letter codes.
#'
#' To see the entire table in the console, run `print(country_code_lookup)`.
#'
#' @seealso [dose_freq_lookup]
#'
#' @export
#'
#' @keywords metadata
#'
#' @family metadata
#'
#' @examples
#' library(tibble)
#' library(dplyr, warn.conflicts = FALSE)
#' library(lubridate)
#'
#' # Create reference dataset for periods
#' adsl <- tribble(#' adsl <- tribble(
#'   ~USUBJID, ~SEX, ~COUNTRY,
#'   "ST01-01", "F",  "AUT",
#'   "ST01-02", "M",  "MWI",
#'   "ST01-03", "F",  "GBR",
#'   "ST01-04", "M",  "CHE",
#'   "ST01-05", "M",  "NOR",
#'   "ST01-06", "F",  "JPN",
#'   "ST01-07", "F",  "USA"
#' )
#'
#' covar <- adsl %>%
#'   derive_vars_merged(
#'   dataset_add = country_code_lookup,
#'   new_vars = exprs(COUNTRYN = country_number, COUNTRYL = country_name),
#'   by_vars = exprs(COUNTRY = country_code)
#'   )
#' covar
#'
#' @rdname country_code_lookup
country_code_lookup <- tibble(
  country_code = c(
    "AFG", "ALA", "ALB", "DZA", "ASM", "AND", "AGO", "AIA", "ATA", "ATG",
    "ARG", "ARM", "ABW", "AUS", "AUT", "AZE", "BHS", "BHR", "BGD", "BRB",
    "BLR", "BEL", "BLZ", "BEN", "BMU", "BTN", "BOL", "BES", "BIH", "BWA",
    "BVT", "BRA", "IOT", "BRN", "BGR", "BFA", "BDI", "CPV", "KHM", "CMR",
    "CAN", "CYM", "CAF", "TCD", "CHL", "CHN", "CXR", "CCK", "COL", "COM",
    "COG", "COD", "COK", "CRI", "CIV", "HRV", "CUB", "CUW", "CYP", "CZE",
    "DNK", "DJI", "DMA", "DOM", "ECU", "EGY", "SLV", "GNQ", "ERI", "EST",
    "SWZ", "ETH", "FLK", "FRO", "FJI", "FIN", "FRA", "GUF", "PYF", "ATF",
    "GAB", "GMB", "GEO", "DEU", "GHA", "GIB", "GRC", "GRL", "GRD", "GLP",
    "GUM", "GTM", "GGY", "GIN", "GNB", "GUY", "HTI", "HMD", "VAT", "HND",
    "HKG", "HUN", "ISL", "IND", "IDN", "IRN", "IRQ", "IRL", "IMN", "ISR",
    "ITA", "JAM", "JPN", "JEY", "JOR", "KAZ", "KEN", "KIR", "PRK", "KOR",
    "KWT", "KGZ", "LAO", "LVA", "LBN", "LSO", "LBR", "LBY", "LIE", "LTU",
    "LUX", "MAC", "MDG", "MWI", "MYS", "MDV", "MLI", "MLT", "MHL", "MTQ",
    "MRT", "MUS", "MYT", "MEX", "FSM", "MDA", "MCO", "MNG", "MNE", "MSR",
    "MAR", "MOZ", "MMR", "NAM", "NRU", "NPL", "NLD", "NCL", "NZL", "NIC",
    "NER", "NGA", "NIU", "NFK", "MKD", "MNP", "NOR", "OMN", "PAK", "PLW",
    "PSE", "PAN", "PNG", "PRY", "PER", "PHL", "PCN", "POL", "PRT", "PRI",
    "QAT", "REU", "ROU", "RUS", "RWA", "BLM", "SHN", "KNA", "LCA", "MAF",
    "SPM", "VCT", "WSM", "SMR", "STP", "SAU", "SEN", "SRB", "SYC", "SLE",
    "SGP", "SXM", "SVK", "SVN", "SLB", "SOM", "ZAF", "SGS", "SSD", "ESP",
    "LKA", "SDN", "SUR", "SJM", "SWE", "CHE", "SYR", "TWN", "TJK", "TZA",
    "THA", "TLS", "TGO", "TKL", "TON", "TTO", "TUN", "TUR", "TKM", "TCA",
    "TUV", "UGA", "UKR", "ARE", "GBR", "USA", "UMI", "URY", "UZB", "VUT",
    "VEN", "VNM", "VGB", "VIR", "WLF", "ESH", "YEM", "ZMB", "ZWE"
  ),
  country_name = c(
    "Afghanistan", "Åland Islands", "Albania", "Algeria", "American Samoa",
    "Andorra", "Angola", "Anguilla", "Antarctica", "Antigua and Barbuda",
    "Argentina", "Armenia", "Aruba", "Australia", "Austria", "Azerbaijan",
    "Bahamas", "Bahrain", "Bangladesh", "Barbados", "Belarus", "Belgium",
    "Belize", "Benin", "Bermuda", "Bhutan", "Bolivia", "Bonaire, Sint Eustatius and Saba",
    "Bosnia and Herzegovina", "Botswana", "Bouvet Island", "Brazil",
    "British Indian Ocean Territory", "Brunei Darussalam", "Bulgaria", "Burkina Faso",
    "Burundi", "Cabo Verde", "Cambodia", "Cameroon", "Canada", "Cayman Islands",
    "Central African Republic", "Chad", "Chile", "China", "Christmas Island",
    "Cocos (Keeling) Islands", "Colombia", "Comoros", "Congo", "Congo, the Democratic Republic of the",
    "Cook Islands", "Costa Rica", "Côte d'Ivoire", "Croatia", "Cuba", "Curaçao",
    "Cyprus", "Czech Republic", "Denmark", "Djibouti", "Dominica", "Dominican Republic",
    "Ecuador", "Egypt", "El Salvador", "Equatorial Guinea", "Eritrea", "Estonia",
    "Eswatini", "Ethiopia", "Falkland Islands (Malvinas)", "Faroe Islands", "Fiji",
    "Finland", "France", "French Guiana", "French Polynesia", "French Southern Territories",
    "Gabon", "Gambia", "Georgia", "Germany", "Ghana", "Gibraltar", "Greece",
    "Greenland", "Grenada", "Guadeloupe", "Guam", "Guatemala", "Guernsey", "Guinea",
    "Guinea-Bissau", "Guyana", "Haiti", "Heard Island and McDonald Islands",
    "Holy See (Vatican City State)", "Honduras", "Hong Kong", "Hungary", "Iceland",
    "India", "Indonesia", "Iran, Islamic Republic of", "Iraq", "Ireland", "Isle of Man",
    "Israel", "Italy", "Jamaica", "Japan", "Jersey", "Jordan", "Kazakhstan", "Kenya",
    "Kiribati", "Korea, Democratic People's Republic of", "Korea, Republic of", "Kuwait",
    "Kyrgyzstan", "Lao People's Democratic Republic", "Latvia", "Lebanon", "Lesotho",
    "Liberia", "Libya", "Liechtenstein", "Lithuania", "Luxembourg", "Macao",
    "Madagascar", "Malawi", "Malaysia", "Maldives", "Mali", "Malta",
    "Marshall Islands", "Martinique", "Mauritania", "Mauritius", "Mayotte",
    "Mexico", "Micronesia, Federated States of", "Moldova, Republic of", "Monaco",
    "Mongolia", "Montenegro", "Montserrat", "Morocco", "Mozambique", "Myanmar",
    "Namibia", "Nauru", "Nepal", "Netherlands", "New Caledonia", "New Zealand",
    "Nicaragua", "Niger", "Nigeria", "Niue", "Norfolk Island", "North Macedonia",
    "Northern Mariana Islands", "Norway", "Oman", "Pakistan", "Palau", "Palestine, State of",
    "Panama", "Papua New Guinea", "Paraguay", "Peru", "Philippines", "Pitcairn",
    "Poland", "Portugal", "Puerto Rico", "Qatar", "Réunion", "Romania", "Russian Federation",
    "Rwanda", "Saint Barthélemy", "Saint Helena, Ascension and Tristan da Cunha",
    "Saint Kitts and Nevis", "Saint Lucia", "Saint Martin (French part)", "Saint Pierre and Miquelon",
    "Saint Vincent and the Grenadines", "Samoa", "San Marino", "Sao Tome and Principe",
    "Saudi Arabia", "Senegal", "Serbia", "Seychelles", "Sierra Leone", "Singapore",
    "Sint Maarten (Dutch part)", "Slovakia", "Slovenia", "Solomon Islands", "Somalia",
    "South Africa", "South Georgia and the South Sandwich Islands", "South Sudan", "Spain",
    "Sri Lanka", "Sudan", "Suriname", "Svalbard and Jan Mayen", "Sweden", "Switzerland",
    "Syrian Arab Republic", "Taiwan, Province of China", "Tajikistan",
    "Tanzania, United Republic of", "Thailand", "Timor-Leste", "Togo", "Tokelau",
    "Tonga", "Trinidad and Tobago", "Tunisia", "Turkey", "Turkmenistan", "Turks and Caicos Islands",
    "Tuvalu", "Uganda", "Ukraine", "United Arab Emirates", "United Kingdom",
    "United States", "United States Minor Outlying Islands", "Uruguay", "Uzbekistan",
    "Vanuatu", "Venezuela, Bolivarian Republic of", "Viet Nam", "Virgin Islands, British",
    "Virgin Islands, U.S.", "Wallis and Futuna", "Western Sahara", "Yemen", "Zambia",
    "Zimbabwe"
  )
) %>%
  arrange(country_code)

# Convert ISO 3166 alpha 3 country codes to numbers 1-249
country_code_lookup$country_number <- seq_len(nrow(country_code_lookup))

