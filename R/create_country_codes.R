#' Country Code Lookup
#'
#' @description
#' These pre-defined country codes are sourced from
#' [ISO Standards](https://www.iso.org/iso-3166-country-codes.html).
#' See also [Wikipedia](https://en.wikipedia.org/wiki/ISO_3166-1_alpha-3).
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
#' adsl <- tribble( #' adsl <- tribble(
#'   ~USUBJID, ~SEX, ~COUNTRY,
#'   "ST01-01", "F", "AUT",
#'   "ST01-02", "M", "MWI",
#'   "ST01-03", "F", "GBR",
#'   "ST01-04", "M", "CHE",
#'   "ST01-05", "M", "NOR",
#'   "ST01-06", "F", "JPN",
#'   "ST01-07", "F", "USA"
#' )
#'
#' covar <- adsl %>%
#'   derive_vars_merged(
#'     dataset_add = country_code_lookup,
#'     new_vars = exprs(COUNTRYN = country_number, COUNTRYL = country_name),
#'     by_vars = exprs(COUNTRY = country_code)
#'   )
#' covar
#'
#' @rdname country_code_lookup
country_code_lookup <- tibble(
  country_code = c(
    "ABW", "AFG", "AGO", "AIA", "ALA", "ALB", "AND", "ARE", "ARG", "ARM",
    "ASM", "ATA", "ATF", "ATG", "AUS", "AUT", "AZE", "BDI", "BEL", "BEN",
    "BES", "BFA", "BGD", "BGR", "BHR", "BHS", "BIH", "BLM", "BLR", "BLZ",
    "BMU", "BOL", "BRA", "BRB", "BRN", "BTN", "BVT", "BWA", "CAF", "CAN",
    "CCK", "CHE", "CHL", "CHN", "CIV", "CMR", "COD", "COG", "COK", "COL",
    "COM", "CPV", "CRI", "CUB", "CUW", "CXR", "CYM", "CYP", "CZE", "DEU",
    "DJI", "DMA", "DNK", "DOM", "DZA", "ECU", "EGY", "ERI", "ESH", "ESP",
    "EST", "ETH", "FIN", "FJI", "FLK", "FRA", "FRO", "FSM", "GAB", "GBR",
    "GEO", "GGY", "GHA", "GIB", "GIN", "GLP", "GMB", "GNB", "GNQ", "GRC",
    "GRD", "GRL", "GTM", "GUF", "GUM", "GUY", "HKG", "HMD", "HND", "HRV",
    "HTI", "HUN", "IDN", "IMN", "IND", "IOT", "IRL", "IRN", "IRQ", "ISL",
    "ISR", "ITA", "JAM", "JEY", "JOR", "JPN", "KAZ", "KEN", "KGZ", "KHM",
    "KIR", "KNA", "KOR", "KWT", "LAO", "LBN", "LBR", "LBY", "LCA", "LIE",
    "LKA", "LSO", "LTU", "LUX", "LVA", "MAC", "MAF", "MAR", "MCO", "MDA",
    "MDG", "MDV", "MEX", "MHL", "MKD", "MLI", "MLT", "MMR", "MNE", "MNG",
    "MNP", "MOZ", "MRT", "MSR", "MTQ", "MUS", "MWI", "MYS", "MYT", "NAM",
    "NCL", "NER", "NFK", "NGA", "NIC", "NIU", "NLD", "NOR", "NPL", "NRU",
    "NZL", "OMN", "PAK", "PAN", "PCN", "PER", "PHL", "PLW", "PNG", "POL",
    "PRI", "PRK", "PRT", "PRY", "PSE", "PYF", "QAT", "REU", "ROU", "RUS",
    "RWA", "SAU", "SDN", "SEN", "SGP", "SGS", "SHN", "SJM", "SLB", "SLE",
    "SLV", "SMR", "SOM", "SPM", "SRB", "SSD", "STP", "SUR", "SVK", "SVN",
    "SWE", "SWZ", "SXM", "SYC", "SYR", "TCA", "TCD", "TGO", "THA", "TJK",
    "TKL", "TKM", "TLS", "TON", "TTO", "TUN", "TUR", "TUV", "TWN", "TZA",
    "UGA", "UKR", "UMI", "URY", "USA", "UZB", "VAT", "VCT", "VEN", "VGB",
    "VIR", "VNM", "VUT", "WLF", "WSM", "YEM", "ZAF", "ZMB", "ZWE"
  ),
  country_name = c(
    "Aruba", "Afghanistan", "Angola", "Anguilla", "Åland Islands",
    "Albania", "Andorra", "United Arab Emirates", "Argentina", "Armenia",
    "American Samoa", "Antarctica", "French Southern Territories",
    "Antigua and Barbuda", "Australia", "Austria", "Azerbaijan", "Burundi",
    "Belgium", "Benin", "Bonaire, Sint Eustatius and Saba", "Burkina Faso",
    "Bangladesh", "Bulgaria", "Bahrain", "Bahamas", "Bosnia and Herzegovina",
    "Saint Barthélemy", "Belarus", "Belize", "Bermuda", "Bolivia (Plurinational State of)",
    "Brazil", "Barbados", "Brunei Darussalam", "Bhutan", "Bouvet Island",
    "Botswana", "Central African Republic", "Canada", "Cocos (Keeling) Islands",
    "Switzerland", "Chile", "China", "Côte d'Ivoire", "Cameroon", "Congo (Democratic Republic of the)",
    "Congo", "Cook Islands", "Colombia", "Comoros", "Cabo Verde",
    "Costa Rica", "Cuba", "Curaçao", "Christmas Island", "Cayman Islands",
    "Cyprus", "Czechia", "Germany", "Djibouti", "Dominica", "Denmark",
    "Dominican Republic", "Algeria", "Ecuador", "Egypt", "Eritrea",
    "Western Sahara", "Spain", "Estonia", "Ethiopia", "Finland", "Fiji",
    "Falkland Islands (Malvinas)", "France", "Faroe Islands", "Micronesia (Federated States of)",
    "Gabon", "United Kingdom of Great Britain and Northern Ireland",
    "Georgia", "Guernsey", "Ghana", "Gibraltar", "Guinea", "Guadeloupe",
    "Gambia", "Guinea-Bissau", "Equatorial Guinea", "Greece", "Grenada",
    "Greenland", "Guatemala", "French Guiana", "Guam", "Guyana",
    "Hong Kong", "Heard Island and McDonald Islands", "Honduras",
    "Croatia", "Haiti", "Hungary", "Indonesia", "Isle of Man",
    "India", "British Indian Ocean Territory", "Ireland", "Iran (Islamic Republic of)",
    "Iraq", "Iceland", "Israel", "Italy", "Jamaica", "Jersey",
    "Jordan", "Japan", "Kazakhstan", "Kenya", "Kyrgyzstan", "Cambodia",
    "Kiribati", "Saint Kitts and Nevis", "Korea (Republic of)",
    "Kuwait", "Lao People's Democratic Republic", "Lebanon", "Liberia",
    "Libya", "Saint Lucia", "Liechtenstein", "Sri Lanka", "Lesotho",
    "Lithuania", "Luxembourg", "Latvia", "Macao", "Saint Martin (French part)",
    "Morocco", "Monaco", "Moldova (Republic of)", "Madagascar",
    "Maldives", "Mexico", "Marshall Islands", "North Macedonia",
    "Mali", "Malta", "Myanmar", "Montenegro", "Mongolia",
    "Northern Mariana Islands", "Mozambique", "Mauritania", "Montserrat",
    "Martinique", "Mauritius", "Malawi", "Malaysia", "Mayotte",
    "Namibia", "New Caledonia", "Niger", "Norfolk Island", "Nigeria",
    "Nicaragua", "Niue", "Netherlands", "Norway", "Nepal", "Nauru",
    "New Zealand", "Oman", "Pakistan", "Panama", "Pitcairn", "Peru",
    "Philippines", "Palau", "Papua New Guinea", "Poland", "Puerto Rico",
    "Korea (Democratic People's Republic of)", "Portugal", "Paraguay",
    "Palestine, State of", "French Polynesia", "Qatar", "Réunion",
    "Romania", "Russian Federation", "Rwanda", "Saudi Arabia",
    "Sudan", "Senegal", "Singapore", "South Georgia and the South Sandwich Islands",
    "Saint Helena, Ascension and Tristan da Cunha", "Svalbard and Jan Mayen",
    "Solomon Islands", "Sierra Leone", "El Salvador", "San Marino",
    "Somalia", "Saint Pierre and Miquelon", "Serbia", "South Sudan",
    "Sao Tome and Principe", "Suriname", "Slovakia", "Slovenia",
    "Sweden", "Eswatini", "Sint Maarten (Dutch part)", "Seychelles",
    "Syrian Arab Republic", "Turks and Caicos Islands", "Chad",
    "Togo", "Thailand", "Tajikistan", "Tokelau", "Turkmenistan",
    "Timor-Leste", "Tonga", "Trinidad and Tobago", "Tunisia",
    "Turkey", "Tuvalu", "Taiwan", "Tanzania, United Republic of",
    "Uganda", "Ukraine", "United States Minor Outlying Islands",
    "Uruguay", "United States of America", "Uzbekistan", "Holy See",
    "Saint Vincent and the Grenadines", "Venezuela (Bolivarian Republic of)",
    "Virgin Islands (British)", "Virgin Islands (U.S.)", "Viet Nam",
    "Vanuatu", "Wallis and Futuna", "Samoa", "Yemen", "South Africa",
    "Zambia", "Zimbabwe"
  )
) %>%
  arrange(country_code)

# Convert ISO 3166 alpha 3 country codes to numbers 1-249
country_code_lookup$country_number <- seq_len(nrow(country_code_lookup))
