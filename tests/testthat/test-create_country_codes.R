## Test 1: Country Codes ----
test_that("create_country_codes Test 1: Country Codes", {
  expected <- tibble::tibble(
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
      "Afghanistan", "Aland Islands", "Albania", "Algeria", "American Samoa",
      "Andorra", "Angola", "Anguilla", "Antarctica", "Antigua and Barbuda",
      "Argentina", "Armenia", "Aruba", "Australia", "Austria", "Azerbaijan",
      "Bahamas", "Bahrain", "Bangladesh", "Barbados", "Belarus", "Belgium",
      "Belize", "Benin", "Bermuda", "Bhutan", "Bolivia, Plurinational State of",
      "Bonaire, Sint Eustatius and Saba", "Bosnia and Herzegovina", "Botswana",
      "Bouvet Island", "Brazil", "British Indian Ocean Territory", "Brunei Darussalam",
      "Bulgaria", "Burkina Faso", "Burundi", "Cabo Verde", "Cambodia", "Cameroon",
      "Canada", "Cayman Islands", "Central African Republic", "Chad", "Chile",
      "China", "Christmas Island", "Cocos (Keeling) Islands", "Colombia",
      "Comoros", "Congo", "Congo, Democratic Republic of the", "Cook Islands",
      "Costa Rica", "Cote d'Ivoire", "Croatia", "Cuba", "Curacao", "Cyprus",
      "Czechia", "Denmark", "Djibouti", "Dominica", "Dominican Republic",
      "Ecuador", "Egypt", "El Salvador", "Equatorial Guinea", "Eritrea", "Estonia",
      "Eswatini", "Ethiopia", "Falkland Islands (Malvinas)", "Faroe Islands", "Fiji",
      "Finland", "France", "French Guiana", "French Polynesia", "French Southern Territories",
      "Gabon", "Gambia", "Georgia", "Germany", "Ghana", "Gibraltar", "Greece",
      "Greenland", "Grenada", "Guadeloupe", "Guam", "Guatemala", "Guernsey", "Guinea",
      "Guinea-Bissau", "Guyana", "Haiti", "Heard Island and McDonald Islands",
      "Holy See", "Honduras", "Hong Kong", "Hungary", "Iceland",
      "India", "Indonesia", "Iran, Islamic Republic of", "Iraq", "Ireland", "Isle of Man",
      "Israel", "Italy", "Jamaica", "Japan", "Jersey", "Jordan", "Kazakhstan", "Kenya",
      "Kiribati", "Korea, Democratic People's Republic of", "Korea, Republic of", "Kuwait",
      "Kyrgyzstan", "Lao People's Democratic Republic", "Latvia", "Lebanon", "Lesotho",
      "Liberia", "Libya", "Liechtenstein", "Lithuania", "Luxembourg", "Macao",
      "Madagascar", "Malawi", "Malaysia", "Maldives", "Mali", "Malta",
      "Marshall Islands", "Martinique", "Mauritania", "Mauritius", "Mayotte",
      "Mexico", "Micronesia, Federated States of", "Moldova, Republic of", "Monaco",
      "Mongolia", "Montenegro", "Montserrat", "Morocco", "Mozambique", "Myanmar",
      "Namibia", "Nauru", "Nepal", "Netherlands, Kingdom of the", "New Caledonia", "New Zealand",
      "Nicaragua", "Niger", "Nigeria", "Niue", "Norfolk Island", "North Macedonia",
      "Northern Mariana Islands", "Norway", "Oman", "Pakistan", "Palau", "Palestine, State of",
      "Panama", "Papua New Guinea", "Paraguay", "Peru", "Philippines", "Pitcairn",
      "Poland", "Portugal", "Puerto Rico", "Qatar", "Reunion", "Romania", "Russian Federation",
      "Rwanda", "Saint Barthelemy", "Saint Helena, Ascension and Tristan da Cunha",
      "Saint Kitts and Nevis", "Saint Lucia", "Saint Martin (French part)",
      "Saint Pierre and Miquelon", "Saint Vincent and the Grenadines", "Samoa",
      "San Marino", "Sao Tome and Principe", "Saudi Arabia", "Senegal", "Serbia",
      "Seychelles", "Sierra Leone", "Singapore", "Sint Maarten (Dutch part)",
      "Slovakia", "Slovenia", "Solomon Islands", "Somalia", "South Africa",
      "South Georgia and the South Sandwich Islands", "South Sudan", "Spain",
      "Sri Lanka", "Sudan", "Suriname", "Svalbard and Jan Mayen", "Sweden", "Switzerland",
      "Syrian Arab Republic", "Taiwan", "Tajikistan",
      "Tanzania, United Republic of", "Thailand", "Timor-Leste", "Togo", "Tokelau",
      "Tonga", "Trinidad and Tobago", "Tunisia", "Turkey", "Turkmenistan",
      "Turks and Caicos Islands", "Tuvalu", "Uganda", "Ukraine", "United Arab Emirates",
      "United Kingdom of Great Britain and Northern Ireland",
      "United States of America", "United States Minor Outlying Islands", "Uruguay", "Uzbekistan",
      "Vanuatu", "Venezuela, Bolivarian Republic of", "Viet Nam", "Virgin Islands (British)",
      "Virgin Islands (U.S.)", "Wallis and Futuna", "Western Sahara", "Yemen", "Zambia",
      "Zimbabwe"
    )
  ) %>%
    arrange(country_code)

  # Convert ISO 3166 alpha 3 country codes to numbers 1-249
  expected$country_number <- as.numeric(seq_len(nrow(expected)))


  actual <- country_code_lookup

  admiraldev::expect_dfs_equal(
    base = expected,
    comp = actual,
    keys = c("country_code")
  )
})
