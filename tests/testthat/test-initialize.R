context("test-initialize")

test_that("deriving predecessor variables from spec works", {
  advs_spec <- tibble::tribble(
    ~Dataset, ~Variable,  ~Source,       ~`Derivation / Comment`,
    "ADVS",   "STUDYID",  "Predecessor", "ADSL.STUDYID",
    "ADVS",   "USUBJID",  "Predecessor", "ADSL.USUBJID",
    "ADVS",   "COUNTRY",  "Predecessor", "ADSL.COUNTRY",
    "ADVS",   "AAGE",     "Predecessor", "ADSL.AAGE",
    "ADVS",   "AAGEU",    "Predecessor", "ADSL.AAGEU",
    "ADVS",   "AGEGR1",   "Predecessor", "ADSL.AGEGR1",
    "ADVS",   "VSTESTCD", "Predecessor", "ADSL.VSTESTCD",
    "ADVS",   "VSSTRESN", "Predecessor", "ADSL.VSSTRESN",
    "ADVS",   "VSSTRESU", "Predecessor", "ADSL.VSSTRESU",
    "ADVS",   "PARAM",    "Assigned",    "See ADVS Parameters tab",
    "ADVS",   "AVAL",     "Derived",     "See ADVS Parameters tab",
    "ADVS",   "AVALU",    "Assigned",    "Set to Standard Units [VS.VSSTRESU]",
  )
  dap_m3 <- structure(list(ADVS = advs_spec), class = "DAP_M3")
  vs <- tibble::tribble(
    ~STUDYID, ~USUBJID, ~DOMAIN, ~VSTESTCD, ~VSSTRESN, ~VSSTRESU,
    "ABC-01", "P01",    "VS",    "HEIGHT",  182,       "cm",
    "ABC-01", "P01",    "VS",    "WEIGHT",  84,        "kg",
    "ABC-01", "P02",    "VS",    "HEIGHT",  174,       "cm",
    "ABC-01", "P02",    "VS",    "WEIGHT",  89,        "kg",
  )
  adsl <- tibble::tribble(
    ~STUDYID, ~USUBJID, ~COUNTRY, ~AAGE, ~AAGEU,  ~AGEGR1, ~VAR,
    "ABC-01", "P01",    "GER",    34,    "YEARS", "<65",   "Y",
    "ABC-01", "P02",    "JPN",    29,    "YEARS", "<65",   "N"
  )
  advs <- vs %>%
    left_join(adsl, by = c("STUDYID", "USUBJID")) %>%
    select(STUDYID, USUBJID, COUNTRY, AAGE, AAGEU, AGEGR1, VSTESTCD, VSSTRESN, VSSTRESU)

  expect_dfs_equal(
    base = advs,
    compare = initialize("ADVS", dap_m3, list(adsl, vs)),
    keys = c("STUDYID", "USUBJID", "VSTESTCD")
  )
})
