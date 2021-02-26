derive_var_chg <- function(bds_dataset) {
  assert_has_variables(bds_dataset, c("AVAL", "BASE"))

  bds_dataset %>%
    mutate(CHG = AVAL - BASE)
}

derive_var_pchg <- function(bds_dataset) {
  assert_has_variables(bds_dataset, c("BASE", "CHG"))

  bds_dataset %>%
    mutate(PCHG = CHG / BASE)
}
