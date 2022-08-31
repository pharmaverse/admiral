atoxgr_criteria_ctcv4 <- system.file("adlb_grading/adlb_grading_spec.xlsx", package = "admiral") %>%
  # Contrary to our usual convention the use of `::` here is explicit. This way we
  # avoid having to list {readxl} in "Imports" and instead get away with just
  # listing it in "Depends".
  readxl::read_excel(sheet = "NCICTCAEv4") %>%
  dplyr::mutate(GRADE_CRITERIA_CODE = gsub("[\r\n]", " ", GRADE_CRITERIA_CODE))

save(atoxgr_criteria_ctcv4, file = "data/atoxgr_criteria_ctcv4.rda")
