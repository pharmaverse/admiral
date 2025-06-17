atoxgr_criteria <- system.file("adlb_grading/adlb_grading_spec.xlsx", package = "admiral")
# Contrary to our usual convention the use of `::` here is explicit. This way we
# avoid having to list {readxl} in "Imports" and instead get away with just
# listing it in "Depends".


# create NCICTCAEv4 metadata for SI units
atoxgr_criteria_ctcv4 <- atoxgr_criteria %>%
  readxl::read_excel(sheet = "NCICTCAEv4") %>%
  dplyr::mutate(GRADE_CRITERIA_CODE = gsub("[\r\n]", " ", GRADE_CRITERIA_CODE))

save(atoxgr_criteria_ctcv4, file = "data/atoxgr_criteria_ctcv4.rda")

# create NCICTCAEv5 metadata for SI units
atoxgr_criteria_ctcv5 <- atoxgr_criteria %>%
  readxl::read_excel(sheet = "NCICTCAEv5") %>%
  dplyr::mutate(GRADE_CRITERIA_CODE = gsub("[\r\n]", " ", GRADE_CRITERIA_CODE))

save(atoxgr_criteria_ctcv5, file = "data/atoxgr_criteria_ctcv5.rda")

# create DAIDs metadata for SI units
atoxgr_criteria_daids <- atoxgr_criteria %>%
  readxl::read_excel(sheet = "DAIDS") %>%
  dplyr::mutate(GRADE_CRITERIA_CODE = gsub("[\r\n]", " ", GRADE_CRITERIA_CODE))

save(atoxgr_criteria_daids, file = "data/atoxgr_criteria_daids.rda")

# create NCICTCAEv4 metadata for USCV units
atoxgr_criteria_ctcv4_uscv <- atoxgr_criteria %>%
  readxl::read_excel(sheet = "NCICTCAEv4_CV") %>%
  dplyr::mutate(GRADE_CRITERIA_CODE = gsub("[\r\n]", " ", GRADE_CRITERIA_CODE))

save(atoxgr_criteria_ctcv4_uscv, file = "data/atoxgr_criteria_ctcv4_uscv.rda")

# create NCICTCAEv5 metadata for USCV units
atoxgr_criteria_ctcv5_uscv <- atoxgr_criteria %>%
  readxl::read_excel(sheet = "NCICTCAEv5_CV") %>%
  dplyr::mutate(GRADE_CRITERIA_CODE = gsub("[\r\n]", " ", GRADE_CRITERIA_CODE))

save(atoxgr_criteria_ctcv5_uscv, file = "data/atoxgr_criteria_ctcv5_uscv.rda")

# create DAIDs metadata for USCV units
atoxgr_criteria_daids_uscv <- atoxgr_criteria %>%
  readxl::read_excel(sheet = "DAIDS_CV") %>%
  dplyr::mutate(GRADE_CRITERIA_CODE = gsub("[\r\n]", " ", GRADE_CRITERIA_CODE))

save(atoxgr_criteria_daids_uscv, file = "data/atoxgr_criteria_daids_uscv.rda")
