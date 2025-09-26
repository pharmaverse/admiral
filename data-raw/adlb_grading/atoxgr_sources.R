atoxgr_criteria <- system.file("adlb_grading/adlb_grading_spec.xlsx", package = "admiral")
# Contrary to our usual convention the use of `::` here is explicit. This way we
# avoid having to list {readxl} in "Imports" and instead get away with just
# listing it in "Depends".


# create NCICTCAEv4 metadata for SI units
atoxgr_criteria_ctcv4 <- atoxgr_criteria %>%
  readxl::read_excel(sheet = "NCICTCAEv4") %>%
  dplyr::mutate(GRADE_CRITERIA_CODE = gsub("[\r\n]", " ", GRADE_CRITERIA_CODE))

save(atoxgr_criteria_ctcv4, file = "data/atoxgr_criteria_ctcv4.rda")

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

# create DAIDs metadata for USCV units
atoxgr_criteria_daids_uscv <- atoxgr_criteria %>%
  readxl::read_excel(sheet = "DAIDS_CV") %>%
  dplyr::mutate(GRADE_CRITERIA_CODE = gsub("[\r\n]", " ", GRADE_CRITERIA_CODE))

save(atoxgr_criteria_daids_uscv, file = "data/atoxgr_criteria_daids_uscv.rda")

# Read in JSON file - for NCICTCAEv5 for CV units
json_file_path <- "data-raw/adlb_grading/ncictcaev5_cv.json"

atoxgr_criteria_ctcv5_uscv <- jsonlite::fromJSON(json_file_path) %>%
  filter(!is.na(GRADE_NA_CODE)) %>%
  dplyr::mutate(
    NEW_GRADE_CODE = "TRUE ~ \"0\"",
    NEW_GRADE_CODE = if_else(
      !is.na(GRADE_1_CODE),
      paste(GRADE_1_CODE, NEW_GRADE_CODE, sep = ", "),
      NEW_GRADE_CODE
    ),
    NEW_GRADE_CODE = if_else(
      !is.na(GRADE_2_CODE),
      paste(GRADE_2_CODE, NEW_GRADE_CODE, sep = ", "),
      NEW_GRADE_CODE
    ),
    NEW_GRADE_CODE = if_else(
      !is.na(GRADE_3_CODE),
      paste(GRADE_3_CODE, NEW_GRADE_CODE, sep = ", "),
      NEW_GRADE_CODE
    ),
    NEW_GRADE_CODE = if_else(
      !is.na(GRADE_4_CODE),
      paste(GRADE_4_CODE, NEW_GRADE_CODE, sep = ", "),
      NEW_GRADE_CODE
    ),
    NEW_GRADE_CODE = if_else(
      !is.na(GRADE_NA_CODE),
      paste(GRADE_NA_CODE, NEW_GRADE_CODE, sep = ", "),
      NEW_GRADE_CODE
    ),
    GRADE_CRITERIA_CODE = paste0("case_when(", NEW_GRADE_CODE, ")")
  )

save(atoxgr_criteria_ctcv5_uscv, file = "data/atoxgr_criteria_ctcv5_uscv.rda")

# Read in JSON file - for NCICTCAEv5 for CV units
json_file_path <- "data-raw/adlb_grading/ncictcaev5.json"

atoxgr_criteria_ctcv5 <- jsonlite::fromJSON(json_file_path) %>%
  filter(!is.na(GRADE_NA_CODE)) %>%
  dplyr::mutate(
    NEW_GRADE_CODE = "TRUE ~ \"0\"",
    NEW_GRADE_CODE = if_else(
      !is.na(GRADE_1_CODE),
      paste(GRADE_1_CODE, NEW_GRADE_CODE, sep = ", "),
      NEW_GRADE_CODE
    ),
    NEW_GRADE_CODE = if_else(
      !is.na(GRADE_2_CODE),
      paste(GRADE_2_CODE, NEW_GRADE_CODE, sep = ", "),
      NEW_GRADE_CODE
    ),
    NEW_GRADE_CODE = if_else(
      !is.na(GRADE_3_CODE),
      paste(GRADE_3_CODE, NEW_GRADE_CODE, sep = ", "),
      NEW_GRADE_CODE
    ),
    NEW_GRADE_CODE = if_else(
      !is.na(GRADE_4_CODE),
      paste(GRADE_4_CODE, NEW_GRADE_CODE, sep = ", "),
      NEW_GRADE_CODE
    ),
    NEW_GRADE_CODE = if_else(
      !is.na(GRADE_NA_CODE),
      paste(GRADE_NA_CODE, NEW_GRADE_CODE, sep = ", "),
      NEW_GRADE_CODE
    ),
    GRADE_CRITERIA_CODE = paste0("case_when(", NEW_GRADE_CODE, ")")
  )

save(atoxgr_criteria_ctcv5, file = "data/atoxgr_criteria_ctcv5.rda")

