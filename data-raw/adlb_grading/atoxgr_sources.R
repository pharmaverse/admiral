atoxgr_criteria <- system.file("adlb_grading/adlb_grading_spec.xlsx", package = "admiral")
# Contrary to our usual convention the use of `::` here is explicit. This way we
# avoid having to list {readxl} in "Imports" and instead get away with just
# listing it in "Depends".


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


atoxgr_json_to_dataframe <- function(dataset, json_file) {
  dataset <- jsonlite::fromJSON(json_file) %>%
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
    ) %>%
    select(-NEW_GRADE_CODE)
}

library(dplyr)

# Read in JSON file - for NCICTCAEv4 for SI units
json_file_path_v4 <- "data-raw/adlb_grading/ncictcaev4.json"

atoxgr_criteria_ctcv4 <- atoxgr_json_to_dataframe(json_file = json_file_path_v4)

save(atoxgr_criteria_ctcv4, file = "data/atoxgr_criteria_ctcv4.rda")

# Read in JSON file - for NCICTCAEv5 for SI units
json_file_path_v5 <- "data-raw/adlb_grading/ncictcaev5.json"

atoxgr_criteria_ctcv5 <- atoxgr_json_to_dataframe(json_file = json_file_path_v5)

save(atoxgr_criteria_ctcv5, file = "data/atoxgr_criteria_ctcv5.rda")

# Read in JSON file - for NCICTCAEv5 for CV units
json_file_path_v5_cv <- "data-raw/adlb_grading/ncictcaev5_cv.json"

atoxgr_criteria_ctcv5_uscv <- atoxgr_json_to_dataframe(json_file = json_file_path_v5_cv)

save(atoxgr_criteria_ctcv5_uscv, file = "data/atoxgr_criteria_ctcv5_uscv.rda")
