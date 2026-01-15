# Contrary to our usual convention the use of `::` here is explicit. This way we
# avoid having to list {fromJSON} in "Imports" and instead get away with just
# listing it in "Depends".

atoxgr_json_to_dataframe <- function(dataset, json_file) {
  dataset <- jsonlite::fromJSON(json_file) %>%
    dplyr::mutate(
      NEW_GRADE_CODE = if_else(
        !is.na(GRADE_NA_CODE),
        "TRUE ~ \"0\"",
        NA_character_
      ),
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
      GRADE_CRITERIA_CODE = if_else(
        !is.na(NEW_GRADE_CODE),
        paste0("case_when(", NEW_GRADE_CODE, ")"),
        "NA_character_"
      )
    ) %>%
    select(
      -NEW_GRADE_CODE, -GRADE_NA_CODE, -GRADE_1_CODE, -GRADE_2_CODE,
      -GRADE_3_CODE, , -GRADE_4_CODE
    )
}

# hold core directory path in an object
json_file_path <- "data-raw/adlb_grading"
save_file_path <- "data"

# Point to JSON files for each criteria

## NCICTCAEv4 for SI units
json_file_path_v4a <- file.path(json_file_path, "ncictcaev4.json")

## NCICTCAEv4 for CV units
json_file_path_v4_uscv <- file.path(json_file_path, "ncictcaev4_uscv.json")

## NCICTCAEv5 for SI units
json_file_path_v5 <- file.path(json_file_path, "ncictcaev5.json")

## NCICTCAEv5 for CV units
json_file_path_v5_uscv <- file.path(json_file_path, "ncictcaev5_uscv.json")

## NCICTCAEv6 for SI units
json_file_path_v6 <- file.path(json_file_path, "ncictcaev6.json")

## NCICTCAEv6 for CV units
json_file_path_v6_uscv <- file.path(json_file_path, "ncictcaev6_uscv.json")

## DAIDs for SI units
json_file_path_daids <- file.path(json_file_path, "DAIDS.json")

## DAIDs for CV units
json_file_path_daids_uscv <- file.path(json_file_path, "DAIDS_uscv.json")

## Read in JSON files and create metadata for each criteria

atoxgr_criteria_ctcv4 <- atoxgr_json_to_dataframe(json_file = json_file_path_v4a)
atoxgr_criteria_ctcv4_uscv <- atoxgr_json_to_dataframe(json_file = json_file_path_v4_uscv)
atoxgr_criteria_ctcv5 <- atoxgr_json_to_dataframe(json_file = json_file_path_v5)
atoxgr_criteria_ctcv5_uscv <- atoxgr_json_to_dataframe(json_file = json_file_path_v5_uscv)
atoxgr_criteria_ctcv6 <- atoxgr_json_to_dataframe(json_file = json_file_path_v6)
atoxgr_criteria_ctcv6_uscv <- atoxgr_json_to_dataframe(json_file = json_file_path_v6_uscv)
atoxgr_criteria_daids <- atoxgr_json_to_dataframe(json_file = json_file_path_daids)
atoxgr_criteria_daids_uscv <- atoxgr_json_to_dataframe(json_file = json_file_path_daids_uscv)

## Save metadata for each criteria

save(atoxgr_criteria_ctcv4,
  file = file.path(save_file_path, "atoxgr_criteria_ctcv4.rda"),
  compress = "bzip2"
)

save(atoxgr_criteria_ctcv4_uscv,
  file = file.path(save_file_path, "atoxgr_criteria_ctcv4_uscv.rda"),
  compress = "bzip2"
)

save(atoxgr_criteria_ctcv5,
  file = file.path(save_file_path, "atoxgr_criteria_ctcv5.rda"),
  compress = "bzip2"
)

save(atoxgr_criteria_ctcv5_uscv,
  file = file.path(save_file_path, "atoxgr_criteria_ctcv5_uscv.rda"),
  compress = "bzip2"
)

save(atoxgr_criteria_ctcv6,
  file = file.path(save_file_path, "atoxgr_criteria_ctcv6.rda"),
  compress = "bzip2"
)

save(atoxgr_criteria_ctcv6_uscv,
  file = file.path(save_file_path, "atoxgr_criteria_ctcv6_uscv.rda"),
  compress = "bzip2"
)

save(atoxgr_criteria_daids,
  file = file.path(save_file_path, "atoxgr_criteria_daids.rda"),
  compress = "bzip2"
)

save(atoxgr_criteria_daids_uscv,
  file = file.path(save_file_path, "atoxgr_criteria_daids_uscv.rda"),
  compress = "bzip2"
)
