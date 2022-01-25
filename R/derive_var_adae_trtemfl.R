#' @title Derive Treatment Emergent Analysis Flag (TRTEMFL) for ADAE
#' 
#' @description
#' This function derives ADAE Treatment Emergent Analysis Flag (ADAE.TRTEMFL). \cr
#' The derivation is based on the default specifications for ADAE.TRTEMFL from the GDSR
#' (Version 2021-Oct.v1: Shenzou12). \cr
#' The following variables have to be available in the input data: 
#' TRTSDTM, ASTDTM, AENDTM, either AESEV and AESEVIN or AETOXGR and AEITOXGR depending on parameter
#' setting, if FAAE is present: AEGRPID
#'
#' @param dataset Input AE dataset
#' 
#' @param dataset_faae Input FAAE dataset (optional). The default is NULL
#' 
#' @param dataset_ex Input EX dataset
#' 
#' @param sevtox Name of variable in which severity or toxicity is stored. Permitted values are
#' "AESEV" or "AETOXGR". The default is "AESEV".
#' 
#' @param subject_keys The key variables to identify a subject provided as a list of variables. 
#' The default is vars(STUDYID, USUBJID).
#' 
#' @author Sandra Ott
#'
#' @details
#' The presence of 'TRTEMFL' is checked and a warning is issued if found. The existing variable
#' will be overwritten.
#' 
#' @return
#' The input AE dataset with the Treatment Emergent Analysis Flag 'TRTEMFL' added.
#'
#' @keywords adam flag derivation adae
#' 
#' @examples
#' \dontrun{
#' # Run function with default parameter setting and faae input
#' derive_var_adae_trtemfl(
#'   dataset = adae, 
#'   dataset_faae = faae, 
#'   dataset_ex = ex) 
#' 
#' # Run function with toxicity grade, custom subject key and no faae input
#' derive_var_adae_trtemfl(
#'   dataset = adae, 
#'   dataset_ex = ex,
#'   sevtox = "AETOXGR",
#'   subject_keys = vars(STUDYID, UNI_ID)
#' ) 
#' }
#' 
#' @export
#' 
derive_var_adae_trtemfl <- 
  function(dataset, dataset_faae = NULL, dataset_ex, sevtox = "AESEV", 
           subject_keys = vars(STUDYID, USUBJID)) {
    
    # Check input parameters
    assert_character_scalar(sevtox, values = c("AESEV", "AETOXGR"), case_sensitive = T)
    assert_vars(subject_keys)
    
    # Check input data
    assert_data_frame(dataset_ex, required_vars = vars(EXDOSE, EXTRT))
    
    assert_data_frame(dataset, required_vars = subject_keys)
    assert_data_frame(dataset, required_vars = vars(TRTSDTM, ASTDTM, AENDTM))
    if (sevtox == "AESEV") {
      assert_data_frame(dataset, required_vars = vars(AESEV, AESEVIN))
    } else {
      assert_data_frame(dataset, required_vars = vars(AETOXGR, AEITOXGR))
    }
    
    if (!is.null(dataset_faae)) {
      assert_data_frame(dataset_faae, required_vars = vars(FAGRPID, FAORRES, FATESTCD, FADTC))
      assert_data_frame(dataset, required_vars = vars(AEGRPID))
    }
    
    # Warn if TRTEMFL variable already exist
    check <- warn_if_vars_exist(dataset, c("TRTEMFL"))
    if (!is.null(check)) {
      warning("Existing variable will be overwritten")
      # Drop variables as they will be re-derived in the next step
      dataset <- select(dataset, -any_of(c("TRTEMFL")))
    }
    
    # Find out if subject had any valid dose 
    filter_ex <- '(EXDOSE > 0 | EXDOSE == 0 & str_detect(EXTRT, "PLACEBO"))'
    # Add filter conditions if needed
    ex_piltk <- names(select(dataset_ex, starts_with("EXPILTK") ))
    if (length(ex_piltk) >= 1) {
      filter_ex_piltk <- paste(ex_piltk[1], "!= '0'")
      if (length(ex_piltk) >= 2) {
        for (i in 2:length(ex_piltk)) {
          filter_ex_piltk <- paste(filter_ex_piltk, "|", ex_piltk[i], "!= '0'")
        }
      }
      filter_ex <-  paste0(filter_ex, " & (", filter_ex_piltk, ")")
    }
    if ("EXOCCUR" %in% names(dataset_ex)) {
        filter_ex <- paste(filter_ex, "& EXOCCUR != 'N'")
    }
    
    # Filter exposure data
    valid_dose <- dataset_ex %>% 
      filter(eval(parse(text = filter_ex))) %>%
      group_by(across(all_of(c(unlist(lapply(subject_keys, as_name)))))) %>% 
      summarise() %>% 
      mutate(VALDOSE="Y") %>% 
      select(unlist(lapply(subject_keys, as_name)), VALDOSE)
    
    # Merge Valid dose flag
    dataset01 <- merge(dataset, valid_dose, by = unlist(lapply(subject_keys, as_name)), all.x = T)
    
    
    # Prepare FAAE data if available
    # (Observations on FAAE are pre-filtered where [FAAE.FATESTCD] either ("SEVERITY" or "SEV") or
    # ("TOXICITY" or "GRADE") depending on user input and
    # [FAAE.FADTC] date parts are complete, and matched on Group ID where
    # [AE.AEGRPID] = [FAAE.FAGRPID]. If [FAAE.FADTC] has missing time parts then impute missing 
    # time with 23:59:59, partially missing times with 23 for missing hours, 59 for missing minutes, 
    # 59 for missing seconds.)

    if (sevtox == "AESEV") {
      filter_sevtox <- '(FATESTCD == "SEVERITY" | FATESTCD == "SEV") & nchar(word(FADTC, 1, sep="T")) == 10'
    } else {
      filter_sevtox <- '(FATESTCD == "TOXICITY" | FATESTCD == "GRADE") & nchar(word(FADTC, 1, sep="T")) == 10'
    }
    
    if (!is.null(dataset_faae)) {
      faae_sevtox <- dataset_faae %>% 
        filter(eval(parse(text = filter_sevtox))) 
      if (dim(faae_sevtox)[1] > 0) {
        faae_sevtox <- faae_sevtox %>% 
          mutate(SEVTOX = case_when(
            FAORRES == "MILD" ~ 1,
            FAORRES == "MODERATE" ~ 2,
            FAORRES == "SEVERE" ~ 3,
            TRUE ~ suppressWarnings(as.numeric(FAORRES))
            )
          ) %>% 
          # Convert FADTC to numeric date/time
          derive_vars_dtm(
            dtc = FADTC,
            new_vars_prefix = "FA",
            date_imputation = NULL,
            time_imputation = "last",
            flag_imputation = "none"
          ) %>% 
          filter(!is.na(FADTM) & !is.na(SEVTOX)) %>% 
          select(unlist(lapply(subject_keys, as_name)), AEGRPID=FAGRPID, FADTM, SEVTOX)
        faae_sevtox$SEVTOX <- as.integer(faae_sevtox$SEVTOX)
        
        
        # Add TRTSDTM 
        trtsdtm <- dataset01 %>% 
          select(unlist(lapply(subject_keys, as_name)), TRTSDTM, AEGRPID) %>% 
          unique()
        faae_sevtox_trt <- merge(faae_sevtox, trtsdtm, by = c(unlist(lapply(subject_keys, as_name)), "AEGRPID"),  
                             all.x = T)
        # Get the last value before TRTSDTM
        faae_sevtox_before <- faae_sevtox_trt %>% 
          filter(FADTM < TRTSDTM) %>% 
          group_by(across(all_of(c(unlist(lapply(subject_keys, as_name)), "AEGRPID")))) %>% 
          arrange("FADTM") %>%  
          slice(n()) %>% 
          select(-TRTSDTM, -FADTM, PRE_SEVTOX = SEVTOX)
        # Get the worst value after or at TRTSDTM
        faae_sevtox_after <- faae_sevtox_trt %>% 
          filter(FADTM >= TRTSDTM) 
        if (dim(faae_sevtox_after)[1] > 0) {
          faae_sevtox_after <- faae_sevtox_after %>% 
            group_by(across(all_of(c(unlist(lapply(subject_keys, as_name)), "AEGRPID")))) %>% 
            summarise(MAXSEVTOX=max(SEVTOX))
        } else {
          faae_sevtox_after$MAXSEVTOX = numeric(nrow(faae_sevtox_after))
        }
        
        # Combine information
        faae_sevtox_comb <- merge(faae_sevtox_before, faae_sevtox_after, 
                                  by = c(unlist(lapply(subject_keys, as_name)), "AEGRPID"),
                                  all = T)
        faae_sevtox_comb <- mutate(faae_sevtox_comb, FAAE = "Y")
  
        # Add on to AE data
        dataset02 <- merge(dataset01, faae_sevtox_comb, 
                           by = c(unlist(lapply(subject_keys, as_name)), "AEGRPID"),
                           all.x = T)
      } 
    } else {
        dataset02 <- dataset01 %>% 
          mutate(PRE_SEVTOX = NA_integer_, MAXSEVTOX = NA_integer_, FAAE = NA_character_)
    }
    
    # Derive TRTEMFL
    # Prepare condition based on user input and variables in data
    if (sevtox == "AESEV" ) {
      cond1 <- '!is.na(ASTDTM) & ASTDTM < TRTSDTM & 
            (AENDTM > TRTSDTM | is.na(AENDTM) ) & 
            ((!is.na(PRE_SEVTOX) & PRE_SEVTOX < MAXSEVTOX | AESEVIN_ < MAXSEVTOX ) | 
             (!is.na(PRE_SEVTOX) & is.na(MAXSEVTOX) & PRE_SEVTOX < AESEV_ ) |
              (is.na(FAAE) & (is.na(AESEV_) | is.na(AESEVIN_) | AESEV_ > AESEVIN_)) )'
    } else if (sevtox == "AETOXGR") {
      cond1 <- '!is.na(ASTDTM) & ASTDTM < TRTSDTM & 
            (AENDTM > TRTSDTM | is.na(AENDTM) ) & 
            ((!is.na(PRE_SEVTOX) & PRE_SEVTOX < MAXSEVTOX | AEITOXGR < MAXSEVTOX) | 
             (!is.na(PRE_SEVTOX) & is.na(MAXSEVTOX) & PRE_SEVTOX < AETOXGR ) |
              (is.na(FAAE) & (is.na(AETOXGR) | is.na(AEITOXGR) | AETOXGR > AEITOXGR)) )'
    } 
    # Convert Severity to numeric for compare
    if (sevtox == "AESEV" ) {
      dataset02 <- dataset02 %>% 
        mutate(
          AESEV_ = case_when(
              AESEV == "MILD" ~ 1L,
              AESEV == "MODERATE" ~ 2L,
              AESEV == "SEVERE" ~ 3L,
              TRUE ~ NA_integer_
            ),
          AESEVIN_ = case_when(
            AESEVIN == "MILD" ~ 1L,
            AESEVIN == "MODERATE" ~ 2L,
            AESEVIN == "SEVERE" ~ 3L,
            TRUE ~ NA_integer_
          )
        )
    } 
    
    dataset <- dataset02 %>% 
      mutate(
        TRTEMFL = case_when(
          # Set to "Y", if the Analysis Start Date/Time [ADAE.ASTDTM] is at or after Datetime of First 
          # Exposure to Treatment [ADSL.TRTSDTM].
          !is.na(TRTSDTM) & TRTSDTM <= ASTDTM ~ "Y",
          
          # Else set to null, if the Analysis End Date/Time [ADAE.AENDTM] is before Datetime of First 
          # Exposure to Treatment [ADSL.TRTSDTM].
          !is.na(AENDTM) & TRTSDTM > AENDTM ~ NA_character_,
          
          # Else set to "Y", if Analysis Start Date/Time [ADAE.ASTDTM] is before Datetime of First 
          # Exposure to Treatment [ADSL.TRTSDTM] but Analysis End Date/Time [ADAE.AENDTM] is either 
          # missing or at or after Datetime of First Exposure to Treatment [ADSL.TRTSDTM] and there is
          # an increase in Severity/Toxicity [FAAE.FAORRES] at a Date/Time of Collection [FAAE.FADTC] 
          # at or after Datetime of First Exposure to Treatment [ADSL.TRTSDTM] from the last 
          # Severity/Toxicity [FAAE.FAORRES] recorded at a Date/Time of Collection [FAAE.FADTC] 
          # before Datetime of First Exposure to Treatment [ADSL.TRTSDTM] or from Initial 
          # Severity/Toxicity [AE.AESEVIN or AE.AEITOXGR]. 
          # Or no related records in FAAE found, but initial Severity/Toxicity [AE.AESEVIN or AE.AEITOXGR]
          # is lower then most extreme Severity/Toxicity [AE.AESEV or AE.AETOXGR] or one of them or 
          # both are missing
           eval(parse(text = cond1))  ~ "Y",
          
          # Else set to null, if the Analysis Start Date/Time [ADAE.ASTDTM] is before Datetime of 
          # First Exposure to Treatment [ADSL.TRTSDTM] or if the subject does not have any valid 
          # dose. 
          !is.na(ASTDTM) & ASTDTM < TRTSDTM | is.na(VALDOSE) ~ NA_character_,
          
          # Set to "Y" in all other cases
          TRUE ~ "Y"
        
          
        )
      ) %>% 
      select(-any_of(c("VALDOSE", "PRE_SEVTOX", "MAXSEVTOX", "FAAE", "AESEVIN_", "AESEV_")))
    
    
    # Return updated dataset
    dataset
    
  }
