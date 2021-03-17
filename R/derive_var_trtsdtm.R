derive_var_trtsdtm <- function(dataset,
                               dataset_ex = ex,
                               filter_ex = exprs(EXDOSE > 0 | (EXDOSE == 0 & str_detect(EXTRT, 'PLACEBO')))){
  derive_merged_vars(dataset,
                     dataset_add = dataset_ex,
                     filter_add = filter_ex,
                     new_vars = exprs(TRTSDTM := ymd(EXSTDTC)),
                     filter_first_order = exprs(EXSTDTC, EXSEQ))
}
