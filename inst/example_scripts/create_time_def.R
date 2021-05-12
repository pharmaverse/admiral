
library(dplyr)
library(lubridate)
data(dm)

# randomly assign number of periods per subject ID
set.seed(391)
max_period_per_usubjid <- sample(1:5, length(dm$USUBJID),
                                 replace = T,
                                 prob = rev(exp(1:5) / sum(exp(1:5))))

# generate period descriptions
period_desc <- paste("period", letters[1:5])

# derive time definitions dataset
set.seed(556)
time_def <- tibble(
  STUDYID = dm$STUDYID[1],
  USUBJID = rep(dm$USUBJID, max_period_per_usubjid),
  # for each subject ID add appropriate number of observations according to number of periods
  APERIOD = unlist(lapply(seq_along(dm$USUBJID), function(i) 1:max_period_per_usubjid[i])),
  APERIODC = period_desc[APERIOD]
) %>%
  left_join(select(dm, USUBJID, start = RFSTDTC, end = RFPENDTC), by = "USUBJID") %>%
  group_by(USUBJID) %>%
  mutate(order_within_subj = row_number(),
         # duration of each subjects participation
         duration = as_date(end) - as_date(start),
         # randomly assign period durations (ends) within the overall subject participation
         period_end_tmp = `if`(is.na(duration[1]),
                               NA,
                               sort(sample(1:as.numeric(duration), n(), replace = F))),
         # end of last duration needs to be the end of participation,
         # otherwise it is the start + duration
         period_end = ifelse(row_number() == n(),
                             end,
                             as.character(as_date(start) + days(period_end_tmp))),
         # start of a period needs to be just after the end of previous period
         period_start = as.character(as_date(dplyr::lag(period_end)) + days(1)),
         # start of the first period needs to be the start of the subjects participation
         period_start = ifelse(row_number() == 1, start, period_start),
         # generate explicit NAs where operations failed
         period_end = ifelse(period_end == "", NA, period_end),
         period_start = ifelse(period_start == "", NA, period_start)) %>%
  ungroup %>%
  left_join(select(dm, USUBJID, TRTP = ARMCD, TRTA = ARM), by = "USUBJID") %>%
  select(STUDYID, USUBJID, APERIOD, APERIODC, TRTP, TRTA,
         STARTDTM = period_start, ENDDTM = period_end)

warnings() # only date parsing -> fine ...

attr(time_def$STUDYID, "label") <- attr(dm$STUDYID, "label")
attr(time_def$USUBJID, "label") <- attr(dm$USUBJID, "label")

save(time_def, file = "data/time_def.rda")
