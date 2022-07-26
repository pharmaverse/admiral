
library(admiral.test)
library(admiral)
library(dplyr)
data(admiral_ex)
ex <- admiral_ex

# check that there is only one start/end date of exposure per subject and visit
check_cond <- ex %>%
  filter(EXDOSE %in% c(0, 54)) %>%
  group_by(USUBJID, VISIT) %>%
  summarise(check1 = n_distinct(EXSTDTC), check2 = n_distinct(EXENDTC)) %>%
  ungroup() %>%
  distinct(check1, check2)
stopifnot(check_cond$check1 == 1 && check_cond$check2 == 1)

dates <- ex %>%
  filter(EXDOSE %in% c(0, 54)) %>%
  mutate(EXSTDTC = as.Date(EXSTDTC), EXENDTC = as.Date(EXENDTC)) %>%
  filter(!is.na(EXSTDTC), !is.na(EXENDTC)) %>%
  group_by(USUBJID, VISIT) %>%
  do(tibble(seq_dates = seq(as.Date(.$EXSTDTC), as.Date(.$EXENDTC), by = "days"))) %>%
  ungroup()

attr(dates$USUBJID, "label") <- "Unique Subject Identifier"
attr(dates$VISIT, "label") <- "Visit Name"

ex_single <- dates %>%
  left_join(filter(ex, EXDOSE %in% c(0, 54)), by = c("USUBJID", "VISIT")) %>%
  # set dates as single doses
  mutate(EXSTDTC = seq_dates,
         EXENDTC = seq_dates) %>%
  group_by(USUBJID) %>%
  # adjust exseq
  arrange(seq_dates) %>%
  mutate(EXSEQ = row_number()) %>%
  # adjust EXSTDY and EXENDY
  mutate(EXSTDY = ifelse(EXSTDTC == min(EXSTDTC), 1, as.numeric(EXSTDTC - min(EXSTDTC)) + 1),
         EXENDY = EXSTDY) %>%
  ungroup() %>%
  mutate(seq_dates = NULL,
         EXENDTC = as.character(EXENDTC),
         EXSTDTC = as.character(EXSTDTC),
         EXDOSFRQ = "ONCE")

attr(ex_single$EXSEQ, "label") <- attr(ex$EXSEQ, "label")
attr(ex_single$EXSTDTC, "label") <- attr(ex$EXSTDTC, "label")
attr(ex_single$EXENDTC, "label") <- attr(ex$EXENDTC, "label")
attr(ex_single$EXDOSFRQ, "label") <- attr(ex$EXDOSFRQ, "label")

save(ex_single, file = file.path("data", "ex_single.rda"), compress = "bzip2")
