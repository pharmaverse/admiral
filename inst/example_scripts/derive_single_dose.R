
library(admiral)
library(dplyr)
library(lubridate)
data(ex)


ex %>% filter(EXDOSE %in% c(0, 54)) %>%
  group_by(USUBJID, VISIT) %>%
  summarise(check1 = n_distinct(EXSTDTC), check2 = n_distinct(EXENDTC)) %>%
  ungroup() %>%
  distinct(check1, check2)
# should only by 1 and 1

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
  mutate(EXSTDTC = seq_dates,
         EXENDTC = seq_dates) %>%
  group_by(USUBJID) %>%
  arrange(seq_dates) %>%
  mutate(EXSEQ = row_number()) %>%
  ungroup() %>%
  mutate(seq_dates = NULL)

attr(ex_single$EXSEQ, "label") <- attr(ex$EXSEQ, "label")
attr(ex_single$EXSTDTC, "label") <- attr(ex$EXSTDTC, "label")
attr(ex_single$EXENDTC, "label") <- attr(ex$EXENDTC, "label")
