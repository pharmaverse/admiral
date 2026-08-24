advs <- tribble(
  ~USUBJID,      ~PARAMCD, ~PARAM,                            ~AVAL, ~AVISIT,
  "01-701-1015", "DIABP",  "Diastolic Blood Pressure (mmHg)",    51, "BASELINE",
  "01-701-1015", "DIABP",  "Diastolic Blood Pressure (mmHg)",    50, "WEEK 2",
  "01-701-1015", "SYSBP",  "Systolic Blood Pressure (mmHg)",    121, "BASELINE",
  "01-701-1015", "SYSBP",  "Systolic Blood Pressure (mmHg)",    121, "WEEK 2",
  "01-701-1028", "DIABP",  "Diastolic Blood Pressure (mmHg)",    79, "BASELINE",
  "01-701-1028", "DIABP",  "Diastolic Blood Pressure (mmHg)",    80, "WEEK 2",
  "01-701-1028", "SYSBP",  "Systolic Blood Pressure (mmHg)",    130, "BASELINE",
  "01-701-1028", "SYSBP",  "Systolic Blood Pressure (mmHg)",     NA, "WEEK 2"
) %>%
  mutate(
    AVALU = "mmHg",
    ADT = case_when(
      AVISIT == "BASELINE" ~ as.Date("2024-01-10"),
      AVISIT == "WEEK 2" ~ as.Date("2024-01-24")
    ),
    ADTF = NA_character_
  )

advs1 <- derive_param_computed(
  advs,
  by_vars = exprs(USUBJID, AVISIT),
  parameters = c("SYSBP", "DIABP"),
  set_values_to = exprs(
    AVAL = (AVAL.SYSBP + 2 * AVAL.DIABP) / 3,
    PARAMCD = "MAP",
    PARAM = "Mean Arterial Pressure (mmHg)",
    AVALU = "mmHg",
    ADT = ADT.SYSBP
  )
)

summary(advs1)
