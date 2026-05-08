# admiral guidelines loaded

# Shared fixtures for `test-derive_var_atoxgr.R`.
# Source this file inside a `test_that()` block to give each test a fresh
# local copy of the expected data objects.

exp_anemia_si <- tibble::tribble(
  ~ATOXDSCL,      ~AVAL,  ~ANRLO, ~ANRHI, ~AVALU, ~ATOXGRL, ~TESTNUM,
  "Not a term",   80,     120,    200,    "G/L",  NA,       1,
  NA_character_,  60,     50,     100,    "G/L",  NA,       2,
  "ANEMIA",       79,     140,    NA,     "G/L",  "3",      3,
  "ANEMIA",       80,     140,    NA,     "G/L",  "2",      4,
  "Anemia",       99,     140,    NA,     "G/L",  "2",      5,
  "Anemia",       100,    140,    NA,     "G/L",  "1",      6,
  # wrong UNIT - GRADE should be missing
  "anemia",       100,    140,    NA,     "G/cL", NA,       7,
  "Anemia",       139,    140,    NA,     "G/L",  "1",      8,
  "ANEMIA",       140,    140,    NA,     "g/L",  "0",      9,
  # ANRLO missing AVAL not grade 2 or 3 - cannot grade
  "Anemia",       140,    NA,     NA,     "G/L",  NA,       10,
  "Anemia",       139,    NA,     NA,     "G/L",  NA,       11,
  "Anemia",       100,    NA,     NA,     "G/L",  NA,       12,
  # ANRLO missing but AVAL satisfies grade 2
  "Anemia",       99,     NA,     NA,     "G/L",  "2",      13,
  # ANRLO missing but AVAL satisfies grade 3
  "Anemia",       79,     NA,     NA,     "G/L",  "3",      14,
  # Unit missing cannot grade
  "Anemia",       140,    140,    NA,     NA,     NA,       15,
  # AVAL missing cannot grade
  "Anemia",       NA,     140,    NA,     "G/L",  NA,       16,
)

exp_anemia_cv <- exp_anemia_si %>%
  mutate(
    AVAL = if_else(is.na(AVAL), NA_real_, AVAL / 10),
    ANRLO = if_else(is.na(ANRLO), NA_real_, ANRLO / 10),
    ANRHI = if_else(is.na(ANRHI), NA_real_, ANRHI / 10),
    AVALU = if_else(toupper(AVALU) == "G/L", "g/dL", AVALU)
  )

exp_leuko_si <- tibble::tribble(
  ~ATOXDSCH,      ~AVAL,  ~ANRLO, ~ANRHI, ~AVALU,    ~ATOXGRH, ~TESTNUM,
  "Not a term",   99,     0,      NA,     "10^9/L",  NA,       1,
  NA,             99,     0,      NA,     "10^9/L",  NA,       2,
  "Leukocytosis", 101,    0,      40,     "10^9/L",  "3",      3,
  "leukocytosis", 100,    0,      40,     "10^9/L",  "0",      4,
  "Leukocytosis", 99,     0,      NA,     "10^9/L",  "0",      5,
  # wrong UNIT - GRADE should be missing
  "Leukocytosis", 99,     0,      40,     "10^9/M",  NA,       6,
  # Unit missing cannot grade
  "Leukocytosis", 99,     0,      40,     NA,        NA,       7,
  # AVAL missing cannot grade
  "Leukocytosis", NA,     0,      40,     "10^9/L",  NA,       8,
)

exp_leuko_cv <- exp_leuko_si %>%
  mutate(
    AVALU = if_else(str_to_upper(AVALU) == "10^9/L", "10^3/uL", AVALU)
  )

exp_leuko_cv2 <- exp_leuko_si %>%
  mutate(
    AVAL = 1000 * AVAL,
    ANRHI = 1000 * ANRHI,
    AVALU = if_else(str_to_upper(AVALU) == "10^9/L", "10^3/mL", AVALU)
  )

exp_aptt <- tibble::tribble(
  ~ATOXDSCH,                                         ~AVAL,  ~ANRHI,  ~ATOXGRH, ~TESTNUM,
  "Not a term",                                      80,     100,     NA,       1,
  NA_character_,                                     60,     100,     NA,       2,
  "Activated partial thromboplastin time prolonged", 251,    100,     "3",      3,
  "Activated Partial thromboplastin time prolonged", 250,    100,     "2",      4,
  "Activated partial Thromboplastin time prolonged", 151,    100,     "2",      5,
  "Activated partial thromboplastin time prolonged", 150,    100,     "1",      6,
  "Activated partial thromboplastin Time prolonged", 101,    100,     "1",      7,
  "Activated partial thromboplastin time prolonged", 100,    100,     "0",      8,
  # ANRHI missing - cannot grade
  "Activated partial thromboplastin time prolonged", 100,    NA,      NA,       9,
  # AVAL missing cannot grade
  "Activated partial thromboplastin time prolonged", NA,     100,     NA,       10,
) %>%
  mutate(AVALU = NA_character_)

exp_alt_ctcv4 <- tibble::tribble(
  ~ATOXDSCH,                            ~AVAL,  ~ANRHI, ~ATOXGRH, ~TESTNUM,
  "Not a term",                         80,     40,     NA,       1,
  NA_character_,                        60,     40,     NA,       2,
  "Alanine aminotransferase Increased", 801,    40,     "4",      3,
  "Alanine aminotransferase Increased", 800,    40,     "3",      4,
  "Alanine aminotransferase Increased", 201,    40,     "3",      5,
  "Alanine aminotransferase Increased", 200,    40,     "2",      6,
  "Alanine aminotransferase Increased", 121,    40,     "2",      7,
  "Alanine aminotransferase Increased", 120,    40,     "1",      8,
  "Alanine aminotransferase Increased", 41,     40,     "1",      9,
  "Alanine aminotransferase Increased", 40,     40,     "0",      10,
  # ANRHI missing - cannot grade
  "Alanine aminotransferase Increased", 100,    NA,     NA,       11,
  # AVAL missing cannot grade
  "Alanine aminotransferase Increased", NA,     40,     NA,       12,
) %>%
  mutate(AVALU = NA_character_)

exp_alkp_ctcv4 <- tibble::tribble(
  ~ATOXDSCH,                         ~AVAL,  ~ANRHI, ~ATOXGRH, ~TESTNUM,
  "Not a term",                      80,     40,     NA,       1,
  NA_character_,                     60,     40,     NA,       2,
  "Alkaline phosphatase increased",  801,    40,     "4",      3,
  "Alkaline phosphatase increased",  800,    40,     "3",      4,
  "Alkaline phosphatase increased",  201,    40,     "3",      5,
  "Alkaline phosphatase increased",  200,    40,     "2",      6,
  "Alkaline phosphatase increased",  101,    40,     "2",      7,
  "Alkaline phosphatase increased",  100,    40,     "1",      8,
  "Alkaline phosphatase increased",  41,     40,     "1",      9,
  "Alkaline phosphatase increased",  40,     40,     "0",      10,
  # ANRHI missing - cannot grade
  "Alkaline phosphatase increased",  100,    NA,     NA,       11,
  # AVAL missing cannot grade
  "Alkaline phosphatase increased",  NA,     40,     NA,       12,
) %>%
  mutate(AVALU = NA_character_)

exp_alkp_ctcv6 <- tibble::tribble(
  ~ATOXDSCH,                         ~AVAL,  ~BASE,  ~ANRHI,  ~ATOXGRH, ~TESTNUM,
  "Not a term",                      80,     40,     40,      NA,       1,
  NA_character_,                     60,     40,     40,      NA,       2,
  "Alkaline phosphatase increased",  41,     40,     40,      "1",      3,
  "Alkaline phosphatase increased",  40,     40,     40,      "0",      4,
  "Alkaline phosphatase increased",  41,     41,     40,      "0",      5,
  "Alkaline phosphatase increased",  41,     40,     41,      "0",      6,
  "Alkaline phosphatase increased",  40,     40,     NA,      "0",      7,
  "Alkaline phosphatase increased",  40,     NA,     40,      "0",      8,
  # ANRHI missing - cannot grade
  "Alkaline phosphatase increased",  100,    40,     NA,      NA,       9,
  # BASE missing - cannot grade
  "Alkaline phosphatase increased",  41,     NA,     40,      NA,       10,
  # AVAL missing cannot grade
  "Alkaline phosphatase increased",  NA,     40,     40,      NA,       11,
) %>%
  mutate(AVALU = NA_character_)

exp_ast_ctcv4 <- tibble::tribble(
  ~ATOXDSCH,                              ~AVAL,  ~ANRHI, ~ATOXGRH, ~TESTNUM,
  "Not a term",                           80,     40,     NA,       1,
  NA_character_,                          60,     40,     NA,       2,
  "Aspartate aminotransferase Increased", 801,    40,     "4",      3,
  "Aspartate aminotransferase Increased", 800,    40,     "3",      4,
  "Aspartate aminotransferase Increased", 201,    40,     "3",      5,
  "Aspartate aminotransferase Increased", 200,    40,     "2",      6,
  "Aspartate aminotransferase Increased", 121,    40,     "2",      7,
  "Aspartate aminotransferase Increased", 120,    40,     "1",      8,
  "Aspartate aminotransferase Increased", 41,     40,     "1",      9,
  "Aspartate aminotransferase Increased", 40,     40,     "0",      10,
  # ANRHI missing - cannot grade
  "Aspartate aminotransferase Increased", 100,    NA,     NA,       11,
  # AVAL missing cannot grade
  "Aspartate aminotransferase Increased", NA,     40,     NA,       12,
) %>%
  mutate(AVALU = NA_character_)

exp_bili_ctcv4 <- tibble::tribble(
  ~ATOXDSCH,                   ~AVAL,  ~ANRHI, ~ATOXGRH, ~TESTNUM,
  "Not a term",                80,     40,     NA,       1,
  NA_character_,               60,     40,     NA,       2,
  "Blood bilirubin increased", 401,    40,     "4",      3,
  "Blood bilirubin increased", 400,    40,     "3",      4,
  "Blood bilirubin increased", 121,    40,     "3",      5,
  "Blood bilirubin increased", 120,    40,     "2",      6,
  "Blood bilirubin increased", 61,     40,     "2",      7,
  "Blood bilirubin increased", 60,     40,     "1",      8,
  "Blood bilirubin increased", 41,     40,     "1",      9,
  "Blood bilirubin increased", 40,     40,     "0",      10,
  # ANRHI missing - cannot grade
  "Blood bilirubin increased", 100,    NA,     NA,       11,
  # AVAL missing cannot grade
  "Blood bilirubin increased", NA,     40,     NA,       12,
) %>%
  mutate(AVALU = NA_character_)

exp_bili_ctcv6_norm <- exp_bili_ctcv4 %>%
  # set BASE to be normal (not HIGH) and create FLAG
  mutate(
    BASE = ANRHI,
    BNRIND = "NORMAL"
  )

exp_bili_ctcv6 <- tibble::tribble(
  ~ATOXDSCH,                   ~AVAL,  ~BASE,  ~ATOXGRH, ~TESTNUM,
  "Not a term",                80,     40,     NA,       13,
  NA_character_,               60,     40,     NA,       14,
  "Blood bilirubin increased", 401,    40,     "4",      15,
  "Blood bilirubin increased", 400,    40,     "3",      16,
  "Blood bilirubin increased", 101,    40,     "3",      17,
  "Blood bilirubin increased", 100,    40,     "2",      18,
  "Blood bilirubin increased", 61,     40,     "2",      19,
  "Blood bilirubin increased", 60,     40,     "1",      20,
  "Blood bilirubin increased", 40,     40,     "1",      21,
  "Blood bilirubin increased", 39,     40,     "0",      22,
  # ANRHI missing - cannot grade
  "Blood bilirubin increased", 100,    NA,     NA,       23,
  # AVAL missing cannot grade
  "Blood bilirubin increased", NA,     40,     NA,       24,
) %>%
  # set BASE to ANRHI then make ANRHI < BASE
  mutate(
    AVALU = NA_character_,
    ANRHI = BASE - 1,
    BNRIND = "HIGH",
    TESTNUM = TESTNUM + 12
  ) %>%
  bind_rows(exp_bili_ctcv6_norm)

exp_cd4_si <- tibble::tribble(
  ~ATOXDSCL,                    ~AVAL,  ~ANRLO, ~AVALU,    ~ATOXGRL, ~TESTNUM,
  "Not a term",                 80,     120,    "10^9/L",  NA,       1,
  NA_character_,                60,     50,     "10^9/L",  NA,       2,
  "CD4 lymphocytes decreased",  0.04,   0.8,    "10^9/L",  "4",      3,
  "CD4 lymphocytes decreased",  0.05,   0.8,    "10^9/l",  "3",      4,
  "CD4 lymphocytes decreased",  0.19,   0.8,    "10^9/L",  "3",      5,
  "CD4 lymphocytes decreased",  0.2,    0.8,    "10^9/L",  "2",      6,
  "CD4 lymphocytes decreased",  0.49,   0.8,    "10^9/L",  "2",      7,
  # wrong unit - grade missing
  "CD4 lymphocytes decreased",  0.49,   0.8,    "10^8/L",  NA,       8,
  "CD4 lymphocytes decreased",  0.5,    0.8,    "10^9/L",  "1",      9,
  "CD4 lymphocytes decreased",  0.79,   0.8,    "10^9/L",  "1",      10,
  "CD4 lymphocytes decreased",  0.8,    0.8,    "10^9/L",  "0",      11,
  # ANRLO missing - AVAL satisfies grade 2 - 4
  "CD4 lymphocytes decreased",  0.04,   NA,     "10^9/L",  "4",      12,
  "CD4 lymphocytes decreased",  0.05,   NA,     "10^9/L",  "3",      13,
  "CD4 lymphocytes decreased",  0.19,   NA,     "10^9/L",  "3",      14,
  "CD4 lymphocytes decreased",  0.2,    NA,     "10^9/L",  "2",      15,
  "CD4 lymphocytes decreased",  0.49,   NA,     "10^9/L",  "2",      16,
  # ANRLO missing - AVAL does NOT satisfies grade 2 - 4
  "CD4 lymphocytes decreased",  0.5,    NA,     "10^9/L",  NA,       17,
  "CD4 lymphocytes decreased",  0.79,   NA,     "10^9/L",  NA,       18,
  "CD4 lymphocytes decreased",  0.8,    NA,     "10^9/L",  NA,       19,
  # Unit missing - cannot grade
  "CD4 lymphocytes decreased",  0.8,    0.8,    NA,        NA,       20,
  # AVAL missing cannot grade
  "CD4 lymphocytes decreased",  NA,     0.8,    "10^9/L",  NA,       21,
)

exp_cd4_cv <- exp_cd4_si %>%
  mutate(
    AVAL = AVAL * 1000,
    ANRLO = ANRLO * 1000,
    AVALU = if_else(str_to_upper(AVALU) == "10^9/L", "1/uL", NA_character_)
  )

exp_choles_si <- tibble::tribble(
  ~ATOXDSCH,          ~AVAL,  ~ANRLO, ~ANRHI, ~AVALU,    ~ATOXGRH, ~TESTNUM,
  "Not a term",       8,      0,      5,      "mmol/L",  NA,       1,
  NA_character_,      10,     0,      5,      "mmol/L",  NA,       2,
  "Cholesterol high", 12.93,  0,      5,      "mmol/L",  "4",      3,
  "Cholesterol High", 12.92,  0,      5,      "mmol/L",  "3",      4,
  "Cholesterol high", 10.35,  0,      5,      "Mmol/L",  "3",      5,
  # wrong unit - grade missing
  "Cholesterol high", 10.35,  0,      5,      "umol/L",  NA,       6,
  "Cholesterol high", 10.34,  0,      5,      "mmol/L",  "2",      7,
  "Cholesterol high", 7.76,   0,      5,      "mmol/L",  "2",      8,
  "Cholesterol high", 7.75,   0,      5,      "mmol/L",  "1",      9,
  "Cholesterol high", 5.1,    0,      5,      "mmol/L",  "1",      10,
  "Cholesterol high", 5,      0,      5,      "mmol/L",  "0",      11,
  # ANRHI missing - AVAL satisfies grade 2 - 4
  "Cholesterol high", 12.93,  0,      NA,     "mmol/L",  "4",      12,
  "Cholesterol High", 12.92,  0,      NA,     "mmol/L",  "3",      13,
  "Cholesterol high", 10.35,  0,      NA,     "Mmol/L",  "3",      14,
  "Cholesterol high", 10.34,  0,      NA,     "mmol/L",  "2",      15,
  "Cholesterol high", 7.76,   0,      NA,     "mmol/L",  "2",      16,
  # ANRHI missing - AVAL does NOT satisfies grade 2 - 4
  "Cholesterol high", 7.75,   0,      NA,     "mmol/L",  NA,       17,
  "Cholesterol high", 5.1,    0,      NA,     "mmol/L",  NA,       18,
  "Cholesterol high", 5,      0,      NA,     "mmol/L",  NA,       19,
  # Unit missing - cannot grade
  "Cholesterol high", 5,      0,      5,      NA,        NA,       20,
  # AVAL missing cannot grade
  "Cholesterol high", NA,     0,      5,      "mmol/L",  NA,       21,
)

exp_choles_cv <- tibble::tribble(
  ~ATOXDSCH,           ~AVAL,  ~ANRLO,  ~ANRHI,   ~AVALU,  ~ATOXGRH, ~TESTNUM,
  "Not a term",        8,      0,       200,     "mg/dL",        NA, 1,
  NA_character_,       10,     0,       200,     "mg/dL",        NA, 2,
  "Cholesterol high",  501,    0,       200,     "mg/dL",       "4", 3,
  "Cholesterol High",  500,    0,       200,     "mg/dL",       "3", 4,
  "Cholesterol high",  401,    0,       200,     "Mg/dL",       "3", 5,
  # wrong unit - grade missing
  "Cholesterol high",  400,    0,       200,     "ug/dL",        NA, 6,
  "Cholesterol high",  400,    0,       200,     "mg/dL",       "2", 7,
  "Cholesterol high",  301,    0,       200,     "mg/dL",       "2", 8,
  "Cholesterol high",  300,    0,       200,     "mg/dL",       "1", 9,
  "Cholesterol high",  201,    0,       200,     "mg/dL",       "1", 10,
  "Cholesterol high",  200,    0,       200,     "mg/dL",       "0", 11,
  # ANRHI missing - AVAL satisfies grade 2 - 4
  "Cholesterol high",  501,    0,       NA,      "mg/dL",       "4", 12,
  "Cholesterol High",  500,    0,       NA,      "mg/dL",       "3", 13,
  "Cholesterol high",  401,    0,       NA,      "MG/dL",       "3", 14,
  "Cholesterol high",  400,    0,       NA,      "mg/dL",       "2", 15,
  "Cholesterol high",  301,    0,       NA,      "mg/dL",       "2", 16,
  # ANRHI missing - AVAL does NOT satisfies grade 2 - 4
  "Cholesterol high",  300,    0,       NA,      "mg/dL",        NA, 17,
  "Cholesterol high",  201,    0,       NA,      "mg/dL",        NA, 18,
  "Cholesterol high",  200,    0,       NA,      "mg/dL",        NA, 19,
  # Unit missing - cannot grade
  "Cholesterol high",  200,    0,       200,          NA,        NA, 20,
  # AVAL missing cannot grade
  "Cholesterol high",  NA,     0,       200,    "mmol/L",        NA, 21,
)

exp_cpk <- tibble::tribble(
  ~ATOXDSCH,        ~AVAL,  ~ANRLO,  ~ANRHI, ~ATOXGRH, ~TESTNUM,
  "Not a term",     80,     0,       40,     NA,       1,
  NA_character_,    60,     0,       40,     NA,       2,
  "CPK increased",  401,    0,       40,     "4",      3,
  "CPK increased",  400,    0,       40,     "3",      4,
  "CPK increased",  201,    0,       40,     "3",      5,
  "CPK increased",  200,    0,       40,     "2",      6,
  "CPK increased",  101,    0,       40,     "2",      7,
  "CPK increased",  100,    0,       40,     "1",      8,
  "CPK increased",  41,     0,       40,     "1",      9,
  "CPK increased",  40,     0,       40,     "0",      10,
  # ANRHI missing - cannot grade
  "CPK increased",  100,    0,       NA,     NA,       11,
  # AVAL missing cannot grade
  "CPK increased",  NA,     0,       40,     NA,       12,
) %>%
  mutate(AVALU = NA_character_)

exp_creatn <- tibble::tribble(
  ~ATOXDSCH,               ~AVAL,  ~BASE, ~ANRHI, ~ATOXGRH, ~V4, ~V5, ~TESTNUM,
  "Not a term",            80,     80,    40,     NA,       "Y", "Y", 1,
  NA_character_,           60,     60,    40,     NA,       "Y", "Y", 2,
  # GRADE derived from AVAL against ANRHI
  "Creatinine increased",  241,    241,   40,     "4",      "Y", "Y", 3,
  "Creatinine increased",  240,    230,   40,     "3",      "Y", "Y", 4,
  "Creatinine increased",  121,    120,   40,     "3",      "Y", "Y", 5,
  "Creatinine increased",  120,    119,   40,     "2",      "Y", "Y", 6,
  "Creatinine increased",  61,     60,    40,     "2",      "Y", "Y", 7,
  "Creatinine increased",  60,     60,    40,     "1",      "Y", "Y", 8,
  "Creatinine increased",  41,     41,    40,     "1",      "Y", "Y", 9,
  "Creatinine increased",  40,     40,    40,     "0",      "Y", "Y", 10,
  # GRADE derived from AVAL against BASE
  "Creatinine increased",  42,     6,     40,     "3",      "Y", "Y", 11,
  "Creatinine increased",  42,     13.9,  40,     "3",      "Y", "Y", 12,
  "Creatinine increased",  42,     14,    40,     "2",      "Y", "Y", 13,
  "Creatinine increased",  42.1,   28,    40,     "2",      "Y", "Y", 14,
  "Creatinine increased",  42,     28,    42,     "1",      "Y", "N", 15,
  "Creatinine increased",  42,     41,    42,     "1",      "Y", "N", 16,
  "Creatinine increased",  42,     42,    42,     "0",      "Y", "N", 17,
  # BASE missing - AVAL <= ANRLO cannot grade as NORMAL
  "Creatinine increased",  42,     NA,    42,     NA,       "Y", "N", 18,
  # ANRHI missing - AVAL <= BASE cannot grade as NORMAL
  "Creatinine increased",  42,     42,    NA,     NA,       "Y", "Y", 19,
  # AVAL missing cannot grade
  "Creatinine increased",  NA,     0,     40,     NA,       "Y", "Y", 20,
) %>%
  mutate(AVALU = NA_character_)

exp_fib <- tibble::tribble(
  ~ATOXDSCL,               ~AVAL,  ~ANRLO,  ~PCHG,  ~AVALU,  ~ATOXGRL, ~TESTNUM,
  "Not a term",            9,      10,      40,      "g/L",        NA, 1,
  NA_character_,           10,     10,      40,      "g/L",        NA, 2,
  # Satisfies < 0.5 for grade 4 - other criteria missing
  "Fibrinogen decreased",  0.49,   NA,      NA,      "g/L",       "4", 3,
  # Satisfies < 0.5 for grade 4 - satisfies grade 3 for other criteria
  "Fibrinogen decreased",  0.49,   1,       -51,     "g/L",       "4", 4,
  # Satisfies < 0.25*LLN for grade 4 - PCHG missing
  "Fibrinogen decreased",  0.5,    2.1,     NA,      "g/L",       "4", 5,
  # Satisfies < 0.25*LLN for grade 4 - PCHG satisfies grade 3
  "Fibrinogen decreased",  0.5,    2.1,     -51,     "g/L",       "4", 6,
  # Satisfies <=75% decrease for grade 4 - LLN  satisfies grade 3
  "Fibrinogen decreased",  0.5,    1.1,     -75,     "g/L",       "4", 7,
  # Satisfies < 0.5*LLN for grade 3 - PCHG missing
  "Fibrinogen decreased",  1,      2.1,     NA,      "g/L",       "3", 8,
  # Satisfies < 0.5*LLN for grade 3 - PCHG satisfies grade 2
  "Fibrinogen decreased",  1,      2.1,     -49,     "g/L",       "3", 9,
  # Satisfies <=50% decrease for grade 3 - LLN  satisfies grade 2
  "Fibrinogen decreased",  1,      2,       -50,     "g/L",       "3", 10,
  # Satisfies < 0.75*LLN for grade 2 - PCHG missing
  "Fibrinogen decreased",  1.5,    2.1,     NA,      "g/L",       "2", 11,
  # Satisfies < 0.75*LLN for grade 2 - PCHG satisfies grade 1
  "Fibrinogen decreased",  1.5,    2.1,     -10,     "g/L",       "2", 12,
  # Satisfies <=25% for grade 2 - LLN satisfies grade 1
  "Fibrinogen decreased",  1.5,    1.6,     -25,     "g/L",       "2", 13,
  # Satisfies < LLN for grade 1 - PCHG missing
  "Fibrinogen decreased",  2,      2.1,     NA,      "g/L",       "1", 14,
  # Satisfies < LLN for grade 1 - PCHG satisfies grade 0
  "Fibrinogen decreased",  2,      2.1,     10,      "g/L",       "1", 15,
  # Satisfies grade 0 - AVAL >= LLN AND no % descrease
  "Fibrinogen decreased",  1.5,    1.5,     0,       "g/L",       "0", 16,
  # Satisfies % decrease for grade 1 - AVAL = LLN so not abnormal
  "Fibrinogen decreased",  1.5,    1.5,     -1,      "g/L",       "0", 17,
  # AVAL >= LLN - PCT missing but its normal so ignore PCT
  "Fibrinogen decreased",  1.5,    1.5,     NA,      "g/L",       "0", 18,
  # Satisfies <=75% decrease for grade 4 - LLN missing do not know its abnormal
  "Fibrinogen decreased",  1,      NA,      -75,     "g/L",        NA, 19,
  # Satisfies <=50% decrease for grade 3 - LLN missing do not know its abnormal
  "Fibrinogen decreased",  1,      NA,      -50,     "g/L",        NA, 20,
  # Satisfies <=25% decrease for grade 2 - LLN missing do not know its abnormal
  "Fibrinogen decreased",  1.5,    NA,      -25,     "g/L",        NA, 21,
  # Satisfies % decrease for grade 1 - LLN missing do not know its abnormal
  "Fibrinogen decreased",  1.5,    NA,      -1,      "g/L",        NA, 22,
  # PCT >= 0 BUT LLN missing cannot grade as NORMAL
  "Fibrinogen decreased",  1.5,    NA,      10,      "g/L",        NA, 23,
  # AVAL missing cannot grade
  "Fibrinogen decreased",  NA,     1.5,     10,      "g/L",        NA, 24,
  # wrong unit cannot grade as it may satisfy grade 4
  "Fibrinogen decreased",  1.5,    1.5,     0,      "g/dL",        NA, 25,
  # missing unit cannot grade as it may satisfy grade 4
  "Fibrinogen decreased",  1.5,    1.5,     0,          NA,        NA, 26,
)

exp_ggt_ctcv4 <- tibble::tribble(
  ~ATOXDSCH,       ~AVAL, ~ANRLO, ~ANRHI, ~ATOXGRH, ~TESTNUM,
  "Not a term",    80,    0,      40,     NA,       1,
  NA_character_,   60,    0,      40,     NA,       2,
  "GGT increased", 801,   0,      40,     "4",      3,
  "GGT increased", 800,   0,      40,     "3",      4,
  "GGT increased", 201,   0,      40,     "3",      5,
  "GGT increased", 200,   0,      40,     "2",      6,
  "GGT increased", 101,   0,      40,     "2",      7,
  "GGT increased", 100,   0,      40,     "1",      8,
  "GGT increased", 41,    0,      40,     "1",      9,
  "GGT increased", 40,    0,      40,     "0",      10,
  # ANRHI missing - cannot grade
  "GGT increased", 100,   0,      NA,     NA,       11,
  # AVAL missing cannot grade
  "GGT increased", NA,    0,      NA,     NA,       12,
) %>%
  mutate(AVALU = NA_character_)

exp_hapt <- tibble::tribble(
  ~ATOXDSCL,               ~AVAL,  ~ANRLO, ~ANRHI, ~ATOXGRL, ~TESTNUM,
  "Not a term",            9,      10,     40,     NA,       1,
  NA_character_,           10,     10,     40,     NA,       2,
  "Haptoglobin decreased", 9,      10,     40,     "1",      3,
  "Haptoglobin decreased", 10,     10,     40,     "0",      4,
  "Haptoglobin decreased", 11,     10,     40,     "0",      5,
  # ANRHI missing - cannot grade
  "Haptoglobin decreased", NA,     10,     40,     NA,       6,
  # AVAL missing cannot grade
  "Haptoglobin decreased", 10,     NA,     NA,     NA,       7,
) %>%
  mutate(AVALU = NA_character_)

exp_hgbi_si <- tibble::tribble(
  ~ATOXDSCH,              ~AVAL, ~BASE, ~ANRHI, ~AVALU, ~ATOXGRH, ~TESTNUM, ~V5,
  "Not a term",           80,    120,   200,     "g/L",       NA,        1, "Y",
  NA_character_,          60,    50,    100,     "g/L",       NA,        2, "Y",
  # BASE greater than ANRHI
  "Hemoglobin increased", 106,   65,    60,      "g/L",      "3",        3, "N",
  "Hemoglobin increased", 105,   65,    60,      "g/L",      "2",        4, "N",
  "Hemoglobin increased", 86,    65,    60,      "g/L",      "2",        5, "N",
  "Hemoglobin increased", 85,    65,    60,      "g/L",      "1",        6, "N",
  "Hemoglobin increased", 66,    65,    60,      "g/L",      "1",        7, "N",
  "Hemoglobin increased", 65,    65,    60,      "g/L",      "0",        8, "N",
  "Hemoglobin increased", NA,    65,    60,      "g/L",       NA,        9, "N",
  # BASE less than or equal to ANRHI
  "Hemoglobin increased", 106,   60,    65,      "g/L",      "3",       10, "Y",
  "Hemoglobin increased", 105,   60,    65,      "g/L",      "2",       11, "Y",
  "Hemoglobin increased", 86,    60,    65,      "g/L",      "2",       12, "Y",
  "Hemoglobin increased", 85,    60,    65,      "g/L",      "1",       13, "Y",
  "Hemoglobin increased", 66,    60,    65,      "g/L",      "1",       14, "Y",
  "Hemoglobin increased", 65,    60,    65,      "g/L",      "0",       15, "Y",
  "Hemoglobin increased", NA,    60,    65,      "g/L",       NA,       16, "Y",
  # BASE missing
  "Hemoglobin increased", 106,   NA,    65,      "g/L",      "3",       17, "N",
  "Hemoglobin increased", 105,   NA,    65,      "g/L",      "2",       18, "N",
  "Hemoglobin increased", 86,    NA,    65,      "g/L",      "2",       19, "N",
  "Hemoglobin increased", 85,    NA,    65,      "g/L",      "1",       20, "N",
  "Hemoglobin increased", 66,    NA,    65,      "g/L",      "1",       21, "N",
  "Hemoglobin increased", 65,    NA,    65,      "g/L",      "0",       22, "N",
  "Hemoglobin increased", NA,    NA,    65,      "g/L",       NA,       23, "N",
  # Unit missing cannot grade
  "Hemoglobin increased", 200,   61,    65,         NA,       NA,       24, "Y",
  # ANRHI missing - cannot grade
  "Hemoglobin increased", 200,   60,    NA,      "g/L",       NA,       25, "Y",
  # AVAL missing cannot grade
  "Hemoglobin increased", NA,    60,    65,      "g/L",       NA,       26, "Y",
)

exp_hgbi_cv <- exp_hgbi_si %>%
  mutate(
    AVAL = AVAL / 10,
    BASE = BASE / 10,
    ANRHI = ANRHI / 10,
    AVALU = if_else(str_to_upper(AVALU) == "G/L", "g/dL", NA_character_)
  )

exp_lip <- tibble::tribble(
  ~ATOXDSCH,          ~AVAL,  ~ANRLO, ~ANRHI, ~ATOXGRH, ~TESTNUM,
  "Not a term",       80,     120,    200,    NA,       1,
  NA_character_,      60,     50,     100,    NA,       2,
  "Lipase IncreaSed", 501,    0,      100,    "4",      3,
  "Lipase Increased", 500,    0,      100,    "3",      4,
  "Lipase Increased", 201,    0,      100,    "3",      5,
  "Lipase Increased", 200,    0,      100,    "2",      6,
  "Lipase Increased", 151,    0,      100,    "2",      7,
  "Lipase Increased", 150,    0,      100,    "1",      8,
  "Lipase Increased", 101,    0,      100,    "1",      9,
  "Lipase Increased", 100,    0,      100,    "0",      10,
  # ANRHI missing cannot grade
  "Lipase Increased", 200,    0,      NA,     NA,       11,
  # AVAL missing cannot grade
  "Lipase Increased", NA,     0,      100,    NA,       12,
) %>%
  mutate(AVALU = NA_character_)

exp_lipv6 <- tibble::tribble(
  ~ATOXDSCH,          ~AVAL,  ~ANRLO, ~ANRHI, ~ATOXGRH, ~TESTNUM,
  "Not a term",       80,     120,    200,    NA,       1,
  NA_character_,      60,     50,     100,    NA,       2,
  "Lipase IncreaSed", 501,    0,      100,    "4",      3,
  "Lipase Increased", 500,    0,      100,    "3",      4,
  "Lipase Increased", 301,    0,      100,    "3",      5,
  "Lipase Increased", 300,    0,      100,    "2",      6,
  "Lipase Increased", 151,    0,      100,    "2",      7,
  "Lipase Increased", 150,    0,      100,    "1",      8,
  "Lipase Increased", 101,    0,      100,    "1",      9,
  "Lipase Increased", 100,    0,      100,    "0",      10,
  # ANRHI missing cannot grade
  "Lipase Increased", 200,    0,      NA,     NA,       11,
  # AVAL missing cannot grade
  "Lipase Increased", NA,     0,      100,    NA,       12,
) %>%
  mutate(AVALU = NA_character_)

exp_lymd_si <- tibble::tribble(
  ~ATOXDSCL,                    ~AVAL,  ~ANRLO, ~ANRHI, ~AVALU,    ~ATOXGRL, ~TESTNUM,
  "Not a term",                 80,     120,    200,    "10^9/L",  NA,       1,
  NA_character_,                60,     50,     100,    "10^9/L",  NA,       2,
  # ANRLO not missing
  "Lymphocyte count decreased", 0.19,   140,    NA,     "10^9/L",  "4",      3,
  "Lymphocyte count decreased", 0.2,    140,    NA,     "10^9/L",  "3",      4,
  "Lymphocyte count decreased", 0.49,   140,    NA,     "10^9/L",  "3",      5,
  "Lymphocyte count decreased", 0.5,    140,    NA,     "10^9/L",  "2",      6,
  "Lymphocyte count decreased", 0.79,   140,    NA,     "10^9/L",  "2",      7,
  "Lymphocyte count decreased", 0.8,    140,    NA,     "10^9/L",  "1",      8,
  "Lymphocyte count decreased", 1.09,   1.1,    NA,     "10^9/L",  "1",      9,
  "Lymphocyte count decreased", 1.1,    1.1,    NA,     "10^9/L",  "0",      10,
  # ANRLO missing - can grade 2-4
  "Lymphocyte count decreased", 0.19,   NA,     NA,     "10^9/L",  "4",      11,
  "Lymphocyte count decreased", 0.2,    NA,     NA,     "10^9/L",  "3",      12,
  "Lymphocyte count decreased", 0.49,   NA,     NA,     "10^9/L",  "3",      13,
  "Lymphocyte count decreased", 0.5,    NA,     NA,     "10^9/L",  "2",      14,
  "Lymphocyte count decreased", 0.79,   NA,     NA,     "10^9/L",  "2",      15,
  # ANRLO missing - can NOT grade 0 or 1
  "Lymphocyte count decreased", 0.8,    NA,     NA,     "10^9/L",  NA,       16,
  "Lymphocyte count decreased", 1.09,   NA,     NA,     "10^9/L",  NA,       17,
  "Lymphocyte count decreased", 1.1,    NA,     NA,     "10^9/L",  NA,       18,
  # Unit missing cannot grade
  "Lymphocyte count decreased", 1.1,    1.1,    NA,     NA,        NA,       19,
  # AVAL missing cannot grade
  "Lymphocyte count decreased", 1.1,    1.1,    NA,     "10^9/L",  "0",      20,
)

exp_lymd_cv <- exp_lymd_si %>%
  mutate(
    AVALU = if_else(str_to_upper(AVALU) == "10^9/L", "10^3/uL", AVALU)
  )

exp_lymd_cv2 <- exp_lymd_si %>%
  mutate(
    AVAL = AVAL * 1000,
    ANRLO = ANRLO * 1000,
    ANRHI = ANRHI * 1000,
    AVALU = if_else(str_to_upper(AVALU) == "10^9/L", "10^3/mL", AVALU)
  )

exp_lymi_si <- tibble::tribble(
  ~ATOXDSCH,                    ~AVAL,  ~ANRLO, ~ANRHI, ~AVALU,    ~ATOXGRH, ~TESTNUM,
  "Not a term",                 80,     120,    200,    "10^9/L",  NA,       1,
  NA_character_,                60,     50,     100,    "10^9/L",  NA,       2,
  "Lymphocyte count increased", 21,     NA,     NA,     "10^9/L",  "3",      3,
  "Lymphocyte count increased", 20,     NA,     NA,     "10^9/L",  "2",      4,
  "Lymphocyte count increased", 4.1,    NA,     NA,     "10^9/L",  "2",      5,
  "Lymphocyte count increased", 4,      NA,     NA,     "10^9/L",  "0",      6,
  # Unit missing cannot grade
  "Lymphocyte count increased", 4,      NA,     NA,     NA,        NA,       7,
  # AVAL missing cannot grade
  "Lymphocyte count increased", NA,     NA,     NA,     "10^9/L",  NA,       8,
)

exp_lymi_cv1 <- exp_lymi_si %>%
  mutate(
    AVALU = if_else(str_to_upper(AVALU) == "10^9/L", "10^3/uL", AVALU)
  )

exp_lymi_cv2 <- exp_lymi_si %>%
  mutate(
    AVAL = AVAL * 1000,
    ANRLO = ANRLO * 1000,
    ANRHI = ANRHI * 1000,
    AVALU = if_else(str_to_upper(AVALU) == "10^9/L", "10^3/mL", AVALU),
    TESTNUM = TESTNUM + 8
  )

exp_lymi_cv <- bind_rows(exp_lymi_cv1, exp_lymi_cv2)

exp_neut_si <- tibble::tribble(
  ~ATOXDSCL,                    ~AVAL,  ~ANRLO, ~ANRHI, ~AVALU,    ~ATOXGRL, ~TESTNUM,
  "Not a term",                 80,     120,    200,    "10^9/L",  NA,       1,
  NA_character_,                60,     50,     100,    "10^9/L",  NA,       2,
  # ANRLO not missing
  "Neutrophil count decreased", 0.49,   2,      NA,     "10^9/L",  "4",      3,
  "Neutrophil count decreased", 0.5,    2,      NA,     "10^9/L",  "3",      4,
  "Neutrophil count decreased", 0.99,   2,      NA,     "10^9/L",  "3",      5,
  "Neutrophil count decreased", 1,      2,      NA,     "10^9/L",  "2",      6,
  "Neutrophil count decreased", 1.49,   2,      NA,     "10^9/L",  "2",      7,
  "Neutrophil count decreased", 1.5,    2,      NA,     "10^9/L",  "1",      8,
  "Neutrophil count decreased", 1.9,    2,      NA,     "10^9/L",  "1",      9,
  "Neutrophil count decreased", 2,      2,      NA,     "10^9/L",  "0",      10,
  # ANRLO missing - can grade 2-4
  "Neutrophil count decreased", 0.49,   NA,     NA,     "10^9/L",  "4",      11,
  "Neutrophil count decreased", 0.5,    NA,     NA,     "10^9/L",  "3",      12,
  "Neutrophil count decreased", 0.99,   NA,     NA,     "10^9/L",  "3",      13,
  "Neutrophil count decreased", 1,      NA,     NA,     "10^9/L",  "2",      14,
  "Neutrophil count decreased", 1.49,   NA,     NA,     "10^9/L",  "2",      15,
  # ANRLO missing - can NOT grade 0 or 1
  "Neutrophil count decreased", 1.5,    NA,     NA,     "10^9/L",  NA,       16,
  "Neutrophil count decreased", 1.9,    NA,     NA,     "10^9/L",  NA,       17,
  "Neutrophil count decreased", 2,      NA,     NA,     "10^9/L",  NA,       18,
  # Unit missing cannot grade
  "Neutrophil count decreased", 2,      2,      NA,     NA,        NA,       19,
  # AVAL missing cannot grade
  "Neutrophil count decreased", NA,     2,      NA,     "10^9/L",  NA,       20,
)

exp_neut_cv1 <- exp_neut_si %>%
  mutate(
    AVALU = if_else(str_to_upper(AVALU) == "10^9/L", "10^3/uL", AVALU)
  )

exp_neut_cv2 <- exp_neut_si %>%
  mutate(
    AVAL = AVAL * 1000,
    ANRLO = ANRLO * 1000,
    ANRHI = ANRHI * 1000,
    AVALU = if_else(str_to_upper(AVALU) == "10^9/L", "10^3/mL", AVALU),
    TESTNUM = TESTNUM + 20
  )

exp_neut_cv <- bind_rows(exp_neut_cv1, exp_neut_cv2)

exp_neut_siv6 <- tibble::tribble(
  ~ATOXDSCL,                    ~AVAL,  ~AVALU,    ~ATOXGRL, ~TESTNUM,
  "Not a term",                 80,     "10^9/L",  NA,       1,
  NA_character_,                60,     "10^9/L",  NA,       2,
  # ANRLO not missing
  "Neutrophil count decreased", 0.09,   "10^9/L",  "4",      3,
  "Neutrophil count decreased", 0.1,    "10^9/L",  "3",      4,
  "Neutrophil count decreased", 0.49,   "10^9/L",  "3",      5,
  "Neutrophil count decreased", 0.5,    "10^9/L",  "2",      6,
  "Neutrophil count decreased", 0.99,   "10^9/L",  "2",      7,
  "Neutrophil count decreased", 1,      "10^9/L",  "1",      8,
  "Neutrophil count decreased", 1.49,   "10^9/L",  "1",      9,
  "Neutrophil count decreased", 1.5,    "10^9/L",  "0",      10,
  # Unit missing cannot grade
  "Neutrophil count decreased", NA,     NA,        NA,       19,
  # AVAL missing cannot grade
  "Neutrophil count decreased", NA,     "10^9/L",  NA,       20,
)

exp_plate_si <- tibble::tribble(
  ~ATOXDSCL,                  ~AVAL, ~ANRLO, ~ANRHI, ~AVALU,   ~ATOXGRL, ~TESTNUM,
  "Not a term",               80,    120,    200,    "10^9/L", NA,       1,
  NA_character_,              60,    50,     100,    "10^9/L", NA,       2,
  # ANRLO not missing
  "Platelet count decreased", 24,    100,    150,    "10^9/L", "4",      3,
  "Platelet count decreased", 25,    100,    150,    "10^9/L", "3",      4,
  "Platelet count decreased", 49,    100,    150,    "10^9/L", "3",      5,
  "Platelet count decreased", 50,    100,    150,    "10^9/L", "2",      6,
  "Platelet count decreased", 74,    100,    150,    "10^9/L", "2",      7,
  "Platelet count decreased", 75,    100,    150,    "10^9/L", "1",      8,
  "Platelet count decreased", 99,    100,    NA,     "10^9/L", "1",      9,
  "Platelet count decreased", 100,   100,    NA,     "10^9/L", "0",      10,
  # ANRLO missing - can grade 2-4      1,
  "Platelet count decreased", 24,    NA,     NA,     "10^9/L", "4",      11,
  "Platelet count decreased", 25,    NA,     NA,     "10^9/L", "3",      12,
  "Platelet count decreased", 49,    NA,     NA,     "10^9/L", "3",      13,
  "Platelet count decreased", 50,    NA,     NA,     "10^9/L", "2",      14,
  "Platelet count decreased", 74,    NA,     NA,     "10^9/L", "2",      15,
  # ANRLO missing - can NOT grade 0 or 1
  "Platelet count decreased", 75,    NA,     NA,     "10^9/L", NA,       16,
  "Platelet count decreased", 99,    NA,     NA,     "10^9/L", NA,       17,
  "Platelet count decreased", 100,   NA,     NA,     "10^9/L", NA,       18,
  # Unit missing cannot grade
  "Platelet count decreased", 100,   100,    NA,     NA,       NA,       19,
  # AVAL missing cannot grade
  "Platelet count decreased", NA,    100,    NA,     "10^9/L", NA,       20,
)

exp_plate_cv <- exp_plate_si %>%
  mutate(
    AVALU = if_else(str_to_upper(AVALU) == "10^9/L", "10^3/uL", AVALU)
  )

exp_plate_cv2 <- exp_plate_si %>%
  mutate(
    AVAL = AVAL * 1000,
    ANRLO = ANRLO * 1000,
    ANRHI = ANRHI * 1000,
    AVALU = if_else(str_to_upper(AVALU) == "10^9/L", "10^3/mL", AVALU)
  )

exp_plate_v6_si <- exp_plate_si %>%
  mutate(
    ATOXDSCL = if_else(
      ATOXDSCL == "Platelet count decreased", "Thrombocytopenia", ATOXDSCL
    ),
    AVAL = case_when(
      AVAL == 24 ~ 9,
      AVAL == 25 ~ 10,
      TRUE ~ AVAL
    )
  )

exp_seri <- tibble::tribble(
  ~ATOXDSCH,                 ~AVAL,  ~ANRLO, ~ANRHI, ~ATOXGRH, ~TESTNUM,
  "Not a term",              80,     120,    200,    NA,       1,
  NA_character_,             60,     50,     100,    NA,       2,
  "Serum amylase increased", 501,    0,      100,    "4",      3,
  "Serum amylase increased", 500,    0,      100,    "3",      4,
  "Serum amylase increased", 201,    0,      100,    "3",      5,
  "Serum amylase increased", 200,    0,      100,    "2",      6,
  "Serum amylase increased", 151,    0,      100,    "2",      7,
  "Serum amylase increased", 150,    0,      100,    "1",      8,
  "Serum amylase increased", 101,    0,      100,    "1",      9,
  "Serum amylase increased", 100,    0,      100,    "0",      10,
  # ANRHI missing cannot grade
  "Serum amylase increased", 200,    0,      NA,     NA,       11,
  # AVAL missing cannot grade
  "Serum amylase increased", NA,     0,      100,    NA,       12,
) %>%
  mutate(AVALU = NA_character_)

exp_wbcd <- tibble::tribble(
  ~ATOXDSCL,                    ~AVAL, ~ANRLO, ~ANRHI, ~AVALU,   ~ATOXGRL, ~TESTNUM,
  "Not a term",                 1,     5,      15,     "10^9/L", NA,       1,
  NA_character_,                2,     5,      15,     "10^9/L", NA,       2,
  # ANRLO not missing
  "White blood cell decreased", 0.9,   5,      15,     "10^9/L", "4",      3,
  "White blood cell decreased", 1,     5,      15,     "10^9/L", "3",      4,
  "White blood cell decreased", 1.9,   5,      15,     "10^9/L", "3",      5,
  "White blood cell decreased", 2,     5,      15,     "10^9/L", "2",      6,
  "White blood cell decreased", 2.9,   5,      15,     "10^9/L", "2",      7,
  "White blood cell decreased", 3,     5,      15,     "10^9/L", "1",      8,
  "White blood cell decreased", 4.9,   5,      15,     "10^9/L", "1",      9,
  "White blood cell decreased", 5,     5,      15,     "10^9/L", "0",      10,
  # ANRLO missing - can grade 2-4
  "White blood cell decreased", 0.9,   NA,     NA,     "10^9/L", "4",      11,
  "White blood cell decreased", 1,     NA,     NA,     "10^9/L", "3",      12,
  "White blood cell decreased", 1.9,   NA,     NA,     "10^9/L", "3",      13,
  "White blood cell decreased", 2,     NA,     NA,     "10^9/L", "2",      14,
  "White blood cell decreased", 2.9,   NA,     NA,     "10^9/L", "2",      15,
  # ANRLO missing - can NOT grade 0 or 1
  "White blood cell decreased", 3,     NA,     NA,     "10^9/L", NA,       16,
  "White blood cell decreased", 3,     NA,     NA,     "10^9/L", NA,       17,
  "White blood cell decreased", 3,     NA,     NA,     "10^9/L", NA,       18,
  # Unit missing cannot grade
  "White blood cell decreased", 3,     100,    NA,     NA,       NA,       19,
  # AVAL missing cannot grade
  "White blood cell decreased", NA,    100,    NA,     "10^9/L", NA,       20,
)

exp_wbcd_uscv1 <- exp_wbcd %>%
  mutate(AVALU = if_else(AVALU == "10^9/L", "10^3/uL", AVALU))

exp_wbcd_uscv2 <- exp_wbcd %>%
  mutate(
    ANRLO = if_else(!is.na(ANRLO), ANRLO * 1000, NA),
    AVAL = if_else(!is.na(AVAL), AVAL * 1000, NA),
    AVALU = if_else(AVALU == "10^9/L", "10^3/mL", AVALU),
    TESTNUM = TESTNUM + 20
  )

exp_wbcd_uscv <- bind_rows(exp_wbcd_uscv1, exp_wbcd_uscv2)

exp_calci_si <- tibble::tribble(
  ~ATOXDSCH,       ~AVAL,  ~ANRLO, ~ANRHI, ~AVALU,    ~ATOXGRH, ~TESTNUM,
  "Not a term",    3.5,    0,      2.5,    "mmol/L",  NA,       1,
  NA_character_,   3.5,    0,      2.5,    "mmol/L",  NA,       2,
  # ANRHI not missing
  "Hypercalcemia", 3.5,    0,      2.5,    "mmol/L",  "4",      3,
  "Hypercalcemia", 3.4,    0,      2.5,    "mmol/L",  "3",      4,
  "Hypercalcemia", 3.2,    0,      2.5,    "mmol/L",  "3",      5,
  "Hypercalcemia", 3.1,    0,      2.5,    "mmol/L",  "2",      6,
  "Hypercalcemia", 3,      0,      2.5,    "mmol/L",  "2",      7,
  "Hypercalcemia", 2.9,    0,      2.5,    "mmol/L",  "1",      8,
  "Hypercalcemia", 2.6,    0,      2.5,    "mmol/L",  "1",      9,
  "Hypercalcemia", 2.5,    0,      2.5,    "mmol/L",  "0",      10,
  # ANRHI missing - can grade 2-4
  "Hypercalcemia", 3.5,    0,      NA,     "mmol/L",  "4",      11,
  "Hypercalcemia", 3.4,    0,      NA,     "mmol/L",  "3",      12,
  "Hypercalcemia", 3.2,    0,      NA,     "mmol/L",  "3",      13,
  "Hypercalcemia", 3.1,    0,      NA,     "mmol/L",  "2",      14,
  "Hypercalcemia", 3,      0,      NA,     "mmol/L",  "2",      15,
  # ANRHI missing - can NOT grade 0 or 1
  "Hypercalcemia", 2.9,    0,      NA,     "mmol/L",  NA,       16,
  "Hypercalcemia", 2.6,    0,      NA,     "mmol/L",  NA,       17,
  "Hypercalcemia", 2.5,    0,      NA,     "mmol/L",  NA,       18,
  # Unit missing cannot grade
  "Hypercalcemia", 2.5,    0,      2.5,    NA,        NA,       19,
  # AVAL missing cannot grade
  "Hypercalcemia", NA,     0,      2.5,    "mmol/L",  NA,       20,
)

exp_calci_cv <- tibble::tribble(
  ~ATOXDSCH,       ~AVAL,  ~ANRLO, ~ANRHI,   ~AVALU,   ~ATOXGRH, ~TESTNUM,
  "Not a term",    13.6,   0,      10.5,    "mg/dL",         NA,        1,
  NA_character_,   13.6,   0,      10.5,    "mg/dL",         NA,        2,
  # ANRHI not missing
  "Hypercalcemia", 13.6,   0,      10.5,    "mg/dL",        "4",        3,
  "Hypercalcemia", 13.5,   0,      10.5,    "mg/dL",        "3",        4,
  "Hypercalcemia", 12.6,   0,      10.5,    "mg/dL",        "3",        5,
  "Hypercalcemia", 12.5,   0,      10.5,    "mg/dL",        "2",        6,
  "Hypercalcemia", 11.6,   0,      10.5,    "mg/dL",        "2",        7,
  "Hypercalcemia", 11.5,   0,      10.5,    "mg/dL",        "1",        8,
  "Hypercalcemia", 10.6,   0,      10.5,    "mg/dL",        "1",        9,
  "Hypercalcemia", 10.5,   0,      10.5,    "mg/dL",        "0",       10,
  # ANRHI missing - can grade 2-4
  "Hypercalcemia", 13.6,   0,      NA,      "mg/dL",        "4",       11,
  "Hypercalcemia", 13.5,   0,      NA,      "mg/dL",        "3",       12,
  "Hypercalcemia", 12.6,   0,      NA,      "mg/dL",        "3",       13,
  "Hypercalcemia", 12.5,   0,      NA,      "mg/dL",        "2",       14,
  "Hypercalcemia", 11.6,   0,      NA,      "mg/dL",        "2",       15,
  # ANRHI missing - can NOT grade 0 or 1 mg/dL
  "Hypercalcemia", 11.5,   0,      NA,      "mg/dL",         NA,       16,
  "Hypercalcemia", 10.6,   0,      NA,      "mg/dL",         NA,       17,
  "Hypercalcemia", 10.5,   0,      NA,      "mg/dL",         NA,       18,
  # Unit missing cannot grade
  "Hypercalcemia", 10.5,   0,      10.5,         NA,         NA,       19,
  # AVAL missing cannot grade
  "Hypercalcemia", NA,     0,      10.5,    "mg/dL",         NA,       20,
)

exp_calioni_si <- tibble::tribble(
  ~ATOXDSCH,                 ~AVAL,  ~ANRLO, ~ANRHI, ~AVALU,    ~ATOXGRH, ~TESTNUM,
  "Not a term",              1.9,    0,      1.3,    "mmol/L",  NA,       1,
  NA_character_,             1.9,    0,      1.3,    "mmol/L",  NA,       2,
  # ANRHI not missing
  "Hypercalcemia (Ionized)", 1.9,    0,      1.3,    "mmol/L",  "4",      3,
  "Hypercalcemia (Ionized)", 1.8,    0,      1.3,    "mmol/L",  "3",      4,
  "Hypercalcemia (Ionized)", 1.7,    0,      1.3,    "mmol/L",  "3",      5,
  "Hypercalcemia (Ionized)", 1.6,    0,      1.3,    "mmol/L",  "2",      6,
  "Hypercalcemia (Ionized)", 1.51,   0,      1.3,    "mmol/L",  "2",      7,
  "Hypercalcemia (Ionized)", 1.5,    0,      1.3,    "mmol/L",  "1",      8,
  "Hypercalcemia (Ionized)", 1.4,    0,      1.3,    "mmol/L",  "1",      9,
  "Hypercalcemia (Ionized)", 1.3,    0,      1.3,    "mmol/L",  "0",      10,
  # ANRHI missing - can grade 2-4
  "Hypercalcemia (Ionized)", 1.9,    0,      NA,     "mmol/L",  "4",      11,
  "Hypercalcemia (Ionized)", 1.8,    0,      NA,     "mmol/L",  "3",      12,
  "Hypercalcemia (Ionized)", 1.7,    0,      NA,     "mmol/L",  "3",      13,
  "Hypercalcemia (Ionized)", 1.6,    0,      NA,     "mmol/L",  "2",      14,
  "Hypercalcemia (Ionized)", 1.51,   0,      NA,     "mmol/L",  "2",      15,
  # ANRHI missing - can NOT grade 0 or 1
  "Hypercalcemia (Ionized)", 1.5,    0,      NA,     "mmol/L",  NA,       16,
  "Hypercalcemia (Ionized)", 1.4,    0,      NA,     "mmol/L",  NA,       17,
  "Hypercalcemia (Ionized)", 1.3,    0,      NA,     "mmol/L",  NA,       18,
  # Unit missing cannot grade
  "Hypercalcemia (Ionized)", 1.3,    0,      1.3,    NA,        NA,       19,
  # AVAL missing cannot grade
  "Hypercalcemia (Ionized)", NA,     0,      1.3,    "mmol/L",  NA,       20,
)

exp_calioni_cv <- tibble::tribble(
  ~ATOXDSCH,                  ~AVAL,  ~ANRLO,  ~ANRHI,   ~AVALU,  ~ATOXGRH,  ~TESTNUM,
  "Not a term",               7.3,    0,       5,       "mg/dL",        NA,  1,
  NA_character_,              7.3,    0,       5,       "mg/dL",        NA,  2,
  # ANRHI not missing
  "Hypercalcemia (Ionized)",  7.3,    0,       5,       "mg/dL",       "4",  3,
  "Hypercalcemia (Ionized)",  7.2,    0,       5,       "mg/dL",       "3",  4,
  "Hypercalcemia (Ionized)",  6.5,    0,       5,       "mg/dL",       "3",  5,
  "Hypercalcemia (Ionized)",  6.4,    0,       5,       "mg/dL",       "2",  6,
  "Hypercalcemia (Ionized)",  6.1,    0,       5,       "mg/dL",       "2",  7,
  "Hypercalcemia (Ionized)",  6.0,    0,       5,       "mg/dL",       "1",  8,
  "Hypercalcemia (Ionized)",  5.1,    0,       5,       "mg/dL",       "1",  9,
  "Hypercalcemia (Ionized)",  5,      0,       5,       "mg/dL",       "0",  10,
  # ANRHI missing - can grade 2-4
  "Hypercalcemia (Ionized)",  7.3,    0,       NA,      "mg/dL",       "4",  11,
  "Hypercalcemia (Ionized)",  7.2,    0,       NA,      "mg/dL",       "3",  12,
  "Hypercalcemia (Ionized)",  6.5,    0,       NA,      "mg/dL",       "3",  13,
  "Hypercalcemia (Ionized)",  6.4,    0,       NA,      "mg/dL",       "2",  14,
  "Hypercalcemia (Ionized)",  6.1,    0,       NA,      "mg/dL",       "2",  15,
  # ANRHI missing - can NOT grade 0 or 1mg/d
  "Hypercalcemia (Ionized)",  6.0,    0,       NA,      "mg/dL",        NA,  16,
  "Hypercalcemia (Ionized)",  5.1,    0,       NA,      "mg/dL",        NA,  17,
  "Hypercalcemia (Ionized)",  5,      0,       NA,      "mg/dL",        NA,  18,
  # Unit missing cannot grade
  "Hypercalcemia (Ionized)",  5,      0,       5,            NA,        NA,  19,
  # AVAL missing cannot grade
  "Hypercalcemia (Ionized)",  NA,     0,       5,       "mg/dL",        NA,  20,
)

exp_glycfi <- tibble::tribble(
  ~ATOXDSCH,                 ~AVAL,  ~ANRLO, ~ANRHI, ~AVALU,    ~ATOXGRH, ~TESTNUM,
  "Not a term",              27.9,   0,      5.3,    "mmol/L",  NA,       1,
  NA_character_,             27.9,   0,      5.3,    "mmol/L",  NA,       2,
  # ANRHI not missing
  "Hyperglycemia (Fasting)", 27.9,   0,      5.3,    "mmol/L",  "4",      3,
  "Hyperglycemia (Fasting)", 27.8,   0,      5.3,    "mmol/L",  "3",      4,
  "Hyperglycemia (Fasting)", 14,     0,      5.3,    "mmol/L",  "3",      5,
  "Hyperglycemia (Fasting)", 13.9,   0,      5.3,    "mmol/L",  "2",      6,
  "Hyperglycemia (Fasting)", 9,      0,      5.3,    "mmol/L",  "2",      7,
  "Hyperglycemia (Fasting)", 8.9,    0,      5.3,    "mmol/L",  "1",      8,
  "Hyperglycemia (Fasting)", 5.4,    0,      5.3,    "mmol/L",  "1",      9,
  "Hyperglycemia (Fasting)", 5.3,    0,      5.3,    "mmol/L",  "0",      10,
  # ANRHI missing - can grade 2-4
  "Hyperglycemia (Fasting)", 27.9,   0,      NA,     "mmol/L",  "4",      11,
  "Hyperglycemia (Fasting)", 27.8,   0,      NA,     "mmol/L",  "3",      12,
  "Hyperglycemia (Fasting)", 14,     0,      NA,     "mmol/L",  "3",      13,
  "Hyperglycemia (Fasting)", 13.9,   0,      NA,     "mmol/L",  "2",      14,
  "Hyperglycemia (Fasting)", 9,      0,      NA,     "mmol/L",  "2",      15,
  # ANRHI missing - can NOT grade 0 or 1
  "Hyperglycemia (Fasting)", 8.9,    0,      NA,     "mmol/L",  NA,       16,
  "Hyperglycemia (Fasting)", 5.4,    0,      NA,     "mmol/L",  NA,       17,
  "Hyperglycemia (Fasting)", 5.3,    0,      NA,     "mmol/L",  NA,       18,
  # Unit missing cannot grade
  "Hyperglycemia (Fasting)", 5.3,    0,      5.3,    NA,        NA,       19,
  # AVAL missing cannot grade
  "Hyperglycemia (Fasting)", NA,     0,      5.3,    "mmol/L",  NA,       20,
)

exp_glyci <- tibble::tribble(
  ~ATOXDSCH,       ~AVAL,  ~ANRLO, ~ANRHI, ~AVALU,    ~ATOXGRH, ~TESTNUM,
  "Not a term",    27.9,   0,      5.3,    "mmol/L",  NA,       1,
  NA_character_,   27.9,   0,      5.3,    "mmol/L",  NA,       2,
  "Hyperglycemia", 27.9,   0,      5.3,    "mmol/L",  "4",      3,
  "Hyperglycemia", 27.8,   0,      5.3,    "mmol/L",  "3",      4,
  "Hyperglycemia", 14,     0,      5.3,    "mmol/L",  "3",      5,
  "Hyperglycemia", 13.9,   0,      5.3,    "mmol/L",  "0",      6,
  "Hyperglycemia", 5.3,    0,      NA,     "mmol/L",  "0",      7,
  # Unit missing cannot grade
  "Hyperglycemia", 5.3,    0,      5.3,    NA,        NA,       8,
  # AVAL missing cannot grade
  "Hyperglycemia", NA,     0,      5.3,    "mmol/L",  NA,       9,
)

exp_glyci_uscv <- tibble::tribble(
  ~ATOXDSCH,        ~AVAL,  ~ANRLO,  ~ANRHI,   ~AVALU,  ~ATOXGRH,  ~TESTNUM,
  "Not a term",     501,    0,       5.3,     "mg/dL",        NA,  1,
  NA_character_,    501,    0,       5.3,     "mg/dL",        NA,  2,
  "Hyperglycemia",  501,    0,       5.3,     "mg/dL",       "4",  3,
  "Hyperglycemia",  500,    0,       5.3,     "mg/dL",       "3",  4,
  "Hyperglycemia",  251,    0,       5.3,     "mg/dL",       "3",  5,
  "Hyperglycemia",  250,    0,       5.3,     "mg/dL",       "0",  6,
  "Hyperglycemia",  250,    0,       NA,      "mg/dL",       "0",  7,
  # Unit missing cannot grade
  "Hyperglycemia",  250,    0,       5.3,          NA,        NA,  8,
  # AVAL missing cannot grade
  "Hyperglycemia",  NA,     0,       5.3,     "mg/dL",        NA,  9,
)

exp_kalei <- tibble::tribble(
  ~ATOXDSCH,       ~AVAL,  ~ANRLO, ~ANRHI, ~AVALU,    ~ATOXGRH, ~TESTNUM,
  "Not a term",    7.1,    0,      5.1,    "mmol/L",  NA,       1,
  NA_character_,   7.1,    0,      5.1,    "mmol/L",  NA,       2,
  # ANRHI not missing
  "Hyperkalemia",  7.1,    0,      5.1,    "mmol/L",  "4",      3,
  "Hyperkalemia",  7,      0,      5.1,    "mmol/L",  "3",      4,
  "Hyperkalemia",  6.1,    0,      5.1,    "mmol/L",  "3",      5,
  "Hyperkalemia",  6,      0,      5.1,    "mmol/L",  "2",      6,
  "Hyperkalemia",  5.6,    0,      5.1,    "mmol/L",  "2",      7,
  "Hyperkalemia",  5.5,    0,      5.1,    "mmol/L",  "1",      8,
  "Hyperkalemia",  5.2,    0,      5.1,    "mmol/L",  "1",      9,
  "Hyperkalemia",  5.1,    0,      5.1,    "mmol/L",  "0",      10,
  # ANRHI missing - can grade 2-4
  "Hyperkalemia",  7.1,    0,      NA,     "mmol/L",  "4",      11,
  "Hyperkalemia",  7,      0,      NA,     "mmol/L",  "3",      12,
  "Hyperkalemia",  6.1,    0,      NA,     "mmol/L",  "3",      13,
  "Hyperkalemia",  6,      0,      NA,     "mmol/L",  "2",      14,
  "Hyperkalemia",  5.6,    0,      NA,     "mmol/L",  "2",      15,
  # ANRHI missing - can NOT grade 0 or 1
  "Hyperkalemia",  5.5,    0,      NA,     "mmol/L",  NA,       16,
  "Hyperkalemia",  5.2,    0,      NA,     "mmol/L",  NA,       17,
  "Hyperkalemia",  5.1,    0,      NA,     "mmol/L",  NA,       18,
  # Unit missing cannot grade
  "Hyperkalemia",  5.1,    0,      5.1,    NA,        NA,       19,
  # AVAL missing cannot grade
  "Hyperkalemia",  NA,     0,      5.1,    "mmol/L",  NA,       20,
)

exp_magni_si <- tibble::tribble(
  ~ATOXDSCH,         ~AVAL, ~ANRLO, ~ANRHI, ~AVALU,   ~ATOXGRH, ~TESTNUM,
  "Not a term",      3.4,   0,      0.8,    "mmol/L", NA,       1,
  NA_character_,     3.4,   0,      0.8,    "mmol/L", NA,       2,
  # ANRHI not missing
  "Hypermagnesemia", 3.4,   0,      0.8,    "mmol/L", "4",      3,
  "Hypermagnesemia", 3.3,   0,      0.8,    "mmol/L", "3",      4,
  "Hypermagnesemia", 1.24,  0,      0.8,    "mmol/L", "3",      5,
  "Hypermagnesemia", 1.23,  0,      0.8,    "mmol/L", "1",      6,
  "Hypermagnesemia", 0.81,  0,      0.8,    "mmol/L", "1",      7,
  "Hypermagnesemia", 0.8,   0,      0.8,    "mmol/L", "0",      8,
  # ANRHI missing - can grade 3-4
  "Hypermagnesemia", 3.4,   0,      NA,     "mmol/L", "4",      9,
  "Hypermagnesemia", 3.3,   0,      NA,     "mmol/L", "3",      10,
  "Hypermagnesemia", 1.24,  0,      NA,     "mmol/L", "3",      11,
  # ANRHI missing - can NOT grade 0 or 1
  "Hypermagnesemia", 1.23,  0,      NA,     "mmol/L", NA,       12,
  "Hypermagnesemia", 0.81,  0,      NA,     "mmol/L", NA,       13,
  "Hypermagnesemia", 0.8,   0,      NA,     "mmol/L", NA,       14,
  # Unit missing cannot grade
  "Hypermagnesemia", 0.8,   0,      0.8,    NA,       NA,       15,
  # AVAL missing cannot grade
  "Hypermagnesemia", NA,    0,      0.8,    "mmol/L", NA,       16,
)

exp_magni_cv <- tibble::tribble(
  ~ATOXDSCH,          ~AVAL,  ~ANRLO,  ~ANRHI,   ~AVALU, ~ATOXGRH,  ~TESTNUM,
  "Not a term",       8.1,    0,       1.2,     "mg/dL",       NA,  1,
  NA_character_,      8.1,    0,       1.2,     "mg/dL",       NA,  2,
  # ANRHI not missing
  "Hypermagnesemia",  8.1,    0,       1.2,     "mg/dL",      "4",  3,
  "Hypermagnesemia",  8,      0,       1.2,     "mg/dL",      "3",  4,
  "Hypermagnesemia",  3.1,    0,       1.2,     "mg/dL",      "3",  5,
  "Hypermagnesemia",  3,      0,       1.2,     "mg/dL",      "1",  6,
  "Hypermagnesemia",  1.21,   0,       1.2,     "mg/dL",      "1",  7,
  "Hypermagnesemia",  1.2,    0,       1.2,     "mg/dL",      "0",  8,
  # ANRHI missing - can grade 3-4
  "Hypermagnesemia",  8.1,    0,       NA,      "mg/dL",      "4",  9,
  "Hypermagnesemia",  8,      0,       NA,      "mg/dL",      "3",  10,
  "Hypermagnesemia",  3.1,    0,       NA,      "mg/dL",      "3",  11,
  # ANRHI missing - can NOT grade 0 or 1
  "Hypermagnesemia",  3,      0,       NA,      "mg/dL",       NA,  12,
  "Hypermagnesemia",  1.21,   0,       NA,      "mg/dL",       NA,  13,
  "Hypermagnesemia",  1.2,    0,       NA,      "mg/dL",       NA,  14,
  # Unit missing cannot grade
  "Hypermagnesemia",  1.2,    0,       1.2,          NA,       NA,  15,
  # AVAL missing cannot grade
  "Hypermagnesemia",  NA,     0,       1.2,     "mg/dL",       NA,  16,
)

exp_natri <- tibble::tribble(
  ~ATOXDSCH,        ~AVAL,  ~ANRLO, ~ANRHI, ~AVALU,    ~ATOXGRH, ~TESTNUM,
  "Not a term",     161,    0,      140,    "mmol/L",  NA,       1,
  NA_character_,    161,    0,      140,    "mmol/L",  NA,       2,
  # ANRHI not missing
  "Hypernatremia",  161,    0,      140,    "mmol/L",  "4",      3,
  "Hypernatremia",  160,    0,      140,    "mmol/L",  "3",      4,
  "Hypernatremia",  156,    0,      140,    "mmol/L",  "3",      5,
  "Hypernatremia",  155,    0,      140,    "mmol/L",  "2",      6,
  "Hypernatremia",  151,    0,      140,    "mmol/L",  "2",      7,
  "Hypernatremia",  150,    0,      140,    "mmol/L",  "1",      8,
  "Hypernatremia",  141,    0,      140,    "mmol/L",  "1",      9,
  "Hypernatremia",  140,    0,      140,    "mmol/L",  "0",      10,
  # ANRHI missing - can grade 3-4
  "Hypernatremia",  161,    0,      NA,     "mmol/L",  "4",      11,
  "Hypernatremia",  160,    0,      NA,     "mmol/L",  "3",      12,
  "Hypernatremia",  156,    0,      NA,     "mmol/L",  "3",      13,
  "Hypernatremia",  155,    0,      NA,     "mmol/L",  "2",      14,
  "Hypernatremia",  151,    0,      NA,     "mmol/L",  "2",      15,
  # ANRHI missing - can NOT grade 0 or 1
  "Hypernatremia",  150,    0,      NA,     "mmol/L",  NA,       16,
  "Hypernatremia",  141,    0,      NA,     "mmol/L",  NA,       17,
  "Hypernatremia",  140,    0,      NA,     "mmol/L",  NA,       18,
  # Unit missing cannot grade
  "Hypernatremia",  140,    0,      140,    NA,        NA,       19,
  # AVAL missing cannot grade
  "Hypernatremia",  NA,     0,      140,    "mmol/L",  NA,       20,
)

exp_trigi_si <- tibble::tribble(
  ~ATOXDSCH,               ~AVAL,  ~ANRLO, ~ANRHI, ~AVALU,    ~ATOXGRH, ~TESTNUM,
  "Not a term",            11.5,   0,      2.1,    "mmol/L",  NA,       1,
  NA_character_,           11.5,   0,      2.1,    "mmol/L",  NA,       2,
  "Hypertriglyceridemia",  11.5,   0,      2.1,    "mmol/L",  "4",      3,
  "Hypertriglyceridemia",  11.4,   0,      2.1,    "mmol/L",  "3",      4,
  "Hypertriglyceridemia",  5.8,    0,      2.1,    "mmol/L",  "3",      5,
  "Hypertriglyceridemia",  5.7,    0,      2.1,    "mmol/L",  "2",      6,
  "Hypertriglyceridemia",  3.43,   0,      2.1,    "mmol/L",  "2",      7,
  "Hypertriglyceridemia",  3.42,   0,      2.1,    "mmol/L",  "1",      8,
  "Hypertriglyceridemia",  1.72,   0,      2.1,    "mmol/L",  "1",      9,
  "Hypertriglyceridemia",  1.71,   0,      2.1,    "mmol/L",  "0",      10,
  # Unit missing cannot grade
  "Hypertriglyceridemia",  1.71,   0,      2.1,    NA,        NA,       11,
  # AVAL missing cannot grade
  "Hypertriglyceridemia",  NA,     0,      2.1,    "mmol/L",  NA,       12,
)

exp_trigi_cv <- tibble::tribble(
  ~ATOXDSCH,               ~AVAL,  ~ANRLO, ~ANRHI,  ~AVALU,   ~ATOXGRH, ~TESTNUM,
  "Not a term",            1001,   0,      50,      "mg/dL",  NA,       1,
  NA_character_,           1001,   0,      50,      "mg/dL",  NA,       2,
  "Hypertriglyceridemia",  1001,   0,      50,      "mg/dL",  "4",      3,
  "Hypertriglyceridemia",  1000,   0,      50,      "mg/dL",  "3",      4,
  "Hypertriglyceridemia",  501,    0,      50,      "mg/dL",  "3",      5,
  "Hypertriglyceridemia",  500,    0,      50,      "mg/dL",  "2",      6,
  "Hypertriglyceridemia",  301,    0,      50,      "mg/dL",  "2",      7,
  "Hypertriglyceridemia",  300,    0,      50,      "mg/dL",  "1",      8,
  "Hypertriglyceridemia",  151,    0,      50,      "mg/dL",  "1",      9,
  "Hypertriglyceridemia",  150,    0,      50,      "mg/dL",  "0",      10,
  # Unit missing cannot grade
  "Hypertriglyceridemia",  150,    0,      50,      NA,       NA,       11,
  # AVAL missing cannot grade
  "Hypertriglyceridemia",  NA,     0,      50,      "mg/dL",  NA,       12,
)

exp_urici_si <- tibble::tribble(
  ~ATOXDSCH,        ~AVAL,  ~ANRLO,  ~ANRHI,    ~AVALU, ~ATOXGRH, ~TESTNUM,
  "Not a term",     591,    0,       200,     "umol/L",       NA, 1,
  NA_character_,    591,    0,       200,     "umol/L",       NA, 2,
  # ANRHI not missing
  "Hyperuricemia",  591,    0,       200,     "umol/L",      "4", 3,
  "Hyperuricemia",  590,    0,       200,     "umol/L",      "3", 4,
  "Hyperuricemia",  201,    0,       200,     "umol/L",      "3", 5,
  "Hyperuricemia",  200,    0,       200,     "umol/L",      "0", 6,
  # ANRHI missing - can grade 4
  "Hyperuricemia",  591,    0,       NA,      "umol/L",      "4", 7,
  # ANRHI missing - can NOT grade 0 or 3
  "Hyperuricemia",  590,    0,       NA,      "umol/L",       NA, 8,
  "Hyperuricemia",  201,    0,       NA,      "umol/L",       NA, 9,
  "Hyperuricemia",  200,    0,       NA,      "umol/L",       NA, 10,
  # Unit missing cannot grade
  "Hyperuricemia",  200,    0,       200,           NA,      "0", 11,
  # AVAL missing cannot grade
  "Hyperuricemia",  NA,     0,       200,     "umol/L",       NA, 12,
)

exp_urici_cv <- tibble::tribble(
  ~ATOXDSCH,        ~AVAL, ~ANRLO, ~ANRHI,  ~AVALU,  ~ATOXGRH, ~TESTNUM,
  "Not a term",     11,    0,      5,      "mg/dL",  NA,       1,
  NA_character_,    11,    0,      5,      "mg/dL",  NA,       2,
  # ANRHI not missing
  "Hyperuricemia",  11,    0,      5,      "mg/dL",  "4",      3,
  "Hyperuricemia",  10,    0,      5,      "mg/dL",  "3",      4,
  "Hyperuricemia",  6,     0,      5,      "mg/dL",  "3",      5,
  "Hyperuricemia",  5,     0,      5,      "mg/dL",  "0",      6,
  # ANRHI missing - can grade 4
  "Hyperuricemia",  11,    0,      NA,     "mg/dL",  "4",      7,
  # ANRHI missing - can NOT grade 0 or 3
  "Hyperuricemia",  10,    0,      NA,     "mg/dL",  NA,       8,
  "Hyperuricemia",  6,     0,      NA,     "mg/dL",  NA,       9,
  "Hyperuricemia",  5,     0,      NA,     "mg/dL",  NA,       10,
  # Unit missing cannot grade
  "Hyperuricemia",  5,     0,      5,           NA,  "0",      11,
  # AVAL missing cannot grade
  "Hyperuricemia",  NA,    0,      5,      "mg/dL",  NA,       12,
)

exp_albd_si <- tibble::tribble(
  ~ATOXDSCL,          ~AVAL,  ~ANRLO, ~ANRHI, ~AVALU, ~ATOXGRL, ~TESTNUM,
  "Not a term",       19,     40,     100,    "G/L",  NA,       1,
  NA_character_,      19,     40,     100,    "G/L",  NA,       2,
  # ANRLO not missing
  "Hypoalbuminemia",  19,     40,     100,    "G/L",  "3",      3,
  "Hypoalbuminemia",  20,     40,     100,    "G/L",  "2",      4,
  "Hypoalbuminemia",  29,     40,     100,    "G/L",  "2",      5,
  "Hypoalbuminemia",  30,     40,     100,    "G/L",  "1",      6,
  "Hypoalbuminemia",  39,     40,     100,    "G/L",  "1",      7,
  "Hypoalbuminemia",  40,     40,     100,    "G/L",  "0",      8,
  "Hypoalbuminemia",  40,     40,     NA,     "G/L",  "0",      9,
  # ANRLO missing - can grade 2-3
  "Hypoalbuminemia",  19,     NA,     100,    "G/L",  "3",      10,
  "Hypoalbuminemia",  20,     NA,     100,    "G/L",  "2",      11,
  "Hypoalbuminemia",  29,     NA,     100,    "G/L",  "2",      12,
  # ANRLO missing - can NOT grade 0 or 1
  "Hypoalbuminemia",  30,     NA,     100,    "G/L",  NA,       13,
  "Hypoalbuminemia",  39,     NA,     100,    "G/L",  NA,       14,
  "Hypoalbuminemia",  40,     NA,     100,    "G/L",  NA,       15,
  # Unit missing cannot grade
  "Hypoalbuminemia",  40,     40,     100,    NA,     NA,       16,
  # AVAL missing cannot grade
  "Hypoalbuminemia",  NA,     40,     100,    "G/L",  NA,       17,
)

exp_albd_cv <- exp_albd_si %>%
  mutate(
    AVAL = AVAL / 10,
    ANRLO = ANRLO / 10,
    ANRHI = ANRHI / 10,
    AVALU = if_else(is.na(AVALU), NA_character_, "g/dL")
  )

exp_calcd_si <- tibble::tribble(
  ~ATOXDSCL,       ~AVAL,  ~ANRLO, ~ANRHI, ~AVALU,    ~ATOXGRL, ~TESTNUM,
  "Not a term",    1.4,    4,      100,    "mmol/L",  NA,       1,
  NA_character_,   1.4,    4,      100,    "mmol/L",  NA,       2,
  # ANRLO not missing
  "Hypocalcemia",  1.4,    4,      100,    "mmol/L",  "4",      3,
  "Hypocalcemia",  1.5,    4,      100,    "mmol/L",  "3",      4,
  "Hypocalcemia",  1.74,   4,      100,    "mmol/L",  "3",      5,
  "Hypocalcemia",  1.75,   4,      100,    "mmol/L",  "2",      6,
  "Hypocalcemia",  1.9,    4,      100,    "mmol/L",  "2",      7,
  "Hypocalcemia",  2,      4,      100,    "mmol/L",  "1",      8,
  "Hypocalcemia",  3,      4,      100,    "mmol/L",  "1",      9,
  "Hypocalcemia",  4,      4,      100,    "mmol/L",  "0",      10,
  # ANRLO missing - can grade 2-4
  "Hypocalcemia",  1.4,    4,      NA,     "mmol/L",  "4",      11,
  "Hypocalcemia",  1.5,    4,      NA,     "mmol/L",  "3",      12,
  "Hypocalcemia",  1.74,   4,      NA,     "mmol/L",  "3",      13,
  "Hypocalcemia",  1.75,   4,      NA,     "mmol/L",  "2",      14,
  "Hypocalcemia",  1.9,    4,      NA,     "mmol/L",  "2",      15,
  # ANRLO missing - can NOT grade 0 or 1
  "Hypocalcemia",  2,      4,      NA,     "mmol/L",  "1",      16,
  "Hypocalcemia",  3,      4,      NA,     "mmol/L",  "1",      17,
  "Hypocalcemia",  4,      4,      NA,     "mmol/L",  "0",      18,
  # Unit missing cannot grade
  "Hypocalcemia",  4,      4,      100,    NA,        NA,       19,
  # AVAL missing cannot grade
  "Hypocalcemia",  NA,     4,      100,    "mmol/L",  NA,       20,
)

exp_calcd_cv <- tibble::tribble(
  ~ATOXDSCL,       ~AVAL,  ~ANRLO,  ~ANRHI,   ~AVALU,  ~ATOXGRL, ~TESTNUM,
  "Not a term",    5.9,    9,       100,     "mg/dL",        NA, 1,
  NA_character_,   5.9,    9,       100,     "mg/dL",        NA, 2,
  # ANRLO not missing
  "Hypocalcemia",  5.9,    9,       100,     "mg/dL",       "4", 3,
  "Hypocalcemia",  6,      9,       100,     "mg/dL",       "3", 4,
  "Hypocalcemia",  6.9,    9,       100,     "mg/dL",       "3", 5,
  "Hypocalcemia",  7,      9,       100,     "mg/dL",       "2", 6,
  "Hypocalcemia",  7.9,    9,       100,     "mg/dL",       "2", 7,
  "Hypocalcemia",  8,      9,       100,     "mg/dL",       "1", 8,
  "Hypocalcemia",  8.9,    9,       100,     "mg/dL",       "1", 9,
  "Hypocalcemia",  9,      9,       100,     "mg/dL",       "0", 10,
  # ANRLO missing - can grade 2-4
  "Hypocalcemia",  5.9,    9,       NA,      "mg/dL",       "4", 11,
  "Hypocalcemia",  6,      9,       NA,      "mg/dL",       "3", 12,
  "Hypocalcemia",  6.9,    9,       NA,      "mg/dL",       "3", 13,
  "Hypocalcemia",  7,      9,       NA,      "mg/dL",       "2", 14,
  "Hypocalcemia",  7.9,    9,       NA,      "mg/dL",       "2", 15,
  # ANRLO missing - can NOT grade 0 or 1
  "Hypocalcemia",  8,      9,       NA,      "mg/dL",       "1", 16,
  "Hypocalcemia",  8.9,    9,       NA,      "mg/dL",       "1", 17,
  "Hypocalcemia",  9,      9,       NA,      "mg/dL",       "0", 18,
  # Unit missing cannot grade
  "Hypocalcemia",  9,      9,       100,          NA,        NA, 19,
  # AVAL missing cannot grade
  "Hypocalcemia",  NA,     9,       100,     "mg/dL",        NA, 20,
)

exp_caliond_si <- tibble::tribble(
  ~ATOXDSCL,                 ~AVAL,  ~ANRLO, ~ANRHI, ~AVALU,    ~ATOXGRL, ~TESTNUM,
  "Not a term",              0.79,   1.3,    100,    "mmol/L",  NA,       1,
  NA_character_,             0.79,   1.3,    100,    "mmol/L",  NA,       2,
  # ANRLO not missing
  "Hypocalcemia (Ionized)",  0.79,   1.3,    100,    "mmol/L",  "4",      3,
  "Hypocalcemia (Ionized)",  0.8,    1.3,    100,    "mmol/L",  "3",      4,
  "Hypocalcemia (Ionized)",  0.89,   1.3,    100,    "mmol/L",  "3",      5,
  "Hypocalcemia (Ionized)",  0.9,    1.3,    100,    "mmol/L",  "2",      6,
  "Hypocalcemia (Ionized)",  0.99,   1.3,    100,    "mmol/L",  "2",      7,
  "Hypocalcemia (Ionized)",  1,      1.3,    100,    "mmol/L",  "1",      8,
  "Hypocalcemia (Ionized)",  1.29,   1.3,    100,    "mmol/L",  "1",      9,
  "Hypocalcemia (Ionized)",  1.3,    1.3,    100,    "mmol/L",  "0",      10,
  # ANRLO missing - can grade 2-4
  "Hypocalcemia (Ionized)",  0.79,   NA,     100,    "mmol/L",  "4",      11,
  "Hypocalcemia (Ionized)",  0.8,    NA,     100,    "mmol/L",  "3",      12,
  "Hypocalcemia (Ionized)",  0.89,   NA,     100,    "mmol/L",  "3",      13,
  "Hypocalcemia (Ionized)",  0.9,    NA,     100,    "mmol/L",  "2",      14,
  "Hypocalcemia (Ionized)",  0.99,   NA,     100,    "mmol/L",  "2",      15,
  # ANRLO missing - can NOT grade 0 or 1
  "Hypocalcemia (Ionized)",  1,      1.3,    NA,     "mmol/L",  "1",      16,
  "Hypocalcemia (Ionized)",  1.29,   1.3,    NA,     "mmol/L",  "1",      17,
  "Hypocalcemia (Ionized)",  1.3,    1.3,    NA,     "mmol/L",  "0",      18,
  # Unit missing cannot grade
  "Hypocalcemia (Ionized)",  1.3,    1.3,    100,    NA,        NA,       19,
  # AVAL missing cannot grade
  "Hypocalcemia (Ionized)",  NA,     1.3,    100,    "mmol/L",  NA,       20,
)

exp_caliond_cv <- exp_caliond_si %>%
  mutate(
    AVAL = AVAL * 4,
    ANRHI = ANRHI * 4,
    ANRLO = ANRLO * 4,
    AVALU = if_else(is.na(AVALU), NA_character_, "mg/dL")
  )

exp_glycd_si <- tibble::tribble(
  ~ATOXDSCL,       ~AVAL,  ~ANRLO,  ~ANRHI,    ~AVALU,  ~ATOXGRL, ~TESTNUM,
  "Not a term",    1.69,   4,       100,     "mmol/L",        NA, 1,
  NA_character_,   1.69,   4,       100,     "mmol/L",        NA, 2,
  # ANRLO not missing
  "Hypoglycemia",  1.69,   4,       100,     "mmol/L",       "4", 3,
  "Hypoglycemia",  1.7,    4,       100,     "mmol/L",       "3", 4,
  "Hypoglycemia",  2.19,   4,       100,     "mmol/L",       "3", 5,
  "Hypoglycemia",  2.2,    4,       100,     "mmol/L",       "2", 6,
  "Hypoglycemia",  2.9,    4,       100,     "mmol/L",       "2", 7,
  "Hypoglycemia",  3,      4,       100,     "mmol/L",       "1", 8,
  "Hypoglycemia",  3.9,    4,       100,     "mmol/L",       "1", 9,
  "Hypoglycemia",  4,      4,       100,     "mmol/L",       "0", 10,
  # ANRLO missing - can grade 2-4
  "Hypoglycemia",  1.69,   NA,      100,     "mmol/L",       "4", 11,
  "Hypoglycemia",  1.7,    NA,      100,     "mmol/L",       "3", 12,
  "Hypoglycemia",  2.19,   NA,      100,     "mmol/L",       "3", 13,
  "Hypoglycemia",  2.2,    NA,      100,     "mmol/L",       "2", 14,
  "Hypoglycemia",  2.9,    NA,      100,     "mmol/L",       "2", 15,
  # ANRLO missing - can NOT grade 0 or 1 mg/dL
  "Hypoglycemia",  3,      NA,      100,     "mmol/L",        NA, 16,
  "Hypoglycemia",  3.9,    NA,      100,     "mmol/L",        NA, 17,
  "Hypoglycemia",  4,      NA,      100,     "mmol/L",        NA, 18,
  # Unit missing cannot grade
  "Hypoglycemia",  4,      4,       100,           NA,        NA, 19,
  # AVAL missing cannot grade
  "Hypoglycemia",  NA,     4,       100,     "mmol/L",        NA, 20,
)

exp_glycd_cv <- tibble::tribble(
  ~ATOXDSCL,       ~AVAL,  ~ANRLO, ~ANRHI,   ~AVALU,  ~ATOXGRL, ~TESTNUM,
  "Not a term",    29,     70,     100,     "mg/dL",        NA, 1,
  NA_character_,   29,     70,     100,      "g/dL",        NA, 2,
  # ANRLO not missing
  "Hypoglycemia",  29,     70,     100,     "mg/dL",       "4", 3,
  "Hypoglycemia",  30,     70,     100,     "mg/dL",       "3", 4,
  "Hypoglycemia",  39,     70,     100,     "mg/dL",       "3", 5,
  "Hypoglycemia",  40,     70,     100,     "mg/dL",       "2", 6,
  "Hypoglycemia",  54,     70,     100,     "mg/dL",       "2", 7,
  "Hypoglycemia",  55,     70,     100,     "mg/dL",       "1", 8,
  "Hypoglycemia",  69,     70,     100,     "mg/dL",       "1", 9,
  "Hypoglycemia",  70,     70,     100,     "mg/dL",       "0", 10,
  # ANRLO missing - can grade 2-4
  "Hypoglycemia",  29,     NA,     100,     "mg/dL",       "4", 11,
  "Hypoglycemia",  30,     NA,     100,     "mg/dL",       "3", 12,
  "Hypoglycemia",  39,     NA,     100,     "mg/dL",       "3", 13,
  "Hypoglycemia",  40,     NA,     100,     "mg/dL",       "2", 14,
  "Hypoglycemia",  54,     NA,     100,     "mg/dL",       "2", 15,
  # ANRLO missing - can NOT grade 0 or 1
  "Hypoglycemia",  55,     NA,     100,     "mg/dL",        NA, 16,
  "Hypoglycemia",  69,     NA,     100,     "mg/dL",        NA, 17,
  "Hypoglycemia",  70,     NA,     100,     "mg/dL",        NA, 18,
  # Unit missing cannot grade
  "Hypoglycemia",  70,     70,     100,          NA,        NA, 19,
  # AVAL missing cannot grade
  "Hypoglycemia",  NA,     70,     100,     "mg/dL",        NA, 20,
)

exp_kaled <- tibble::tribble(
  ~ATOXDSCL,      ~AVAL,  ~ANRLO, ~ANRHI, ~AVALU,    ~ATOXGRL, ~TESTNUM,
  "Not a term",   2.49,   4,      100,    "mmol/L",  NA,       1,
  NA_character_,  2.49,   4,      100,    "mmol/L",  NA,       2,
  # ANRLO not missing
  "Hypokalemia",  2.49,   4,      100,    "mmol/L",  "4",      3,
  "Hypokalemia",  2.5,    4,      100,    "mmol/L",  "3",      4,
  "Hypokalemia",  2.9,    4,      100,    "mmol/L",  "3",      5,
  "Hypokalemia",  3,      4,      100,    "mmol/L",  "2",      6,
  "Hypokalemia",  3.9,    4,      100,    "mmol/L",  "2",      7,
  "Hypokalemia",  4,      4,      100,    "mmol/L",  "0",      8,
  # ANRLO missing - can grade 3-4
  "Hypokalemia",  2.49,   NA,     100,    "mmol/L",  "4",      9,
  "Hypokalemia",  2.5,    NA,     100,    "mmol/L",  "3",      10,
  "Hypokalemia",  2.9,    NA,     100,    "mmol/L",  "3",      11,
  # ANRLO missing - can NOT grade 0 or 2
  "Hypokalemia",  3,      NA,     100,    "mmol/L",  NA,       12,
  "Hypokalemia",  3.9,    NA,     100,    "mmol/L",  NA,       13,
  "Hypokalemia",  4,      NA,     NA,     "mmol/L",  NA,       14,
  # Unit missing cannot grade
  "Hypokalemia",  4,      4,      100,    NA,        NA,       15,
  # AVAL missing cannot grade
  "Hypokalemia",  NA,     4,      100,    "mmol/L",  NA,       16,
)

exp_magnd_si <- tibble::tribble(
  ~ATOXDSCL,         ~AVAL,  ~ANRLO, ~ANRHI, ~AVALU,    ~ATOXGRL, ~TESTNUM,
  "Not a term",      0.29,   1,      100,    "mmol/L",  NA,       1,
  NA_character_,     0.29,   1,      100,    "mmol/L",  NA,       2,
  # ANRLO not missing
  "Hypomagnesemia",  0.29,   1,      100,    "mmol/L",  "4",      3,
  "Hypomagnesemia",  0.3,    1,      100,    "mmol/L",  "3",      4,
  "Hypomagnesemia",  0.39,   1,      100,    "mmol/L",  "3",      5,
  "Hypomagnesemia",  0.4,    1,      100,    "mmol/L",  "2",      6,
  "Hypomagnesemia",  0.49,   1,      100,    "mmol/L",  "2",      7,
  "Hypomagnesemia",  0.5,    1,      100,    "mmol/L",  "1",      8,
  "Hypomagnesemia",  0.9,    1,      100,    "mmol/L",  "1",      9,
  "Hypomagnesemia",  1,      1,      100,    "mmol/L",  "0",      10,
  # ANRLO missing - can grade 2-4
  "Hypomagnesemia",  0.29,   NA,     100,    "mmol/L",  "4",      11,
  "Hypomagnesemia",  0.3,    NA,     100,    "mmol/L",  "3",      12,
  "Hypomagnesemia",  0.39,   NA,     100,    "mmol/L",  "3",      13,
  "Hypomagnesemia",  0.4,    NA,     100,    "mmol/L",  "2",      14,
  "Hypomagnesemia",  0.49,   NA,     100,    "mmol/L",  "2",      15,
  # ANRLO missing - can NOT grade 0 or 1
  "Hypomagnesemia",  0.5,    NA,     100,    "mmol/L",  NA,       16,
  "Hypomagnesemia",  0.9,    NA,     100,    "mmol/L",  NA,       17,
  "Hypomagnesemia",  1,      NA,     100,    "mmol/L",  NA,       18,
  # Unit missing cannot grade
  "Hypomagnesemia",  1,      1,      100,    NA,        NA,       19,
  # AVAL missing cannot grade
  "Hypomagnesemia",  NA,     1,      100,    "mmol/L",  NA,       20,
)

exp_magnd_cv <- tibble::tribble(
  ~ATOXDSCL,         ~AVAL,  ~ANRLO,  ~ANRHI,   ~AVALU,  ~ATOXGRL, ~TESTNUM,
  "Not a term",      0.69,   1.5,     10,      "mg/dL",        NA, 1,
  NA_character_,     0.69,   1.5,     10,      "mg/dL",        NA, 2,
  # ANRLO not missing
  "Hypomagnesemia",  0.69,   1.5,     10,      "mg/dL",       "4", 3,
  "Hypomagnesemia",  0.7,    1.5,     10,      "mg/dL",       "3", 4,
  "Hypomagnesemia",  0.89,   1.5,     10,      "mg/dL",       "3", 5,
  "Hypomagnesemia",  0.9,    1.5,     10,      "mg/dL",       "2", 6,
  "Hypomagnesemia",  1.19,   1.5,     10,      "mg/dL",       "2", 7,
  "Hypomagnesemia",  1.2,    1.5,     10,      "mg/dL",       "1", 8,
  "Hypomagnesemia",  1.49,   1.5,     10,      "mg/dL",       "1", 9,
  "Hypomagnesemia",  1.5,    1.5,     10,      "mg/dL",       "0", 10,
  # ANRLO missing - can grade 2-4
  "Hypomagnesemia",  0.69,   NA,      10,      "mg/dL",       "4", 11,
  "Hypomagnesemia",  0.7,    NA,      10,      "mg/dL",       "3", 12,
  "Hypomagnesemia",  0.89,   NA,      10,      "mg/dL",       "3", 13,
  "Hypomagnesemia",  0.9,    NA,      10,      "mg/dL",       "2", 14,
  "Hypomagnesemia",  1.19,   NA,      10,      "mg/dL",       "2", 15,
  # ANRLO missing - can NOT grade 0 or 1
  "Hypomagnesemia",  1.2,    NA,      10,      "mg/dL",        NA, 16,
  "Hypomagnesemia",  1.49,   NA,      10,      "mg/dL",        NA, 17,
  "Hypomagnesemia",  1.5,    NA,      10,      "mg/dL",        NA, 18,
  # Unit missing cannot grade
  "Hypomagnesemia",  1.5,    1.5,     10,           NA,        NA, 19,
  # AVAL missing cannot grade
  "Hypomagnesemia",  NA,     1.5,     10,      "mg/dL",        NA, 20,
)

exp_natrd <- tibble::tribble(
  ~ATOXDSCL,       ~AVAL,  ~ANRLO,   ~ANRHI, ~AVALU,    ~ATOXGRL, ~TESTNUM,
  "Not a term",    119,    140,      100,    "mmol/L",  NA,       1,
  NA_character_,   119,    140,      100,    "mmol/L",  NA,       2,
  # ANRLO not missing
  "Hyponatremia",  119,    140,      100,    "mmol/L",  "4",      3,
  "Hyponatremia",  120,    140,      100,    "mmol/L",  "3",      4,
  "Hyponatremia",  129,    140,      100,    "mmol/L",  "3",      5,
  "Hyponatremia",  130,    140,      100,    "mmol/L",  "1",      6,
  "Hyponatremia",  139,    140,      100,    "mmol/L",  "1",      7,
  "Hyponatremia",  140,    140,      100,    "mmol/L",  "0",      8,
  # ANRLO missing - can grade 3-4
  "Hyponatremia",  119,    NA,       100,    "mmol/L",  "4",      9,
  "Hyponatremia",  120,    NA,       100,    "mmol/L",  "3",      10,
  "Hyponatremia",  129,    NA,       100,    "mmol/L",  "3",      11,
  # ANRLO missing - can NOT grade 0 or 1
  "Hyponatremia",  130,    NA,       100,    "mmol/L",  NA,       12,
  "Hyponatremia",  139,    NA,       100,    "mmol/L",  NA,       13,
  "Hyponatremia",  140,    NA,       100,    "mmol/L",  NA,       14,
  # Unit missing cannot grade
  "Hyponatremia",  140,    140,      100,    NA,        NA,       15,
  # AVAL missing cannot grade
  "Hyponatremia",  NA,     140,      100,    "mmol/L",  NA,       16,
)

exp_albl_daids_si <- tibble::tribble(
  ~ATOXDSCL,      ~AVAL,  ~ANRLO, ~AVALU, ~ATOXGRL, ~TESTNUM,
  "Not a term",   35,     40,     "g/L",  NA,       1,
  NA_character_,  35,     40,     "g/L",  NA,       2,
  # ANRLO not missing
  "Albumin, Low", 19,     40,     "g/L",  "3",      3,
  "Albumin, Low", 20,     40,     "g/L",  "2",      4,
  "Albumin, Low", 29,     40,     "g/L",  "2",      5,
  "Albumin, Low", 30,     40,     "g/L",  "1",      6,
  "Albumin, Low", 39,     40,     "g/L",  "1",      7,
  "Albumin, Low", 40,     40,     "g/L",  "0",      8,
  # ANRLO missing - can grade 2 and 3
  "Albumin, Low", 19,     NA,     "g/L",  "3",      9,
  "Albumin, Low", 20,     NA,     "g/L",  "2",      10,
  "Albumin, Low", 29,     NA,     "g/L",  "2",      11,
  # ANRLO missing - can NOT grade 0 or 1
  "Albumin, Low", 30,     NA,     "g/L",  NA,       12,
  "Albumin, Low", 39,     NA,     "g/L",  NA,       13,
  "Albumin, Low", 40,     NA,     "g/L",  NA,       14,
  # AVALU missing cannot grade
  "Albumin, Low", 40,     40,     NA,     NA,       15,
  # AVAL missing cannot grade
  "Albumin, Low", NA,     40,     "g/L",  NA,       16,
)

exp_albl_daids_cv <- exp_albl_daids_si %>%
  mutate(
    AVAL = AVAL / 10,
    ANRLO = ANRLO / 10,
    AVALU = if_else(is.na(AVALU), NA_character_, "g/dL")
  )

exp_bicad_daids <- tibble::tribble(
  ~ATOXDSCL,           ~AVAL,  ~ANRLO, ~AVALU,    ~ATOXGRL, ~TESTNUM,
  "Not a term",        22,     20,     "mmol/L",  NA,       1,
  NA_character_,       22,     20,     "mmol/L",  NA,       2,
  # ANRLO not missing
  "Bicarbonate, Low",  7.9,    20,     "mmol/L",  "4",      3,
  "Bicarbonate, Low",  8,      20,     "mmol/L",  "3",      4,
  "Bicarbonate, Low",  10.9,   20,     "mmol/L",  "3",      5,
  "Bicarbonate, Low",  11,     20,     "mmol/L",  "2",      6,
  "Bicarbonate, Low",  15.9,   20,     "mmol/L",  "2",      7,
  "Bicarbonate, Low",  16,     20,     "mmol/L",  "1",      8,
  "Bicarbonate, Low",  19,     20,     "mmol/L",  "1",      9,
  "Bicarbonate, Low",  20,     20,     "mmol/L",  "0",      10,
  # ANRLO missing - can grade 2-4
  "Bicarbonate, Low",  7.9,    NA,     "mmol/L",  "4",      11,
  "Bicarbonate, Low",  8,      NA,     "mmol/L",  "3",      12,
  "Bicarbonate, Low",  10.9,   NA,     "mmol/L",  "3",      13,
  "Bicarbonate, Low",  11,     NA,     "mmol/L",  "2",      14,
  "Bicarbonate, Low",  15.9,   NA,     "mmol/L",  "2",      15,
  # ANRLO missing - can NOT grade 0 or 1
  "Bicarbonate, Low",  16,     NA,     "mmol/L",  NA,       16,
  "Bicarbonate, Low",  19,     NA,     "mmol/L",  NA,       17,
  "Bicarbonate, Low",  20,     NA,     "mmol/L",  NA,       18,
  # Unit missing cannot grade
  "Bicarbonate, Low",  20,     20,     NA,        NA,       19,
  # AVAL missing cannot grade
  "Bicarbonate, Low",  NA,     20,     "mmol/L",  NA,       20,
)

exp_dbiligt28d_daids_si <- tibble::tribble(
  ~ATOXDSCH,                  ~AVAL,  ~ANRHI, ~AVALU,   ~ATOXGRH, ~TESTNUM,
  "Not a term",               7,      8,      "umol/L", NA,       1,
  NA_character_,              7,      8,      "umol/L", NA,       2,
  # ANRHI not missing
  "Direct Bilirubin, High",   8.1,    8,      "umol/L", "4",      3,
  "Direct Bilirubin, High",   8,      8,      "umol/L", "0",      4,
  "Direct Bilirubin, High",   7.9,    8,      "umol/L", "0",      5,
  # ANRHI missing cannot still grade 0 or 4
  "Direct Bilirubin, High",   8.1,    NA,     "umol/L", NA,       6,
  "Direct Bilirubin, High",   8,      NA,     "umol/L", NA,       7,
  "Direct Bilirubin, High",   7.9,    NA,     "umol/L", NA,       8,
  # AVAL missing cannot grade
  "Direct Bilirubin, High",   NA,     8,      "umol/L", NA,       9,
) %>%
  mutate(
    BRTHDT = lubridate::ymd("2023-01-01"),
    ADT = lubridate::ymd("2023-01-30")
  )

exp_dbiligt28d_daids_cv <- exp_dbiligt28d_daids_si %>%
  mutate(
    AVALU = "mg/dL",
    TESTNUM = TESTNUM + 30
  )

exp_dbilile28d_daids_si <- tibble::tribble(
  ~ATOXDSCH,                  ~AVAL,  ~ANRHI, ~AVALU,   ~ATOXGRH, ~TESTNUM,
  "Not a term",               7,      8,      "umol/L", NA,       10,
  NA_character_,              7,      8,      "umol/L", NA,       11,
  # ANRHI not missing
  "Direct Bilirubin, High",   34.3,   8,      "umol/L", "4",      12,
  "Direct Bilirubin, High",   34.2,   8,      "umol/L", "3",      13,
  "Direct Bilirubin, High",   25.66,  8,      "umol/L", "3",      14,
  "Direct Bilirubin, High",   25.65,  8,      "umol/L", "2",      15,
  "Direct Bilirubin, High",   17.19,  8,      "umol/L", "2",      16,
  "Direct Bilirubin, High",   17.1,   8,      "umol/L", "1",      17,
  "Direct Bilirubin, High",   8,      8,      "umol/L", "1",      18,
  "Direct Bilirubin, High",   7.9,    8,      "umol/L", "0",      19,
  # ANRHI missing can still grade 2 - 4
  "Direct Bilirubin, High",   34.3,   NA,     "umol/L", "4",      20,
  "Direct Bilirubin, High",   34.2,   NA,     "umol/L", "3",      21,
  "Direct Bilirubin, High",   25.66,  NA,     "umol/L", "3",      22,
  "Direct Bilirubin, High",   25.65,  NA,     "umol/L", "2",      23,
  "Direct Bilirubin, High",   17.19,  NA,     "umol/L", "2",      24,
  # ANRHI missing cannot still grade 0 - 1
  "Direct Bilirubin, High",   17.1,   NA,     "umol/L", NA,       25,
  "Direct Bilirubin, High",   8,      NA,     "umol/L", NA,       26,
  "Direct Bilirubin, High",   7.9,    NA,     "umol/L", NA,       27,
  # AVAL missing cannot grade
  "Direct Bilirubin, High",   NA,     8,      "umol/L", NA,       28,
) %>%
  mutate(
    BRTHDT = lubridate::ymd("2023-01-01"),
    ADT = lubridate::ymd("2023-01-29")
  )

exp_dbilile28d_daids_cv <- tibble::tribble(
  ~ATOXDSCH,                  ~AVAL,  ~ANRHI, ~AVALU,  ~ATOXGRH, ~TESTNUM,
  "Not a term",               2.1,    0.5,    "mg/dL", NA,       10,
  NA_character_,              2.1,    0.5,    "mg/dL", NA,       11,
  # ANRHI not missing
  "Direct Bilirubin, High",   2.1,    0.5,    "mg/dL", "4",      12,
  "Direct Bilirubin, High",   2,      0.5,    "mg/dL", "3",      13,
  "Direct Bilirubin, High",   1.51,   0.5,    "mg/dL", "3",      14,
  "Direct Bilirubin, High",   1.5,    0.5,    "mg/dL", "2",      15,
  "Direct Bilirubin, High",   1.1,    0.5,    "mg/dL", "2",      16,
  "Direct Bilirubin, High",   1,      0.5,    "mg/dL", "1",      17,
  "Direct Bilirubin, High",   0.5,    0.5,    "mg/dL", "1",      18,
  "Direct Bilirubin, High",   0.49,   0.5,    "mg/dL", "0",      19,
  # ANRHI missing can still grade 2 - 4
  "Direct Bilirubin, High",   2.1,    NA,     "mg/dL", "4",      20,
  "Direct Bilirubin, High",   2,      NA,     "mg/dL", "3",      21,
  "Direct Bilirubin, High",   1.51,   NA,     "mg/dL", "3",      22,
  "Direct Bilirubin, High",   1.5,    NA,     "mg/dL", "2",      23,
  "Direct Bilirubin, High",   1.1,    NA,     "mg/dL", "2",      24,
  # ANRHI missing cannot still grade 0 - 1
  "Direct Bilirubin, High",   1,      NA,     "mg/dL", NA,       25,
  "Direct Bilirubin, High",   0.5,    NA,     "mg/dL", NA,       26,
  "Direct Bilirubin, High",   0.49,   NA,     "mg/dL", NA,       27,
  # AVAL missing cannot grade
  "Direct Bilirubin, High",   NA,     0.5,    "mg/dL", NA,       28,
) %>%
  mutate(
    BRTHDT = lubridate::ymd("2023-01-01"),
    ADT = lubridate::ymd("2023-01-29")
  )

exp_dbilinoage_daids_si <- exp_dbilile28d_daids_si %>%
  filter(TESTNUM %in% c(18, 19)) %>%
  mutate(
    ADT = if_else(TESTNUM == 18, NA, ADT),
    BRTHDT = if_else(TESTNUM == 19, NA, BRTHDT),
    ATOXGRH = NA_character_,
    TESTNUM = if_else(TESTNUM == 18, 29, 30)
  )

exp_dbilinoage_daids_cv <- exp_dbilile28d_daids_cv %>%
  filter(TESTNUM %in% c(48, 49)) %>%
  mutate(
    ADT = if_else(TESTNUM == 48, NA, ADT),
    BRTHDT = if_else(TESTNUM == 49, NA, BRTHDT),
    ATOXGRH = NA_character_,
    TESTNUM = if_else(TESTNUM == 48, 59, 60)
  )

exp_dbili_daids_si <- exp_dbilinoage_daids_si %>%
  bind_rows(
    exp_dbilile28d_daids_si,
    exp_dbiligt28d_daids_si
  )

exp_dbili_daids_cv <- exp_dbilinoage_daids_cv %>%
  bind_rows(
    exp_dbilile28d_daids_cv,
    exp_dbiligt28d_daids_cv
  )

exp_tbiligt28d_daids <- tibble::tribble(
  ~ATOXDSCH,                ~AVAL,  ~ANRHI, ~AVALU,    ~ATOXGRH, ~TESTNUM,
  "Not a term",             9,      10,     "umol/L",  NA,       1,
  NA_character_,            9,      10,     "umol/L",  NA,       2,
  # ANRHI not missing
  "Total Bilirubin, High",  50,     10,     "umol/L",  "4",      3,
  "Total Bilirubin, High",  49,     10,     "umol/L",  "3",      4,
  "Total Bilirubin, High",  26,     10,     "umol/L",  "3",      5,
  "Total Bilirubin, High",  25.9,   10,     "umol/L",  "2",      6,
  "Total Bilirubin, High",  16,     10,     "umol/L",  "2",      7,
  "Total Bilirubin, High",  15.9,   10,     "umol/L",  "1",      8,
  "Total Bilirubin, High",  11,     10,     "umol/L",  "1",      9,
  "Total Bilirubin, High",  10.9,   10,     "umol/L",  "0",      10,
  # Unit missing can grade - grade based on comparison of AVAL with ANRHI
  "Total Bilirubin, High",  10.9,   10,     NA,        "0",      11,
  # ANRHI missing - cannot grade
  "Total Bilirubin, High",  10.9,   NA,     "umol/L",  NA,       12,
  # AVAL missing cannot grade
  "Total Bilirubin, High",  NA,     10,     "umol/L",  NA,       13,
) %>%
  mutate(
    BRTHDT = lubridate::ymd("2023-01-01"),
    ADT = lubridate::ymd("2023-01-30")
  )

exp_tbilile28d_daids <- exp_tbiligt28d_daids %>%
  mutate(
    ADT = lubridate::ymd("2023-01-29"),
    ATOXGRH = NA_character_,
    TESTNUM = TESTNUM + 13
  )

exp_tbilinoage_daids <- exp_tbiligt28d_daids %>%
  filter(TESTNUM %in% c(10, 11)) %>%
  mutate(
    ADT = if_else(TESTNUM == 10, NA, ADT),
    BRTHDT = if_else(TESTNUM == 11, NA, BRTHDT),
    ATOXGRH = NA_character_,
    TESTNUM = if_else(TESTNUM == 10, 27, 28)
  )

exp_tbili_daids <- exp_tbilinoage_daids %>%
  bind_rows(
    exp_tbiligt28d_daids,
    exp_tbilile28d_daids
  )

exp_calcige7d_daids_si <- tibble::tribble(
  ~ATOXDSCH,       ~AVAL,  ~AVALU,    ~ATOXGRH, ~TESTNUM,
  "Not a term",    3.5,    "mmol/L",  NA,       1,
  NA_character_,   3.5,    "mmol/L",  NA,       2,
  # ANRHI not missing
  "Calcium, High", 3.38,   "mmol/L",  "4",      3,
  "Calcium, High", 3.37,   "mmol/L",  "3",      4,
  "Calcium, High", 3.13,   "mmol/L",  "3",      5,
  "Calcium, High", 3.12,   "mmol/L",  "2",      6,
  "Calcium, High", 2.88,   "mmol/L",  "2",      7,
  "Calcium, High", 2.87,   "mmol/L",  "1",      8,
  "Calcium, High", 2.65,   "mmol/L",  "1",      9,
  "Calcium, High", 2.64,   "mmol/L",  "0",      10,
  # Unit missing cannot grade
  "Calcium, High", 2.5,    NA,        NA,       11,
  # AVAL missing cannot grade
  "Calcium, High", NA,     "mmol/L",  NA,       12,
) %>%
  mutate(
    BRTHDT = lubridate::ymd("2023-01-01"),
    ADT = lubridate::ymd("2023-01-08")
  )

exp_calcige7d_daids_cv <- tibble::tribble(
  ~ATOXDSCH,        ~AVAL,   ~AVALU,  ~ATOXGRH,  ~TESTNUM,
  "Not a term",     13.5,   "mg/dL",        NA,  1,
  NA_character_,    13.5,   "mg/dL",        NA,  2,
  # ANRHI not missing
  "Calcium, High",  13.5,   "mg/dL",       "4",  3,
  "Calcium, High",  13.49,  "mg/dL",       "3",  4,
  "Calcium, High",  12.5,   "mg/dL",       "3",  5,
  "Calcium, High",  12.49,  "mg/dL",       "2",  6,
  "Calcium, High",  11.5,   "mg/dL",       "2",  7,
  "Calcium, High",  11.49,  "mg/dL",       "1",  8,
  "Calcium, High",  10.6,   "mg/dL",       "1",  9,
  "Calcium, High",  10.59,  "mg/dL",       "0",  10,
  # Unit missing cannot grade
  "Calcium, High",  10.5,        NA,        NA,  11,
  # AVAL missing cannot grade
  "Calcium, High",  NA,     "mg/dL",        NA,  12,
) %>%
  mutate(
    BRTHDT = lubridate::ymd("2023-01-01"),
    ADT = lubridate::ymd("2023-01-08")
  )

exp_calcilt7d_daids_si <- tibble::tribble(
  ~ATOXDSCH,       ~AVAL,  ~AVALU,    ~ATOXGRH, ~TESTNUM,
  "Not a term",    3.5,    "mmol/L",  NA,       13,
  NA_character_,   3.5,    "mmol/L",  NA,       14,
  # ANRHI not missing
  "Calcium, High", 3.38,   "mmol/L",  "4",      15,
  "Calcium, High", 3.37,   "mmol/L",  "3",      16,
  "Calcium, High", 3.23,   "mmol/L",  "3",      17,
  "Calcium, High", 3.22,   "mmol/L",  "2",      18,
  "Calcium, High", 3.1,    "mmol/L",  "2",      19,
  "Calcium, High", 3.09,   "mmol/L",  "1",      20,
  "Calcium, High", 2.88,   "mmol/L",  "1",      21,
  "Calcium, High", 2.87,   "mmol/L",  "0",      22,
  # Unit missing cannot grade
  "Calcium, High", 3.5,    NA,        NA,       23,
  # AVAL missing cannot grade
  "Calcium, High", NA,     "mmol/L",  NA,       24,
) %>%
  mutate(
    BRTHDT = lubridate::ymd("2023-01-01"),
    ADT = lubridate::ymd("2023-01-07")
  )

exp_calcilt7d_daids_cv <- tibble::tribble(
  ~ATOXDSCH,       ~AVAL,  ~AVALU,   ~ATOXGRH, ~TESTNUM,
  "Not a term",    13.5,   "mg/dL",  NA,       13,
  NA_character_,   13.5,   "mg/dL",  NA,       14,
  # ANRHI not missing
  "Calcium, High", 13.5,   "mg/dL",  "4",      15,
  "Calcium, High", 13.49,  "mg/dL",  "3",      16,
  "Calcium, High", 12.9,   "mg/dL",  "3",      17,
  "Calcium, High", 12.89,  "mg/dL",  "2",      18,
  "Calcium, High", 12.4,   "mg/dL",  "2",      19,
  "Calcium, High", 12.39,  "mg/dL",  "1",      20,
  "Calcium, High", 11.5,   "mg/dL",  "1",      21,
  "Calcium, High", 11.49,  "mg/dL",  "0",      22,
  # Unit missing cannot grade
  "Calcium, High", 11.49,  NA,       NA,       23,
  # AVAL missing cannot grade
  "Calcium, High", NA,     "mg/dL",  NA,       24,
) %>%
  mutate(
    BRTHDT = lubridate::ymd("2023-01-01"),
    ADT = lubridate::ymd("2023-01-07")
  )

exp_calcinoage_daids_si <- exp_calcige7d_daids_si %>%
  filter(TESTNUM %in% c(9, 10)) %>%
  mutate(
    ADT = if_else(TESTNUM == 9, NA, ADT),
    BRTHDT = if_else(TESTNUM == 10, NA, BRTHDT),
    ATOXGRH = NA_character_,
    TESTNUM = if_else(TESTNUM == 9, 25, 26)
  )

exp_calcinoage_daids_cv <- exp_calcige7d_daids_cv %>%
  filter(TESTNUM %in% c(9, 10)) %>%
  mutate(
    ADT = if_else(TESTNUM == 9, NA, ADT),
    BRTHDT = if_else(TESTNUM == 10, NA, BRTHDT),
    ATOXGRH = NA_character_,
    TESTNUM = if_else(TESTNUM == 9, 25, 26)
  )

exp_calci_daids_si <- exp_calcinoage_daids_si %>%
  bind_rows(
    exp_calcige7d_daids_si,
    exp_calcilt7d_daids_si
  )

exp_calci_daids_cv <- exp_calcinoage_daids_cv %>%
  bind_rows(
    exp_calcige7d_daids_cv,
    exp_calcilt7d_daids_cv
  )

exp_calioni_daids_si <- tibble::tribble(
  ~ATOXDSCH,                 ~AVAL, ~ANRLO, ~ANRHI, ~AVALU,   ~ATOXGRH, ~TESTNUM,
  "Not a term",              1.8,   1.1,    1.4,    "mmol/L", NA,       1,
  NA_character_,             1.79,  1.1,    1.4,    "mmol/L", NA,       2,
  # ANRHI not missing
  "Calcium (Ionized), High", 1.8,   1.1,    1.4,    "mmol/L", "4",      3,
  "Calcium (Ionized), High", 1.79,  1.1,    1.4,    "mmol/L", "3",      4,
  "Calcium (Ionized), High", 1.6,   1.1,    1.4,    "mmol/L", "3",      5,
  "Calcium (Ionized), High", 1.59,  1.1,    1.4,    "mmol/L", "2",      6,
  "Calcium (Ionized), High", 1.5,   1.1,    1.4,    "mmol/L", "2",      7,
  "Calcium (Ionized), High", 1.49,  1.1,    1.4,    "mmol/L", "1",      8,
  "Calcium (Ionized), High", 1.41,  1.1,    1.4,    "mmol/L", "1",      9,
  "Calcium (Ionized), High", 1.4,   1.1,    1.4,    "mmol/L", "0",      10,
  # ANRHI missing - can grade 2-4
  "Calcium (Ionized), High", 1.8,   1.1,    NA,     "mmol/L", "4",      11,
  "Calcium (Ionized), High", 1.79,  1.1,    NA,     "mmol/L", "3",      12,
  "Calcium (Ionized), High", 1.6,   1.1,    NA,     "mmol/L", "3",      13,
  "Calcium (Ionized), High", 1.59,  1.1,    NA,     "mmol/L", "2",      14,
  "Calcium (Ionized), High", 1.5,   1.1,    NA,     "mmol/L", "2",      15,
  # ANRHI missing - can NOT grade 0 or 1
  "Calcium (Ionized), High", 1.49,  1.1,    NA,     "mmol/L", NA,       16,
  "Calcium (Ionized), High", 1.41,  1.1,    NA,     "mmol/L", NA,       17,
  "Calcium (Ionized), High", 1.4,   1.1,    NA,     "mmol/L", NA,       18,
  # Unit missing cannot grade
  "Calcium (Ionized), High", 1.3,   1.1,    1.4,    NA,       NA,       19,
  # AVAL missing cannot grade
  "Calcium (Ionized), High", NA,    1.1,    1.4,    "mmol/L", NA,       20,
)

exp_calioni_daids_cv <- exp_calioni_daids_si %>%
  mutate(
    AVAL = 4 * AVAL,
    ANRLO = 4 * ANRLO,
    ANRHI = 4 * ANRHI,
    AVALU = if_else(is.na(AVALU), NA_character_, "mg/dL")
  )

exp_calcdge7d_daids_si <- tibble::tribble(
  ~ATOXDSCL,       ~AVAL,  ~AVALU,    ~ATOXGRL, ~TESTNUM,
  "Not a term",    2.2,    "mmol/L",  NA,       1,
  NA_character_,   2.2,    "mmol/L",  NA,       2,
  # ANRLO not missing
  "Calcium, Low",  1.52,   "mmol/L",  "4",      3,
  "Calcium, Low",  1.53,   "mmol/L",  "3",      4,
  "Calcium, Low",  1.74,   "mmol/L",  "3",      5,
  "Calcium, Low",  1.75,   "mmol/L",  "2",      6,
  "Calcium, Low",  1.94,   "mmol/L",  "2",      7,
  "Calcium, Low",  1.95,   "mmol/L",  "1",      8,
  "Calcium, Low",  2.09,   "mmol/L",  "1",      9,
  "Calcium, Low",  2.1,    "mmol/L",  "0",      10,
  # Unit missing cannot grade
  "Calcium, Low",  2.5,    NA,        NA,       11,
  # AVAL missing cannot grade
  "Calcium, Low",  NA,     "mmol/L",  NA,       12,
) %>%
  mutate(
    BRTHDT = lubridate::ymd("2023-01-01"),
    ADT = lubridate::ymd("2023-01-08")
  )

exp_calcdge7d_daids_cv <- tibble::tribble(
  ~ATOXDSCL,       ~AVAL,   ~AVALU,  ~ATOXGRL,  ~TESTNUM,
  "Not a term",    6.09,   "mg/dL",        NA,  1,
  NA_character_,   6.09,   "mg/dL",        NA,  2,
  # ANRLO not missing
  "Calcium, Low",  6.09,   "mg/dL",       "4",  3,
  "Calcium, Low",  6.1,    "mg/dL",       "3",  4,
  "Calcium, Low",  6.9,    "mg/dL",       "3",  5,
  "Calcium, Low",  7,      "mg/dL",       "2",  6,
  "Calcium, Low",  7.7,    "mg/dL",       "2",  7,
  "Calcium, Low",  7.8,    "mg/dL",       "1",  8,
  "Calcium, Low",  8.3,    "mg/dL",       "1",  9,
  "Calcium, Low",  8.4,    "mg/dL",       "0",  10,
  # Unit missing cannot grade
  "Calcium, Low",  8.4,         NA,        NA,  11,
  # AVAL missing cannot grade
  "Calcium, Low",  NA,     "mg/dL",        NA,  12,
) %>%
  mutate(
    BRTHDT = lubridate::ymd("2023-01-01"),
    ADT = lubridate::ymd("2023-01-08")
  )

exp_calcdlt7d_daids_si <- tibble::tribble(
  ~ATOXDSCL,      ~AVAL, ~AVALU,   ~ATOXGRL, ~TESTNUM,
  "Not a term",   2.2,   "mmol/L", NA,       13,
  NA_character_,  2.2,   "mmol/L", NA,       14,
  # ANRLO not missing
  "Calcium, Low", 1.37,  "mmol/L", "4",      15,
  "Calcium, Low", 1.38,  "mmol/L", "3",      16,
  "Calcium, Low", 1.49,  "mmol/L", "3",      17,
  "Calcium, Low", 1.5,   "mmol/L", "2",      18,
  "Calcium, Low", 1.62,  "mmol/L", "2",      19,
  "Calcium, Low", 1.63,  "mmol/L", "1",      20,
  "Calcium, Low", 1.87,  "mmol/L", "1",      21,
  "Calcium, Low", 1.88,  "mmol/L", "0",      22,
  # Unit missing cannot grade
  "Calcium, Low", 2.2,   NA,       NA,       23,
  # AVAL missing cannot grade
  "Calcium, Low", NA,    "mmol/L", NA,       24,
) %>%
  mutate(
    BRTHDT = lubridate::ymd("2023-01-01"),
    ADT = lubridate::ymd("2023-01-07")
  )

exp_calcdlt7d_daids_cv <- tibble::tribble(
  ~ATOXDSCL,       ~AVAL,   ~AVALU,  ~ATOXGRL,  ~TESTNUM,
  "Not a term",    5.49,   "mg/dL",        NA,  13,
  NA_character_,   5.49,   "mg/dL",        NA,  14,
  # ANRLO not missing
  "Calcium, Low",  5.49,   "mg/dL",       "4",  15,
  "Calcium, Low",  5.5,    "mg/dL",       "3",  16,
  "Calcium, Low",  5.99,   "mg/dL",       "3",  17,
  "Calcium, Low",  6,      "mg/dL",       "2",  18,
  "Calcium, Low",  6.49,   "mg/dL",       "2",  19,
  "Calcium, Low",  6.5,    "mg/dL",       "1",  20,
  "Calcium, Low",  7.49,   "mg/dL",       "1",  21,
  "Calcium, Low",  7.5,    "mg/dL",       "0",  22,
  # Unit missing cannot grade
  "Calcium, Low",  7.5,         NA,        NA,  23,
  # AVAL missing cannot grade
  "Calcium, Low",  NA,     "mg/dL",        NA,  24,
) %>%
  mutate(
    BRTHDT = lubridate::ymd("2023-01-01"),
    ADT = lubridate::ymd("2023-01-07")
  )

exp_calcdnoage_daids_si <- exp_calcdge7d_daids_si %>%
  filter(TESTNUM %in% c(9, 10)) %>%
  mutate(
    ADT = if_else(TESTNUM == 9, NA, ADT),
    BRTHDT = if_else(TESTNUM == 10, NA, BRTHDT),
    ATOXGRL = NA_character_,
    TESTNUM = if_else(TESTNUM == 9, 25, 26)
  )

exp_calcdnoage_daids_cv <- exp_calcdge7d_daids_cv %>%
  filter(TESTNUM %in% c(9, 10)) %>%
  mutate(
    ADT = if_else(TESTNUM == 9, NA, ADT),
    BRTHDT = if_else(TESTNUM == 10, NA, BRTHDT),
    ATOXGRL = NA_character_,
    TESTNUM = if_else(TESTNUM == 9, 25, 26)
  )

exp_calcd_daids_si <- exp_calcdnoage_daids_si %>%
  bind_rows(
    exp_calcdge7d_daids_si,
    exp_calcdlt7d_daids_si
  )

exp_calcd_daids_cv <- exp_calcdnoage_daids_cv %>%
  bind_rows(
    exp_calcdge7d_daids_cv,
    exp_calcdlt7d_daids_cv
  )

exp_caliond_daids_si <- tibble::tribble(
  ~ATOXDSCL,                ~AVAL, ~ANRLO, ~ANRHI, ~AVALU,   ~ATOXGRL, ~TESTNUM,
  "Not a term",             0.79,  1.1,    100,    "mmol/L", NA,       1,
  NA_character_,            0.79,  1.1,    100,    "mmol/L", NA,       2,
  # ANRLO not missing
  "Calcium (Ionized), Low", 0.79,  1.1,    100,    "mmol/L", "4",      3,
  "Calcium (Ionized), Low", 0.8,   1.1,    100,    "mmol/L", "3",      4,
  "Calcium (Ionized), Low", 0.89,  1.1,    100,    "mmol/L", "3",      5,
  "Calcium (Ionized), Low", 0.9,   1.1,    100,    "mmol/L", "2",      6,
  "Calcium (Ionized), Low", 0.99,  1.1,    100,    "mmol/L", "2",      7,
  "Calcium (Ionized), Low", 1,     1.1,    100,    "mmol/L", "1",      8,
  "Calcium (Ionized), Low", 1.09,  1.1,    100,    "mmol/L", "1",      9,
  "Calcium (Ionized), Low", 1.1,   1.1,    100,    "mmol/L", "0",      10,
  # ANRLO missing - can grade 2-4
  "Calcium (Ionized), Low", 0.79,  NA,     100,    "mmol/L", "4",      11,
  "Calcium (Ionized), Low", 0.8,   NA,     100,    "mmol/L", "3",      12,
  "Calcium (Ionized), Low", 0.89,  NA,     100,    "mmol/L", "3",      13,
  "Calcium (Ionized), Low", 0.9,   NA,     100,    "mmol/L", "2",      14,
  "Calcium (Ionized), Low", 0.99,  NA,     100,    "mmol/L", "2",      15,
  # ANRLO missing - can NOT grade 0 or 1
  "Calcium (Ionized), Low", 1,     NA,     100,    "mmol/L", NA,       16,
  "Calcium (Ionized), Low", 1.09,  NA,     100,    "mmol/L", NA,       17,
  "Calcium (Ionized), Low", 1.1,   NA,     100,    "mmol/L", NA,       18,
  # Unit missing cannot grade
  "Calcium (Ionized), Low", 1.1,   1.1,    100,    NA,       NA,       19,
  # AVAL missing cannot grade
  "Calcium (Ionized), Low", NA,    1.1,    100,    "mmol/L", NA,       20,
)

exp_caliond_daids_cv <- exp_caliond_daids_si %>%
  mutate(
    AVAL = AVAL * 4,
    ANRLO = ANRLO * 4,
    AVALU = if_else(is.na(AVALU), NA_character_, "mg/dL")
  )

exp_glucdge1m_daids_si <- tibble::tribble(
  ~ATOXDSCL,      ~AVAL,  ~AVALU,   ~ATOXGRL, ~TESTNUM,
  "Not a term",   9.5,    "mg/L",   NA,       1,
  NA_character_,  4.1,    "mmol/L", NA,       2,
  "Glucose, Low", 1.66,   "mmol/L", "4",      3,
  "Glucose, Low", 1.67,   "mmol/L", "3",      4,
  "Glucose, Low", 2.21,   "mmol/L", "3",      5,
  "Glucose, Low", 2.22,   "mmol/L", "2",      6,
  "Glucose, Low", 3.04,   "mmol/L", "2",      7,
  "Glucose, Low", 3.05,   "mmol/L", "1",      8,
  "Glucose, Low", 3.54,   "mmol/L", "1",      9,
  "Glucose, Low", 3.55,   "mmol/L", "0",      10,
  # AVALU missing cannot grade
  "Glucose, Low", 4,      NA,       NA,       11,
  # AVAL missing cannot grade
  "Glucose, Low", NA,     "mmol/L", NA,       12,
) %>%
  mutate(
    BRTHDT = lubridate::ymd("2022-11-30"),
    ADT = lubridate::ymd("2022-12-30"),
  )

exp_glucdge1m_daids_cv <- tibble::tribble(
  ~ATOXDSCL, ~AVAL, ~AVALU, ~ATOXGRL, ~TESTNUM,
  "Not a term", 29, "mg/L", NA, 1,
  NA_character_, 29, "mg/dL", NA, 2,
  "Glucose, Low", 29, "mg/dL", "4", 3,
  "Glucose, Low", 30, "mg/dL", "3", 4,
  "Glucose, Low", 39, "mg/dL", "3", 5,
  "Glucose, Low", 40, "mg/dL", "2", 6,
  "Glucose, Low", 54, "mg/dL", "2", 7,
  "Glucose, Low", 55, "mg/dL", "1", 8,
  "Glucose, Low", 64, "mg/dL", "1", 9,
  "Glucose, Low", 65, "mg/dL", "0", 10,
  # AVALU missing cannot grade
  "Glucose, Low", 65, NA, NA, 11,
  # AVAL missing cannot grade
  "Glucose, Low", NA, "mg/dL", NA, 12,
) %>%
  mutate(
    BRTHDT = lubridate::ymd("2022-11-30"),
    ADT = lubridate::ymd("2022-12-30"),
  )

exp_glucdlt1m_daids_si <- tibble::tribble(
  ~ATOXDSCL,      ~AVAL,  ~AVALU,   ~ATOXGRL, ~TESTNUM,
  "Not a term",   9.5,    "mg/L",   NA,       13,
  NA_character_,  4.1,    "mmol/L", NA,       14,
  "Glucose, Low", 1.66,   "mmol/L", "4",      15,
  "Glucose, Low", 1.67,   "mmol/L", "3",      16,
  "Glucose, Low", 2.21,   "mmol/L", "3",      17,
  "Glucose, Low", 2.22,   "mmol/L", "2",      18,
  "Glucose, Low", 2.77,   "mmol/L", "2",      19,
  "Glucose, Low", 2.78,   "mmol/L", "1",      20,
  "Glucose, Low", 2.99,   "mmol/L", "1",      21,
  "Glucose, Low", 3,      "mmol/L", "0",      22,
  # AVALU missing cannot grade
  "Glucose, Low", 4,      NA,       NA,       23,
  # AVAL missing cannot grade
  "Glucose, Low", NA,     "mmol/L", NA,       24,
) %>%
  mutate(
    BRTHDT = lubridate::ymd("2022-11-30"),
    ADT = lubridate::ymd("2022-12-29"),
  )

exp_glucdlt1m_daids_cv <- tibble::tribble(
  ~ATOXDSCL,      ~AVAL,  ~AVALU,   ~ATOXGRL, ~TESTNUM,
  "Not a term",   29,     "mg/L",   NA,       13,
  NA_character_,  29,     "mg/dL",  NA,       14,
  "Glucose, Low", 29,     "mg/dL",  "4",      15,
  "Glucose, Low", 30,     "mg/dL",  "3",      16,
  "Glucose, Low", 39,     "mg/dL",  "3",      17,
  "Glucose, Low", 40,     "mg/dL",  "2",      18,
  "Glucose, Low", 49,     "mg/dL",  "2",      19,
  "Glucose, Low", 50,     "mg/dL",  "1",      20,
  "Glucose, Low", 54,     "mg/dL",  "1",      21,
  "Glucose, Low", 55,     "mg/dL",  "0",      22,
  # AVALU missing cannot grade
  "Glucose, Low", 55,     NA,       NA,       23,
  # AVAL missing cannot grade
  "Glucose, Low", NA,     "mg/dL",  NA,       24,
) %>%
  mutate(
    BRTHDT = lubridate::ymd("2022-11-30"),
    ADT = lubridate::ymd("2022-12-29"),
  )

exp_glucdnoage_daids_si <- exp_glucdge1m_daids_si %>%
  filter(TESTNUM %in% c(9, 10)) %>%
  mutate(
    ATOXGRL = NA_character_,
    ADT = if_else(TESTNUM == 9, NA, ADT),
    BRTHDT = if_else(TESTNUM == 10, NA, BRTHDT),
    TESTNUM = if_else(TESTNUM == 9, 25, 26)
  )

exp_glucdnoage_daids_cv <- exp_glucdge1m_daids_cv %>%
  filter(TESTNUM %in% c(9, 10)) %>%
  mutate(
    ATOXGRL = NA_character_,
    ADT = if_else(TESTNUM == 9, NA, ADT),
    BRTHDT = if_else(TESTNUM == 10, NA, BRTHDT),
    TESTNUM = if_else(TESTNUM == 9, 25, 26)
  )

exp_glucd_daids_si <- exp_glucdnoage_daids_si %>%
  bind_rows(
    exp_glucdge1m_daids_si,
    exp_glucdlt1m_daids_si
  )

exp_glucd_daids_cv <- exp_glucdnoage_daids_cv %>%
  bind_rows(
    exp_glucdge1m_daids_cv,
    exp_glucdlt1m_daids_cv
  )

exp_cholfige18y_daids_si <- tibble::tribble(
  ~ATOXDSCH,                    ~AVAL,  ~AVALU,    ~ATOXGRH, ~TESTNUM,
  "Not a term",                 3.5,    "mmol/L",  NA,       1,
  NA_character_,                3.5,    "mmol/L",  NA,       2,
  "Cholesterol, Fasting, High", 7.77,   "mmol/L",  "3",      3,
  "Cholesterol, Fasting, High", 7.76,   "mmol/L",  "2",      4,
  "Cholesterol, Fasting, High", 6.19,   "mmol/L",  "2",      5,
  "Cholesterol, Fasting, High", 6.18,   "mmol/L",  "1",      6,
  "Cholesterol, Fasting, High", 5.18,   "mmol/L",  "1",      7,
  "Cholesterol, Fasting, High", 5.17,   "mmol/L",  "0",      8,
  # Unit missing cannot grade
  "Cholesterol, Fasting, High", 3.5,    NA,        NA,       9,
  # AVAL missing cannot grade
  "Cholesterol, Fasting, High", NA,     "mmol/L",  NA,       10,
) %>%
  mutate(
    BRTHDT = lubridate::ymd("2005-01-08"),
    ADT = lubridate::ymd("2023-01-08")
  )

exp_cholfige18y_daids_cv <- tibble::tribble(
  ~ATOXDSCH,                    ~AVAL,  ~AVALU,    ~ATOXGRH, ~TESTNUM,
  "Not a term",                 300,    "mg/dL",   NA,       1,
  NA_character_,                300,    "mg/dL",   NA,       2,
  "Cholesterol, Fasting, High", 300,    "mg/dL",   "3",      3,
  "Cholesterol, Fasting, High", 299,    "mg/dL",   "2",      4,
  "Cholesterol, Fasting, High", 240,    "mg/dL",   "2",      5,
  "Cholesterol, Fasting, High", 239,    "mg/dL",   "1",      6,
  "Cholesterol, Fasting, High", 200,    "mg/dL",   "1",      7,
  "Cholesterol, Fasting, High", 199,    "mg/dL",   "0",      8,
  # Unit missing cannot grade
  "Cholesterol, Fasting, High", 200,    NA,        NA,       9,
  # AVAL missing cannot grade
  "Cholesterol, Fasting, High", NA,     "mg/dL",   NA,       10,
) %>%
  mutate(
    BRTHDT = lubridate::ymd("2005-01-08"),
    ADT = lubridate::ymd("2023-01-08")
  )

exp_cholfilt18y_daids_si <- tibble::tribble(
  ~ATOXDSCH,                    ~AVAL,  ~AVALU,    ~ATOXGRH, ~TESTNUM,
  "Not a term",                 3.5,    "mmol/L",  NA,       13,
  NA_character_,                3.5,    "mmol/L",  NA,       14,
  "Cholesterol, Fasting, High", 7.77,   "mmol/L",  "3",      15,
  "Cholesterol, Fasting, High", 7.76,   "mmol/L",  "2",      16,
  "Cholesterol, Fasting, High", 5.15,   "mmol/L",  "2",      17,
  "Cholesterol, Fasting, High", 5.14,   "mmol/L",  "1",      18,
  "Cholesterol, Fasting, High", 4.4,    "mmol/L",  "1",      19,
  "Cholesterol, Fasting, High", 4.39,   "mmol/L",  "0",      20,
  # Unit missing cannot grade
  "Cholesterol, Fasting, High", 3.5,    NA,        NA,       21,
  # AVAL missing cannot grade
  "Cholesterol, Fasting, High", NA,     "mmol/L",  NA,       22,
) %>%
  mutate(
    BRTHDT = lubridate::ymd("2005-01-08"),
    ADT = lubridate::ymd("2023-01-07")
  )

exp_cholfilt18y_daids_cv <- tibble::tribble(
  ~ATOXDSCH,                    ~AVAL,  ~AVALU,   ~ATOXGRH, ~TESTNUM,
  "Not a term",                 300,    "mg/dL",  NA,       13,
  NA_character_,                300,    "mg/dL",  NA,       14,
  "Cholesterol, Fasting, High", 300,    "mg/dL",  "3",      15,
  "Cholesterol, Fasting, High", 299,    "mg/dL",  "2",      16,
  "Cholesterol, Fasting, High", 200,    "mg/dL",  "2",      17,
  "Cholesterol, Fasting, High", 199,    "mg/dL",  "1",      18,
  "Cholesterol, Fasting, High", 170,    "mg/dL",  "1",      19,
  "Cholesterol, Fasting, High", 169,    "mg/dL",  "0",      20,
  # Unit missing cannot grade
  "Cholesterol, Fasting, High", 169,    NA,       NA,       21,
  # AVAL missing cannot grade
  "Cholesterol, Fasting, High", NA,     "mg/dL",  NA,       22,
) %>%
  mutate(
    BRTHDT = lubridate::ymd("2005-01-08"),
    ADT = lubridate::ymd("2023-01-07")
  )

exp_cholfinoage_daids_si <- exp_cholfige18y_daids_si %>%
  filter(TESTNUM %in% c(9, 10)) %>%
  mutate(
    ADT = if_else(TESTNUM == 9, NA, ADT),
    BRTHDT = if_else(TESTNUM == 10, NA, BRTHDT),
    ATOXGRH = NA_character_,
    TESTNUM = if_else(TESTNUM == 9, 25, 26)
  )

exp_cholfinoage_daids_cv <- exp_cholfige18y_daids_cv %>%
  filter(TESTNUM %in% c(9, 10)) %>%
  mutate(
    ADT = if_else(TESTNUM == 9, NA, ADT),
    BRTHDT = if_else(TESTNUM == 10, NA, BRTHDT),
    ATOXGRH = NA_character_,
    TESTNUM = if_else(TESTNUM == 9, 25, 26)
  )

exp_cholfi_daids_si <- exp_cholfinoage_daids_si %>%
  bind_rows(
    exp_cholfige18y_daids_si,
    exp_cholfilt18y_daids_si
  )

exp_cholfi_daids_cv <- exp_cholfinoage_daids_cv %>%
  bind_rows(
    exp_cholfige18y_daids_cv,
    exp_cholfilt18y_daids_cv
  )

exp_ldlfige18y_daids_si <- tibble::tribble(
  ~ATOXDSCH,            ~AVAL,  ~AVALU,    ~ATOXGRH, ~TESTNUM,
  "Not a term",         3.1,    "mmol/L",  NA,       1,
  NA_character_,        3.1,    "mmol/L",  NA,       2,
  "LDL, Fasting, High", 4.9,    "mmol/L",  "3",      3,
  "LDL, Fasting, High", 4.89,   "mmol/L",  "2",      4,
  "LDL, Fasting, High", 4.12,   "mmol/L",  "2",      5,
  "LDL, Fasting, High", 4.11,   "mmol/L",  "1",      6,
  "LDL, Fasting, High", 3.37,   "mmol/L",  "1",      7,
  "LDL, Fasting, High", 3.36,   "mmol/L",  "0",      8,
  # Unit missing cannot grade
  "LDL, Fasting, High", 3.1,    NA,        NA,       9,
  # AVAL missing cannot grade
  "LDL, Fasting, High", NA,     "mmol/L",  NA,       10,
) %>%
  mutate(
    BRTHDT = lubridate::ymd("2005-01-08"),
    ADT = lubridate::ymd("2023-01-08")
  )

exp_ldlfige18y_daids_cv <- tibble::tribble(
  ~ATOXDSCH,             ~AVAL,   ~AVALU,  ~ATOXGRH,  ~TESTNUM,
  "Not a term",          190,    "mg/dL",        NA,  1,
  NA_character_,         190,    "mg/dL",        NA,  2,
  "LDL, Fasting, High",  190,    "mg/dL",       "3",  3,
  "LDL, Fasting, High",  189,    "mg/dL",       "2",  4,
  "LDL, Fasting, High",  160,    "mg/dL",       "2",  5,
  "LDL, Fasting, High",  159,    "mg/dL",       "1",  6,
  "LDL, Fasting, High",  130,    "mg/dL",       "1",  7,
  "LDL, Fasting, High",  129,    "mg/dL",       "0",  8,
  # Unit missing cannot grade
  "LDL, Fasting, High",  129,         NA,        NA,  9,
  # AVAL missing cannot grade
  "LDL, Fasting, High",  NA,     "mg/dL",        NA,  10,
) %>%
  mutate(
    BRTHDT = lubridate::ymd("2005-01-08"),
    ADT = lubridate::ymd("2023-01-08")
  )

exp_ldlfilt18y_daids_si <- tibble::tribble(
  ~ATOXDSCH,            ~AVAL,  ~AVALU,    ~ATOXGRH, ~TESTNUM,
  "Not a term",         2.8,    "mmol/L",  NA,       11,
  NA_character_,        2.8,    "mmol/L",  NA,       12,
  "LDL, Fasting, High", 4.9,    "mmol/L",  "3",      13,
  "LDL, Fasting, High", 4.89,   "mmol/L",  "2",      14,
  "LDL, Fasting, High", 3.34,   "mmol/L",  "2",      15,
  "LDL, Fasting, High", 3.33,   "mmol/L",  "1",      16,
  "LDL, Fasting, High", 2.85,   "mmol/L",  "1",      17,
  "LDL, Fasting, High", 2.84,   "mmol/L",  "0",      18,
  # Unit missing cannot grade
  "LDL, Fasting, High", 3.5,    NA,        NA,       19,
  # AVAL missing cannot grade
  "LDL, Fasting, High", NA,     "mmol/L",  NA,       20,
) %>%
  mutate(
    BRTHDT = lubridate::ymd("2020-01-07"),
    ADT = lubridate::ymd("2023-01-07")
  )

exp_ldlfilt18y_daids_cv <- tibble::tribble(
  ~ATOXDSCH,            ~AVAL,  ~AVALU,   ~ATOXGRH, ~TESTNUM,
  "Not a term",         190,    "mg/dL",  NA,       11,
  NA_character_,        190,    "mg/dL",  NA,       12,
  "LDL, Fasting, High", 190,    "mg/dL",  "3",      13,
  "LDL, Fasting, High", 189,    "mg/dL",  "2",      14,
  "LDL, Fasting, High", 130,    "mg/dL",  "2",      15,
  "LDL, Fasting, High", 129,    "mg/dL",  "1",      16,
  "LDL, Fasting, High", 110,    "mg/dL",  "1",      17,
  "LDL, Fasting, High", 109,    "mg/dL",  "0",      18,
  # Unit missing cannot grade
  "LDL, Fasting, High", 110,    NA,       NA,       19,
  # AVAL missing cannot grade
  "LDL, Fasting, High", NA,     "mg/dL",  NA,       20,
) %>%
  mutate(
    BRTHDT = lubridate::ymd("2020-01-07"),
    ADT = lubridate::ymd("2023-01-07")
  )

exp_ldlfinoage_daids_si <- exp_ldlfige18y_daids_si %>%
  filter(TESTNUM %in% c(7, 8)) %>%
  mutate(
    ADT = if_else(TESTNUM == 7, NA, ADT),
    BRTHDT = if_else(TESTNUM == 8, NA, BRTHDT),
    ATOXGRH = NA_character_,
    TESTNUM = if_else(TESTNUM == 7, 25, 26)
  )

exp_ldlfinoage_daids_cv <- exp_ldlfige18y_daids_cv %>%
  filter(TESTNUM %in% c(7, 8)) %>%
  mutate(
    ADT = if_else(TESTNUM == 7, NA, ADT),
    BRTHDT = if_else(TESTNUM == 8, NA, BRTHDT),
    ATOXGRH = NA_character_,
    TESTNUM = if_else(TESTNUM == 7, 25, 26)
  )

exp_ldlfile2y_daids_si <- exp_ldlfige18y_daids_si %>%
  filter(TESTNUM %in% c(7, 8)) %>%
  mutate(
    BRTHDT = if_else(TESTNUM == 7, lubridate::ymd("2021-01-07"), lubridate::ymd("2022-01-07")),
    ATOXGRH = NA_character_,
    TESTNUM = if_else(TESTNUM == 7, 27, 28)
  )

exp_ldlfile2y_daids_cv <- exp_ldlfige18y_daids_cv %>%
  filter(TESTNUM %in% c(7, 8)) %>%
  mutate(
    BRTHDT = if_else(TESTNUM == 7, lubridate::ymd("2021-01-07"), lubridate::ymd("2022-01-07")),
    ATOXGRH = NA_character_,
    TESTNUM = if_else(TESTNUM == 7, 27, 28)
  )

exp_ldlfi_daids_si <- exp_ldlfile2y_daids_si %>%
  bind_rows(
    exp_ldlfinoage_daids_si,
    exp_ldlfige18y_daids_si,
    exp_ldlfilt18y_daids_si
  )

exp_ldlfi_daids_cv <- exp_ldlfile2y_daids_cv %>%
  bind_rows(
    exp_ldlfinoage_daids_cv,
    exp_ldlfige18y_daids_cv,
    exp_ldlfilt18y_daids_cv
  )

exp_magd_daids_si <- tibble::tribble(
  ~ATOXDSCL,         ~AVAL,  ~AVALU,    ~ATOXGRL, ~TESTNUM,
  "Not a term",      0.5,    "mmol/L",  NA,       1,
  NA_character_,     0.5,    "mmol/L",  NA,       2,
  "Magnesium, Low",  0.29,   "mmol/L",  "4",      3,
  "Magnesium, Low",  0.3,    "mmol/L",  "3",      4,
  "Magnesium, Low",  0.44,   "mmol/L",  "3",      5,
  "Magnesium, Low",  0.45,   "mmol/L",  "2",      6,
  "Magnesium, Low",  0.59,   "mmol/L",  "2",      7,
  "Magnesium, Low",  0.6,    "mmol/L",  "1",      8,
  "Magnesium, Low",  0.69,   "mmol/L",  "1",      9,
  "Magnesium, Low",  0.7,    "mmol/L",  "0",      10,
  # Unit missing cannot grade
  "Magnesium, Low",  0.5,    NA,        NA,       11,
  # AVAL missing cannot grade
  "Magnesium, Low",  NA,     "mmol/L",  NA,       12,
)

exp_magd_daids_cv <- exp_magd_daids_si %>%
  mutate(
    AVAL = AVAL * 2.431,
    AVALU = if_else(is.na(AVALU), NA_character_, "mg/dL")
  )

exp_phosd_daids_gt14y_si <- tibble::tribble(
  ~AVAL, ~ANRLO, ~AVALU,   ~ATOXGRL, ~TESTNUM,
  0.9,   0.8,    "MM3",    NA,       1,
  0.31,  0.8,    "mmol/L", "4",      2,
  0.32,  0.8,    "mmol/L", "3",      3,
  0.44,  0.8,    "mmol/L", "3",      4,
  0.45,  0.8,    "mmol/L", "2",      5,
  0.64,  0.8,    "mmol/L", "2",      6,
  0.65,  0.8,    "mmol/L", "1",      7,
  0.79,  0.8,    "mmol/L", "1",      8,
  0.8,   0.8,    "mmol/L", "0",      9,
  # missing ANRLO - can grade 2 - 4
  0.31,  0.8,    "mmol/L", "4",      10,
  0.32,  0.8,    "mmol/L", "3",      11,
  0.44,  0.8,    "mmol/L", "3",      12,
  0.45,  0.8,    "mmol/L", "2",      13,
  0.64,  0.8,    "mmol/L", "2",      14,
  # missing ANRLO - can grade 0 - 1
  0.65,  0.8,    "mmol/L", "1",      15,
  0.79,  0.8,    "mmol/L", "1",      16,
  0.8,   0.8,    "mmol/L", "0",      17,
  # missing AVAL
  NA,    0.8,    "mmol/L", NA,       18,
  # missing UNIT
  1,     0.8,    NA,       NA,       19,
) %>%
  mutate(
    ATOXDSCL = "Phosphate, Low",
    BRTHDT = lubridate::ymd("2008-07-01"),
    ADT = lubridate::ymd("2023-07-01")
  )

exp_phosd_daids_gt14y_cv <- tibble::tribble(
  ~AVAL,  ~ANRLO,   ~AVALU,  ~ATOXGRL,  ~TESTNUM,
  0.9,    2.5,       "MM3",        NA,  1,
  0.9,    2.5,     "mg/dL",       "4",  2,
  1,      2.5,     "mg/dL",       "3",  3,
  1.3,    2.5,     "mg/dL",       "3",  4,
  1.4,    2.5,     "mg/dL",       "2",  5,
  1.9,    2.5,     "mg/dL",       "2",  6,
  2,      2.5,     "mg/dL",       "1",  7,
  2.49,   2.5,     "mg/dL",       "1",  8,
  2.5,    2.5,     "mg/dL",       "0",  9,
  # missing ANRLO - can grade 2 - 4
  0.9,    2.5,     "mg/dL",       "4",  10,
  1,      2.5,     "mg/dL",       "3",  11,
  1.3,    2.5,     "mg/dL",       "3",  12,
  1.4,    2.5,     "mg/dL",       "2",  13,
  1.9,    2.5,     "mg/dL",       "2",  14,
  # missing ANRLO - can grade 0 - 1
  2,      2.5,     "mg/dL",       "1",  15,
  2.49,   2.5,     "mg/dL",       "1",  16,
  2.5,    2.5,     "mg/dL",       "0",  17,
  # missing AVAL
  NA,     2.5,     "mg/dL",        NA,  18,
  # missing UNIT
  1,      2.5,          NA,        NA,  19,
) %>%
  mutate(
    ATOXDSCL = "Phosphate, Low",
    BRTHDT = lubridate::ymd("2008-07-01"),
    ADT = lubridate::ymd("2023-07-01")
  )

exp_phosd_daids_le14y_si <- tibble::tribble(
  ~AVAL,  ~AVALU,    ~ATOXGRL, ~TESTNUM,
  1.2,    "MM3",     NA,       20,
  0.47,   "mmol/L",  "4",      21,
  0.48,   "mmol/L",  "3",      22,
  0.8,    "mmol/L",  "3",      23,
  0.81,   "mmol/L",  "2",      24,
  0.96,   "mmol/L",  "2",      25,
  0.97,   "mmol/L",  "1",      26,
  1.12,   "mmol/L",  "1",      27,
  1.13,   "mmol/L",  "0",      28,
  NA,     "mmol/L",  NA,       29,
  1,      NA,        NA,       30,
) %>%
  mutate(
    ATOXDSCL = "Phosphate, Low",
    BRTHDT = lubridate::ymd("2022-07-01"),
    ADT = lubridate::ymd("2023-07-01")
  )

exp_phosd_daids_le14y_cv <- tibble::tribble(
  ~AVAL,  ~AVALU,    ~ATOXGRL, ~TESTNUM,
  1.49,   "MM3",     NA,       20,
  1.49,   "mg/dL",   "4",      21,
  1.5,    "mg/dL",   "3",      22,
  2.4,    "mg/dL",   "3",      23,
  2.5,    "mg/dL",   "2",      24,
  2.9,    "mg/dL",   "2",      25,
  3,      "mg/dL",   "1",      26,
  3.49,   "mg/dL",   "1",      27,
  3.5,    "mg/dL",   "0",      28,
  NA,     "mg/dL",   NA,       29,
  3.5,    NA,        NA,       30,
) %>%
  mutate(
    ATOXDSCL = "Phosphate, Low",
    BRTHDT = lubridate::ymd("2022-07-01"),
    ADT = lubridate::ymd("2023-07-01")
  )

exp_phosd_daids_lt1y_si <- tibble::tribble(
  ~AVAL,  ~AVALU,    ~ATOXGRL, ~TESTNUM,
  1.5,    "MM3",     NA,       31,
  0.47,   "mmol/L",  "4",      32,
  0.48,   "mmol/L",  "3",      33,
  0.8,    "mmol/L",  "3",      34,
  0.81,   "mmol/L",  "2",      35,
  1.12,   "mmol/L",  "2",      36,
  1.13,   "mmol/L",  "1",      37,
  1.44,   "mmol/L",  "1",      38,
  1.45,   "mmol/L",  "0",      39,
  NA,     "mmol/L",  NA,       40,
  1.5,    NA,        NA,       41,
) %>%
  mutate(
    ATOXDSCL = "Phosphate, Low",
    BRTHDT = lubridate::ymd("2023-07-01"),
    ADT = lubridate::ymd("2023-07-02")
  )

exp_phosd_daids_lt1y_cv <- tibble::tribble(
  ~AVAL,  ~AVALU,    ~ATOXGRL, ~TESTNUM,
  1.49,   "MM3",     NA,       31,
  1.49,   "mg/dL",   "4",      32,
  1.5,    "mg/dL",   "3",      33,
  2.4,    "mg/dL",   "3",      34,
  2.5,    "mg/dL",   "2",      35,
  3.4,    "mg/dL",   "2",      36,
  3.5,    "mg/dL",   "1",      37,
  4.49,   "mg/dL",   "1",      38,
  4.5,    "mg/dL",   "0",      39,
  NA,     "mg/dL",   NA,       40,
  4.5,    NA,        NA,       41,
) %>%
  mutate(
    ATOXDSCL = "Phosphate, Low",
    BRTHDT = lubridate::ymd("2023-07-01"),
    ADT = lubridate::ymd("2023-07-02")
  )

exp_phosd_daids_noage_si <- exp_phosd_daids_gt14y_si %>%
  filter(TESTNUM %in% c(8, 9)) %>%
  mutate(
    ADT = if_else(TESTNUM == 8, NA, ADT),
    BRTHDT = if_else(TESTNUM == 9, NA, BRTHDT),
    ATOXGRL = NA,
    TESTNUM = if_else(TESTNUM == 8, 42, 43)
  )

exp_phosd_daids_noage_cv <- exp_phosd_daids_gt14y_cv %>%
  filter(TESTNUM %in% c(8, 9)) %>%
  mutate(
    ADT = if_else(TESTNUM == 8, NA, ADT),
    BRTHDT = if_else(TESTNUM == 9, NA, BRTHDT),
    ATOXGRL = NA,
    TESTNUM = if_else(TESTNUM == 8, 42, 43)
  )

exp_phosd_daids_si <- exp_phosd_daids_gt14y_si %>%
  bind_rows(
    exp_phosd_daids_le14y_si,
    exp_phosd_daids_lt1y_si,
    exp_phosd_daids_noage_si
  )

exp_phosd_daids_cv <- exp_phosd_daids_gt14y_cv %>%
  bind_rows(
    exp_phosd_daids_le14y_cv,
    exp_phosd_daids_lt1y_cv,
    exp_phosd_daids_noage_cv
  )

exp_poti_daids <- tibble::tribble(
  ~ATOXDSCH,          ~AVAL,  ~AVALU,    ~ATOXGRH, ~TESTNUM,
  "Not a term",       5,      "mmol/L",  NA,       1,
  NA_character_,      5,      "mmol/L",  NA,       2,
  "Potassium, High",  7,      "mmol/L",  "4",      3,
  "Potassium, High",  6.9,    "mmol/L",  "3",      4,
  "Potassium, High",  6.5,    "mmol/L",  "3",      5,
  "Potassium, High",  6.4,    "mmol/L",  "2",      6,
  "Potassium, High",  6,      "mmol/L",  "2",      7,
  "Potassium, High",  5.9,    "mmol/L",  "1",      8,
  "Potassium, High",  5.6,    "mmol/L",  "1",      9,
  "Potassium, High",  5.5,    "mmol/L",  "0",      10,
  # Unit missing cannot grade
  "Potassium, High",  5,      NA,        NA,       11,
  # AVAL missing cannot grade
  "Potassium, High",  NA,     "mmol/L",  NA,       12,
)

exp_potd_daids <- tibble::tribble(
  ~ATOXDSCL,         ~AVAL,  ~AVALU,    ~ATOXGRL, ~TESTNUM,
  "Not a term",      3,      "mmol/L",  NA,       1,
  NA_character_,     3,      "mmol/L",  NA,       2,
  "Potassium, Low",  1.9,    "mmol/L",  "4",      3,
  "Potassium, Low",  2,      "mmol/L",  "3",      4,
  "Potassium, Low",  2.4,    "mmol/L",  "3",      5,
  "Potassium, Low",  2.5,    "mmol/L",  "2",      6,
  "Potassium, Low",  2.9,    "mmol/L",  "2",      7,
  "Potassium, Low",  3,      "mmol/L",  "1",      8,
  "Potassium, Low",  3.3,    "mmol/L",  "1",      9,
  "Potassium, Low",  3.4,    "mmol/L",  "0",      10,
  # Unit missing cannot grade
  "Potassium, Low",  3,      NA,        NA,       11,
  # AVAL missing cannot grade
  "Potassium, Low",  NA,     "mmol/L",  NA,       12,
)

exp_sodi_daids <- tibble::tribble(
  ~ATOXDSCH,       ~AVAL,  ~AVALU,    ~ATOXGRH, ~TESTNUM,
  "Not a term",    146,    "mmol/L",  NA,       1,
  NA_character_,   146,    "mmol/L",  NA,       2,
  "Sodium, High",  160,    "mmol/L",  "4",      3,
  "Sodium, High",  159,    "mmol/L",  "3",      4,
  "Sodium, High",  154,    "mmol/L",  "3",      5,
  "Sodium, High",  153,    "mmol/L",  "2",      6,
  "Sodium, High",  150,    "mmol/L",  "2",      7,
  "Sodium, High",  149,    "mmol/L",  "1",      8,
  "Sodium, High",  146,    "mmol/L",  "1",      9,
  "Sodium, High",  145,    "mmol/L",  "0",      10,
  # Unit missing cannot grade
  "Sodium, High",  140,    NA,        NA,       11,
  # AVAL missing cannot grade
  "Sodium, High",  NA,     "mmol/L",  NA,       12,
)

exp_sodd_daids <- tibble::tribble(
  ~ATOXDSCL,      ~AVAL,  ~AVALU,    ~ATOXGRL, ~TESTNUM,
  "Not a term",   119,    "mmol/L",  NA,       1,
  NA_character_,  119,    "mmol/L",  NA,       2,
  "Sodium, Low",  120,    "mmol/L",  "4",      3,
  "Sodium, Low",  121,    "mmol/L",  "3",      4,
  "Sodium, Low",  124,    "mmol/L",  "3",      5,
  "Sodium, Low",  125,    "mmol/L",  "2",      6,
  "Sodium, Low",  129,    "mmol/L",  "2",      7,
  "Sodium, Low",  130,    "mmol/L",  "1",      8,
  "Sodium, Low",  134,    "mmol/L",  "1",      9,
  "Sodium, Low",  135,    "mmol/L",  "0",      10,
  # Unit missing cannot grade
  "Sodium, Low",  140,    NA,        NA,       11,
  # AVAL missing cannot grade
  "Sodium, Low",  NA,     "mmol/L",  NA,       12,
)

exp_cd4d_daids_gt5y <- tibble::tribble(
  ~AVAL,  ~AVALU,    ~ATOXGRL, ~TESTNUM,
  0.1,    "MM3",     NA,       1,
  0.09,   "10^9/L",  "4",      2,
  0.1,    "10^9/L",  "3",      3,
  0.19,   "10^9/L",  "3",      4,
  0.2,    "10^9/L",  "2",      5,
  0.29,   "10^9/L",  "2",      6,
  0.3,    "10^9/L",  "1",      7,
  0.39,   "10^9/L",  "1",      8,
  0.4,    "10^9/L",  "0",      9,
  NA,     "10^9/L",  NA,       10,
  1,      NA,        NA,       11,
) %>%
  mutate(
    ATOXDSCL = "Absolute CD4+ Count, Low",
    BRTHDT = lubridate::ymd("2023-07-01"),
    ADT = lubridate::ymd("2029-07-01")
  )

exp_cd4d_daids_le5y <- exp_cd4d_daids_gt5y %>%
  mutate(
    ADT = lubridate::ymd("2028-07-01"),
    ATOXGRL = NA_character_,
    TESTNUM = TESTNUM + 11
  )

exp_cd4d_daids_noage <- exp_cd4d_daids_gt5y %>%
  filter(TESTNUM %in% c(5, 6)) %>%
  mutate(
    ADT = if_else(TESTNUM == 5, NA, ADT),
    BRTHDT = if_else(TESTNUM == 6, NA, BRTHDT),
    ATOXGRL = NA_character_,
    TESTNUM = if_else(TESTNUM == 5, 23, 24)
  )

exp_cd4d_daids_si <- exp_cd4d_daids_gt5y %>%
  bind_rows(
    exp_cd4d_daids_le5y,
    exp_cd4d_daids_noage
  )

exp_cd4d_daids_cv <- exp_cd4d_daids_si %>%
  mutate(
    AVAL = AVAL * 1000,
    AVALU = if_else(str_to_upper(AVALU) == "10^9/L", "1/uL", AVALU)
  )

exp_lymphd_daids_gt5y <- tibble::tribble(
  ~AVAL,  ~AVALU,    ~ATOXGRL, ~TESTNUM,
  0.35,   "MM3",     NA,       1,
  0.34,   "10^9/L",  "4",      2,
  0.35,   "10^9/L",  "3",      3,
  0.49,   "10^9/L",  "3",      4,
  0.5,    "10^9/L",  "2",      5,
  0.59,   "10^9/L",  "2",      6,
  0.6,    "10^9/L",  "1",      7,
  0.64,   "10^9/L",  "1",      8,
  0.65,   "10^9/L",  "0",      9,
  0.66,   "10^9/L",  "0",      10,
  NA,     "10^9/L",  NA,       11,
  1,      NA,        NA,       12,
) %>%
  mutate(
    ATOXDSCL = "Absolute Lymphocyte Count, Low",
    BRTHDT = lubridate::ymd("2023-07-01"),
    ADT = lubridate::ymd("2029-07-01")
  )

exp_lymphd_daids_le5y <- exp_lymphd_daids_gt5y %>%
  mutate(
    ADT = lubridate::ymd("2028-07-01"),
    ATOXGRL = NA_character_,
    TESTNUM = TESTNUM + 12
  )

exp_lymphd_daids_noage <- exp_lymphd_daids_gt5y %>%
  filter(TESTNUM %in% c(5, 6)) %>%
  mutate(
    ADT = if_else(TESTNUM == 5, NA, ADT),
    BRTHDT = if_else(TESTNUM == 6, NA, BRTHDT),
    ATOXGRL = NA_character_,
    TESTNUM = if_else(TESTNUM == 5, 25, 26)
  )

exp_lymphd_daids_si <- exp_lymphd_daids_gt5y %>%
  bind_rows(
    exp_lymphd_daids_le5y,
    exp_lymphd_daids_noage
  )

exp_lymphd_daids_cv <- exp_lymphd_daids_si %>%
  mutate(
    AVALU = if_else(str_to_upper(AVALU) == "10^9/L", "10^3/uL", AVALU)
  )

exp_lymphd_daids_cv2 <- exp_lymphd_daids_si %>%
  mutate(
    AVAL = AVAL * 1000,
    AVALU = if_else(str_to_upper(AVALU) == "10^9/L", "10^3/mL", AVALU)
  )

exp_ancd_daids_gt7d <- tibble::tribble(
  ~AVAL,  ~AVALU,    ~ATOXGRL, ~TESTNUM,
  0.3,    "MM3",     NA,       1,
  0.399,  "10^9/L",  "4",      2,
  0.4,    "10^9/L",  "3",      3,
  0.599,  "10^9/L",  "3",      4,
  0.6,    "10^9/L",  "2",      5,
  0.799,  "10^9/L",  "2",      6,
  0.8,    "10^9/L",  "1",      7,
  1,      "10^9/L",  "1",      8,
  1.01,   "10^9/L",  "0",      9,
  NA,     "10^9/L",  NA,       10,
  1,      NA,        NA,       11,
) %>%
  mutate(
    ATOXDSCL = "Absolute Neutrophil Count (ANC), Low",
    BRTHDT = lubridate::ymd("2023-07-01"),
    ADT = lubridate::ymd("2023-07-09")
  )

exp_ancd_daids_ge2d <- tibble::tribble(
  ~AVAL,  ~AVALU,    ~ATOXGRL, ~TESTNUM,
  0.7,    "MM3",     NA,       12,
  0.749,  "10^9/L",  "4",      13,
  0.75,   "10^9/L",  "3",      14,
  0.999,  "10^9/L",  "3",      15,
  1,      "10^9/L",  "2",      16,
  1.249,  "10^9/L",  "2",      17,
  1.25,   "10^9/L",  "1",      18,
  1.5,    "10^9/L",  "1",      19,
  1.51,   "10^9/L",  "0",      20,
  NA,     "10^9/L",  NA,       21,
  1,      NA,        NA,       22,
) %>%
  mutate(
    ATOXDSCL = "Absolute Neutrophil Count (ANC), Low",
    BRTHDT = lubridate::ymd("2023-07-01"),
    ADT = lubridate::ymd("2023-07-03")
  )

exp_ancd_daids_le1d <- tibble::tribble(
  ~AVAL,  ~AVALU,    ~ATOXGRL, ~TESTNUM,
  1.5,    "MM3",     NA,       23,
  1.499,  "10^9/L",  "4",      24,
  1.5,    "10^9/L",  "3",      25,
  2.999,  "10^9/L",  "3",      26,
  3,      "10^9/L",  "2",      27,
  3.999,  "10^9/L",  "2",      28,
  4,      "10^9/L",  "1",      29,
  5,      "10^9/L",  "1",      30,
  5.01,   "10^9/L",  "0",      31,
  NA,     "10^9/L",  NA,       32,
  5,      NA,        NA,       33,
) %>%
  mutate(
    ATOXDSCL = "Absolute Neutrophil Count (ANC), Low",
    BRTHDT = lubridate::ymd("2023-07-01"),
    ADT = lubridate::ymd("2023-07-02")
  )

exp_ancd_daids <- exp_ancd_daids_gt7d %>%
  bind_rows(
    exp_ancd_daids_ge2d,
    exp_ancd_daids_le1d
  )

exp_ancd_daids_noage <- exp_ancd_daids %>%
  filter(TESTNUM %in% c(20, 31)) %>%
  mutate(
    ADT = if_else(TESTNUM == 20, NA, ADT),
    BRTHDT = if_else(TESTNUM == 31, NA, BRTHDT),
    ATOXGRL = NA,
    TESTNUM = case_when(
      TESTNUM == 20 ~ 34,
      TESTNUM == 31 ~ 35
    )
  )

exp_ancd_daids_si <- exp_ancd_daids %>%
  bind_rows(exp_ancd_daids_noage)

exp_ancd_daids_cv <- exp_ancd_daids_si %>%
  mutate(
    AVALU = if_else(str_to_upper(AVALU) == "10^9/L", "10^3/uL", AVALU)
  )

exp_ancd_daids_cv2 <- exp_ancd_daids_si %>%
  mutate(
    AVAL = AVAL * 1000,
    AVALU = if_else(str_to_upper(AVALU) == "10^9/L", "10^3/mL", AVALU)
  )

exp_fibd_daids_si <- tibble::tribble(
  ~ATOXDSCL,              ~AVAL, ~ANRLO, ~AVALU, ~ATOXGRL, ~TESTNUM,
  "Not a term",           2,     1,      "g/L",  NA,       1,
  NA_character_,          2,     1,      "g/L",  NA,       2,
  "Fibrinogen Decreased", 2,     1,      "g/dL", NA,       3,
  # test first half of criteria
  "Fibrinogen Decreased", 0.49,  1,      "g/L",  "4",      4,
  "Fibrinogen Decreased", 0.5,   1,      "g/L",  "3",      5,
  "Fibrinogen Decreased", 0.74,  1,      "g/L",  "3",      6,
  "Fibrinogen Decreased", 0.75,  1,      "g/L",  "2",      7,
  "Fibrinogen Decreased", 0.99,  1,      "g/L",  "2",      8,
  "Fibrinogen Decreased", 1,     1,      "g/L",  "1",      9,
  "Fibrinogen Decreased", 1.99,  1,      "g/L",  "1",      10,
  "Fibrinogen Decreased", 2,     1,      "g/L",  "0",      11,
  # test second half of criteria
  "Fibrinogen Decreased", 0.74,  3,      "g/L",  "4",      12,
  "Fibrinogen Decreased", 0.75,  3,      "g/L",  "3",      13,
  "Fibrinogen Decreased", 1.49,  3,      "g/L",  "3",      14,
  "Fibrinogen Decreased", 1.5,   3,      "g/L",  "2",      15,
  "Fibrinogen Decreased", 2.24,  3,      "g/L",  "2",      16,
  "Fibrinogen Decreased", 2.25,  3,      "g/L",  "1",      17,
  "Fibrinogen Decreased", 2.99,  3,      "g/L",  "1",      18,
  "Fibrinogen Decreased", 3,     3,      "g/L",  "0",      19,
  # TEST for missing values
  "Fibrinogen Decreased", 0.49,  NA,     "g/L",  "4",      20,
  "Fibrinogen Decreased", 0.5,   NA,     "g/L",  NA,       21,
  "Fibrinogen Decreased", 2,     1,      NA,     NA,       22,
  "Fibrinogen Decreased", NA,    1,      "g/L",  NA,       23,
)

exp_fibd_daids_cv <- exp_fibd_daids_si %>%
  mutate(
    AVAL = 100 * AVAL,
    ANRLO = 100 * ANRLO,
    AVALU = if_else(str_to_upper(AVALU) == "G/L", "mg/dL", NA_character_)
  )

exp_hgbd_daids_ge13ym <- tibble::tribble(
  ~ATOXDSCL,         ~AVAL, ~AVALU, ~ATOXGRL, ~SEX, ~TESTNUM,
  "Hemoglobin, Low", 69,    "MM3",  NA,       "M",  1,
  "Hemoglobin, Low", 69,    "g/L",  "4",      "M",  2,
  "Hemoglobin, Low", 70,    "g/L",  "3",      "M",  3,
  "Hemoglobin, Low", 89,    "g/L",  "3",      "M",  4,
  "Hemoglobin, Low", 90,    "g/L",  "2",      "M",  5,
  "Hemoglobin, Low", 99,    "g/L",  "2",      "M",  6,
  "Hemoglobin, Low", 100,   "g/L",  "1",      "M",  7,
  "Hemoglobin, Low", 109,   "g/L",  "1",      "M",  8,
  "Hemoglobin, Low", 110,   "g/L",  "0",      "M",  9,
  "Hemoglobin, Low", NA,    "g/L",  NA,       "M",  10,
  "Hemoglobin, Low", 110,   NA,     NA,       "M",  11,
  "Hemoglobin, Low", 110,   "g/L",  NA,       NA,   12,
) %>%
  mutate(
    BRTHDT = lubridate::ymd("2010-07-01"),
    ADT = lubridate::ymd("2023-07-01")
  )

exp_hgbd_daids_ge13yf <- tibble::tribble(
  ~ATOXDSCL,         ~AVAL, ~AVALU, ~ATOXGRL, ~SEX, ~TESTNUM,
  "Hemoglobin, Low", 64,    "MM3",  NA,       "F",  13,
  "Hemoglobin, Low", 64,    "g/L",  "4",      "F",  14,
  "Hemoglobin, Low", 65,    "g/L",  "3",      "F",  15,
  "Hemoglobin, Low", 84,    "g/L",  "3",      "F",  16,
  "Hemoglobin, Low", 85,    "g/L",  "2",      "F",  17,
  "Hemoglobin, Low", 94,    "g/L",  "2",      "F",  18,
  "Hemoglobin, Low", 95,    "g/L",  "1",      "F",  19,
  "Hemoglobin, Low", 104,   "g/L",  "1",      "F",  20,
  "Hemoglobin, Low", 105,   "g/L",  "0",      "F",  21,
  "Hemoglobin, Low", NA,    "g/L",  NA,       "F",  22,
  "Hemoglobin, Low", 110,   NA,     NA,       "F",  23,
  "Hemoglobin, Low", 110,   "g/L",  NA,       NA,   24,
) %>%
  mutate(
    BRTHDT = lubridate::ymd("2010-07-01"),
    ADT = lubridate::ymd("2023-07-01")
  )

exp_hgbd_daids_lt13y <- tibble::tribble(
  ~ATOXDSCL,         ~AVAL, ~AVALU, ~ATOXGRL, ~TESTNUM,
  "Hemoglobin, Low", 64,    "MM3",  NA,       25,
  "Hemoglobin, Low", 64,    "g/L",  "4",      26,
  "Hemoglobin, Low", 65,    "g/L",  "3",      27,
  "Hemoglobin, Low", 84,    "g/L",  "3",      28,
  "Hemoglobin, Low", 85,    "g/L",  "2",      29,
  "Hemoglobin, Low", 94,    "g/L",  "2",      30,
  "Hemoglobin, Low", 95,    "g/L",  "1",      31,
  "Hemoglobin, Low", 104,   "g/L",  "1",      32,
  "Hemoglobin, Low", 105,   "g/L",  "0",      33,
  "Hemoglobin, Low", NA,    "g/L",  NA,       34,
  "Hemoglobin, Low", 110,   NA,     NA,       35,
) %>%
  mutate(
    BRTHDT = lubridate::ymd("2010-07-01"),
    ADT = lubridate::ymd("2023-06-30")
  )

exp_hgbd_daids_le56d <- tibble::tribble(
  ~ATOXDSCL,         ~AVAL, ~AVALU, ~ATOXGRL, ~TESTNUM,
  "Hemoglobin, Low", 59,    "MM3",  NA,       36,
  "Hemoglobin, Low", 59,    "g/L",  "4",      37,
  "Hemoglobin, Low", 60,    "g/L",  "3",      38,
  "Hemoglobin, Low", 69,    "g/L",  "3",      39,
  "Hemoglobin, Low", 70,    "g/L",  "2",      40,
  "Hemoglobin, Low", 84,    "g/L",  "2",      41,
  "Hemoglobin, Low", 85,    "g/L",  "1",      42,
  "Hemoglobin, Low", 96,    "g/L",  "1",      43,
  "Hemoglobin, Low", 97,    "g/L",  "0",      44,
  "Hemoglobin, Low", NA,    "g/L",  NA,       45,
  "Hemoglobin, Low", 110,   NA,     NA,       46,
) %>%
  mutate(
    BRTHDT = lubridate::ymd("2023-07-01"),
    ADT = lubridate::ymd("2023-08-26")
  )

exp_hgbd_daids_le35d <- tibble::tribble(
  ~ATOXDSCL,         ~AVAL, ~AVALU, ~ATOXGRL, ~TESTNUM,
  "Hemoglobin, Low", 66,    "MM3",  NA,       47,
  "Hemoglobin, Low", 66,    "g/L",  "4",      48,
  "Hemoglobin, Low", 67,    "g/L",  "3",      49,
  "Hemoglobin, Low", 79,    "g/L",  "3",      50,
  "Hemoglobin, Low", 80,    "g/L",  "2",      51,
  "Hemoglobin, Low", 94,    "g/L",  "2",      52,
  "Hemoglobin, Low", 95,    "g/L",  "1",      53,
  "Hemoglobin, Low", 110,   "g/L",  "1",      54,
  "Hemoglobin, Low", 111,   "g/L",  "0",      55,
  "Hemoglobin, Low", NA,    "g/L",  NA,       56,
  "Hemoglobin, Low", 110,   NA,     NA,       57,
) %>%
  mutate(
    BRTHDT = lubridate::ymd("2023-07-01"),
    ADT = lubridate::ymd("2023-08-05")
  )

exp_hgbd_daids_le21d <- tibble::tribble(
  ~ATOXDSCL,         ~AVAL, ~AVALU, ~ATOXGRL, ~TESTNUM,
  "Hemoglobin, Low", 79,    "MM3",  NA,       58,
  "Hemoglobin, Low", 79,    "g/L",  "4",      59,
  "Hemoglobin, Low", 80,    "g/L",  "3",      60,
  "Hemoglobin, Low", 89,    "g/L",  "3",      61,
  "Hemoglobin, Low", 90,    "g/L",  "2",      62,
  "Hemoglobin, Low", 109,   "g/L",  "2",      63,
  "Hemoglobin, Low", 110,   "g/L",  "1",      64,
  "Hemoglobin, Low", 130,   "g/L",  "1",      65,
  "Hemoglobin, Low", 131,   "g/L",  "0",      66,
  "Hemoglobin, Low", NA,    "g/L",  NA,       67,
  "Hemoglobin, Low", 110,   NA,     NA,       68,
) %>%
  mutate(
    BRTHDT = lubridate::ymd("2023-07-01"),
    ADT = lubridate::ymd("2023-07-22")
  )

exp_hgbd_daids_le7d <- tibble::tribble(
  ~ATOXDSCL,         ~AVAL, ~AVALU, ~ATOXGRL, ~TESTNUM,
  "Hemoglobin, Low", 89,    "MM3",  NA,       69,
  "Hemoglobin, Low", 89,    "g/L",  "4",      70,
  "Hemoglobin, Low", 90,    "g/L",  "3",      71,
  "Hemoglobin, Low", 99,    "g/L",  "3",      72,
  "Hemoglobin, Low", 100,   "g/L",  "2",      73,
  "Hemoglobin, Low", 129,   "g/L",  "2",      74,
  "Hemoglobin, Low", 130,   "g/L",  "1",      75,
  "Hemoglobin, Low", 140,   "g/L",  "1",      76,
  "Hemoglobin, Low", 141,   "g/L",  "0",      77,
  "Hemoglobin, Low", NA,    "g/L",  NA,       78,
  "Hemoglobin, Low", 110,   NA,     NA,       79,
) %>%
  mutate(
    BRTHDT = lubridate::ymd("2023-07-01"),
    ADT = lubridate::ymd("2023-07-08")
  )

exp_hgbd_daids <- exp_hgbd_daids_ge13ym %>%
  bind_rows(
    exp_hgbd_daids_ge13yf,
    exp_hgbd_daids_lt13y,
    exp_hgbd_daids_le56d,
    exp_hgbd_daids_le35d,
    exp_hgbd_daids_le21d,
    exp_hgbd_daids_le7d
  )

exp_hgbd_daids_noage <- exp_hgbd_daids %>%
  filter(TESTNUM %in% c(5, 17, 29)) %>%
  mutate(
    ADT = NA,
    ATOXGRL = NA,
    TESTNUM = case_when(
      TESTNUM == 5 ~ 80,
      TESTNUM == 17 ~ 81,
      TESTNUM == 29 ~ 82
    )
  )

exp_hgbd_daids <- exp_hgbd_daids %>%
  bind_rows(exp_hgbd_daids_noage)

exp_methi_daids <- tibble::tribble(
  ~ATOXDSCH,         ~AVAL, ~AVALU,  ~ATOXGRH, ~TESTNUM,
  "Not a term",      20,    "%",     NA,       1,
  NA_character_,     20,    "%",     NA,       2,
  "Methemoglobin",   20,    "%",     "4",      3,
  "Methemoglobin",   19,    "%",     "3",      4,
  "Methemoglobin",   15,    "%",     "3",      5,
  "Methemoglobin",   14,    "%",     "2",      6,
  "Methemoglobin",   10,    "%",     "2",      7,
  "Methemoglobin",   9,     "%",     "1",      8,
  "Methemoglobin",   5,     "%",     "1",      9,
  "Methemoglobin",   4,     "%",     "0",      10,
  # Unit wrong - cannot grade
  "Methemoglobin",   100,   NA,      NA,       11,
  # AVAL missing cannot grade
  "Methemoglobin",   NA,    "%",     NA,       12,
)

exp_plated_daids_si <- tibble::tribble(
  ~ATOXDSCL,              ~AVAL, ~AVALU,    ~ATOXGRL, ~TESTNUM,
  "Not a term",           126,   "10^9/L",  NA,       1,
  NA_character_,          120,   "10^9/L",  NA,       2,
  "Platelets, Decreased", 115,   "MM3",     NA,       3,
  "Platelets, Decreased", 24.9,  "10^9/L",  "4",      4,
  "Platelets, Decreased", 25,    "10^9/L",  "3",      5,
  "Platelets, Decreased", 49.9,  "10^9/L",  "3",      6,
  "Platelets, Decreased", 50,    "10^9/L",  "2",      7,
  "Platelets, Decreased", 99.9,  "10^9/L",  "2",      8,
  "Platelets, Decreased", 100,   "10^9/L",  "1",      9,
  "Platelets, Decreased", 124.9, "10^9/L",  "1",      10,
  "Platelets, Decreased", 125,   "10^9/L",  "0",      11,
)

exp_plated_daids_cv <- exp_plated_daids_si %>%
  mutate(
    AVALU = if_else(str_to_upper(AVALU) == "10^9/L", "10^3/uL", AVALU)
  )

exp_plated_daids_cv2 <- exp_plated_daids_si %>%
  mutate(
    AVAL = AVAL * 1000,
    AVALU = if_else(str_to_upper(AVALU) == "10^9/L", "10^3/mL", AVALU)
  )

exp_wbcd_daids_gt7d <- tibble::tribble(
  ~ATOXDSCL,        ~AVAL, ~AVALU,    ~ATOXGRL, ~TESTNUM,
  "Not a term",     1,     "10^9/L",  NA,       1,
  NA_character_,    2,     "10^9/L",  NA,       2,
  "WBC, Decreased", 0.9,   "MM3",     NA,       3,
  "WBC, Decreased", 0.9,   "10^9/L",  "4",      4,
  "WBC, Decreased", 1,     "10^9/L",  "3",      5,
  "WBC, Decreased", 1.49,  "10^9/L",  "3",      6,
  "WBC, Decreased", 1.5,   "10^9/L",  "2",      7,
  "WBC, Decreased", 1.99,  "10^9/L",  "2",      8,
  "WBC, Decreased", 2,     "10^9/L",  "1",      9,
  "WBC, Decreased", 2.49,  "10^9/L",  "1",      10,
  "WBC, Decreased", 2.5,   "10^9/L",  "0",      11,
) %>%
  mutate(
    BRTHDT = lubridate::ymd("2023-07-01"),
    ADT = lubridate::ymd("2023-07-09")
  )

exp_wbcd_daids_le7d <- tibble::tribble(
  ~ATOXDSCL,        ~AVAL, ~AVALU,   ~ATOXGRL, ~TESTNUM,
  "WBC, Decreased", 2.49,  "MM3",    NA,       12,
  "WBC, Decreased", 2.49,  "10^9/L", "4",      13,
  "WBC, Decreased", 2.5,   "10^9/L", "3",      14,
  "WBC, Decreased", 3.99,  "10^9/L", "3",      15,
  "WBC, Decreased", 4,     "10^9/L", "2",      16,
  "WBC, Decreased", 5.49,  "10^9/L", "2",      17,
  "WBC, Decreased", 5.5,   "10^9/L", "1",      18,
  "WBC, Decreased", 6.99,  "10^9/L", "1",      19,
  "WBC, Decreased", 7,     "10^9/L", "0",      20,
) %>%
  mutate(
    BRTHDT = lubridate::ymd("2023-07-01"),
    ADT = lubridate::ymd("2023-07-08")
  )

exp_wbcd_daids_noage <- exp_wbcd_daids_gt7d %>%
  filter(TESTNUM %in% c(10, 11)) %>%
  mutate(
    BRTHDT = if_else(TESTNUM == 10, NA, BRTHDT),
    ADT = if_else(TESTNUM == 11, NA, ADT),
    ATOXGRL = NA_character_,
    TESTNUM = TESTNUM + 11
  )

exp_wbcd_daids_si <- exp_wbcd_daids_gt7d %>%
  bind_rows(
    exp_wbcd_daids_le7d,
    exp_wbcd_daids_noage
  )

exp_wbcd_daids_cv <- exp_wbcd_daids_si %>%
  mutate(
    AVALU = if_else(str_to_upper(AVALU) == "10^9/L", "10^3/uL", AVALU)
  )

exp_wbcd_daids_cv2 <- exp_wbcd_daids_si %>%
  mutate(
    AVAL = AVAL * 1000,
    AVALU = if_else(str_to_upper(AVALU) == "10^9/L", "10^3/mL", AVALU)
  )
