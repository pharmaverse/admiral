
# derive_var_atoxgr ----

test_that("derive_var_atoxgr: Test 1 ATOXGR cannot be graded", {
  exp_out_1 <- tibble::tribble(
    ~ATOXDSCL,          ~ATOXDSCH,        ~ATOXGRL,      ~ATOXGRH,       ~ATOXGR,
    NA_character_,      NA_character_,    NA_character_, NA_character_,  NA_character_,
    "Hypophosphatemia", NA_character_,    NA_character_, NA_character_,  NA_character_,
    NA_character_,      "Hyperglycemia",  NA_character_, NA_character_,  NA_character_,
    "Hypoglycemia",     "Hyperglycemia",  NA_character_, NA_character_,  NA_character_,
    # ATOXGRL is ungradable so cannot say ATOXGR is normal
    "Hypoglycemia",     "Hyperglycemia",  NA_character_, "0",            NA_character_,
    # ATOXGRH is ungradable so cannot say ATOXGR is normal
    "Hypoglycemia",     "Hyperglycemia",  "0",           NA_character_,  NA_character_,
  )

  input_1 <- exp_out_1 %>%
    select(-ATOXGR)

  expect_equal(
    derive_var_atoxgr(
      input_1,
      lotox_description_var = ATOXDSCL,
      hitox_description_var = ATOXDSCH
    ),
    exp_out_1
  )
})

test_that("derive_var_atoxgr: Test 2 ATOXGR = 0 (normal)", {
  exp_out_2 <- tibble::tribble(
    ~ATOXDSCL,          ~ATOXDSCH,        ~ATOXGRL,      ~ATOXGRH,       ~ATOXGR,
    "Hypoglycemia",     "Hyperglycemia",  "0",           "0",            "0",
    NA_character_,      "INR Increased",  NA_character_, "0",            "0",
    "Hypophosphatemia", NA_character_,    "0",           NA_character_,  "0",
  )

  input_2 <- exp_out_2 %>%
    select(-ATOXGR)

  expect_equal(
    derive_var_atoxgr(
      input_2,
      lotox_description_var = ATOXDSCL,
      hitox_description_var = ATOXDSCH
    ),
    exp_out_2
  )
})

test_that("derive_var_atoxgr: Test 3 ATOXGR > 0 (HYPER)", {
  exp_out_3 <- tibble::tribble(
    ~ATOXDSCL,          ~ATOXDSCH,        ~ATOXGRL,      ~ATOXGRH,       ~ATOXGR,
    "Hypoglycemia",     "Hyperglycemia",  NA_character_, "1",            "1",
    "Hypoglycemia",     "Hyperglycemia",  "0",           "4",            "4",
    NA_character_,      "INR Increased",  NA_character_, "2",            "2",
  )

  input_3 <- exp_out_3 %>%
    select(-ATOXGR)

  expect_equal(
    derive_var_atoxgr(
      input_3,
      lotox_description_var = ATOXDSCL,
      hitox_description_var = ATOXDSCH
    ),
    exp_out_3
  )
})

test_that("derive_var_atoxgr: Test 4 ATOXGR < 0 (HYPO)", {
  exp_out_4 <- tibble::tribble(
    ~ATOXDSCL,          ~ATOXDSCH,        ~ATOXGRL,      ~ATOXGRH,       ~ATOXGR,
    "Hypoglycemia",     "Hyperglycemia",  "3",           NA_character_,  "-3",
    "Hypoglycemia",     "Hyperglycemia",  "1",           "0",            "-1",
    "Hypophosphatemia", NA_character_,    "4",           NA_character_,  "-4",
  )

  input_4 <- exp_out_4 %>%
    select(-ATOXGR)

  expect_equal(
    derive_var_atoxgr(
      input_4,
      lotox_description_var = ATOXDSCL,
      hitox_description_var = ATOXDSCH
    ),
    exp_out_4
  )
})


## Blood and lymphatic system disorders ----
## Grade 3: <80 g/L
## Grade 2: <100 - 80g/L
## Grade 1: <LLN - 100 g/L

### 1. Anemia ----
test_that("derive_var_atoxgr_dir: Test 1 NCICTCAEv4 Anemia", {
  exp_out_ctcv4_1 <- tibble::tribble(
    ~ATOXDSCL,      ~AVAL,  ~ANRLO, ~ANRHI, ~AVALU, ~ATOXGRL,
    "Not a term",   80,     120,    200,    "G/L",  NA,
    NA_character_,  60,     50,     100,    "G/L",  NA,
    "ANEMIA",       79,     140,    NA,     "G/L",  "3",
    "ANEMIA",       80,     140,    NA,     "G/L",  "2",
    "Anemia",       99,     140,    NA,     "G/L",  "2",
    "Anemia",       100,    140,    NA,     "G/L",  "1",
    # wrong UNIT - GRADE should be missing
    "anemia",       100,    140,    NA,     "G/dL", NA,
    "Anemia",       139,    140,    NA,     "G/L",  "1",
    "ANEMIA",       140,    140,    NA,     "g/L",  "0",
    # ANRLO missing AVAL not grade 2 or 3 - cannot grade
    "Anemia",       140,    NA,     NA,     "G/L",  NA,
    "Anemia",       139,    NA,     NA,     "G/L",  NA,
    "Anemia",       100,    NA,     NA,     "G/L",  NA,
    # ANRLO missing but AVAL satisfies grade 2
    "Anemia",       99,     NA,     NA,     "G/L",  "2",
    # ANRLO missing but AVAL satisfies grade 3
    "Anemia",       79,     NA,     NA,     "G/L",  "3",
    # Unit missing cannot grade
    "Anemia",       140,    140,    NA,     NA,     NA,
    # AVAL missing cannot grade
    "Anemia",       NA,     140,    NA,     "G/L",  NA,
  )
  input_ctcv4_1 <- exp_out_ctcv4_1 %>%
    select(-ATOXGRL)

  actual_output_ctcv4_1 <- derive_var_atoxgr_dir(
    input_ctcv4_1,
    new_var = ATOXGRL,
    meta_criteria = atoxgr_criteria_ctcv4,
    tox_description_var = ATOXDSCL,
    criteria_direction = "L",
    get_unit_expr = AVALU
  )

  expect_dfs_equal(
    base = exp_out_ctcv4_1,
    compare = actual_output_ctcv4_1,
    keys = c("ATOXDSCL", "AVAL", "ANRLO", "ANRHI", "AVALU")
  )
})

### 2. Leukocytosis ----
### Grade 3: >100,000/mm3

test_that("derive_var_atoxgr_dir: Test 2 NCICTCAEv4 Leukocytosis", {
  exp_out_ctcv4_2 <- tibble::tribble(
    ~ATOXDSCL,      ~AVAL,  ~ANRLO, ~ANRHI, ~AVALU,    ~ATOXGRL,
    "Not a term",   99,     0,      NA,     "10^9/L",  NA,
    NA,             99,     0,      NA,     "10^9/L",  NA,
    "Leukocytosis", 101,    0,      40,     "10^9/L",  "3",
    "leukocytosis", 100,    0,      40,     "10^9/L",  "0",
    "Leukocytosis", 99,     0,      NA,     "10^9/L",  "0",
    # wrong UNIT - GRADE should be missing
    "Leukocytosis", 99,     0,      40,     "10^9/M",  NA,
    # Unit missing cannot grade
    "Leukocytosis", 99,     0,      40,     NA,        NA,
    # AVAL missing cannot grade
    "Leukocytosis", NA,     0,      40,     "10^9/L",  NA,
  )
  input_ctcv4_2 <- exp_out_ctcv4_2 %>%
    select(-ATOXGRL)

  actual_output_ctcv4_2 <- derive_var_atoxgr_dir(
    input_ctcv4_2,
    new_var = ATOXGRL,
    meta_criteria = atoxgr_criteria_ctcv4,
    tox_description_var = ATOXDSCL,
    criteria_direction = "L",
    get_unit_expr = AVALU
  )

  expect_dfs_equal(
    base = exp_out_ctcv4_2,
    compare = actual_output_ctcv4_2,
    keys = c("ATOXDSCL", "AVAL", "ANRLO", "ANRHI", "AVALU")
  )
})

## Investigations ----

### 3. Activated partial thromboplastin time prolonged ----
### Grade 3: >2.5 x ULN
### Grade 2: >1.5 - 2.5 x ULN
### Grade 1: >ULN - 1.5 x ULN

test_that("derive_var_atoxgr_dir: Test 3 CTCAEv4 Activated partial thromboplastin time prolonged", {
  exp_out_ctcv4_3 <- tibble::tribble(
    ~ATOXDSCH,                                         ~AVAL,  ~ANRHI,  ~AVALU,         ~ATOXGRH,
    "Not a term",                                      80,     100,     NA_character_,  NA,
    NA_character_,                                     60,     100,     NA_character_,  NA,
    "Activated partial thromboplastin time prolonged", 251,    100,     NA_character_,  "3",
    "Activated Partial thromboplastin time prolonged", 250,    100,     NA_character_,  "2",
    "Activated partial Thromboplastin time prolonged", 151,    100,     NA_character_,  "2",
    "Activated partial thromboplastin time prolonged", 150,    100,     NA_character_,  "1",
    "Activated partial thromboplastin Time prolonged", 101,    100,     NA_character_,  "1",
    "Activated partial thromboplastin time prolonged", 100,    100,     NA_character_,  "0",
    # ANRHI missing - cannot grade
    "Activated partial thromboplastin time prolonged", 100,    NA,      NA_character_,  NA,
    # AVAL missing cannot grade
    "Activated partial thromboplastin time prolonged", NA,     100,     NA_character_,  NA,
  )
  input_ctcv4_3 <- exp_out_ctcv4_3 %>%
    select(-ATOXGRH)

  actual_output_ctcv4_3 <- derive_var_atoxgr_dir(
    input_ctcv4_3,
    new_var = ATOXGRH,
    meta_criteria = atoxgr_criteria_ctcv4,
    tox_description_var = ATOXDSCH,
    criteria_direction = "H",
    get_unit_expr = AVALU
  )

  expect_dfs_equal(
    base = exp_out_ctcv4_3,
    compare = actual_output_ctcv4_3,
    keys = c("ATOXDSCH", "AVAL", "ANRHI", "AVALU")
  )
})

### 4. Alanine aminotransferase increased ----
### Grade 4: >20.0 x ULN
### Grade 3: >5.0 - 20.0 x ULN
### Grade 2: >3.0 - 5.0 x ULN
### Grade 1: >ULN - 3.0 x ULN

test_that("derive_var_atoxgr_dir: Test 4 NCICTCAEv4 Alanine aminotransferase increased", {
  exp_out_ctcv4_4 <- tibble::tribble(
    ~ATOXDSCH,                            ~AVAL,  ~ANRHI, ~AVALU,         ~ATOXGRH,
    "Not a term",                         80,     40,     NA_character_,  NA,
    NA_character_,                        60,     40,     NA_character_,  NA,
    "Alanine aminotransferase Increased", 801,    40,     NA_character_,  "4",
    "Alanine aminotransferase Increased", 800,    40,     NA_character_,  "3",
    "Alanine aminotransferase Increased", 201,    40,     NA_character_,  "3",
    "Alanine aminotransferase Increased", 200,    40,     NA_character_,  "2",
    "Alanine aminotransferase Increased", 121,    40,     NA_character_,  "2",
    "Alanine aminotransferase Increased", 120,    40,     NA_character_,  "1",
    "Alanine aminotransferase Increased", 41,     40,     NA_character_,  "1",
    "Alanine aminotransferase Increased", 40,     40,     NA_character_,  "0",
    # ANRHI missing - cannot grade
    "Alanine aminotransferase Increased", 100,    NA,     NA_character_,  NA,
    # AVAL missing cannot grade
    "Alanine aminotransferase Increased", NA,     40,     NA_character_,  NA,
  )
  input_ctcv4_4 <- exp_out_ctcv4_4 %>%
    select(-ATOXGRH)

  actual_output_ctcv4_4 <- derive_var_atoxgr_dir(
    input_ctcv4_4,
    new_var = ATOXGRH,
    meta_criteria = atoxgr_criteria_ctcv4,
    tox_description_var = ATOXDSCH,
    criteria_direction = "H",
    get_unit_expr = AVALU
  )

  expect_dfs_equal(
    base = exp_out_ctcv4_4,
    compare = actual_output_ctcv4_4,
    keys = c("ATOXDSCH", "AVAL", "ANRHI", "AVALU")
  )
})

### 5. Alkaline phosphatase increased ----
### Grade 4: >20.0 x ULN
### Grade 3: >5.0 - 20.0 x ULN
### Grade 2: >2.5 - 5.0 x ULN
### Grade 1: >ULN - 2.5 x ULN

test_that("derive_var_atoxgr_dir: Test 5 NCICTCAEv4 Alkaline phosphatase increased", {
  exp_out_ctcv4_5 <- tibble::tribble(
    ~ATOXDSCH,                         ~AVAL,  ~ANRHI, ~AVALU,         ~ATOXGRH,
    "Not a term",                      80,     40,     NA_character_,  NA,
    NA_character_,                     60,     40,     NA_character_,  NA,
    "Alkaline phosphatase increased",  801,    40,     NA_character_,  "4",
    "Alkaline phosphatase increased",  800,    40,     NA_character_,  "3",
    "Alkaline phosphatase increased",  201,    40,     NA_character_,  "3",
    "Alkaline phosphatase increased",  200,    40,     NA_character_,  "2",
    "Alkaline phosphatase increased",  101,    40,     NA_character_,  "2",
    "Alkaline phosphatase increased",  100,    40,     NA_character_,  "1",
    "Alkaline phosphatase increased",  41,     40,     NA_character_,  "1",
    "Alkaline phosphatase increased",  40,     40,     NA_character_,  "0",
    # ANRHI missing - cannot grade
    "Alkaline phosphatase increased",  100,    NA,     NA_character_,  NA,
    # AVAL missing cannot grade
    "Alkaline phosphatase increased",  NA,     40,     NA_character_,  NA,
  )
  input_ctcv4_5 <- exp_out_ctcv4_5 %>%
    select(-ATOXGRH)

  actual_output_ctcv4_5 <- derive_var_atoxgr_dir(
    input_ctcv4_5,
    new_var = ATOXGRH,
    meta_criteria = atoxgr_criteria_ctcv4,
    tox_description_var = ATOXDSCH,
    criteria_direction = "H",
    get_unit_expr = AVALU
  )

  expect_dfs_equal(
    base = exp_out_ctcv4_5,
    compare = actual_output_ctcv4_5,
    keys = c("ATOXDSCH", "AVAL", "ANRHI", "AVALU")
  )
})

### 6. Aspartate aminotransferase increased ----
### Grade 4: >20.0 x ULN
### Grade 3: >5.0 - 20.0 x ULN
### Grade 2: >3.0 - 5.0 x ULN
### Grade 1: >ULN - 3.0 x ULN

test_that("derive_var_atoxgr_dir: Test 6 NCICTCAEv4 Aspartate aminotransferase increased", {
  exp_out_ctcv4_6 <- tibble::tribble(
    ~ATOXDSCH,                              ~AVAL,  ~ANRHI, ~AVALU,         ~ATOXGRH,
    "Not a term",                           80,     40,     NA_character_,  NA,
    NA_character_,                          60,     40,     NA_character_,  NA,
    "Aspartate aminotransferase Increased", 801,    40,     NA_character_,  "4",
    "Aspartate aminotransferase Increased", 800,    40,     NA_character_,  "3",
    "Aspartate aminotransferase Increased", 201,    40,     NA_character_,  "3",
    "Aspartate aminotransferase Increased", 200,    40,     NA_character_,  "2",
    "Aspartate aminotransferase Increased", 121,    40,     NA_character_,  "2",
    "Aspartate aminotransferase Increased", 120,    40,     NA_character_,  "1",
    "Aspartate aminotransferase Increased", 41,     40,     NA_character_,  "1",
    "Aspartate aminotransferase Increased", 40,     40,     NA_character_,  "0",
    # ANRHI missing - cannot grade
    "Aspartate aminotransferase Increased", 100,    NA,     NA_character_,  NA,
    # AVAL missing cannot grade
    "Aspartate aminotransferase Increased", NA,     40,     NA_character_,  NA,
  )
  input_ctcv4_6 <- exp_out_ctcv4_6 %>%
    select(-ATOXGRH)

  actual_output_ctcv4_6 <- derive_var_atoxgr_dir(
    input_ctcv4_6,
    new_var = ATOXGRH,
    meta_criteria = atoxgr_criteria_ctcv4,
    tox_description_var = ATOXDSCH,
    criteria_direction = "H",
    get_unit_expr = AVALU
  )

  expect_dfs_equal(
    base = exp_out_ctcv4_6,
    compare = actual_output_ctcv4_6,
    keys = c("ATOXDSCH", "AVAL", "ANRHI", "AVALU")
  )
})

### 7. Blood bilirubin increased ----
### Grade 4: >10.0 x ULN
### Grade 3: >3.0 - 10.0 x ULN
### Grade 2: >3.0 - 1.5 x ULN
### Grade 1: >ULN - 1.5 x ULN

test_that("derive_var_atoxgr_dir: Test 7 NCICTCAEv4 Blood bilirubin increased", {
  exp_out_ctcv4_7 <- tibble::tribble(
    ~ATOXDSCH,                   ~AVAL,  ~ANRHI, ~AVALU,         ~ATOXGRH,
    "Not a term",                80,     40,     NA_character_,  NA,
    NA_character_,               60,     40,     NA_character_,  NA,
    "Blood bilirubin increased", 401,    40,     NA_character_,  "4",
    "Blood bilirubin increased", 400,    40,     NA_character_,  "3",
    "Blood bilirubin increased", 121,    40,     NA_character_,  "3",
    "Blood bilirubin increased", 120,    40,     NA_character_,  "2",
    "Blood bilirubin increased", 61,     40,     NA_character_,  "2",
    "Blood bilirubin increased", 60,     40,     NA_character_,  "1",
    "Blood bilirubin increased", 41,     40,     NA_character_,  "1",
    "Blood bilirubin increased", 40,     40,     NA_character_,  "0",
    # ANRHI missing - cannot grade
    "Blood bilirubin increased", 100,    NA,     NA_character_,  NA,
    # AVAL missing cannot grade
    "Blood bilirubin increased", NA,     40,     NA_character_,  NA,
  )
  input_ctcv4_7 <- exp_out_ctcv4_7 %>%
    select(-ATOXGRH)

  actual_output_ctcv4_7 <- derive_var_atoxgr_dir(
    input_ctcv4_7,
    new_var = ATOXGRH,
    meta_criteria = atoxgr_criteria_ctcv4,
    tox_description_var = ATOXDSCH,
    criteria_direction = "H",
    get_unit_expr = AVALU
  )

  expect_dfs_equal(
    base = exp_out_ctcv4_7,
    compare = actual_output_ctcv4_7,
    keys = c("ATOXDSCH", "AVAL", "ANRHI", "AVALU")
  )
})

### 8. CD4 Lymphocytes decreased ----
### Grade 4: <0.05 x 10e9 /L
### Grade 3: <0.2 x 0.05 - 10e9 /L
### Grade 2: <0.5 - 0.2 x 10e9 /L
### Grade 1: <LLN - 0.5 x 10e9 /L

test_that("derive_var_atoxgr_dir: Test 8 NCICTCAEv4 CD4 Lymphocytes decreased", {
  exp_out_ctcv4_8 <- tibble::tribble(
    ~ATOXDSCL,                    ~AVAL,  ~ANRLO, ~AVALU,    ~ATOXGRL,
    "Not a term",                 80,     120,    "10^9/L",  NA,
    NA_character_,                60,     50,     "10^9/L",  NA,
    "CD4 lymphocytes decreased",  0.04,   0.8,    "10^9/L",  "4",
    "CD4 lymphocytes decreased",  0.05,   0.8,    "10^9/l",  "3",
    "CD4 lymphocytes decreased",  0.19,   0.8,    "10^9/L",  "3",
    "CD4 lymphocytes decreased",  0.2,    0.8,    "10^9/L",  "2",
    "CD4 lymphocytes decreased",  0.49,   0.8,    "10^9/L",  "2",
    # wrong unit - grade missing
    "CD4 lymphocytes decreased",  0.49,   0.8,    "10^8/L",  NA,
    "CD4 lymphocytes decreased",  0.5,    0.8,    "10^9/L",  "1",
    "CD4 lymphocytes decreased",  0.79,   0.8,    "10^9/L",  "1",
    "CD4 lymphocytes decreased",  0.8,    0.8,    "10^9/L",  "0",
    # ANRLO missing - AVAL satisfies grade 2 - 4
    "CD4 lymphocytes decreased",  0.04,   NA,     "10^9/L",  "4",
    "CD4 lymphocytes decreased",  0.05,   NA,     "10^9/L",  "3",
    "CD4 lymphocytes decreased",  0.19,   NA,     "10^9/L",  "3",
    "CD4 lymphocytes decreased",  0.2,    NA,     "10^9/L",  "2",
    "CD4 lymphocytes decreased",  0.49,   NA,     "10^9/L",  "2",
    # ANRLO missing - AVAL does NOT satisfies grade 2 - 4
    "CD4 lymphocytes decreased",  0.5,    NA,     "10^9/L",  NA,
    "CD4 lymphocytes decreased",  0.79,   NA,     "10^9/L",  NA,
    "CD4 lymphocytes decreased",  0.8,    NA,     "10^9/L",  NA,
    # Unit missing - cannot grade
    "CD4 lymphocytes decreased",  0.8,    0.8,    NA,        NA,
    # AVAL missing cannot grade
    "CD4 lymphocytes decreased",  NA,     0.8,    "10^9/L",  NA,
  )

  input_ctcv4_8 <- exp_out_ctcv4_8 %>%
    select(-ATOXGRL)

  actual_output_ctcv4_8 <- derive_var_atoxgr_dir(
    input_ctcv4_8,
    new_var = ATOXGRL,
    meta_criteria = atoxgr_criteria_ctcv4,
    tox_description_var = ATOXDSCL,
    criteria_direction = "L",
    get_unit_expr = AVALU
  )

  expect_dfs_equal(
    base = exp_out_ctcv4_8,
    compare = actual_output_ctcv4_8,
    keys = c("ATOXDSCL", "AVAL", "ANRLO", "AVALU")
  )
})

### 9. Cholesterol high ----
### Grade 4: >12.92 mmol/L
### Grade 3: >10.34 - 12.92 mmol/L
### Grade 2: >7.75 - 10.34 mmol/L
### Grade 1: >ULN - 7.75 mmol/L

test_that("derive_var_atoxgr_dir: Test 9 NCICTCAEv4 Cholesterol high", {
  exp_out_ctcv4_9 <- tibble::tribble(
    ~ATOXDSCH,          ~AVAL,  ~ANRLO, ~ANRHI, ~AVALU,    ~ATOXGRH,
    "Not a term",       8,      0,      5,      "mmol/L",  NA,
    NA_character_,      10,     0,      5,      "mmol/L",  NA,
    "Cholesterol high", 12.93,  0,      5,      "mmol/L",  "4",
    "Cholesterol High", 12.92,  0,      5,      "mmol/L",  "3",
    "Cholesterol high", 10.35,  0,      5,      "Mmol/L",  "3",
    # wrong unit - grade missing
    "Cholesterol high", 10.35,  0,      5,      "umol/L",  NA,
    "Cholesterol high", 10.34,  0,      5,      "mmol/L",  "2",
    "Cholesterol high", 7.76,   0,      5,      "mmol/L",  "2",
    "Cholesterol high", 7.75,   0,      5,      "mmol/L",  "1",
    "Cholesterol high", 5.1,    0,      5,      "mmol/L",  "1",
    "Cholesterol high", 5,      0,      5,      "mmol/L",  "0",
    # ANRHI missing - AVAL satisfies grade 2 - 4
    "Cholesterol high", 12.93,  0,      NA,     "mmol/L",  "4",
    "Cholesterol High", 12.92,  0,      NA,     "mmol/L",  "3",
    "Cholesterol high", 10.35,  0,      NA,     "Mmol/L",  "3",
    "Cholesterol high", 10.34,  0,      NA,     "mmol/L",  "2",
    "Cholesterol high", 7.76,   0,      NA,     "mmol/L",  "2",
    # ANRHI missing - AVAL does NOT satisfies grade 2 - 4
    "Cholesterol high", 7.75,   0,      NA,     "mmol/L",  NA,
    "Cholesterol high", 5.1,    0,      NA,     "mmol/L",  NA,
    "Cholesterol high", 5,      0,      NA,     "mmol/L",  NA,
    # Unit missing - cannot grade
    "Cholesterol high", 5,      0,      5,      NA,        NA,
    # AVAL missing cannot grade
    "Cholesterol high", NA,     0,      5,      "mmol/L",  NA,
  )
  input_ctcv4_9 <- exp_out_ctcv4_9 %>%
    select(-ATOXGRH)

  actual_output_ctcv4_9 <- derive_var_atoxgr_dir(
    input_ctcv4_9,
    new_var = ATOXGRH,
    meta_criteria = atoxgr_criteria_ctcv4,
    tox_description_var = ATOXDSCH,
    criteria_direction = "H",
    get_unit_expr = AVALU
  )

  expect_dfs_equal(
    base = exp_out_ctcv4_9,
    compare = actual_output_ctcv4_9,
    keys = c("ATOXDSCH", "AVAL", "ANRLO", "ANRHI", "AVALU")
  )
})

### 10. CPK increased ----
### Grade 4: >10.0 x ULN
### Grade 3: >5.0 - 10.0 x ULN
### Grade 2: >2.5 - 5.0 x ULN
### Grade 1: >ULN - 2.5 x ULN

test_that("derive_var_atoxgr_dir: Test 10 NCICTCAEv4 CPK increased", {
  exp_out_ctcv4_10 <- tibble::tribble(
    ~ATOXDSCH,        ~AVAL,  ~ANRLO,  ~ANRHI,  ~AVALU,         ~ATOXGRH,
    "Not a term",     80,     0,       40,      NA_character_,  NA,
    NA_character_,    60,     0,       40,      NA_character_,  NA,
    "CPK increased",  401,    0,       40,      NA_character_,  "4",
    "CPK increased",  400,    0,       40,      NA_character_,  "3",
    "CPK increased",  201,    0,       40,      NA_character_,  "3",
    "CPK increased",  200,    0,       40,      NA_character_,  "2",
    "CPK increased",  101,    0,       40,      NA_character_,  "2",
    "CPK increased",  100,    0,       40,      NA_character_,  "1",
    "CPK increased",  41,     0,       40,      NA_character_,  "1",
    "CPK increased",  40,     0,       40,      NA_character_,  "0",
    # ANRHI missing - cannot grade
    "CPK increased",  100,    0,       NA,      NA_character_,  NA,
    # AVAL missing cannot grade
    "CPK increased",  NA,     0,       40,      NA_character_,  NA,
  )
  input_ctcv4_10 <- exp_out_ctcv4_10 %>%
    select(-ATOXGRH)

  actual_output_ctcv4_10 <- derive_var_atoxgr_dir(
    input_ctcv4_10,
    new_var = ATOXGRH,
    meta_criteria = atoxgr_criteria_ctcv4,
    tox_description_var = ATOXDSCH,
    criteria_direction = "H",
    get_unit_expr = AVALU
  )

  expect_dfs_equal(
    base = exp_out_ctcv4_10,
    compare = actual_output_ctcv4_10,
    keys = c("ATOXDSCH", "AVAL", "ANRLO", "ANRHI", "AVALU")
  )
})

### 11. Creatinine increased ----
### Grade 4: >6.0 x ULN
### Grade 3: >3.0 baseline; >3.0 - 6.0 x ULN
### Grade 2: >1.5 - 3.0 x baseline; >1.5 - 3.0 x ULN
### Grade 1: >1 - 1.5 x baseline; >ULN - 1.5 x ULN

test_that("derive_var_atoxgr_dir: Test 10 NCICTCAEv4 Creatinine increased", {
  exp_out_ctcv4_11 <- tibble::tribble(
    ~ATOXDSCH,               ~AVAL,  ~BASE, ~ANRHI, ~AVALU,         ~ATOXGRH,
    "Not a term",            80,     80,    40,     NA_character_,  NA,
    NA_character_,           60,     60,    40,     NA_character_,  NA,
    # GRADE derived from AVAL against ANRHI
    "Creatinine increased",  241,    241,   40,     NA_character_,  "4",
    "Creatinine increased",  240,    230,   40,     NA_character_,  "3",
    "Creatinine increased",  121,    120,   40,     NA_character_,  "3",
    "Creatinine increased",  120,    119,   40,     NA_character_,  "2",
    "Creatinine increased",  61,     60,    40,     NA_character_,  "2",
    "Creatinine increased",  60,     60,    40,     NA_character_,  "1",
    "Creatinine increased",  41,     41,    40,     NA_character_,  "1",
    "Creatinine increased",  40,     40,    40,     NA_character_,  "0",
    # GRADE derived from AVAL against BASE
    "Creatinine increased",  42,     6,     40,     NA_character_,  "3",
    "Creatinine increased",  42,     13.9,  40,     NA_character_,  "3",
    "Creatinine increased",  42,     14,    40,     NA_character_,  "2",
    "Creatinine increased",  42.1,   28,    40,     NA_character_,  "2",
    "Creatinine increased",  42,     28,    42,     NA_character_,  "1",
    "Creatinine increased",  42,     41,    42,     NA_character_,  "1",
    "Creatinine increased",  42,     42,    42,     NA_character_,  "0",
    # BASE missing - AVAL <= ANRLO cannot grade as NORMAL
    "Creatinine increased",  42,     NA,    42,     NA_character_,  NA,
    # ANRHI missing - AVAL <= BASE cannot grade as NORMAL
    "Creatinine increased",  42,     42,    NA,     NA_character_,  NA,
    # AVAL missing cannot grade
    "Creatinine increased",  NA,     0,     40,     NA_character_,  NA,
  )
  input_ctcv4_11 <- exp_out_ctcv4_11 %>%
    select(-ATOXGRH)

  actual_output_ctcv4_11 <- derive_var_atoxgr_dir(
    input_ctcv4_11,
    new_var = ATOXGRH,
    meta_criteria = atoxgr_criteria_ctcv4,
    tox_description_var = ATOXDSCH,
    criteria_direction = "H",
    get_unit_expr = AVALU
  )

  expect_dfs_equal(
    base = exp_out_ctcv4_11,
    compare = actual_output_ctcv4_11,
    keys = c("ATOXDSCH", "AVAL", "BASE", "ANRHI", "AVALU")
  )
})

### 12. Fibrinogen decreased decreased ----
### Grade 4: <0.25 x LLN or 75% decrease from baseline or absolute value <50 mg/dL
### Grade 3: <0.5 - 0.25 x LLN or 50 - <75% decrease from baseline
### Grade 2: <0.75 - 0.5 x LLN or 25 - <50% decrease from baseline
### Grade 1: <1.0 - 0.75 x LLN or <25% decrease from baseline

test_that("derive_var_atoxgr_dir: Test 12 NCICTCAEv4 Fibrinogen decreased", {
  exp_out_ctcv4_12 <- tibble::tribble(
    ~ATOXDSCL,               ~AVAL,  ~ANRLO, ~PCHG,  ~AVALU,  ~ATOXGRL,
    "Not a term",            9,      10,     40,     "g/L",   NA,
    NA_character_,           10,     10,     40,     "g/L",   NA,
    # Satisfies < 0.5 for grade 4 - other criteria missing
    "Fibrinogen decreased",  0.49,   NA,     NA,     "g/L",   "4",
    # Satisfies < 0.5 for grade 4 - satisfies grade 3 for other criteria
    "Fibrinogen decreased",  0.49,   1,      -51,    "g/L",   "4",
    # Satisfies < 0.25*LLN for grade 4 - PCHG missing
    "Fibrinogen decreased",  0.5,    2.1,    NA,     "g/L",   "4",
    # Satisfies < 0.25*LLN for grade 4 - PCHG satisfies grade 3
    "Fibrinogen decreased",  0.5,    2.1,    -51,    "g/L",   "4",
    # Satisfies <=75% decrease for grade 4 - LLN  missing
    "Fibrinogen decreased",  1,      NA,     -75,    "g/L",   "4",
    # Satisfies <=75% decrease for grade 4 - LLN  satisfies grade 3
    "Fibrinogen decreased",  1,      0.49,   -75,    "g/L",   "4",
    # Satisfies < 0.5*LLN for grade 3 - PCHG missing
    "Fibrinogen decreased",  1,      2.1,    NA,     "g/L",   "3",
    # Satisfies < 0.5*LLN for grade 3 - PCHG satisfies grade 2
    "Fibrinogen decreased",  1,      2.1,    -49,    "g/L",   "3",
    # Satisfies <=50% decrease for grade 3 - LLN  missing
    "Fibrinogen decreased",  1,      NA,     -50,    "g/L",   "3",
    # Satisfies <=50% decrease for grade 3 - LLN  satisfies grade 2
    "Fibrinogen decreased",  1,      2,      -50,    "g/L",   "3",
    # Satisfies < 0.75*LLN for grade 2 - PCHG missing
    "Fibrinogen decreased",  1.5,    2.1,    NA,     "g/L",   "2",
    # Satisfies < 0.75*LLN for grade 2 - PCHG satisfies grade 1
    "Fibrinogen decreased",  1.5,    2.1,    -10,    "g/L",   "2",
    # Satisfies <=25% for grade 2 - LLN missing
    "Fibrinogen decreased",  1.5,    NA,     -25,    "g/L",   "2",
    # Satisfies <=25% for grade 2 - LLN satisfies grade 1
    "Fibrinogen decreased",  1.5,    1.6,    -25,    "g/L",   "2",
    # Satisfies < LLN for grade 1 - PCHG missing
    "Fibrinogen decreased",  2,      2.1,    NA,     "g/L",   "1",
    # Satisfies < LLN for grade 1 - PCHG satisfies grade 0
    "Fibrinogen decreased",  2,      2.1,    10,     "g/L",   "1",
    # Satisfies % decrease for grade 1 - LLN missing
    "Fibrinogen decreased",  1.5,    NA,     -1,     "g/L",   "1",
    # Satisfies % decrease for grade 1 - AVAL = LLN
    "Fibrinogen decreased",  1.5,    1.5,    -1,     "g/L",   "1",
    # Satisfies grade 0 - AVAL >= LLN AND no % descrease
    "Fibrinogen decreased",  1.5,    1.5,    0,      "g/L",   "0",
    # AVAL >= LLN BUT PCT missing cannot grade as NORMAL
    "Fibrinogen decreased",  1.5,    1.5,    NA,     "g/L",   NA,
    # PCT >= 0 BUT LLN missing cannot grade as NORMAL
    "Fibrinogen decreased",  1.5,    NA,     10,     "g/L",   NA,
    # AVAL missing cannot grade
    "Fibrinogen decreased",  NA,     1.5,    10,     "g/L",   NA,
    # wrong unit cannot grade as it may satisfy grade 4
    "Fibrinogen decreased",  1.5,    1.5,    0,      "g/dL",  NA,
    # missing unit cannot grade as it may satisfy grade 4
    "Fibrinogen decreased",  1.5,    1.5,    0,      NA,      NA,
  )
  input_ctcv4_12 <- exp_out_ctcv4_12 %>%
    select(-ATOXGRL)

  actual_output_ctcv4_12 <- derive_var_atoxgr_dir(
    input_ctcv4_12,
    new_var = ATOXGRL,
    meta_criteria = atoxgr_criteria_ctcv4,
    tox_description_var = ATOXDSCL,
    criteria_direction = "L",
    get_unit_expr = AVALU
  )

  expect_dfs_equal(
    base = exp_out_ctcv4_12,
    compare = actual_output_ctcv4_12,
    keys = c("ATOXDSCL", "AVAL", "ANRLO", "PCHG", "AVALU")
  )
})

### 13. GGT increased ----
### Grade 4: >20.0 x ULN
### Grade 3: >5.0 - 20.0 x ULN
### Grade 2: >2.5 - 5.0 x ULN
### Grade 1: >ULN - 2.5 x ULN

test_that("derive_var_atoxgr_dir: Test 13 NCICTCAEv4 GGT increased", {
  exp_out_ctcv4_13 <- tibble::tribble(
    ~ATOXDSCH,       ~AVAL, ~ANRLO, ~ANRHI, ~AVALU,        ~ATOXGRH,
    "Not a term",    80,    0,      40,     NA_character_, NA,
    NA_character_,   60,    0,      40,     NA_character_, NA,
    "GGT increased", 801,   0,      40,     NA_character_, "4",
    "GGT increased", 800,   0,      40,     NA_character_, "3",
    "GGT increased", 201,   0,      40,     NA_character_, "3",
    "GGT increased", 200,   0,      40,     NA_character_, "2",
    "GGT increased", 101,   0,      40,     NA_character_, "2",
    "GGT increased", 100,   0,      40,     NA_character_, "1",
    "GGT increased", 41,    0,      40,     NA_character_, "1",
    "GGT increased", 40,    0,      40,     NA_character_, "0",
    # ANRHI missing - cannot grade
    "GGT increased", 100,   0,      NA,     NA_character_, NA,
    # AVAL missing cannot grade
    "GGT increased", NA,    0,      NA,     NA_character_, NA,
  )
  input_ctcv4_13 <- exp_out_ctcv4_13 %>%
    select(-ATOXGRH)

  actual_output_ctcv4_13 <- derive_var_atoxgr_dir(
    input_ctcv4_13,
    new_var = ATOXGRH,
    meta_criteria = atoxgr_criteria_ctcv4,
    tox_description_var = ATOXDSCH,
    criteria_direction = "H",
    get_unit_expr = AVALU
  )

  expect_dfs_equal(
    base = exp_out_ctcv4_13,
    compare = actual_output_ctcv4_13,
    keys = c("ATOXDSCH", "AVAL", "ANRLO", "ANRHI", "AVALU")
  )
})

### 14. Haptoglobin decreased ----
### Grade 1: <LLN

test_that("derive_var_atoxgr_dir: Test 14 NCICTCAEv4 Haptoglobin decreased", {
  exp_out_ctcv4_14 <- tibble::tribble(
    ~ATOXDSCL,               ~AVAL,  ~ANRLO, ~ANRHI, ~AVALU,         ~ATOXGRL,
    "Not a term",            9,      10,     40,     NA_character_,  NA,
    NA_character_,           10,     10,     40,     NA_character_,  NA,
    "Haptoglobin decreased", 9,      10,     40,     NA_character_,  "1",
    "Haptoglobin decreased", 10,     10,     40,     NA_character_,  "0",
    "Haptoglobin decreased", 11,     10,     40,     NA_character_,  "0",
    # ANRHI missing - cannot grade
    "Haptoglobin decreased", NA,     10,     40,     NA_character_,  NA,
    # AVAL missing cannot grade
    "Haptoglobin decreased", 10,     NA,     NA,     NA_character_,  NA,
  )
  input_ctcv4_14 <- exp_out_ctcv4_14 %>%
    select(-ATOXGRL)

  actual_output_ctcv4_14 <- derive_var_atoxgr_dir(
    input_ctcv4_14,
    new_var = ATOXGRL,
    meta_criteria = atoxgr_criteria_ctcv4,
    tox_description_var = ATOXDSCL,
    criteria_direction = "L",
    get_unit_expr = AVALU
  )

  expect_dfs_equal(
    base = exp_out_ctcv4_14,
    compare = actual_output_ctcv4_14,
    keys = c("ATOXDSCL", "AVAL", "ANRLO", "ANRHI", "AVALU")
  )
})

### 15. Hemoglobin increased ----
### Grade 3: Increase in >4 gm/dL above ULN or above baseline if baseline is above ULN
### Grade 2: Increase in >2 - 4 gm/dL above ULN or above baseline if baseline is above ULN
### Grade 1: Increase in >0 - 2 gm/dL above ULN or above baseline if baseline is above ULN

test_that("derive_var_atoxgr_dir: Test 15 NCICTCAEv4 Hemoglobin increased", {
  exp_out_ctcv4_15 <- tibble::tribble(
    ~ATOXDSCH,              ~AVAL,  ~BASE,  ~ANRHI, ~AVALU, ~ATOXGRH, ~TESTNUM,
    "Not a term",           80,     120,    200,    "g/L",  NA,       1,
    NA_character_,          60,     50,     100,    "g/L",  NA,       2,
    # BASE greater than ANRHI
    "Hemoglobin increased", 106,    65,     60,     "g/L",  "3",      3,
    "Hemoglobin increased", 105,    65,     60,     "g/L",  "2",      4,
    "Hemoglobin increased", 86,     65,     60,     "g/L",  "2",      5,
    "Hemoglobin increased", 85,     65,     60,     "g/L",  "1",      6,
    "Hemoglobin increased", 66,     65,     60,     "g/L",  "1",      7,
    "Hemoglobin increased", 65,     65,     60,     "g/L",  "0",      8,
    "Hemoglobin increased", NA,     65,     60,     "g/L",  NA,       9,
    # BASE less than or equal to ANRHI
    "Hemoglobin increased", 106,    60,     65,     "g/L",  "3",      10,
    "Hemoglobin increased", 105,    60,     65,     "g/L",  "2",      11,
    "Hemoglobin increased", 86,     60,     65,     "g/L",  "2",      12,
    "Hemoglobin increased", 85,     60,     65,     "g/L",  "1",      13,
    "Hemoglobin increased", 66,     60,     65,     "g/L",  "1",      14,
    "Hemoglobin increased", 65,     60,     65,     "g/L",  "0",      15,
    "Hemoglobin increased", NA,     60,     65,     "g/L",  NA,       16,
    # BASE missing
    "Hemoglobin increased", 106,    NA,     65,     "g/L",  "3",      17,
    "Hemoglobin increased", 105,    NA,     65,     "g/L",  "2",      18,
    "Hemoglobin increased", 86,     NA,     65,     "g/L",  "2",      19,
    "Hemoglobin increased", 85,     NA,     65,     "g/L",  "1",      20,
    "Hemoglobin increased", 66,     NA,     65,     "g/L",  "1",      21,
    "Hemoglobin increased", 65,     NA,     65,     "g/L",  "0",      22,
    "Hemoglobin increased", NA,     NA,     65,     "g/L",  NA,       23,
    # Unit missing cannot grade
    "Hemoglobin increased", 200,    61,     65,     NA,     NA,       24,
    # ANRHI missing - cannot grade
    "Hemoglobin increased", 200,    60,     NA,     "g/L",  NA,       25,
    # AVAL missing cannot grade
    "Hemoglobin increased", NA,     60,     65,     "g/L",  NA,       26,
  )
  input_ctcv4_15 <- exp_out_ctcv4_15 %>%
    select(-ATOXGRH)

  actual_output_ctcv4_15 <- derive_var_atoxgr_dir(
    input_ctcv4_15,
    new_var = ATOXGRH,
    meta_criteria = atoxgr_criteria_ctcv4,
    tox_description_var = ATOXDSCH,
    criteria_direction = "H",
    get_unit_expr = AVALU
  )

  expect_dfs_equal(
    base = exp_out_ctcv4_15,
    compare = actual_output_ctcv4_15,
    keys = c("TESTNUM")
  )
})

### 16. INR increased ----
### Grade 3: >2.5 x ULN; >2.5 times above baseline if on anticoagulation
### Grade 2: >1.5 - 2.5 x ULN; >1.5 - 2.5 times above baseline if on anticoagulation
### Grade 1: >1 - 1.5 x ULN; >1 - 1.5 times above baseline if on anticoagulation

test_that("derive_var_atoxgr_dir: Test 16 NCICTCAEv4 INR increased", {
  exp_out_ctcv4_16 <- tibble::tribble(
    ~ATOXDSCH,       ~AVAL, ~BASE, ~ANRHI, ~AVALU,        ~ATOXGRH, ~TESTNUM,
    "Not a term",    80,    120,   200,    NA_character_, NA,       1,
    NA_character_,   60,    50,    100,    NA_character_, NA,       2,
    # GRADE derived from AVAL against ANRHI
    "INR IncreaSed", 251,   200,   100,    NA_character_, "3",      3,
    "INR Increased", 250,   199,   100,    NA_character_, "2",      4,
    "INR Increased", 151,   150,   100,    NA_character_, "2",      5,
    "INR Increased", 150,   150,   100,    NA_character_, "1",      6,
    "INR Increased", 101,   150,   100,    NA_character_, "1",      7,
    "INR Increased", 100,   100,   100,    NA_character_, "0",      8,
    # GRADE derived from AVAL against BASE
    "INR IncreaSed", 251,   100,   200,    NA_character_, "3",      9,
    "INR Increased", 250,   100,   199,    NA_character_, "2",      10,
    "INR Increased", 151,   100,   150,    NA_character_, "2",      11,
    "INR Increased", 150,   100,   150,    NA_character_, "1",      12,
    "INR Increased", 101,   100,   150,    NA_character_, "1",      13,
    "INR Increased", 100,   100,   100,    NA_character_, "0",      14,
    # BASE missing - AVAL <= ANRLO cannot grade as NORMAL
    "INR Increased", 100,   NA,    100,    NA_character_, NA,       15,
    # ANRHI missing - AVAL <= BASE cannot grade as NORMAL
    "INR Increased", 100,   100,   NA,     NA_character_, NA,       16,
    # AVAL missing cannot grade
    "INR Increased", NA,    100,   100,    NA_character_, NA,       17,
  )
  input_ctcv4_16 <- exp_out_ctcv4_16 %>%
    select(-ATOXGRH)

  actual_output_ctcv4_16 <- derive_var_atoxgr_dir(
    input_ctcv4_16,
    new_var = ATOXGRH,
    meta_criteria = atoxgr_criteria_ctcv4,
    tox_description_var = ATOXDSCH,
    criteria_direction = "H",
    get_unit_expr = AVALU
  )

  expect_dfs_equal(
    base = exp_out_ctcv4_16,
    compare = actual_output_ctcv4_16,
    keys = c("TESTNUM")
  )
})

### 17. Lipase increased ----
### Grade 4: >5.0 x ULN
### Grade 3: >2.0 - 5.0 x ULN
### Grade 2: >1.5 - 2.0 x ULN
### Grade 1: >ULN - 1.5 x ULN

test_that("derive_var_atoxgr_dir: Test 17 NCICTCAEv4 Lipase increased", {
  exp_out_ctcv4_17 <- tibble::tribble(
    ~ATOXDSCH,          ~AVAL,  ~ANRLO, ~ANRHI, ~AVALU,         ~ATOXGRH,
    "Not a term",       80,     120,    200,    NA_character_,  NA,
    NA_character_,      60,     50,     100,    NA_character_,  NA,
    "Lipase IncreaSed", 501,    0,      100,    NA_character_,  "4",
    "Lipase Increased", 500,    0,      100,    NA_character_,  "3",
    "Lipase Increased", 201,    0,      100,    NA_character_,  "3",
    "Lipase Increased", 200,    0,      100,    NA_character_,  "2",
    "Lipase Increased", 151,    0,      100,    NA_character_,  "2",
    "Lipase Increased", 150,    0,      100,    NA_character_,  "1",
    "Lipase Increased", 101,    0,      100,    NA_character_,  "1",
    "Lipase Increased", 100,    0,      100,    NA_character_,  "0",
    # ANRHI missing cannot grade
    "Lipase Increased", 200,    0,      NA,     NA_character_,  NA,
    # AVAL missing cannot grade
    "Lipase Increased", NA,     0,      100,    NA_character_,  NA,
  )
  input_ctcv4_17 <- exp_out_ctcv4_17 %>%
    select(-ATOXGRH)

  actual_output_ctcv4_17 <- derive_var_atoxgr_dir(
    input_ctcv4_17,
    new_var = ATOXGRH,
    meta_criteria = atoxgr_criteria_ctcv4,
    tox_description_var = ATOXDSCH,
    criteria_direction = "H",
    get_unit_expr = AVALU
  )

  expect_dfs_equal(
    base = exp_out_ctcv4_17,
    compare = actual_output_ctcv4_17,
    keys = c("AVAL", "ANRLO", "ANRHI", "AVALU")
  )
})

### 18. Lymphocyte count decreased ----
### Grade 4: <0.2 x 10e9 /L
### Grade 3: <0.5 - 0.2 x 10e9 /L
### Grade 2: <0.8 - 0.5 x 10e9 /L
### Grade 1: <LLN - 0.8 x 10e9/L

test_that("derive_var_atoxgr_dir: Test 18 NCICTCAEv4 Lymphocyte count decreased", {
  exp_out_ctcv4_18 <- tibble::tribble(
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
  input_ctcv4_18 <- exp_out_ctcv4_18 %>%
    select(-ATOXGRL)

  actual_output_ctcv4_18 <- derive_var_atoxgr_dir(
    input_ctcv4_18,
    new_var = ATOXGRL,
    meta_criteria = atoxgr_criteria_ctcv4,
    tox_description_var = ATOXDSCL,
    criteria_direction = "L",
    get_unit_expr = AVALU
  )

  expect_dfs_equal(
    base = exp_out_ctcv4_18,
    compare = actual_output_ctcv4_18,
    keys = c("TESTNUM")
  )
})

### 19. Lymphocyte count increased ----
### Grade 3: >20,000/mm3
### Grade 2: >4000/mm3 - 20,000/mm3
test_that("derive_var_atoxgr_dir: Test 19 NCICTCAEv4 Lymphocyte count increased", {
  exp_out_ctcv4_19 <- tibble::tribble(
    ~ATOXDSCH,                    ~AVAL,  ~ANRLO, ~ANRHI, ~AVALU,    ~ATOXGRH,
    "Not a term",                 80,     120,    200,    "10^9/L",  NA,
    NA_character_,                60,     50,     100,    "10^9/L",  NA,
    "Lymphocyte count increased", 21,     NA,     NA,     "10^9/L",  "3",
    "Lymphocyte count increased", 20,     NA,     NA,     "10^9/L",  "2",
    "Lymphocyte count increased", 4.1,    NA,     NA,     "10^9/L",  "2",
    "Lymphocyte count increased", 4,      NA,     NA,     "10^9/L",  "0",
    # Unit missing cannot grade
    "Lymphocyte count increased", 4,      NA,     NA,     NA,        NA,
    # AVAL missing cannot grade
    "Lymphocyte count increased", NA,     NA,     NA,     "10^9/L",  NA,
  )
  input_ctcv4_19 <- exp_out_ctcv4_19 %>%
    select(-ATOXGRH)

  actual_output_ctcv4_19 <- derive_var_atoxgr_dir(
    input_ctcv4_19,
    new_var = ATOXGRH,
    meta_criteria = atoxgr_criteria_ctcv4,
    tox_description_var = ATOXDSCH,
    criteria_direction = "H",
    get_unit_expr = AVALU
  )

  expect_dfs_equal(
    base = exp_out_ctcv4_19,
    compare = actual_output_ctcv4_19,
    keys = c("AVAL", "ANRLO", "ANRHI", "AVALU")
  )
})

### 20. Neutrophil count decreased ----
### Grade 4: <25.0 x 10e9 /L
### Grade 3: <1.0 - 0.5 x 10e9 /L
### Grade 2: <1.5 - 1.0 x 10e9 /L
### Grade 1: <LLN - 1.5 x 10e9 /L

test_that("derive_var_atoxgr_dir: Test 20 NCICTCAEv4 Neutrophil count decreased", {
  exp_out_ctcv4_20 <- tibble::tribble(
    ~ATOXDSCL,                    ~AVAL,  ~ANRLO, ~ANRHI, ~AVALU,    ~ATOXGRL,
    "Not a term",                 80,     120,    200,    "10^9/L",  NA,
    NA_character_,                60,     50,     100,    "10^9/L",  NA,
    # ANRLO not missing
    "Neutrophil count decreased", 0.49,   2,      NA,     "10^9/L",  "4",
    "Neutrophil count decreased", 0.5,    2,      NA,     "10^9/L",  "3",
    "Neutrophil count decreased", 0.99,   2,      NA,     "10^9/L",  "3",
    "Neutrophil count decreased", 1,      2,      NA,     "10^9/L",  "2",
    "Neutrophil count decreased", 1.49,   2,      NA,     "10^9/L",  "2",
    "Neutrophil count decreased", 1.5,    2,      NA,     "10^9/L",  "1",
    "Neutrophil count decreased", 1.9,    2,      NA,     "10^9/L",  "1",
    "Neutrophil count decreased", 2,      2,      NA,     "10^9/L",  "0",
    # ANRLO missing - can grade 2-4
    "Neutrophil count decreased", 0.49,   NA,     NA,     "10^9/L",  "4",
    "Neutrophil count decreased", 0.5,    NA,     NA,     "10^9/L",  "3",
    "Neutrophil count decreased", 0.99,   NA,     NA,     "10^9/L",  "3",
    "Neutrophil count decreased", 1,      NA,     NA,     "10^9/L",  "2",
    "Neutrophil count decreased", 1.49,   NA,     NA,     "10^9/L",  "2",
    # ANRLO missing - can NOT grade 0 or 1
    "Neutrophil count decreased", 1.5,    NA,     NA,     "10^9/L",  NA,
    "Neutrophil count decreased", 1.9,    NA,     NA,     "10^9/L",  NA,
    "Neutrophil count decreased", 2,      NA,     NA,     "10^9/L",  NA,
    # Unit missing cannot grade
    "Neutrophil count decreased", 2,      2,      NA,     NA,        NA,
    # AVAL missing cannot grade
    "Neutrophil count decreased", NA,     2,      NA,     "10^9/L",  NA,
  )
  input_ctcv4_20 <- exp_out_ctcv4_20 %>%
    select(-ATOXGRL)

  actual_output_ctcv4_20 <- derive_var_atoxgr_dir(
    input_ctcv4_20,
    new_var = ATOXGRL,
    meta_criteria = atoxgr_criteria_ctcv4,
    tox_description_var = ATOXDSCL,
    criteria_direction = "L",
    get_unit_expr = AVALU
  )

  expect_dfs_equal(
    base = exp_out_ctcv4_20,
    compare = actual_output_ctcv4_20,
    keys = c("AVAL", "ANRLO", "ANRHI", "AVALU")
  )
})

### 21. Platelet count decreased ----
### Grade 4: <25.0 x 10e9 /L
### Grade 3: <50.0 - 25.0 x 10e9 /L
### Grade 2: <75.0 - 50.0 x 10e9 /L
### Grade 1: <LLN - 75.0 x 10e9 /L

test_that("derive_var_atoxgr_dir: Test 21 NCICTCAEv4 Platelet count decreased", {
  exp_out_ctcv4_21 <- tibble::tribble(
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
  input_ctcv4_21 <- exp_out_ctcv4_21 %>%
    select(-ATOXGRL)

  actual_output_ctcv4_21 <- derive_var_atoxgr_dir(
    input_ctcv4_21,
    new_var = ATOXGRL,
    meta_criteria = atoxgr_criteria_ctcv4,
    tox_description_var = ATOXDSCL,
    criteria_direction = "L",
    get_unit_expr = AVALU
  )

  expect_dfs_equal(
    base = exp_out_ctcv4_21,
    compare = actual_output_ctcv4_21,
    keys = c("TESTNUM")
  )
})

### 22. Serum amylase increased ----
### Grade 4: >5.0 x ULN
### Grade 3: >2.0 - 5.0 x ULN
### Grade 2: >1.5 - 2.0 x ULN
### Grade 1: >ULN - 1.5 x ULN

test_that("derive_var_atoxgr_dir: Test 22 NCICTCAEv4 Serum amylase increased", {
  exp_out_ctcv4_22 <- tibble::tribble(
    ~ATOXDSCH,                 ~AVAL,  ~ANRLO, ~ANRHI, ~AVALU,         ~ATOXGRH,
    "Not a term",              80,     120,    200,    NA_character_,  NA,
    NA_character_,             60,     50,     100,    NA_character_,  NA,
    "Serum amylase increased", 501,    0,      100,    NA_character_,  "4",
    "Serum amylase increased", 500,    0,      100,    NA_character_,  "3",
    "Serum amylase increased", 201,    0,      100,    NA_character_,  "3",
    "Serum amylase increased", 200,    0,      100,    NA_character_,  "2",
    "Serum amylase increased", 151,    0,      100,    NA_character_,  "2",
    "Serum amylase increased", 150,    0,      100,    NA_character_,  "1",
    "Serum amylase increased", 101,    0,      100,    NA_character_,  "1",
    "Serum amylase increased", 100,    0,      100,    NA_character_,  "0",
    # ANRHI missing cannot grade
    "Serum amylase increased", 200,    0,      NA,     NA_character_,  NA,
    # AVAL missing cannot grade
    "Serum amylase increased", NA,     0,      100,    NA_character_,  NA,
  )
  input_ctcv4_22 <- exp_out_ctcv4_22 %>%
    select(-ATOXGRH)

  actual_output_ctcv4_22 <- derive_var_atoxgr_dir(
    input_ctcv4_22,
    new_var = ATOXGRH,
    meta_criteria = atoxgr_criteria_ctcv4,
    tox_description_var = ATOXDSCH,
    criteria_direction = "H",
    get_unit_expr = AVALU
  )

  expect_dfs_equal(
    base = exp_out_ctcv4_22,
    compare = actual_output_ctcv4_22,
    keys = c("AVAL", "ANRLO", "ANRHI", "AVALU")
  )
})

### 23. White blood cell decreased ----
### Grade 4: <1.0 x 10e9 /L
### Grade 3: <2.0 - 1.0 x 10e9 /L
### Grade 2: <3.0 - 2.0 x 10e9 /L
### Grade 1: <LLN - 3.0 x 10e9 /L

test_that("derive_var_atoxgr_dir: Test 23 NCICTCAEv4 White blood cell decreased", {
  exp_out_ctcv4_23 <- tibble::tribble(
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
  input_ctcv4_23 <- exp_out_ctcv4_23 %>%
    select(-ATOXGRL)

  actual_output_ctcv4_23 <- derive_var_atoxgr_dir(
    input_ctcv4_23,
    new_var = ATOXGRL,
    meta_criteria = atoxgr_criteria_ctcv4,
    tox_description_var = ATOXDSCL,
    criteria_direction = "L",
    get_unit_expr = AVALU
  )

  expect_dfs_equal(
    base = exp_out_ctcv4_23,
    compare = actual_output_ctcv4_23,
    keys = c("TESTNUM")
  )
})

# Metabolism and nutrition disorders ----

### 24. Hypercalcemia ----
### Grade 4: >3.4 mmol/L
### Grade 3: >3.1 - 3.4 mmol/L
### Grade 2: >2.9 - 3.1 mmol/L
### Grade 1: >ULN - 2.9 mmol/L

test_that("derive_var_atoxgr_dir: Test 24 NCICTCAEv4 Hypercalcemia", {
  exp_out_ctcv4_24 <- tibble::tribble(
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
  input_ctcv4_24 <- exp_out_ctcv4_24 %>%
    select(-ATOXGRH)

  actual_output_ctcv4_24 <- derive_var_atoxgr_dir(
    input_ctcv4_24,
    new_var = ATOXGRH,
    meta_criteria = atoxgr_criteria_ctcv4,
    tox_description_var = ATOXDSCH,
    criteria_direction = "H",
    get_unit_expr = AVALU
  )

  expect_dfs_equal(
    base = exp_out_ctcv4_24,
    compare = actual_output_ctcv4_24,
    keys = c("TESTNUM")
  )
})

### 25. Hypercalcemia (Ionized) ----
### Grade 4: >1.8 mmol/L
### Grade 3: >1.6 - 1.8 mmol/L
### Grade 2: >1.5 - 1.6 mmol/L
### Grade 1: >ULN - 1.5 mmol/L

test_that("derive_var_atoxgr_dir: Test 25 NCICTCAEv4 Hypercalcemia (Ionized)", {
  exp_out_ctcv4_25 <- tibble::tribble(
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
    # Unit missing cannot grade       1,
    "Hypercalcemia (Ionized)", 1.3,    0,      1.3,    NA,        NA,       19,
    # AVAL missing cannot grade
    "Hypercalcemia (Ionized)", NA,     0,      1.3,    "mmol/L",  NA,       20,
  )
  input_ctcv4_25 <- exp_out_ctcv4_25 %>%
    select(-ATOXGRH)

  actual_output_ctcv4_25 <- derive_var_atoxgr_dir(
    input_ctcv4_25,
    new_var = ATOXGRH,
    meta_criteria = atoxgr_criteria_ctcv4,
    tox_description_var = ATOXDSCH,
    criteria_direction = "H",
    get_unit_expr = AVALU
  )

  expect_dfs_equal(
    base = exp_out_ctcv4_25,
    compare = actual_output_ctcv4_25,
    keys = c("TESTNUM")
  )
})

### 26. Hyperglycemia (Fasting) ----
### Grade 4: >27.8 mmol/L
### Grade 3: >13.9 - 27.8 mmol/L
### Grade 2: >8.9 - 13.9 mmol/L
### Grade 1: >ULN - 8.9 mmol/L

test_that("derive_var_atoxgr_dir: Test 26 NCICTCAEv4 Hyperglycemia (Fasting)", {
  exp_out_ctcv4_26 <- tibble::tribble(
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
  input_ctcv4_26 <- exp_out_ctcv4_26 %>%
    select(-ATOXGRH)

  actual_output_ctcv4_26 <- derive_var_atoxgr_dir(
    input_ctcv4_26,
    new_var = ATOXGRH,
    meta_criteria = atoxgr_criteria_ctcv4,
    tox_description_var = ATOXDSCH,
    criteria_direction = "H",
    get_unit_expr = AVALU
  )

  expect_dfs_equal(
    base = exp_out_ctcv4_26,
    compare = actual_output_ctcv4_26,
    keys = c("TESTNUM")
  )
})

### 27. Hyperglycemia ----
### Grade 4: >27.8 mmol/L
### Grade 3: >13.9 - 27.8 mmol/L

test_that("derive_var_atoxgr_dir: Test 27 NCICTCAEv4 Hyperglycemia", {
  exp_out_ctcv4_27 <- tibble::tribble(
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
  input_ctcv4_27 <- exp_out_ctcv4_27 %>%
    select(-ATOXGRH)

  actual_output_ctcv4_27 <- derive_var_atoxgr_dir(
    input_ctcv4_27,
    new_var = ATOXGRH,
    meta_criteria = atoxgr_criteria_ctcv4,
    tox_description_var = ATOXDSCH,
    criteria_direction = "H",
    get_unit_expr = AVALU
  )

  expect_dfs_equal(
    base = exp_out_ctcv4_27,
    compare = actual_output_ctcv4_27,
    keys = c("TESTNUM")
  )
})

### 28. Hyperkalemia ----
### Grade 4: >7.0 mmol/L
### Grade 3: >6.0 - 7.0 mmol/L
### Grade 2: >5.5 - 6.0 mmol/L
### Grade 1: >ULN - 5.5 mmol/L

test_that("derive_var_atoxgr_dir: Test 28 NCICTCAEv4 Hyperkalemia", {
  exp_out_ctcv4_28 <- tibble::tribble(
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
  input_ctcv4_28 <- exp_out_ctcv4_28 %>%
    select(-ATOXGRH)

  actual_output_ctcv4_28 <- derive_var_atoxgr_dir(
    input_ctcv4_28,
    new_var = ATOXGRH,
    meta_criteria = atoxgr_criteria_ctcv4,
    tox_description_var = ATOXDSCH,
    criteria_direction = "H",
    get_unit_expr = AVALU
  )

  expect_dfs_equal(
    base = exp_out_ctcv4_28,
    compare = actual_output_ctcv4_28,
    keys = c("TESTNUM")
  )
})

### 29. Hypermagnesemia ----
### Grade 4: >3.30 mmol/L
### Grade 3: >1.23 - 3.30 mmol/L
### Grade 1: >ULN - 1.23 mmol/L

test_that("derive_var_atoxgr_dir: Test 29 NCICTCAEv4 Hypermagnesemia", {
  exp_out_ctcv4_29 <- tibble::tribble(
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
  input_ctcv4_29 <- exp_out_ctcv4_29 %>%
    select(-ATOXGRH)

  actual_output_ctcv4_29 <- derive_var_atoxgr_dir(
    input_ctcv4_29,
    new_var = ATOXGRH,
    meta_criteria = atoxgr_criteria_ctcv4,
    tox_description_var = ATOXDSCH,
    criteria_direction = "H",
    get_unit_expr = AVALU
  )

  expect_dfs_equal(
    base = exp_out_ctcv4_29,
    compare = actual_output_ctcv4_29,
    keys = c("TESTNUM")
  )
})

### 30. Hypernatremia ----
### Grade 4: >160 mmol/L
### Grade 3: >155 - 160 mmol/L
### Grade 2: >150 - 155 mmol/L
### Grade 1: >ULN - 150 mmol/L

test_that("derive_var_atoxgr_dir: Test 30 NCICTCAEv4 Hypernatremia", {
  exp_out_ctcv4_30 <- tibble::tribble(
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
  input_ctcv4_30 <- exp_out_ctcv4_30 %>%
    select(-ATOXGRH)

  actual_output_ctcv4_30 <- derive_var_atoxgr_dir(
    input_ctcv4_30,
    new_var = ATOXGRH,
    meta_criteria = atoxgr_criteria_ctcv4,
    tox_description_var = ATOXDSCH,
    criteria_direction = "H",
    get_unit_expr = AVALU
  )

  expect_dfs_equal(
    base = exp_out_ctcv4_30,
    compare = actual_output_ctcv4_30,
    keys = c("TESTNUM")
  )
})

### 31. Hypertriglyceridemia ----
### Grade 4: >11.4 mmol/L
### Grade 3: >5.7 mmol/L - 11.4 mmol/L
### Grade 2: >3.42 mmol/L - 5.7 mmol/L
### Grade 1: 1.71 mmol/L - 3.42 mmol/L

test_that("derive_var_atoxgr_dir: Test 31 NCICTCAEv4 Hypertriglyceridemia", {
  exp_out_ctcv4_31 <- tibble::tribble(
    ~ATOXDSCH,               ~AVAL,  ~ANRLO, ~ANRHI, ~AVALU,    ~ATOXGRH,
    "Not a term",            11.5,   0,      2.1,    "mmol/L",  NA,
    NA_character_,           11.5,   0,      2.1,    "mmol/L",  NA,
    "Hypertriglyceridemia",  11.5,   0,      2.1,    "mmol/L",  "4",
    "Hypertriglyceridemia",  11.4,   0,      2.1,    "mmol/L",  "3",
    "Hypertriglyceridemia",  5.8,    0,      2.1,    "mmol/L",  "3",
    "Hypertriglyceridemia",  5.7,    0,      2.1,    "mmol/L",  "2",
    "Hypertriglyceridemia",  3.43,   0,      2.1,    "mmol/L",  "2",
    "Hypertriglyceridemia",  3.42,   0,      2.1,    "mmol/L",  "1",
    "Hypertriglyceridemia",  1.72,   0,      2.1,    "mmol/L",  "1",
    "Hypertriglyceridemia",  1.71,   0,      2.1,    "mmol/L",  "0",
    # Unit missing cannot grade
    "Hypertriglyceridemia",  1.71,   0,      2.1,    NA,        NA,
    # AVAL missing cannot grade
    "Hypertriglyceridemia",  NA,     0,      2.1,    "mmol/L",  NA,
  )
  input_ctcv4_31 <- exp_out_ctcv4_31 %>%
    select(-ATOXGRH)

  actual_output_ctcv4_31 <- derive_var_atoxgr_dir(
    input_ctcv4_31,
    new_var = ATOXGRH,
    meta_criteria = atoxgr_criteria_ctcv4,
    tox_description_var = ATOXDSCH,
    criteria_direction = "H",
    get_unit_expr = AVALU
  )

  expect_dfs_equal(
    base = exp_out_ctcv4_31,
    compare = actual_output_ctcv4_31,
    keys = c("ATOXDSCH", "AVAL", "AVALU")
  )
})

### 32. Hyperuricemia ----
### Grade 4: >0.59 mmol/L;
### Grade 3: >ULN - 10 mg/dL (0.59 mmol/L)

test_that("derive_var_atoxgr_dir: Test 32 NCICTCAEv4 Hyperuricemia", {
  exp_out_ctcv4_32 <- tibble::tribble(
    ~ATOXDSCH,        ~AVAL,  ~ANRLO, ~ANRHI, ~AVALU,    ~ATOXGRH,
    "Not a term",     591,    0,      200,    "umol/L",  NA,
    NA_character_,    591,    0,      200,    "umol/L",  NA,
    # ANRHI not missing
    "Hyperuricemia",  591,    0,      200,    "umol/L",  "4",
    "Hyperuricemia",  590,    0,      200,    "umol/L",  "3",
    "Hyperuricemia",  201,    0,      200,    "umol/L",  "3",
    "Hyperuricemia",  200,    0,      200,    "umol/L",  "0",
    # ANRHI missing - can grade 4
    "Hyperuricemia",  591,    0,      NA,     "umol/L",  "4",
    # ANRHI missing - can NOT grade 0 or 3
    "Hyperuricemia",  590,    0,      NA,     "umol/L",  NA,
    "Hyperuricemia",  201,    0,      NA,     "umol/L",  NA,
    "Hyperuricemia",  200,    0,      NA,     "umol/L",  NA,
    # Unit missing cannot grade
    "Hyperuricemia",  200,    0,      200,    NA,        NA,
    # AVAL missing cannot grade
    "Hyperuricemia",  NA,     0,      200,    "umol/L",  NA,
  )
  input_ctcv4_32 <- exp_out_ctcv4_32 %>%
    select(-ATOXGRH)


  actual_output_ctcv4_32 <- derive_var_atoxgr_dir(
    input_ctcv4_32,
    new_var = ATOXGRH,
    meta_criteria = atoxgr_criteria_ctcv4,
    tox_description_var = ATOXDSCH,
    criteria_direction = "H",
    get_unit_expr = AVALU
  )

  expect_dfs_equal(
    base = exp_out_ctcv4_32,
    compare = actual_output_ctcv4_32,
    keys = c("ATOXDSCH", "AVAL", "ANRHI", "AVALU")
  )
})

### 33. Hypoalbuminemia ----
### Grade 3: <20 g/L
### Grade 2: <30 - 20 g/L
### Grade 1: <LLN - 30 g/L

test_that("derive_var_atoxgr_dir: Test 33 NCICTCAEv4 Hypoalbuminemia", {
  exp_out_ctcv4_33 <- tibble::tribble(
    ~ATOXDSCL,          ~AVAL,  ~ANRLO, ~ANRHI, ~AVALU, ~ATOXGRL,
    "Not a term",       19,     40,     100,    "G/L",  NA,
    NA_character_,      19,     40,     100,    "G/L",  NA,
    # ANRLO not missing
    "Hypoalbuminemia",  19,     40,     100,    "G/L",  "3",
    "Hypoalbuminemia",  20,     40,     100,    "G/L",  "2",
    "Hypoalbuminemia",  29,     40,     100,    "G/L",  "2",
    "Hypoalbuminemia",  30,     40,     100,    "G/L",  "1",
    "Hypoalbuminemia",  39,     40,     100,    "G/L",  "1",
    "Hypoalbuminemia",  40,     40,     100,    "G/L",  "0",
    "Hypoalbuminemia",  40,     40,     NA,     "G/L",  "0",
    # ANRLO missing - can grade 2-3
    "Hypoalbuminemia",  19,     NA,     100,    "G/L",  "3",
    "Hypoalbuminemia",  20,     NA,     100,    "G/L",  "2",
    "Hypoalbuminemia",  29,     NA,     100,    "G/L",  "2",
    # ANRLO missing - can NOT grade 0 or 1
    "Hypoalbuminemia",  30,     NA,     100,    "G/L",  NA,
    "Hypoalbuminemia",  39,     NA,     100,    "G/L",  NA,
    "Hypoalbuminemia",  40,     NA,     100,    "G/L",  NA,
    # Unit missing cannot grade
    "Hypoalbuminemia",  40,     40,     100,    NA,     NA,
    # AVAL missing cannot grade
    "Hypoalbuminemia",  NA,     40,     100,    "G/L",  NA,
  )
  input_ctcv4_33 <- exp_out_ctcv4_33 %>%
    select(-ATOXGRL)

  actual_output_ctcv4_33 <- derive_var_atoxgr_dir(
    input_ctcv4_33,
    new_var = ATOXGRL,
    meta_criteria = atoxgr_criteria_ctcv4,
    tox_description_var = ATOXDSCL,
    criteria_direction = "L",
    get_unit_expr = AVALU
  )

  expect_dfs_equal(
    base = exp_out_ctcv4_33,
    compare = actual_output_ctcv4_33,
    keys = c("ATOXDSCL", "AVAL", "ANRLO", "ANRHI", "AVALU")
  )
})


### 34. Hypocalcemia ----
### Grade 4: <1.5 mmol/L
### Grade 3: <1.75 - 1.5 mmol/L
### Grade 2: <2.0 - 1.75 mmol/L
### Grade 1: <LLN - 2.0 mmol/L
test_that("derive_var_atoxgr_dir: Test 34 NCICTCAEv4 Hypocalcemia", {
  exp_out_ctcv4_34 <- tibble::tribble(
    ~ATOXDSCL,       ~AVAL,  ~ANRLO, ~ANRHI, ~AVALU,    ~ATOXGRL,
    "Not a term",    1.4,    4,      100,    "mmol/L",  NA,
    NA_character_,   1.4,    4,      100,    "mmol/L",  NA,
    # ANRLO not missing
    "Hypocalcemia",  1.4,    4,      100,    "mmol/L",  "4",
    "Hypocalcemia",  1.5,    4,      100,    "mmol/L",  "3",
    "Hypocalcemia",  1.74,   4,      100,    "mmol/L",  "3",
    "Hypocalcemia",  1.75,   4,      100,    "mmol/L",  "2",
    "Hypocalcemia",  1.9,    4,      100,    "mmol/L",  "2",
    "Hypocalcemia",  2,      4,      100,    "mmol/L",  "1",
    "Hypocalcemia",  3,      4,      100,    "mmol/L",  "1",
    "Hypocalcemia",  4,      4,      100,    "mmol/L",  "0",
    # ANRLO missing - can grade 2-4
    "Hypocalcemia",  1.4,    4,      NA,     "mmol/L",  "4",
    "Hypocalcemia",  1.5,    4,      NA,     "mmol/L",  "3",
    "Hypocalcemia",  1.74,   4,      NA,     "mmol/L",  "3",
    "Hypocalcemia",  1.75,   4,      NA,     "mmol/L",  "2",
    "Hypocalcemia",  1.9,    4,      NA,     "mmol/L",  "2",
    # ANRLO missing - can NOT grade 0 or 1
    "Hypocalcemia",  2,      4,      NA,     "mmol/L",  "1",
    "Hypocalcemia",  3,      4,      NA,     "mmol/L",  "1",
    "Hypocalcemia",  4,      4,      NA,     "mmol/L",  "0",
    # Unit missing cannot grade
    "Hypocalcemia",  4,      4,      100,    NA,        NA,
    # AVAL missing cannot grade
    "Hypocalcemia",  NA,     4,      100,    "mmol/L",  NA,
  )
  input_ctcv4_34 <- exp_out_ctcv4_34 %>%
    select(-ATOXGRL)

  actual_output_ctcv4_34 <- derive_var_atoxgr_dir(
    input_ctcv4_34,
    new_var = ATOXGRL,
    meta_criteria = atoxgr_criteria_ctcv4,
    tox_description_var = ATOXDSCL,
    criteria_direction = "L",
    get_unit_expr = AVALU
  )

  expect_dfs_equal(
    base = exp_out_ctcv4_34,
    compare = actual_output_ctcv4_34,
    keys = c("ATOXDSCL", "AVAL", "ANRLO", "ANRHI", "AVALU")
  )
})

### 35. Hypocalcemia (Ionized) ----
### Grade 4: <0.8 mmol/L
### Grade 3: <0.9 - 0.8 mmol/L
### Grade 2: <1.0 - 0.9 mmol/L
### Grade 1: <LLN - 1.0 mmol/L
test_that("derive_var_atoxgr_dir: Test 35 NCICTCAEv4 Hypocalcemia (Ionized)", {
  exp_out_ctcv4_35 <- tibble::tribble(
    ~ATOXDSCL,                 ~AVAL,  ~ANRLO, ~ANRHI, ~AVALU,    ~ATOXGRL,
    "Not a term",              0.79,   1.3,    100,    "mmol/L",  NA,
    NA_character_,             0.79,   1.3,    100,    "mmol/L",  NA,
    # ANRLO not missing
    "Hypocalcemia (Ionized)",  0.79,   1.3,    100,    "mmol/L",  "4",
    "Hypocalcemia (Ionized)",  0.8,    1.3,    100,    "mmol/L",  "3",
    "Hypocalcemia (Ionized)",  0.89,   1.3,    100,    "mmol/L",  "3",
    "Hypocalcemia (Ionized)",  0.9,    1.3,    100,    "mmol/L",  "2",
    "Hypocalcemia (Ionized)",  0.99,   1.3,    100,    "mmol/L",  "2",
    "Hypocalcemia (Ionized)",  1,      1.3,    100,    "mmol/L",  "1",
    "Hypocalcemia (Ionized)",  1.29,   1.3,    100,    "mmol/L",  "1",
    "Hypocalcemia (Ionized)",  1.3,    1.3,    100,    "mmol/L",  "0",
    # ANRLO missing - can grade 2-4
    "Hypocalcemia (Ionized)",  0.79,   NA,     100,    "mmol/L",  "4",
    "Hypocalcemia (Ionized)",  0.8,    NA,     100,    "mmol/L",  "3",
    "Hypocalcemia (Ionized)",  0.89,   NA,     100,    "mmol/L",  "3",
    "Hypocalcemia (Ionized)",  0.9,    NA,     100,    "mmol/L",  "2",
    "Hypocalcemia (Ionized)",  0.99,   NA,     100,    "mmol/L",  "2",
    # ANRLO missing - can NOT grade 0 or 1
    "Hypocalcemia (Ionized)",  1,      1.3,    NA,     "mmol/L",  "1",
    "Hypocalcemia (Ionized)",  1.29,   1.3,    NA,     "mmol/L",  "1",
    "Hypocalcemia (Ionized)",  1.3,    1.3,    NA,     "mmol/L",  "0",
    # Unit missing cannot grade
    "Hypocalcemia (Ionized)",  1.3,    1.3,    100,    NA,        NA,
    # AVAL missing cannot grade
    "Hypocalcemia (Ionized)",  NA,     1.3,    100,    "mmol/L",  NA,
  )
  input_ctcv4_35 <- exp_out_ctcv4_35 %>%
    select(-ATOXGRL)

  actual_output_ctcv4_35 <- derive_var_atoxgr_dir(
    input_ctcv4_35,
    new_var = ATOXGRL,
    meta_criteria = atoxgr_criteria_ctcv4,
    tox_description_var = ATOXDSCL,
    criteria_direction = "L",
    get_unit_expr = AVALU
  )

  expect_dfs_equal(
    base = exp_out_ctcv4_35,
    compare = actual_output_ctcv4_35,
    keys = c("ATOXDSCL", "AVAL", "ANRLO", "ANRHI", "AVALU")
  )
})

### 36. Hypoglycemia ----
### Grade 4: <1.7 mmol/L
### Grade 3: <2.2 - 1.7 mmol/L
### Grade 2: <3.0 - 2.2 mmol/L
### Grade 1: <LLN - 3.0 mmol/L
test_that("derive_var_atoxgr_dir: Test 36 NCICTCAEv4 Hypoglycemia", {
  exp_out_ctcv4_36 <- tibble::tribble(
    ~ATOXDSCL,       ~AVAL,  ~ANRLO, ~ANRHI, ~AVALU,    ~ATOXGRL,
    "Not a term",    1.69,   4,      100,    "mmol/L",  NA,
    NA_character_,   1.69,   4,      100,    "mmol/L",  NA,
    # ANRLO not missing
    "Hypoglycemia",  1.69,   4,      100,    "mmol/L",  "4",
    "Hypoglycemia",  1.7,    4,      100,    "mmol/L",  "3",
    "Hypoglycemia",  2.19,   4,      100,    "mmol/L",  "3",
    "Hypoglycemia",  2.2,    4,      100,    "mmol/L",  "2",
    "Hypoglycemia",  2.9,    4,      100,    "mmol/L",  "2",
    "Hypoglycemia",  3,      4,      100,    "mmol/L",  "1",
    "Hypoglycemia",  3.9,    4,      100,    "mmol/L",  "1",
    "Hypoglycemia",  4,      4,      100,    "mmol/L",  "0",
    # ANRLO missing - can grade 2-4
    "Hypoglycemia",  1.69,   NA,     100,    "mmol/L",  "4",
    "Hypoglycemia",  1.7,    NA,     100,    "mmol/L",  "3",
    "Hypoglycemia",  2.19,   NA,     100,    "mmol/L",  "3",
    "Hypoglycemia",  2.2,    NA,     100,    "mmol/L",  "2",
    "Hypoglycemia",  2.9,    NA,     100,    "mmol/L",  "2",
    # ANRLO missing - can NOT grade 0 or 1
    "Hypoglycemia",  3,      NA,     100,    "mmol/L",  NA,
    "Hypoglycemia",  3.9,    NA,     100,    "mmol/L",  NA,
    "Hypoglycemia",  4,      NA,     100,    "mmol/L",  NA,
    # Unit missing cannot grade
    "Hypoglycemia",  4,      4,      100,    NA,        NA,
    # AVAL missing cannot grade
    "Hypoglycemia",  NA,     4,      100,    "mmol/L",  NA,
  )
  input_ctcv4_36 <- exp_out_ctcv4_36 %>%
    select(-ATOXGRL)

  actual_output_ctcv4_36 <- derive_var_atoxgr_dir(
    input_ctcv4_36,
    new_var = ATOXGRL,
    meta_criteria = atoxgr_criteria_ctcv4,
    tox_description_var = ATOXDSCL,
    criteria_direction = "L",
    get_unit_expr = AVALU
  )

  expect_dfs_equal(
    base = exp_out_ctcv4_36,
    compare = actual_output_ctcv4_36,
    keys = c("ATOXDSCL", "AVAL", "ANRLO", "ANRHI", "AVALU")
  )
})

### 37. Hypokalemia ----
### Grade 4: <2.5 mmol/L
### Grade 3: <3.0 - 2.5 mmol/L
### Grade 2: <LLN - 3.0 mmol/L
test_that("derive_var_atoxgr_dir: Test 37 NCICTCAEv4 Hypokalemia", {
  exp_out_ctcv4_37 <- tibble::tribble(
    ~ATOXDSCL,      ~AVAL,  ~ANRLO, ~ANRHI, ~AVALU,    ~ATOXGRL,
    "Not a term",   2.49,   4,      100,    "mmol/L",  NA,
    NA_character_,  2.49,   4,      100,    "mmol/L",  NA,
    # ANRLO not missing
    "Hypokalemia",  2.49,   4,      100,    "mmol/L",  "4",
    "Hypokalemia",  2.5,    4,      100,    "mmol/L",  "3",
    "Hypokalemia",  2.9,    4,      100,    "mmol/L",  "3",
    "Hypokalemia",  3,      4,      100,    "mmol/L",  "2",
    "Hypokalemia",  3.9,    4,      100,    "mmol/L",  "2",
    "Hypokalemia",  4,      4,      100,    "mmol/L",  "0",
    # ANRLO missing - can grade 3-4
    "Hypokalemia",  2.49,   NA,     100,    "mmol/L",  "4",
    "Hypokalemia",  2.5,    NA,     100,    "mmol/L",  "3",
    "Hypokalemia",  2.9,    NA,     100,    "mmol/L",  "3",
    # ANRLO missing - can NOT grade 0 or 2
    "Hypokalemia",  3,      NA,     100,    "mmol/L",  NA,
    "Hypokalemia",  3.9,    NA,     100,    "mmol/L",  NA,
    "Hypokalemia",  4,      NA,     NA,     "mmol/L",  NA,
    # Unit missing cannot grade
    "Hypokalemia",  4,      4,      100,    NA,        NA,
    # AVAL missing cannot grade
    "Hypokalemia",  NA,     4,      100,    "mmol/L",  NA,
  )
  input_ctcv4_37 <- exp_out_ctcv4_37 %>%
    select(-ATOXGRL)

  actual_output_ctcv4_37 <- derive_var_atoxgr_dir(
    input_ctcv4_37,
    new_var = ATOXGRL,
    meta_criteria = atoxgr_criteria_ctcv4,
    tox_description_var = ATOXDSCL,
    criteria_direction = "L",
    get_unit_expr = AVALU
  )

  expect_dfs_equal(
    base = exp_out_ctcv4_37,
    compare = actual_output_ctcv4_37,
    keys = c("ATOXDSCL", "AVAL", "ANRLO", "ANRHI", "AVALU")
  )
})

### 38. Hypomagnesemia ----
### Grade 4: <0.3 mmol/L
### Grade 3: <0.4 - 0.3 mmol/L
### Grade 2: <0.5 - 0.4 mmol/L
### Grade 1: <LLN - 0.5 mmol/L

test_that("derive_var_atoxgr_dir: Test 38 NCICTCAEv4 Hypomagnesemia", {
  exp_out_ctcv4_38 <- tibble::tribble(
    ~ATOXDSCL,         ~AVAL,  ~ANRLO, ~ANRHI, ~AVALU,    ~ATOXGRL,
    "Not a term",      0.29,   1,      100,    "mmol/L",  NA,
    NA_character_,     0.29,   1,      100,    "mmol/L",  NA,
    # ANRLO not missing
    "Hypomagnesemia",  0.29,   1,      100,    "mmol/L",  "4",
    "Hypomagnesemia",  0.3,    1,      100,    "mmol/L",  "3",
    "Hypomagnesemia",  0.39,   1,      100,    "mmol/L",  "3",
    "Hypomagnesemia",  0.4,    1,      100,    "mmol/L",  "2",
    "Hypomagnesemia",  0.49,   1,      100,    "mmol/L",  "2",
    "Hypomagnesemia",  0.5,    1,      100,    "mmol/L",  "1",
    "Hypomagnesemia",  0.9,    1,      100,    "mmol/L",  "1",
    "Hypomagnesemia",  1,      1,      100,    "mmol/L",  "0",
    # ANRLO missing - can grade 2-4
    "Hypomagnesemia",  0.29,   NA,     100,    "mmol/L",  "4",
    "Hypomagnesemia",  0.3,    NA,     100,    "mmol/L",  "3",
    "Hypomagnesemia",  0.39,   NA,     100,    "mmol/L",  "3",
    "Hypomagnesemia",  0.4,    NA,     100,    "mmol/L",  "2",
    "Hypomagnesemia",  0.49,   NA,     100,    "mmol/L",  "2",
    # ANRLO missing - can NOT grade 0 or 1
    "Hypomagnesemia",  0.5,    NA,     100,    "mmol/L",  NA,
    "Hypomagnesemia",  0.9,    NA,     100,    "mmol/L",  NA,
    "Hypomagnesemia",  1,      NA,     100,    "mmol/L",  NA,
    # Unit missing cannot grade
    "Hypomagnesemia",  1,      1,      100,    NA,        NA,
    # AVAL missing cannot grade
    "Hypomagnesemia",  NA,     1,      100,    "mmol/L",  NA,
  )
  input_ctcv4_38 <- exp_out_ctcv4_38 %>%
    select(-ATOXGRL)

  actual_output_ctcv4_38 <- derive_var_atoxgr_dir(
    input_ctcv4_38,
    new_var = ATOXGRL,
    meta_criteria = atoxgr_criteria_ctcv4,
    tox_description_var = ATOXDSCL,
    criteria_direction = "L",
    get_unit_expr = AVALU
  )

  expect_dfs_equal(
    base = exp_out_ctcv4_38,
    compare = actual_output_ctcv4_38,
    keys = c("ATOXDSCL", "AVAL", "ANRLO", "ANRHI", "AVALU")
  )
})

### 39. Hyponatremia ----
### Grade 4: <120 mmol/L
### Grade 3: <130 - 120 mmol/L
### Grade 1: <LLN - 130 mmol/L

test_that("derive_var_atoxgr_dir: Test 39 NCICTCAEv4 Hyponatremia", {
  exp_out_ctcv4_39 <- tibble::tribble(
    ~ATOXDSCL,       ~AVAL,  ~ANRLO,   ~ANRHI, ~AVALU,    ~ATOXGRL,
    "Not a term",    119,    140,      100,    "mmol/L",  NA,
    NA_character_,   119,    140,      100,    "mmol/L",  NA,
    # ANRLO not missing
    "Hyponatremia",  119,    140,      100,    "mmol/L",  "4",
    "Hyponatremia",  120,    140,      100,    "mmol/L",  "3",
    "Hyponatremia",  129,    140,      100,    "mmol/L",  "3",
    "Hyponatremia",  130,    140,      100,    "mmol/L",  "1",
    "Hyponatremia",  139,    140,      100,    "mmol/L",  "1",
    "Hyponatremia",  140,    140,      100,    "mmol/L",  "0",
    # ANRLO missing - can grade 3-4
    "Hyponatremia",  119,    NA,       100,    "mmol/L",  "4",
    "Hyponatremia",  120,    NA,       100,    "mmol/L",  "3",
    "Hyponatremia",  129,    NA,       100,    "mmol/L",  "3",
    # ANRLO missing - can NOT grade 0 or 1
    "Hyponatremia",  130,    NA,       100,    "mmol/L",  NA,
    "Hyponatremia",  139,    NA,       100,    "mmol/L",  NA,
    "Hyponatremia",  140,    NA,       100,    "mmol/L",  NA,
    # Unit missing cannot grade
    "Hyponatremia",  140,    140,      100,    NA,        NA,
    # AVAL missing cannot grade
    "Hyponatremia",  NA,     140,      100,    "mmol/L",  NA,
  )
  input_ctcv4_39 <- exp_out_ctcv4_39 %>%
    select(-ATOXGRL)

  actual_output_ctcv4_39 <- derive_var_atoxgr_dir(
    input_ctcv4_39,
    new_var = ATOXGRL,
    meta_criteria = atoxgr_criteria_ctcv4,
    tox_description_var = ATOXDSCL,
    criteria_direction = "L",
    get_unit_expr = AVALU
  )

  expect_dfs_equal(
    base = exp_out_ctcv4_39,
    compare = actual_output_ctcv4_39,
    keys = c("ATOXDSCL", "AVAL", "ANRLO", "ANRHI", "AVALU")
  )
})

### 40. Hypophosphatemia ----
### Grade 4: <0.3 mmol/L
### Grade 3: <0.6 - 0.3 mmol/L
### Grade 2: <0.8 - 0.6 mmol/L
### Grade 1: <LLN - 0.8 mmol/L

test_that("derive_var_atoxgr_dir: Test 40 NCICTCAEv4 Hypophosphatemia", {
  exp_out_ctcv4_40 <- tibble::tribble(
    ~ATOXDSCL,           ~AVAL,  ~ANRLO, ~ANRHI, ~AVALU,    ~ATOXGRL,
    "Not a term",        0.29,   1,      100,    "mmol/L",  NA,
    NA_character_,       0.29,   1,      100,    "mmol/L",  NA,
    # ANRLO not missing
    "Hypophosphatemia",  0.29,   1,      100,    "mmol/L",  "4",
    "Hypophosphatemia",  0.3,    1,      100,    "mmol/L",  "3",
    "Hypophosphatemia",  0.59,   1,      100,    "mmol/L",  "3",
    "Hypophosphatemia",  0.6,    1,      100,    "mmol/L",  "2",
    "Hypophosphatemia",  0.79,   1,      100,    "mmol/L",  "2",
    "Hypophosphatemia",  0.8,    1,      100,    "mmol/L",  "1",
    "Hypophosphatemia",  0.9,    1,      100,    "mmol/L",  "1",
    "Hypophosphatemia",  1,      1,      100,    "mmol/L",  "0",
    # ANRLO missing - can grade 3-4
    "Hypophosphatemia",  0.29,   NA,     100,    "mmol/L",  "4",
    "Hypophosphatemia",  0.3,    NA,     100,    "mmol/L",  "3",
    "Hypophosphatemia",  0.59,   NA,     100,    "mmol/L",  "3",
    "Hypophosphatemia",  0.6,    NA,     100,    "mmol/L",  "2",
    "Hypophosphatemia",  0.79,   NA,     100,    "mmol/L",  "2",
    # ANRLO missing - can NOT grade 0 or 1
    "Hypophosphatemia",  0.8,    NA,     100,    "mmol/L",  NA,
    "Hypophosphatemia",  0.9,    NA,     100,    "mmol/L",  NA,
    "Hypophosphatemia",  1,      NA,     100,    "mmol/L",  NA,
    # Unit missing cannot grade
    "Hypophosphatemia",  1,      1,      100,    NA,        NA,
    # AVAL missing cannot grade
    "Hypophosphatemia",  NA,     1,      100,    "mmol/L",  NA,
  )
  input_ctcv4_40 <- exp_out_ctcv4_40 %>%
    select(-ATOXGRL)

  actual_output_ctcv4_40 <- derive_var_atoxgr_dir(
    input_ctcv4_40,
    new_var = ATOXGRL,
    meta_criteria = atoxgr_criteria_ctcv4,
    tox_description_var = ATOXDSCL,
    criteria_direction = "L",
    get_unit_expr = AVALU
  )

  expect_dfs_equal(
    base = exp_out_ctcv4_40,
    compare = actual_output_ctcv4_40,
    keys = c("ATOXDSCL", "AVAL", "ANRLO", "ANRHI", "AVALU")
  )
})
