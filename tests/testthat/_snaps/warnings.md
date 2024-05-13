# warn_if_vars_exist Test 1: warning if a variable already exists in the input dataset

    Code
      warn_if_vars_exist(dm, "AGE")
    Condition
      Warning:
      Variable "AGE" already exists in the dataset.

---

    Code
      warn_if_vars_exist(dm, c("AGE", "AGEU", "ARM"))
    Condition
      Warning:
      Variables "AGE", "AGEU", and "ARM" already exist in the dataset.

---

    Code
      warn_if_vars_exist(dm, c("AAGE", "AGEU", "ARM"))
    Condition
      Warning:
      Variables "AGEU" and "ARM" already exist in the dataset.

# warn_if_invalid_dtc Test 2: Warning if vector contains unknown datetime format

    Code
      warn_if_invalid_dtc(dtc = "20210406T12:30:30")
    Condition
      Warning:
      Dataset contains incorrect datetime format: --DTC may be incorrectly imputed on row(s)
      * Row 1 : --DTC = 20210406T12:30:30
      * ISO representations of the form YYYY-MM-DDThh:mm:ss.ddd are expected, e.g., 2003-12-15T13:15:17.123. Missing parts at the end can be omitted. Missing parts in the middle must be represented by a dash, e.g., 2003---15.

