# derive_vars_query Test 5: Derive decides between TERMCHAR and TERMNUM based on type

    Code
      derive_vars_query(mutate(my_ae, AELLTCD = as.logical(AELLTCD)), query)
    Condition
      Error in `value[[3L]]()`:
      ! Error in derive_vars_query call of get_vars_query: Error in `get_vars_query()`:
      ! The source variables (values of `SRCVAR`) must be numeric or character.
      i AELLTCD is of type logical

# derive_vars_query Test 6: Error is given if both TERMCHAR/TERMNUM are NA/empty

    Code
      derive_vars_query(my_ae, query)
    Condition
      Error in `value[[3L]]()`:
      ! Error in derive_vars_query call of get_vars_query: Error in `get_vars_query()`:
      ! Either `TERMCHAR` or `TERMNUM` need to be specified in `dataset_queries`. They both cannot be NA or empty.

# assert_valid_queries Test 9: assert_valid_queries checks

    Code
      assert_valid_queries(mutate(query, PREFIX = c("30", "55")), "test")
    Condition
      Error in `assert_valid_queries()`:
      ! `PREFIX` in `test` must start with 2-3 letters. Problem with "30" and "55".

---

    Code
      assert_valid_queries(mutate(query, PREFIX = c("AA", "BB")), "test")
    Condition
      Error in `assert_valid_queries()`:
      ! `PREFIX` in `test` must end with 2-digit numbers. Issue with "AA" and "BB".

---

    Code
      assert_valid_queries(mutate(query, GRPNAME = c("", "A")), "test")
    Condition
      Error in `assert_valid_queries()`:
      ! `GRPNAME` in `test` cannot be empty string or NA.

---

    Code
      assert_valid_queries(mutate(query, GRPID = as.character(GRPID)), "test")
    Condition
      Error in `assert_valid_queries()`:
      ! `GRPID` in `test` must be numeric.

---

    Code
      assert_valid_queries(mutate(query, SCOPE = letters[1:2]), "test")
    Condition
      Error in `assert_valid_queries()`:
      ! `SCOPE` in `test` can only be 'BROAD', 'NARROW' or `NA`.

---

    Code
      assert_valid_queries(mutate(query, SCOPEN = 10:11), "test")
    Condition
      Error in `assert_valid_queries()`:
      ! `SCOPEN` in `test` must be one of 1, 2, or NA. Issue with 10 and 11.

---

    Code
      assert_valid_queries(mutate(query, PREFIX = c("CQ40", "CQ40")), "test")
    Condition
      Error in `assert_valid_queries()`:
      ! In `test` `GRPNAME` of "CQ40" is not unique.

---

    Code
      assert_valid_queries(mutate(query, PREFIX = c("CQ40", "CQ40"), GRPNAME = c(
        "My Query 1", "My Query 1")), "test")
    Condition
      Error in `assert_valid_queries()`:
      ! In `test` `GRPID` of "CQ40" is not unique.

