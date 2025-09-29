# create_query_data Test 5: issues error if SMQs without meddra_version are requested

    Code
      create_query_data(queries = list(pregsmq), get_terms_fun = get_smq)
    Condition
      Error in `assert_db_requirements()`:
      ! `version` is not specified. This is expected for baskets.
      i A basket is requested by query 1:
      <query> object
      prefix: "SMQ02"
      name: auto
      add_scope_num: FALSE
      definition:
        <basket_select> object
        name: "Pregnancy and neonatal topics (SMQ)"
        id: NULL
        scope: "NARROW"
        type: "smq"

# create_query_data Test 7: issues error if SDGs without meddra_version are requested

    Code
      create_query_data(queries = list(sdg), get_terms_fun = get_sdg)
    Condition
      Error in `assert_db_requirements()`:
      ! `version` is not specified. This is expected for baskets.
      i A basket is requested by query 1:
      <query> object
      prefix: "SDG01"
      name: auto
      add_scope_num: FALSE
      definition:
        <basket_select> object
        name: "5-aminosalicylates for ulcerative colitis"
        id: NULL
        scope: "NA"
        type: "sdg"

# create_query_data Test 8: error if no `get_terms_fun` provided

    Code
      create_query_data(queries = list(query(prefix = "SMQ02", id = auto, definition = basket_select(
        name = "Pregnancy and neonatal topics (SMQ)", scope = "NARROW", type = "smq"))))
    Condition
      Error in `assert_db_requirements()`:
      ! `get_terms_fun` is not specified. This is expected for baskets.
      i A basket is requested by query 1:
      <query> object
      prefix: "SMQ02"
      name: auto
      id: auto
      add_scope_num: FALSE
      definition:
        <basket_select> object
        name: "Pregnancy and neonatal topics (SMQ)"
        id: NULL
        scope: "NARROW"
        type: "smq"

# create_query_data Test 9: catching error from user function (get_terms_fun)

    Code
      create_query_data(queries = list(pregsmq), version = "20.0", get_terms_fun = faulty_fun)
    Condition
      Error in `value[[3L]]()`:
      ! An error occurred while calling the function `faulty_fun()` provided to the `get_terms_fun` argument.
      This could be due to incorrect handling of input parameters inside `faulty_fun()`.
      Current arguments passed to `faulty_fun()`:
      - version: 20.0
      - basket_select: Pregnancy and neonatal topics (SMQ), NULL, NARROW, and smq
      - keep_id: TRUE
      Error message: Intentional error for testing

# query Test 10: error if name = auto for non SMQs/SDGs

    Code
      sdg <- query(prefix = "CQ01", definition = cqterms)
    Condition
      Error in `validate_query()`:
      ! The auto keyword can be used for baskets only.
      i It was provided for the `name` element.

# query Test 11: error if id = auto for non SMQs/SDGs

    Code
      sdg <- query(name = "My CQ", id = auto, prefix = "CQ01", definition = cqterms)
    Condition
      Error in `validate_query()`:
      ! The auto keyword can be used for baskets only.
      i It was provided for the `id` element.

# query Test 12: error if invalid definition

    Code
      sdg <- query(name = "My CQ", prefix = "CQ01", definition = 1)
    Condition
      Error in `validate_query()`:
      ! `definition` expects a <basket_select> object, a data frame, or a list of data frames and <basket_select> objects.
      i An object of the following class was provided: <numeric>

# validate_query Test 13: error if definition is not a data frame or basket_select

    Code
      validate_query(obj)
    Condition
      Error in `validate_query()`:
      ! Each element of the list in the `definition` field must be a data frame or an object of class <basket_select> but the following are not:
      Element 1 is a string.
      Element 2 is a string.
      Element 3 is a string.

# assert_terms Test 14: error if SRCVAR missing

    Code
      assert_terms(terms = select(cqterms, -SRCVAR), source_text = "my test data")
    Condition
      Error in `assert_terms()`:
      ! Required variable `SRCVAR` is missing in my test data.

# assert_terms Test 15: error if SRCVAR and GRPNAME missing

    Code
      assert_terms(terms = select(cqterms, -SRCVAR), source_text = "my test data",
      expect_grpname = TRUE)
    Condition
      Error in `assert_terms()`:
      ! Required variables `SRCVAR` and `GRPNAME` are missing in my test data.

# assert_terms Test 16: error if TERMCHAR and TERMNUM missing

    Code
      assert_terms(terms = select(cqterms, SRCVAR), source_text = "my test data")
    Condition
      Error in `assert_terms()`:
      ! Variable `TERMCHAR` or `TERMNUM` is required.
      None of them is in my test data.
      i Provided variables: `SRCVAR`

# assert_terms Test 17: error if no data frame

    Code
      assert_terms(terms = 42, source_text = "object returned by calling get_mysmq()")
    Condition
      Error in `assert_terms()`:
      ! object returned by calling get_mysmq() is not a data frame but a function.

# assert_terms Test 18: error if no observations

    Code
      assert_terms(terms = filter(cqterms, TERMNUM == 42), source_text = "object returned by calling get_my_smq")
    Condition
      Error in `assert_terms()`:
      ! object returned by calling get_my_smq does not contain any observations.

# assert_terms Test 19: error if GRPNAME is missing

    Code
      assert_terms(terms = cqterms, expect_grpname = TRUE, source_text = "object returned by calling get_my_smq")
    Condition
      Error in `assert_terms()`:
      ! Required variable `GRPNAME` is missing in object returned by calling get_my_smq.

# assert_terms Test 20: error if GRPID is missing

    Code
      assert_terms(terms = cqterms, expect_grpid = TRUE, source_text = "object returned by calling get_my_smq")
    Condition
      Error in `assert_terms()`:
      ! Required variable `GRPID` is missing in object returned by calling get_my_smq.

# basket_select Test 21: error if name and id specified

    Code
      basket_select(name = "My SMQ", id = 42, scope = "NARROW", type = "smq")
    Condition
      Error in `validate_basket_select()`:
      ! Either `id` or `name` has to be null.

# basket_select Test 22: error if neither name nor id specified

    Code
      basket_select(scope = "NARROW", type = "smq")
    Condition
      Error in `validate_basket_select()`:
      ! Either `id` or `name` has to be non null.

# basket_select Test 23: error if type is not specified

    Code
      basket_select(id = 42, scope = "NARROW")
    Condition
      Error in `basket_select()`:
      ! argument "type" is missing, with no default

# basket_select Test 24: error if arguments inside ... are not named

    Code
      basket_select(name = "Noninfectious meningitis", scope = "NARROW", type = "smq",
        "CHECK 1", "CHECK 3")
    Condition
      Error in `basket_select()`:
      ! All arguments inside `...` must be named

