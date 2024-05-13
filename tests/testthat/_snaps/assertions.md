# assert_filter_cond Test 3: `assert_filter_cond` works as expected

    Code
      assert_filter_cond(arg = fc)
    Condition
      Error:
      ! Argument `fc` must be a filter condition, but is a string

# assert_data_frame Test 4: error if not a dataframe

    Code
      example_fun(c(1, 2, 3))
    Condition
      Error in `example_fun()`:
      ! Argument `dataset` must be class <data.frame>, but is a double vector.

# assert_data_frame Test 7: error if dataframe is grouped

    Code
      example_fun(data)
    Condition
      Error in `example_fun()`:
      ! Argument `dataset` must not be a grouped dataset, please `ungroup()` it.

# assert_data_frame Test 8: error if an expected variable is missing

    Code
      example_fun(data)
    Condition
      Error in `example_fun()`:
      ! Required variable `USUBJID` is missing in `dataset`

# assert_data_frame Test 9: error if expected variables are missing

    Code
      example_fun(data)
    Condition
      Error in `example_fun()`:
      ! Required variables `STUDYID` and `USUBJID` are missing in `dataset`

# assert_character_scalar Test 15: error if `arg` not in values

    Code
      check_unit("month")
    Condition
      Error in `check_unit()`:
      ! Argument `duration_unit` must be equal to one of "years", "months", "weeks", "days", "hours", "minutes", or "seconds".

---

    Code
      check_unit2("month")
    Condition
      Error in `check_unit2()`:
      ! Argument `duration_unit` must be equal to one of "YEARS", "MONTHS", "WEEKS", "DAYS", "HOURS", "MINUTES", or "SECONDS".

# assert_character_scalar Test 16: error if not character

    Code
      example_fun2(2)
    Condition
      Error in `example_fun2()`:
      ! Argument `msg_type` must be a scalar of class <character>, but is a number.

# assert_character_scalar Test 17: error if input is a vector

    Code
      example_fun2(c("admiral", "admiralonco"))
    Condition
      Error in `example_fun2()`:
      ! Argument `msg_type` must be a scalar of class <character>, but is length 2

# assert_character_vector Test 18: error if `arg` not a character vector

    Code
      assert_character_vector(arg)
    Condition
      Error:
      ! Argument `arg` must be <character>, but is a double vector.

# assert_character_vector Test 19: error if `arg` is not in values

    Code
      example_fun(character = c("oak", "mint"))
    Condition
      Error in `example_fun()`:
      ! Argument `character` must be <character> with values "test" and "oak".

# assert_character_vector Test 20: arg_name correctly displayed in name check

    Code
      example_fun(character = c(tree = "oak", "test"))
    Condition
      Error in `example_fun()`:
      ! All elements of `character` argument must be named.
      i The indices of the unnamed elements are 2

# assert_logical_scalar Test 22: error if `arg` is not TRUE or FALSE

    Code
      example_fun("test")
    Condition
      Error in `example_fun()`:
      ! Argument `arg` must be either TRUE or FALSE, but is a string.

# assert_symbol Test 24: `assert_symbol` throws an error if `arg` is missing

    Code
      example_fun(f())
    Condition
      Error in `example_fun()`:
      ! Argument `arg` must be a <symbol>, but is missing.

# assert_symbol Test 25: `assert_symbol` throws an error if `arg` is not a symbol

    Code
      example_fun(f(NULL))
    Condition
      Error in `example_fun()`:
      ! Argument `arg` must be a <symbol>, but is NULL.

# assert_expr Test 29: `assert_expr` throws an error if `arg` is missing

    Code
      assert_expr()
    Condition
      Error:
      ! Argument `arg` cannot be missing.

---

    Code
      example_fun()
    Condition
      Error in `example_fun()`:
      ! Argument `data` cannot be missing.

# assert_expr Test 30: `assert_expr` throws an error if `arg` is not an expression

    Code
      var <- c(1, 2)
      assert_expr(var)
    Condition
      Error:
      ! Argument `var` must be an expression, but is a double vector

# assert_vars Test 32: error if unexpected input

    Code
      assert_vars(AVAL + 1)
    Condition
      Error:
      ! Argument `AVAL + 1` must be a list of <symbol>, e.g., `exprs(USUBJID, VISIT)`.

---

    Code
      assert_vars(rlang::quos(USUBJID, PARAMCD))
    Condition
      Error:
      ! Each element of the list in argument `rlang::quos(USUBJID, PARAMCD)` must be class/type <symbol>.
      i But, element 1 is a <quosure> object, and element 2 is a <quosure> object

---

    Code
      assert_vars(c("USUBJID", "PARAMCD", "VISIT"))
    Condition
      Error:
      ! Argument `c("USUBJID", "PARAMCD", "VISIT")` must be class <list>, but is a character vector.

---

    Code
      assert_vars(exprs(USUBJID, AVAL + 2))
    Condition
      Error:
      ! Each element of the list in argument `exprs(USUBJID, AVAL + 2)` must be class/type <symbol>.
      i But, element 2 is a call

---

    Code
      assert_vars(exprs(APERSDT = APxxSDT, APxxEDT), expect_names = TRUE)
    Condition
      Error:
      ! All elements of `exprs(APERSDT = APxxSDT, APxxEDT)` argument must be named.
      i The indices of the unnamed elements are 2

# assert_vars Test 33: error if some elements of `arg` are not unquoted variable names

    Code
      example_fun(exprs(USUBJID, PARAMCD, NULL))
    Condition
      Error in `example_fun()`:
      ! Each element of the list in argument `arg` must be class/type <symbol>.
      i But, element 3 is NULL

# assert_integer_scalar Test 35: error if chosen subset not in subsets

    Code
      example_fun(1)
    Condition
      Error in `assert_integer_scalar()`:
      ! Argument `subset` must be equal to one of "positive", "non-negative", "negative", or "none".

# assert_integer_scalar Test 37: error if `arg` is not an integer scalar

    Code
      example_fun(1.5)
    Condition
      Error in `example_fun()`:
      ! Argument `arg` must be an integer scalar.

# assert_numeric_vector Test 39: error if `arg` is not a numeric vector

    Code
      example_fun(TRUE)
    Condition
      Error in `example_fun()`:
      ! Argument `arg` must be a numeric vector, but it is `TRUE`.

---

    Code
      example_fun(arg)
    Condition
      Error in `example_fun()`:
      ! Argument `arg` must be a numeric vector, but it is NULL.

---

    Code
      example_fun("1.5")
    Condition
      Error in `example_fun()`:
      ! Argument `arg` must be a numeric vector, but it is a string.

# assert_s3_class Test 40: error if `arg` is not an object of a specific class S3

    Code
      example_fun("test")
    Condition
      Error in `example_fun()`:
      ! Argument `arg` must be class <factor>, but is a string.

# assert_s3_class Test 42: error if `arg` is NULL and optional is FALSE

    Code
      example_fun(NULL)
    Condition
      Error in `example_fun()`:
      ! Argument `arg` must be class <factor>, but is NULL.

# assert_list_of Test 44: error if `arg` is not a list of specific class S3 objects

    Code
      example_fun(list("test"))
    Condition
      Error in `example_fun()`:
      ! Each element of the list in argument `arg` must be class/type <factor>.
      i But, element 1 is a string

# assert_list_of Test 46: error if `arg` is NULL and optional is FALSE

    Code
      example_fun(NULL)
    Condition
      Error in `example_fun()`:
      ! Argument `arg` must be class <list>, but is NULL.

# assert_list_of Test 48: error if `arg` is not a named list (no elements named)

    Code
      mylist <- list(1, 2, 3)
      assert_list_of(mylist, cls = "numeric", named = TRUE)
    Condition
      Error:
      ! All elements of `mylist` argument must be named.
      i The indices of the unnamed elements are 1, 2, and 3

# assert_list_of Test 49: error if `arg` is not a named list (some elements named)

    Code
      mylist <- list(1, 2, 3, d = 4)
      assert_list_of(mylist, cls = "numeric", named = TRUE)
    Condition
      Error:
      ! All elements of `mylist` argument must be named.
      i The indices of the unnamed elements are 1, 2, and 3

# assert_named Test 52: error if no elements are named

    Code
      arg <- c(1, 2)
      assert_named(arg)
    Condition
      Error:
      ! All elements of `arg` argument must be named.
      i The indices of the unnamed elements are 1 and 2

# assert_function Test 56: error if `arg` is not a function

    Code
      example_fun(5)
    Condition
      Error in `example_fun()`:
      ! Argument `arg` must be a function, but is a number.

# assert_function Test 59: error if  `params`  is missing with no default

    Code
      example_fun(sum)
    Condition
      Error in `example_fun()`:
      ! "x" is not an argument of the function specified for `arg`.

---

    Code
      example_fun(sum)
    Condition
      Error in `example_fun()`:
      ! "x" and "y" are not arguments of the function specified for `arg`.

# assert_unit Test 64: error if there are multiple units in the input dataset

    Code
      assert_unit(advs, param = "WEIGHT", required_unit = "kg", get_unit_expr = VSSTRESU)
    Condition
      Error:
      ! Multiple units "kg" and "lb" found for "WEIGHT". Please review and update the units.

# assert_unit Test 65: error if unexpected unit in the input dataset

    Code
      assert_unit(advs, param = "WEIGHT", required_unit = "lb", get_unit_expr = VSSTRESU)
    Condition
      Error:
      ! It is expected that "WEIGHT" has unit of "lb". In the input dataset the unit is "kg".

# assert_param_does_not_exist Test 66: error if parameter exists in the input dataset

    Code
      assert_param_does_not_exist(advs, param = "WEIGHT")
    Condition
      Error:
      ! The parameter code "WEIGHT" already exists in dataset `advs`.

# assert_varval_list Test 68: error if `arg` is not a list of var-value expressions

    Code
      example_fun(c("USUBJID", "PARAMCD", "VISIT"))
    Condition
      Error in `example_fun()`:
      ! Argument `arg` must be a named list of expressions where each element is a symbol, character scalar, numeric scalar, an expression, or NA, but is a character vector.
      i To create a list of expressions use `exprs()`.

# assert_varval_list Test 69: error if `arg` is not a list of var-value expressions

    Code
      example_fun(exprs(USUBJID, PARAMCD, NULL))
    Condition
      Error in `example_fun()`:
      ! Argument `arg` must be a list of expressions where each element is a symbol, character scalar, numeric scalar, an expression, or NA, but is a list.
      i To create a list of expressions use `exprs()`.

# assert_varval_list Test 70: error if `required_elements` are missing from `arg`

    Code
      example_fun(exprs(DTHSEQ = AESEQ))
    Condition
      Error in `example_fun()`:
      ! The following required elements are missing from argument `arg`: "DTHDOM".

# assert_varval_list Test 72: error if `accept_expr` is TRUE and value is invalid

    Code
      example_fun(exprs(DTHSEQ = TRUE))
    Condition
      Error in `example_fun()`:
      ! The elements of the list in argument `arg` must be a symbol, character scalar, numeric scalar, an expression, or NA.
      i "DTHSEQ" = `TRUE` is of type <logical>

# assert_varval_list Test 73: error if `accept_expr` is FALSE and value is invalid

    Code
      example_fun(exprs(DTHSEQ = exprs()))
    Condition
      Error in `example_fun()`:
      ! The elements of the list in argument `arg` must be a symbol, character scalar, numeric scalar, or NA.
      i "DTHSEQ" = `exprs()` is of type <language>

# assert_list_element Test 82: error if the elements do not fulfill the condition

    Code
      assert_list_element(list(list(var = expr(DTHDT), val = 1), list(var = expr(
        EOSDT), val = -1), list(var = expr(EOSDT), val = -2)), element = "val",
      condition = val >= 0, message_text = "List element {.val val} must be `>=0` in argument {.arg {arg_name}}:",
      arg_name = "input")
    Condition
      Error:
      ! List element "val" must be `>=0` in argument `input`:
      i  But, `input[[2]]$val = -1`, and `input[[3]]$val = -2`

# assert_one_to_one Test 83: error if there is a one to many mapping

    Code
      assert_one_to_one(dm, exprs(DOMAIN), exprs(USUBJID))
    Condition
      Error:
      ! For some values of "DOMAIN" there is more than one value of "USUBJID"
      i Call `admiral::get_one_to_many_dataset()` to get all one-to-many values.

# assert_date_var Test 86: error if variable is not a date or datetime variable

    Code
      example_fun(dataset = my_data, var = USUBJID)
    Condition
      Error in `example_fun()`:
      ! Column "USUBJID" in dataset `dataset` must be a date or datetime, but is a character vector.

# assert_date_vector Test 90: error if `arg` is NULL and optional is FALSE

    Code
      example_fun(NULL)
    Condition
      Error in `example_fun()`:
      ! Argument `arg` must be a date or datetime, but is NULL.

# assert_atomic_vector Test 91: error if input is not atomic vector

    Code
      assert_atomic_vector(x)
    Condition
      Error:
      ! Argument `x` must be an atomic vector, but is a list.

# assert_same_type Test 93: error if different type

    Code
      assert_same_type(true_value, false_value, missing_value)
    Condition
      Error:
      ! Arguments `true_value`, `false_value`, and `missing_value` must be the same type.
      i Argument types are `true_value` <character>, `false_value` <character>, `missing_value` <double>

