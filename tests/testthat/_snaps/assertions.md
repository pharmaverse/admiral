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

# assert_vars Test 13: error if unexpected input

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

# assert_character_scalar Test 19: error if `arg` not in values

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

# assert_character_vector Test 20: error if `arg` not a character vector

    Code
      assert_character_vector(arg)
    Condition
      Error:
      ! Argument `arg` must be <character>, but is a double vector.

# assert_character_vector Test 21: error if `arg` is not in values

    Code
      example_fun(character = c("oak", "mint"))
    Condition
      Error in `example_fun()`:
      ! Argument `character` must be <character> with values "test" and "oak".

# assert_character_vector Test 22: arg_name correctly displayed in name check

    Code
      example_fun(character = c(tree = "oak", "test"))
    Condition
      Error in `example_fun()`:
      ! All elements of `character` argument must be named.
      i The indices of the unnamed elements are 2

# assert_logical_scalar Test 24: error if `arg` is not TRUE or FALSE

    Code
      example_fun("test")
    Condition
      Error in `example_fun()`:
      ! Argument `arg` must be either TRUE or FALSE, but is a string.

# assert_symbol Test 26: `assert_symbol` throws an error if `arg` is missing

    Code
      example_fun(f())
    Condition
      Error in `example_fun()`:
      ! Argument `arg` must be a <symbol>, but is missing.

# assert_symbol Test 27: `assert_symbol` throws an error if `arg` is not a symbol

    Code
      example_fun(f(NULL))
    Condition
      Error in `example_fun()`:
      ! Argument `arg` must be a <symbol>, but is NULL.

# assert_expr Test 31: `assert_expr` throws an error if `arg` is missing

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

# assert_expr Test 32: `assert_expr` throws an error if `arg` is not an expression

    Code
      var <- c(1, 2)
      assert_expr(var)
    Condition
      Error:
      ! Argument `var` must be an expression, but is a double vector

# assert_s3_class Test 43: error if `arg` is NULL and optional is FALSE

    Code
      example_fun(NULL)
    Condition
      Error in `example_fun()`:
      ! Argument `arg` must be class <factor>, but is NULL.

# assert_list_of Test 47: error if `arg` is NULL and optional is FALSE

    Code
      example_fun(NULL)
    Condition
      Error in `example_fun()`:
      ! Argument `arg` must be class <list>, but is NULL.

# assert_list_of Test 49: error if `arg` is not a named list (no elements named)

    Code
      mylist <- list(1, 2, 3)
      assert_list_of(mylist, cls = "numeric", named = TRUE)
    Condition
      Error:
      ! All elements of `mylist` argument must be named.
      i The indices of the unnamed elements are 1, 2, and 3

# assert_list_of Test 50: error if `arg` is not a named list (some elements named)

    Code
      mylist <- list(1, 2, 3, d = 4)
      assert_list_of(mylist, cls = "numeric", named = TRUE)
    Condition
      Error:
      ! All elements of `mylist` argument must be named.
      i The indices of the unnamed elements are 1, 2, and 3

# assert_named Test 53: error if no elements are named

    Code
      arg <- c(1, 2)
      assert_named(arg)
    Condition
      Error:
      ! All elements of `arg` argument must be named.
      i The indices of the unnamed elements are 1 and 2

# assert_function Test 57: error if `arg` is not a function

    Code
      example_fun(5)
    Condition
      Error in `example_fun()`:
      ! Argument `arg` must be a function, but is a number.

# assert_unit Test 65: error if there are multiple units in the input dataset

    Code
      assert_unit(advs, param = "WEIGHT", required_unit = "kg", get_unit_expr = VSSTRESU)
    Condition
      Error:
      ! Multiple units "kg" and "lb" found for "WEIGHT". Please review and update the units.

# assert_unit Test 66: error if unexpected unit in the input dataset

    Code
      assert_unit(advs, param = "WEIGHT", required_unit = "lb", get_unit_expr = VSSTRESU)
    Condition
      Error:
      ! It is expected that "WEIGHT" has unit of "lb". In the input dataset the unit is "kg".

# assert_param_does_not_exist Test 67: error if parameter exists in the input dataset

    Code
      assert_param_does_not_exist(advs, param = "WEIGHT")
    Condition
      Error:
      ! The parameter code "WEIGHT" already exists in dataset `advs`.

# assert_list_element Test 83: error if the elements do not fulfill the condition

    Code
      assert_list_element(list(list(var = expr(DTHDT), val = 1), list(var = expr(
        EOSDT), val = -1), list(var = expr(EOSDT), val = -2)), element = "val",
      condition = val >= 0, message_text = "List element {.val val} must be `>=0` in argument {.arg {arg_name}}:",
      arg_name = "input")
    Condition
      Error:
      ! List element "val" must be `>=0` in argument `input`:
      i  But, `input[[2]]$val = -1`, and `input[[3]]$val = -2`

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

