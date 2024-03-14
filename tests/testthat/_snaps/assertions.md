# assert_character_scalar Test 17: error if `arg` not in values

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

# assert_s3_class Test 41: error if `arg` is NULL and optional is FALSE

    Code
      example_fun(NULL)
    Condition
      Error in `example_fun()`:
      ! Argument `arg` must be class <factor>, but is NULL.

# assert_list_of Test 45: error if `arg` is NULL and optional is FALSE

    Code
      example_fun(NULL)
    Condition
      Error in `example_fun()`:
      ! Argument `arg` must be class <list>, but is NULL.

# assert_list_of Test 47: error if `arg` is not a named list (no elements named)

    Code
      mylist <- list(1, 2, 3)
      assert_list_of(mylist, cls = "numeric", named = TRUE)
    Condition
      Error:
      ! All elements of `mylist` argument must be named.
      i The indices of the unnamed elements are 1, 2, and 3

# assert_list_of Test 48: error if `arg` is not a named list (some elements named)

    Code
      mylist <- list(1, 2, 3, d = 4)
      assert_list_of(mylist, cls = "numeric", named = TRUE)
    Condition
      Error:
      ! All elements of `mylist` argument must be named.
      i The indices of the unnamed elements are 1, 2, and 3

# assert_named Test 51: error if no elements are named

    Code
      arg <- c(1, 2)
      assert_named(arg)
    Condition
      Error:
      ! All elements of `arg` argument must be named.
      i The indices of the unnamed elements are 1 and 2

# assert_same_type Test 91: error if different type

    Code
      assert_same_type(true_value, false_value, missing_value)
    Condition
      Error:
      ! Arguments `true_value`, `false_value`, and `missing_value` must be the same type.
      i Argument types are `true_value` <character>, `false_value` <character>, `missing_value` <double>

