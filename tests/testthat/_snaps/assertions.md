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

