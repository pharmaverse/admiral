# friendly_type_of ----
## Test 1: friendly_type_of() supports objects ----
test_that("friendly_type_of Test 1: friendly_type_of() supports objects", {
  expect_equal(friendly_type_of(mtcars), "a <data.frame> object")
  expect_equal(friendly_type_of(quo(1)), "a <quosure> object")
})

## Test 2: friendly_type_of() supports matrices and arrays (#141) ----
test_that("friendly_type_of Test 2: friendly_type_of() supports matrices and arrays (#141)", {
  expect_equal(friendly_type_of(list()), "an empty list")
  expect_equal(friendly_type_of(matrix(list(1, 2))), "a list matrix")
  expect_equal(friendly_type_of(array(list(1, 2, 3), dim = 1:3)), "a list array")

  expect_equal(friendly_type_of(matrix(1:3)), "an integer matrix")
  expect_equal(friendly_type_of(array(1:3, dim = 1:3)), "an integer array")

  expect_equal(friendly_type_of(matrix(letters)), "a character matrix")
  expect_equal(friendly_type_of(array(letters[1:3], dim = 1:3)), "a character array")
})

## Test 3: friendly_type_of() handles scalars ----
test_that("friendly_type_of Test 3: friendly_type_of() handles scalars", {
  expect_equal(friendly_type_of(NA), "`NA`")

  expect_equal(friendly_type_of(TRUE), "`TRUE`")
  expect_equal(friendly_type_of(FALSE), "`FALSE`")

  expect_equal(friendly_type_of(1L), "an integer")
  expect_equal(friendly_type_of(1.0), "a number")
  expect_equal(friendly_type_of(1i), "a complex number")
  expect_equal(friendly_type_of(as.raw(1)), "a raw value")

  expect_equal(friendly_type_of("foo"), "a string")
  expect_equal(friendly_type_of(""), "`\"\"`")

  expect_equal(friendly_type_of(list(1)), "a list")

  expect_equal(friendly_type_of(matrix(NA)), "a logical matrix")
  expect_equal(friendly_type_of(matrix(1)), "a double matrix")
})

## Test 4: friendly_type_of() edge cases ----
test_that("friendly_type_of Test 4: friendly_type_of() edge cases", {
  expect_equal(friendly_type_of(), "absent")
  expect_equal(friendly_type_of(NULL), "NULL")
  expect_equal(friendly_type_of(1:2, length = TRUE), "an integer vector of length 2")
})

# .rlang_as_friendly_type ----
## Test 5: .rlang_as_friendly_type() works as expected ----
test_that(".rlang_as_friendly_type Test 5: .rlang_as_friendly_type() works as expected", {
  expect_equal(.rlang_as_friendly_type(typeof(list(test = 1:3))), "a list")
  expect_equal(.rlang_as_friendly_type(typeof(NULL)), "NULL")
  expect_equal(.rlang_as_friendly_type(typeof(new.env(parent = emptyenv()))), "an environment")
  expect_equal(
    .rlang_as_friendly_type(
      typeof(xml2::read_xml("<foo><bar /></foo>")$node)
    ),
    "a pointer"
  )

  test_weakref <- new_weakref(new.env(parent = emptyenv()),
    finalizer = function(e) message("finalized")
  )
  expect_equal(.rlang_as_friendly_type(typeof(test_weakref)), "a weak reference")

  # For S4 classes -- think this triggers messages but the expect_equal comes fine
  # setClass("Person", slots = c(name = "character", age = "numeric"))
  # john <- new("Person", name = "John Smith", age = NA_real_)
  # expect_equal(.rlang_as_friendly_type(typeof(john)), "an S4 object")

  # Skip name
  expect_equal(.rlang_as_friendly_type(typeof(sym("test symbol"))), "a symbol")
  expect_equal(.rlang_as_friendly_type(typeof(expr(1 + 1))), "a call")
  expect_equal(.rlang_as_friendly_type(typeof(pairlist(x = 1, y = 2))), "a pairlist node")
  expect_equal(.rlang_as_friendly_type(typeof(expression(x <- 4, x))), "an expression vector")

  # Unsure what char is in compat_friendly_type.R line 118
  # Promise seems impossible because it stops being a promise at evaluation?
  # Unsure how to check for `...`
  # Unsure how to check for `any`
  expect_equal(.rlang_as_friendly_type(typeof(compiler::compile(quote(1 + 3)))), "an internal bytecode object")

  # Skip primitive
  # Skip builtin
  expect_equal(.rlang_as_friendly_type(typeof(switch)), "a primitive function")
  expect_equal(.rlang_as_friendly_type(typeof(mean)), "a function")
})

# .rlang_stop_unexpected_typeof ----
## Test 6: .rlang_stop_unexpected_typeof() works ----
test_that(".rlang_stop_unexpected_typeof Test 6: .rlang_stop_unexpected_typeof() works", {
  expect_error(.rlang_stop_unexpected_typeof("test"), "Unexpected type <character>.")
})

# stop_input_type ----
## Test 7: stop_input_type() works ----
test_that("stop_input_type Test 7: stop_input_type() works", {
  expect_error(stop_input_type(1, what = "character"), "`1` must be character, not a number")
})
