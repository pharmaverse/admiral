

#  hello_admiral -------
## Test 1: hello_admiral works ----
test_that("hello_admiral() Test #1  works", {
    expect_no_error(hello_admiral(T))
    expect_no_error(hello_admiral(F))
    expect_no_error(hello_admiral(0))
    expect_no_error(hello_admiral())
    }
)

## Test2: hello_admiral  throws Error if input not acceptable ----
test_that("hello_admiral() Test #2 Error if input not acceptable type", {
   expect_error(hello_admiral(""))
   expect_error(hello_admiral("junk"))
   expect_error(hello_admiral(junk))
 }
)

## Test3: hello_admiral ..... -----
test_that("hello_admiral() Test #3: Return value is `interpreted string literal`", {
   S ="Welcome to Admiral family"
   expect_identical(hello_admiral(), cat(S, "\n"))
}
)


