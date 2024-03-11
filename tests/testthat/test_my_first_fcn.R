test_that("hello admiral without hw",{
  expect_message(
    my_first_fcn(),
    "Welcome to the admiral family!"
  )
})

test_that("hello admiral with hw",{
  expect_message(
    my_first_fcn(TRUE),
    "Welcome to the admiral family!"
  )
})
