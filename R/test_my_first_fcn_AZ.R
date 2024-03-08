library(admiraldev)

test_that("hello admiral without hw",{
  expect_message(
    my_first_fcn(),
    "^welcome to the admiral family!\\n"
  )
})

test_that("hello admiral with hw",{
  expect_message(
    my_first_fcn(TRUE),
    "^welcome to the admiral family!\\n"
  )
})
