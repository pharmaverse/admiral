test_that("hello admiral greets without hw", {
  expect_message(
    wlcm_admiral(),
    "^Welcome to the admiral family!\\n"
  )
})
test_that("hello admiral greets with hw", {
  expect_message(
    wlcm_admiral(TRUE),
    "^Welcome to the admiral family!\\n"
  )
})







