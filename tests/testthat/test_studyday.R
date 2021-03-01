context('Study Day Derivation')

test_that('Derivation works' ,
          {expect_identical(studyday(ymd("20000101"), ymd("20000101")), 1)
           expect_identical(studyday(ymd("20000101"), ymd_hms("1999-12-01T23:00:00")), -31)
           expect_identical(studyday(ymd_hms("20000101T22:45:11"), ymd("20000106")), 6)})
test_that('Checks work',
          {expect_error(studyday(), 'Argument refdate is not specified.')
           expect_error(studyday(ymd('20200202')), 'Argument date is not specified.')
           expect_error(studyday(ymd('20200202'), 1), 'Argument date=1 is not a lubridate date.')
           expect_error(studyday('hello world', ymd('20200202')), 'Argument refdate=hello world is not a lubridate date.')})