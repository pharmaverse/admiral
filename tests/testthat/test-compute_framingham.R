test_that("Framingham Equation - Male, Treated for Hypertension = N", {
  expect_equal(
    compute_framingham(
      sysbp = 133,
      chol = 216.16,
      cholhdl = 54.91,
      age=53,
      sex = "M",
      smokefl = "N",
      diabetfl = "N",
      trthypfl = "N"
    ),
    ((1 - (0.88936 ^ exp((3.06117 * log(53))
                       + (1.12370 * log(216.16))
                       - (0.93263 * log(54.91))
                       + (1.93303 * log(133))
                       - 23.9802)
           )) * 100)
    )
})

test_that("Framingham Equation - Male, Treated for Hypertension = Y", {
  expect_equal(
    compute_framingham(
      sysbp = 133,
      chol = 216.16,
      cholhdl = 54.91,
      age=53,
      sex = "M",
      smokefl = "N",
      diabetfl = "N",
      trthypfl = "Y"
    ),
    ((1 - (0.88936 ^ exp((3.06117 * log(53))
                         + (1.12370 * log(216.16))
                         - (0.93263 * log(54.91))
                         + (1.99881 * log(133))
                         - 23.9802)
    )) * 100)
  )
})

test_that("Framingham Equation - Male, Treated for Hypertension = N, Diabetic = Y", {
expect_equal(
  compute_framingham(
    sysbp = 133,
    chol = 216.16,
    cholhdl = 54.91,
    age=53,
    sex = "M",
    smokefl = "N",
    diabetfl = "Y",
    trthypfl = "N"
  ),
  ((1 - (0.88936 ^ exp((3.06117 * log(53))
                       + (1.12370 * log(216.16))
                       - (0.93263 * log(54.91))
                       + (1.93303 * log(133))
                       + 0.57367
                       - 23.9802)
  )) * 100)
)
})

test_that("Framingham Equation - Male, Treated for Hypertension = N, Smoker = Y", {
  expect_equal(
    compute_framingham(
      sysbp = 133,
      chol = 216.16,
      cholhdl = 54.91,
      age=53,
      sex = "M",
      smokefl = "Y",
      diabetfl = "N",
      trthypfl = "N"
    ),
    ((1 - (0.88936 ^ exp((3.06117 * log(53))
                         + (1.12370 * log(216.16))
                         - (0.93263 * log(54.91))
                         + (1.93303 * log(133))
                         + 0.65451
                         - 23.9802)
    )) * 100)
  )
})

test_that("Framingham Equation - FeMale, Treated for Hypertension = N", {
  expect_equal(
    compute_framingham(
      sysbp = 133,
      chol = 216.16,
      cholhdl = 54.91,
      age=53,
      sex = "F",
      smokefl = "N",
      diabetfl = "N",
      trthypfl = "N"
    ),
    ((1 - (0.95012 ^ exp((2.32888 * log(53))
                         + (1.20904 * log(216.16))
                         - (.70833 * log(54.91))
                         + (2.76157 * log(133))
                         - 26.1931)
    )) * 100)
  )
})

test_that("Framingham Equation - FeMale, Treated for Hypertension = Y", {
  expect_equal(
    compute_framingham(
      sysbp = 133,
      chol = 216.16,
      cholhdl = 54.91,
      age=53,
      sex = "F",
      smokefl = "N",
      diabetfl = "N",
      trthypfl = "Y"
    ),
    ((1 - (0.95012 ^ exp((2.32888 * log(53))
                         + (1.20904 * log(216.16))
                         - (.70833 * log(54.91))
                         + (2.82263 * log(133))
                         - 26.1931)
    )) * 100)
  )
})

test_that("Framingham Equation - FeMale, Treated for Hypertension = N, Diabetic = Y", {
  expect_equal(
    compute_framingham(
      sysbp = 133,
      chol = 216.16,
      cholhdl = 54.91,
      age=53,
      sex = "F",
      smokefl = "N",
      diabetfl = "Y",
      trthypfl = "N"
    ),
    ((1 - (0.95012 ^ exp((2.32888 * log(53))
                         + (1.20904 * log(216.16))
                         - (.70833 * log(54.91))
                         + (2.76157 * log(133))
                         + 0.69154
                         - 26.1931)
    )) * 100)
  )
})

test_that("Framingham Equation - FMale, Treated for Hypertension = N, Diabetic = Y", {
  expect_equal(
    compute_framingham(
      sysbp = 133,
      chol = 216.16,
      cholhdl = 54.91,
      age=53,
      sex = "F",
      smokefl = "Y",
      diabetfl = "N",
      trthypfl = "N"
    ),
    ((1 - (0.95012 ^ exp((2.32888 * log(53))
                         + (1.20904 * log(216.16))
                         - (.70833 * log(54.91))
                         + (2.76157 * log(133))
                         + 0.52873
                         - 26.1931)
    )) * 100)
  )
})

