context("test additional functions")

test_that("test dampen function",{
  skip_on_cran()
  
  expect_error(dampen_vec(
    vec = 20 * (sin(pi/(365)))^2 + rnorm(300, mean = 0, sd = 5),
    tau = 40 ))
  
  dampen <- dampen_vec(
    vec = 20 * (sin(pi/(365)))^2 + rnorm(365, mean = 0, sd = 5),
    tau = 40 )
  
  # output must be numeric
  expect_type(dampen, "double")
})

test_that("ftemp_kphio",{
  skip_on_cran()
  
  # output must be numeric
  temp <- ftemp_kphio(20)
  
  # value is double (single value not numeric = multiple values)
  expect_type(temp , "double")
  
  temp <- ftemp_kphio(20, c4 = TRUE)
  expect_type(temp, "double")
})

test_that("soilm",{
  skip_on_cran()

  
  moisture <- soilmstress(
    soilm = 1,
    meanalpha = 1.0,
    apar_soilm = 0.0,
    bpar_soilm = 0.685)
  
  # output must be numeric
  expect_type(moisture, "double")
})

test_that("quadratic functions",{
  skip_on_cran()
  
  expect_warning(
    QUADP(6, 2, 1)
    )
  
  expect_equal(
    QUADP(0, 0, 0),
    0
  )
  
  expect_equal(
    QUADP(NA, 0, 0),
    NA
  )
  
  expect_warning(
    QUADM(6, 2, 1)
    )
  
  expect_equal(
    QUADM(0, 0, 0),
    0
  )
  
  expect_equal(
    QUADM(NA, 0, 0),
    NA
  )
  
  expect_equal(
    QUADM(0, 0, 6),
    0
  )
  
  expect_equal(
    QUADP(0, 0, 6),
    0
  )
  
  expect_equal(
    QUADM(0, 5, 5),
    -1
  )
  
  expect_equal(
    QUADP(0, 5, 5),
    -1
  )
})