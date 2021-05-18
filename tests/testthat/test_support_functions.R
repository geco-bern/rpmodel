context("test additional functions")

test_that("test dampen function",{
  skip_on_cran()
  
  dampen <- dampen_vec(
    vec = 20 * (sin(pi/(365)))^2 + rnorm(365, mean = 0, sd = 5),
    tau = 40 )
  
  # output must be numeric
  expect_type(dampen, "numeric")
})

test_that("calc_ftemp_kphio",{
  skip_on_cran()
  
  # output must be numeric
  expect_type(calc_ftemp_kphio(20), "numeric")
  expect_type(calc_ftemp_kphio(20, c4 = TRUE), "numeric")
})

test_that("calc_soilm",{
  skip_on_cran()
  
  # output must be numeric
  expect_type(
    calc_soilmstress(
      soilm = 1,
      meanalpha = 1.0,
      apar_soilm = 0.0,
      bpar_soilm = 0.685)
    , "numeric")
})