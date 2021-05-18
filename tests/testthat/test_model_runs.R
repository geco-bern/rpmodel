context("test model runs")

test_that("default model run",{
  skip_on_cran()
  
  # no atmospheric pressure given
  expect_warning(rpmodel( 
    tc             = 20,
    vpd            = 1000,
    co2            = 400,
    fapar          = 1,
    ppfd           = 300,
    elv            = 0,
    kphio          = 0.049977,
    beta           = 146,
    c4             = FALSE,
    method_optci   = "prentice14",
    method_jmaxlim = "none",
    do_ftemp_kphio = FALSE,
    do_soilmstress = FALSE,
    verbose        = TRUE
  ))
  
  out_pmodel <- rpmodel( 
    tc             = 20,
    vpd            = 1000,
    co2            = 400,
    fapar          = 1,
    ppfd           = 300,
    elv            = 0,
    kphio          = 0.049977,
    beta           = 146,
    patm = 1024,
    c4             = FALSE,
    method_optci   = "prentice14",
    method_jmaxlim = "none",
    do_ftemp_kphio = FALSE,
    do_soilmstress = FALSE,
    verbose        = TRUE
  )
  
  # output must be a list (no atmosphere warning)
  expect_type(out_pmodel, "list")
})

test_that("jmax wang17",{
  skip_on_cran()
  
  out_pmodel <- rpmodel( 
    tc             = 20,
    vpd            = 1000,
    co2            = 400,
    fapar          = 1,
    ppfd           = 300,
    elv            = 0,
    kphio          = 0.049977,
    beta           = 146,
    patm = 1024,
    c4             = FALSE,
    method_optci   = "prentice14",
    method_jmaxlim = "wang17",
    do_ftemp_kphio = FALSE,
    do_soilmstress = FALSE,
    verbose        = TRUE
  )
  
  # output must be a list (no atmosphere warning)
  expect_type(out_pmodel, "list")
})

test_that("jmax smith19",{
  skip_on_cran()
  
  out_pmodel <- rpmodel( 
    tc             = 20,
    vpd            = 1000,
    co2            = 400,
    fapar          = 1,
    ppfd           = 300,
    elv            = 0,
    kphio          = 0.049977,
    beta           = 146,
    patm = 1024,
    c4             = FALSE,
    method_optci   = "prentice14",
    method_jmaxlim = "smith19",
    do_ftemp_kphio = FALSE,
    do_soilmstress = FALSE,
    verbose        = TRUE
  )
  
  # output must be a list (no atmosphere warning)
  expect_type(out_pmodel, "list")
})
