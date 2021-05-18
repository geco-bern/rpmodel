context("test inst model runs")

test_that("inst model run",{
  skip_on_cran()
  
  # no atmospheric pressure given
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
  
  # inst call
  inst_out <- inst_rpmodel(
    out_pmodel,
    tc = 20,
    vpd = 1000,
    co2 = 400,
    fapar = 1,
    patm = 1024,
    ppfd = 300,
    kphio = 1)
  
  # output must be a list (no atmosphere warning)
  expect_type(inst_out, "list")
})
