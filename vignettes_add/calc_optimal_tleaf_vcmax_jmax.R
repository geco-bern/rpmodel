# .........................................................
# Remove block after debugging
library(dplyr)
library(purrr)
library(tidyr)
library(ggplot2)
source("vignettes_add/eval_environment_chi.R")
source("vignettes_add/eval_environment_jmax.R")
source("vignettes_add/eval_environment_gs_vcmax.R")
source("vignettes_add/QUADP.R")
library(rpmodel)

library(tealeaves)
# .........................................................

calc_optimal_tcleaf_vcmax_jmax <- function(tc_leaf = 25,
                                           patm = 101325,
                                           co2 = 400,
                                           vpd = 1000,
                                           ppfd = 130,
                                           fapar = 1,
                                           kphio = 0.05,
                                           beta = 146,
                                           c_cost = 0.41,
                                           vcmax_start = 20,
                                           gs_start = 0.5,
                                           jmax_start = 40,
                                           method_jmaxlim_inst = "smith37") {
  

  optimise_this_tcleaf_vcmax_jmax <- function( par, args, iabs, kphio, beta, c_cost, method_jmaxlim_inst, maximize=FALSE, return_all=FALSE ){
    
    ## Parameters to be optimized
    vcmax <- par[1]
    gs    <- par[2]
    jmax  <- par[3]
    
    ## Arguments to calculate variables
    tc_leaf <- args[1]
    patm    <- args[2]
    co2     <- args[3]
    vpd     <- args[4]
    
    ## Local variables based on arguments -> TODO: maybe move this to outside?
    kmm       <- calc_kmm(tc_leaf, patm)
    gammastar <- calc_gammastar(tc_leaf, patm)
    ns_star   <- calc_viscosity_h2o(tc_leaf, patm) / calc_viscosity_h2o(25, 101325)
    ca        <- ( 1.0e-6 ) * co2 * patm         # Pa, atms. CO2; done like this here because function missing in current rpmodel build
    vpd       <- vpd
    
    ## Electron transport is limiting
    ## Solve quadratic equation system using: A(Fick's Law) = A(Jmax Limitation)
    ## This leads to a quadratic equation:
    ## A * ci^2 + B * ci + C  = 0
    ## 0 = a + b*x + c*x^2
    
    ## Jmax Limitation following Smith (1937):
    if (method_jmaxlim_inst == "smith37") {
      ## A = gs (ca - ci)
      ## A = kphio * iabs (ci-gammastar)/ci+2*gammastar) * L
      ## L = 1 / sqrt(1 + ((4 * kphio * iabs)/jmax)^2)
      
      ## with
      L <- 1.0 / sqrt(1.0 + ((4.0 * kphio * iabs)/jmax)^2)
      A <- -gs
      B <- gs * ca - 2 * gammastar * gs - L * kphio * iabs
      C <- 2 * gammastar * gs * ca + L * kphio * iabs * gammastar
      
      ci_j <- QUADM(A, B, C)
      a_j  <- kphio * iabs * (ci_j - gammastar)/(ci_j + 2 * gammastar) * L  
    }
    
    ## Jmax Limitation following Farquhar (1989):
    if (method_jmaxlim_inst == "farquhar89") {
      ## A = gs (ca - ci)
      ## A = j/4 * (ci-gammastar)/ci+2*gammastar)
      ## j = (kphio * iabs + jmax - sqrt(( kphio * iabs + jmax)^2 - (4 * kphio * theta * iabs * jmax))) / (2*theta)
      
      ## with
      theta <- 0.85
      j <- (kphio * iabs + jmax - sqrt(( kphio * iabs + jmax)^2 - (4 * kphio * theta * iabs * jmax))) / (2 * theta)
      A <- -gs
      B <- gs * ca - 2 * gammastar * gs - j/4
      C <- 2 * gammastar * gs * ca + gammastar * j/4
      
      ci_j <- ci_j <- QUADM(A, B, C)
      a_j <- j/4 * (ci_j - gammastar)/(ci_j + 2 * gammastar)
    }
    
    ## Rubisco is limiting
    ## Solve Eq. system
    ## A = gs (ca- ci)
    ## A = Vcmax * (ci - gammastar)/(ci + Kmm)
    
    ## This leads to a quadratic equation:
    ## A * ci^2 + B * ci + C  = 0
    ## 0 = a + b*x + c*x^2
    
    ## with
    A <- -1.0 * gs
    B <- gs * ca - gs * kmm - vcmax
    C <- gs * ca * kmm + vcmax * gammastar
    
    ci_c <- QUADM(A, B, C)
    a_c <- vcmax * (ci_c - gammastar) / (ci_c + kmm)
    
    ## Take minimum of the two assimilation rates and maximum of the two ci
    assim <- min( a_j, a_c )
    ci <- max(ci_c, ci_j)
    
    ## only cost ratio is defined. for this here we need absolute values. Set randomly
    cost_transp <- 1.6 * ns_star * gs * vpd
    cost_vcmax  <- beta * vcmax
    cost_jmax   <- c_cost * jmax
    
    ## Option B: This is equivalent to the P-model with its optimization of ci:ca.
    if (assim<=0) {
      net_assim <- -(999999999.9)
    } else {
      net_assim <- -(cost_transp + cost_vcmax + cost_jmax) / assim
    }
    
    if (maximize) net_assim <- -net_assim

    # print(par)
    # print(net_assim)

    if (return_all){
      return(tibble(vcmax_mine = vcmax, jmax_mine = jmax, gs_mine = gs, ci_mine = ci, chi_mine = ci / ca, a_c_mine = a_c, a_j_mine = a_j, assim = assim, ci_c_mine = ci_c, ci_j_mine = ci_j, cost_transp = cost_transp, cost_vcmax = cost_vcmax, cost_jmax = cost_jmax, net_assim = net_assim, method_jmaxlim_inst = method_jmaxlim_inst))
    } else {
      return( net_assim )
    }
  }
  
  out_optim <- optimr::optimr(
    par        = c( vcmax_start, gs_start, jmax_start ), # starting values
    lower      = c( vcmax_start*0.001, gs_start*0.001, jmax_start*0.001 ),
    upper      = c( vcmax_start*1000, gs_start*1000, jmax_start*1000 ),
    fn         = optimise_this_tcleaf_vcmax_jmax,
    args       = c(tc_leaf, patm, co2, vpd),
    iabs       = (ppfd * fapar),
    kphio      = kphio,
    beta       = beta,
    c_cost     = c_cost/4,
    method_jmaxlim_inst = method_jmaxlim_inst,
    method     = "L-BFGS-B",
    maximize   = TRUE,
    control    = list( maxit = 1000 )
  )
  
  varlist <- optimise_this_tcleaf_vcmax_jmax(
    par = out_optim$par,
    args = c(tc_leaf, patm, co2, vpd),
    iabs = (fapar * ppfd),
    kphio,
    beta,
    c_cost / 4,
    method_jmaxlim_inst,
    maximize = FALSE,
    return_all = TRUE
  )

  return(varlist)
}

# .............................................................................
# Additional functions #### 

LeafEnergyBalance <- function(Tleaf = 21.5,  # Input in degC
                              Tair = 20,     # Input in degC
                              gs = 0.30,     # Input in mol/m2/d/Pa
                              PPFD = 130,    # Input in mol/m2/d
                              VPD = 1000,    # Input in Pa
                              Patm = 101325, # Input in Pa
                              Wind = 2,      # Input in m/s
                              Wleaf = 0.02,
                              StomatalRatio = 1, # 2 for amphistomatous
                              LeafAbs = 0.5, # in shortwave range, much less than PAR
                              returnwhat = c("balance", "fluxes")
) {
  
  
  ## Edits to make function runnable within rpmodel
  ## Functions needed in here and not available directly from plantecophy package
  Tk <- function(x){x+273.15}
  
  esat <- function(TdegC, Pa=101){  
    a <- 611.21
    b <- 17.502
    c <- 240.97
    f <- 1.0007 + 3.46 * 10^-8 * Pa * 1000
    esatval <- f * a * (exp(b * TdegC/(c + TdegC)))
    return(esatval)
  }
  
  ## Correct input units to fit calculations below
  PPFD <- PPFD * 10^6 / (3600*24) # mol/m2/d to umol/m2/s
  gs   <- gs * 10^-6 * Patm       # mol/m2/d/Pa to ppm mol/m2/s
  Patm <- Patm / 1000             # Pa to kPa
  VPD  <- VPD  / 1000             # Pa to kPa
  
  
  ## Original function:
  returnwhat <- match.arg(returnwhat)
  
  # Constants
  Boltz <- 5.67 * 10^-8     # w M-2 K-4
  Emissivity <- 0.95        # -
  LatEvap <- 2.54           # MJ kg-1
  CPAIR <- 1010.0           # J kg-1 K-1
  
  H2OLV0 <- 2.501e6         # J kg-1
  H2OMW <- 18e-3            # J kg-1
  AIRMA <- 29.e-3           # mol mass air (kg/mol)
  AIRDENS <- 1.204          # kg m-3
  UMOLPERJ <- 4.57
  DHEAT <- 21.5e-6          # molecular diffusivity for heat
  
  # Density of dry air
  AIRDENS <- Patm*1000/(287.058 * Tk(Tair))
  
  # Latent heat of water vapour at air temperature (J mol-1)
  LHV <- (H2OLV0 - 2.365E3 * Tair) * H2OMW
  
  # Const s in Penman-Monteith equation  (Pa K-1)
  SLOPE <- (esat(Tair + 0.1) - esat(Tair)) / 0.1
  
  # Radiation conductance (mol m-2 s-1)
  Gradiation <- 4.*Boltz*Tk(Tair)^3 * Emissivity / (CPAIR * AIRMA)
  
  # See Leuning et al (1995) PC&E 18:1183-1200 Appendix E
  # Boundary layer conductance for heat - single sided, forced convection
  CMOLAR <- Patm*1000 / (8.314 * Tk(Tair))   # .Rgas() in package...
  Gbhforced <- 0.003 * sqrt(Wind/Wleaf) * CMOLAR
  
  # Free convection
  GRASHOF <- 1.6E8 * abs(Tleaf-Tair) * (Wleaf^3) # Grashof number
  Gbhfree <- 0.5 * DHEAT * (GRASHOF^0.25) / Wleaf * CMOLAR
  
  # Total conductance to heat (both leaf sides)
  Gbh <- 2*(Gbhfree + Gbhforced)
  
  # Heat and radiative conductance
  Gbhr <- Gbh + 2*Gradiation
  
  # Boundary layer conductance for water (mol m-2 s-1)
  Gbw <- StomatalRatio * 1.075 * Gbh  # Leuning 1995
  gw <- gs*Gbw/(gs + Gbw)
  
  # Longwave radiation
  # (positive flux is heat loss from leaf)
  Rlongup <- Emissivity*Boltz*Tk(Tleaf)^4
  
  # Rnet
  Rsol <- 2*PPFD/UMOLPERJ   # W m-2
  Rnet <- LeafAbs*Rsol - Rlongup   # full
  
  # Isothermal net radiation (Leuning et al. 1995, Appendix)
  ea <- esat(Tair) - 1000*VPD
  ema <- 0.642*(ea/Tk(Tair))^(1/7)
  Rnetiso <- LeafAbs*Rsol - (1 - ema)*Boltz*Tk(Tair)^4 # isothermal net radiation
  
  # Isothermal version of the Penmon-Monteith equation
  GAMMA <- CPAIR*AIRMA*Patm*1000/LHV
  ET <- (1/LHV) * (SLOPE * Rnetiso + 1000*VPD * Gbh * CPAIR * AIRMA) / (SLOPE + GAMMA * Gbhr/gw)
  
  # Latent heat loss
  lambdaET <- LHV * ET
  
  # Heat flux calculated using Gradiation (Leuning 1995, Eq. 11)
  Y <- 1/(1 + Gradiation/Gbh)
  H2 <- Y*(Rnetiso - lambdaET)
  
  # Heat flux calculated from leaf-air T difference.
  # (positive flux is heat loss from leaf)
  H <- -CPAIR * AIRDENS * (Gbh/CMOLAR) * (Tair - Tleaf)
  
  # Leaf-air temperature difference recalculated from energy balance.
  # (same equation as above!)
  Tleaf2 <- Tair + H2/(CPAIR * AIRDENS * (Gbh/CMOLAR))
  
  # Difference between input Tleaf and calculated, this will be minimized.
  # EnergyBal <- (Tleaf - Tleaf2)           # OLD, needed to work with uniroot()
  EnergyBal <- (Tleaf - Tleaf2)^2         # NEW, needed to work with optimr()
  # EnergyBal <- abs(Tleaf - Tleaf2)        # NEW, needs more iterations than ()^2
  
  if(returnwhat == "balance"){
    
    return(EnergyBal)                      # OLD
    
    # out <- list(tc_leaf       = Tleaf,     # NEW
    # 						tc_leaf_star  = Tleaf2,    # NEW
    # 						eps           = EnergyBal) # NEW
    # return(out)                            # NEW
  }
  
  if(returnwhat == "fluxes"){
    
    l <- data.frame(ELEAFeb=1000*ET, Gradiation=Gradiation, Rsol=Rsol, Rnetiso=Rnetiso, Rlongup=Rlongup, H=H, lambdaET=lambdaET, gw=gw, Gbh=Gbh, H2=H2, Tleaf2=Tleaf2)
    return(l)
  }
}

calc_tc_leaf_from_tc_air <- function(tc_air   = 25,   # input in degC
                                     gs       = 0.30, # input in mol/m2/d/Pa
                                     vpd      = 1000, # input in Pa
                                     patm     = 101325, # input in Pa
                                     ppfd     = 130,  # input in mol/m2/d
                                     wind     = 2,    # input in m/s
                                     wleaf    = 0.02,
                                     stoma_r  = 1, 
                                     leaf_abs = 0.5) {
  
  
  # LeafEnergyBalance is equivalent to "maximize_this_tc_leaf" with its Tleaf getting optimized
  
    sol_optimr <-	optimr::optimr(
    # Parameter boundaries to optimize within:
    par       = 15,
    lower     = 15 - tc_air, # OLD: 0 
    upper     = 15 + tc_air, # OLD: 40 
    
    # Function to optimize and its inputs:
    fn        = LeafEnergyBalance,
    Tair      = tc_air,    # input in degC
    gs        = gs,        # input in mol/m2/d/Pa
    VPD       = vpd,       # input in Pa
    Patm      = patm,      # input in Pa
    PPFD      = ppfd,      # input in mol/m2/d
    Wind      = wind,      # input in m/s
    Wleaf     = wleaf,
    StomatalRatio = stoma_r,        # 2 for amphistomatous
    LeafAbs   = leaf_abs,
    
    # Optimr settings:
    method    = "L-BFGS-B",
    control   = list( maxit = 100, maximize = TRUE )
    )
    
    return(sol_optimr$par)
}


maximize_this_tc_leaf <- function(tc_leaf   = 25, # This gets optimized in optimr()
                                  tc_air    = 25,
                                  tc_growth = 25,
                                  tc_home   = 25,
                                  ppfd      = 130, # mol/m2/d
                                  fapar     = 1,
                                  co2       = 400,
                                  patm      = 101325,
                                  vpd       = 1000,
                                  kphio     = 0.05,
                                  toggles   = NA,
                                  wind      = 2,
                                  wleaf     = 0.02,
                                  stoma_r   = 1,
                                  leaf_abs  = 0.5,
                                  method_jmaxlim_inst = "smith37",
                                  beta      = 146,
                                  c_cost    = 0.41) {
  
  varlist_opt_tcleaf <- calc_optimal_tcleaf_vcmax_jmax(tc_leaf = tc_leaf,
                                                       patm = patm,
                                                       co2 = co2,
                                                       vpd = vpd,
                                                       ppfd = ppfd,
                                                       fapar = fapar,
                                                       kphio = kphio,
                                                       beta = beta,
                                                       c_cost = c_cost,
                                                       method_jmaxlim_inst = method_jmaxlim_inst,
                                                       vcmax_start = 20,
                                                       gs_start = 0.5,
                                                       jmax_start = 40)
  gs <- varlist_opt_tcleaf$gs_mine
  
  if (T) {
  ## Via plantecophys energy balance
  tc_leaf_x <- calc_tc_leaf_from_tc_air(tc_air = tc_air,
                                        # gs       = 0.15 / 101325,
                                        gs       = gs,
                                        wind     = wind,
                                        wleaf    = wleaf,
                                        stoma_r  = stoma_r,
                                        leaf_abs = leaf_abs)
  
  } else {
    ## Via tealeaves energy balance
    
    # Get leaf parameters
    leaf_par <- make_leafpar(
      replace = list(
        g_sw = set_units(gs, "mol/m^2/d/Pa")))
    
    # Get environmental parameters Next, we'll change the air temperature to 25 degree C (= 298.15 K)
    enviro_par <- make_enviropar(
      replace = list(
        T_air = set_units(tc_air + 273.15, "K")))
    # Get physical constants
    constants  <- make_constants()
    
    # Get tc_leaf
    tc_leaf_x <- tleaf(leaf_par, enviro_par, constants, quiet = TRUE)$T_leaf %>% 
      set_units("degree_Celsius") %>% drop_units()
  }
  
  # Get difference between tc_leaf and tc_leaf_x
  eps <- (tc_leaf_x - tc_leaf)^2
  
  return(eps)
}


calc_tc_leaf_final <- function(tc_air    = 25,
                               tc_growth = 25,
                               tc_home   = 25,
                               ppfd      = 130,
                               fapar     = 1,
                               patm      = 101325,
                               co2       = 400,
                               vpd       = 1000,
                               kphio     = 0.05,
                               toggles   = NA,
                               wind      = 2,
                               wleaf     = 0.02,
                               stoma_r   = 1,
                               leaf_abs  = 0.5,
                               method_jmaxlim_inst = "smith37",
                               beta      = 146, 
                               c_cost    = 0.41
){
  
  # Get optimized tc_leaf

  # Call optimize()

  
  sol_optimize <- tryCatch(
    {
      sol_optimize <- optimize(maximize_this_tc_leaf,
                               interval  = c(1, 40),
                               tc_air    = tc_air,
                               tc_growth = tc_growth,
                               tc_home   = tc_home,
                               ppfd      = ppfd,
                               fapar     = fapar,
                               co2       = co2,
                               patm      = patm,
                               vpd       = vpd,
                               kphio		 = kphio,
                               toggles   = toggles,
                               wind      = wind,
                               wleaf     = wleaf,
                               stoma_r   = stoma_r,
                               leaf_abs  = leaf_abs,
                               method_jmaxlim_inst = method_jmaxlim_inst,
                               beta      = beta,
                               c_cost    = c_cost)
    },
    warning = function(cond){
      # print("There was a warning")
      return(NA)
    },
    error = function(cond){
      # print("This message will not be printed.")
      return(NA)
    },
    finally = {
      #pass
    })
  
  if (length(sol_optimize) == 1) {
    return (tc_air)
  }
  
  return(sol_optimize$minimum)
}

# TEST ZONE -------------------------------------------------------------------

if (F) { # F to avoid calculation when sourcing
  
  # Parameters:
  tc_air    <- 40
  ppfd      <- 130
  co2       <- 400
  patm      <- 101325
  vpd       <- 1000
  kphio		  <- 0.05
  fapar     <- 1
  method_jmaxlim_inst <- "smith37"
  beta      <- 146
  c_cost    <- 0.41
  
  vcmax_start <- 20
  jmax_start  <- 40
  gs_start    <- 0.5

  tc_leaf <- calc_tc_leaf_final(tc_air = tc_air,
                                ppfd = ppfd,
                                fapar = fapar,
                                co2 = co2,
                                patm = patm,
                                vpd = vpd,
                                kphio = kphio,
                                method_jmaxlim_inst = method_jmaxlim_inst,
                                beta = beta,
                                c_cost = c_cost)
  
  varlist_final <- calc_optimal_tcleaf_vcmax_jmax(tc_leaf = tc_leaf,
                                 patm = patm,
                                 co2 = co2,
                                 vpd = vpd,
                                 ppfd = ppfd,
                                 fapar = fapar,
                                 kphio = kphio,
                                 beta = beta,
                                 c_cost = c_cost,
                                 vcmax_start = vcmax_start,
                                 gs_start = gs_start,
                                 jmax_start = jmax_start,
                                 method_jmaxlim_inst = method_jmaxlim_inst)
}
# calc_tc_leaf_final()
