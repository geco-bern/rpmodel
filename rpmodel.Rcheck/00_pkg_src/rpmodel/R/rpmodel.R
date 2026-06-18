#' Invokes a P-model function call
#'
#' R implementation of the P-model and its 
#' corollary predictions (Prentice et al., 2014; Han et al., 2017).
#'
#' @param tc Temperature, relevant for photosynthesis (deg C)
#' @param vpd Vapour pressure deficit (Pa)
#' @param co2 Atmospheric CO2 concentration (ppm)
#' @param fapar (Optional) Fraction of absorbed photosynthetically active
#'  radiation (unitless, defaults to \code{NA})
#' @param ppfd Incident photosynthetic photon flux density 
#'  (mol m-2 d-1, defaults to \code{NA}). Note that the units of 
#'  \code{ppfd} (per area and per time) determine the units of outputs 
#'  \code{lue}, \code{gpp}, \code{vcmax}, and \code{rd}. For example, 
#'  if \code{ppfd} is provided in units of mol m-2 month-1, then
#'  respective output variables are returned as per unit months.
#' @param patm Atmospheric pressure (Pa). When provided, overrides
#'  \code{elv}, otherwise \code{patm} is calculated using standard
#'  atmosphere (101325 Pa), corrected for elevation (argument \code{elv}),
#'  using the function \link{calc_patm}.
#' @param elv Elevation above sea-level (m.a.s.l.). Is used only for 
#'  calculating atmospheric pressure (using standard atmosphere (101325 Pa),
#'  corrected for elevation (argument \code{elv}), using the function
#' \link{calc_patm}), if argument \code{patm} is not provided. If argument
#' \code{patm} is provided, \code{elv} is overridden.
#' @param kphio Apparent quantum yield efficiency (unitless). Defaults to
#'  0.081785 for \code{method_jmaxlim="wang17", do_ftemp_kphio=TRUE, 
#'  do_soilmstress=FALSE}, 0.087182 for \code{method_jmaxlim="wang17",
#'  do_ftemp_kphio=TRUE, do_soilmstress=TRUE}, and 0.049977 for 
#'  \code{method_jmaxlim="wang17", do_ftemp_kphio=FALSE, do_soilmstress=FALSE},
#'  corresponding to the empirically fitted value as presented in Stocker et al.
#'  (2019) Geosci. Model Dev. for model setup 'BRC', 'FULL', and 'ORG' 
#'  respectively, corresponding to \eqn{(a_L b_L)/4} in 
#'  Eq.20 in Stocker et al. (2020) for C3 photosynthesis. For C4 photosynthesis
#'  (\code{c4 = TRUE}), \code{kphio} defaults to 1.0, corresponding to the 
#'   parametrisation by  Cai & Prentice (2020).
#' @param beta Unit cost ratio. Defaults to 146.0 (see Stocker et al., 2019) for
#'   C3 plants and 146/9 for C4 plants.
#' @param soilm (Optional, used only if \code{do_soilmstress==TRUE}) Relative 
#'  soil moisture as a fraction of field capacity (unitless). Defaults to 1.0 
#'  (no soil moisture stress). This information is used to calculate
#'  an empirical soil moisture stress factor (\link{calc_soilmstress}) whereby
#'  the sensitivity is determined by average aridity, defined by the local 
#'  annual mean ratio of actual over potential evapotranspiration, supplied by
#'  argument \code{meanalpha}.
#' @param meanalpha (Optional, used only if \code{do_soilmstress==TRUE}) Local 
#'  annual mean ratio of actual over potential evapotranspiration, measure for 
#'  average aridity. Defaults to 1.0. Only scalar numbers are accepted. If 
#'  a vector is provided, only the first element will be used.
#' @param apar_soilm (Optional, used only if \code{do_soilmstress==TRUE}) 
#'  Parameter determining the sensitivity of the empirical soil moisture stress 
#'  function. Defaults to 0.0, the empirically fitted value as presented in 
#'  Stocker et al. (2019) Geosci. Model Dev. for model setup 'FULL' 
#'  (corresponding to a setup with \code{method_jmaxlim="wang17", 
#'  do_ftemp_kphio=TRUE, do_soilmstress=TRUE}).
#' @param bpar_soilm (Optional, used only if \code{do_soilmstress==TRUE}) 
#'  Parameter determining the sensitivity of the empirical soil moisture stress
#'  function. Defaults to 0.7330, the empirically fitted value as presented in 
#'  Stocker et al. (2019) Geosci. Model Dev. for model setup 'FULL' 
#'  (corresponding to a setup with \code{method_jmaxlim="wang17", 
#'  do_ftemp_kphio=TRUE, do_soilmstress=TRUE}).
#' @param c4 (Optional) A logical value specifying whether the C3 or C4
#'  photosynthetic pathway is followed.Defaults to \code{FALSE}. If \code{TRUE},
#'  the leaf-internal CO2 concentration is still estimated using beta but
#'  \eqn{m} (returned variable \code{mj}) tends to 1, and \eqn{m'} tends to
#'  0.669 (with \code{c = 0.41}) to represent CO2 concentrations within the leaf.
#'  With \code{do_ftemp_kphio = TRUE}, a C4-specific temperature dependence of
#'  the quantum yield efficiency is used (see \link{ftemp_kphio}).
#' @param method_jmaxlim (Optional) A character string specifying which method 
#'  is to be used for factoring in Jmax limitation. Defaults to \code{"wang17"},
#'  based on Wang Han et al. 2017 Nature Plants and (Smith 1937). Available is 
#'  also \code{"smith19"}, following the method by Smith et al., 2019 Ecology 
#'  Letters, and \code{"none"} for ignoring effects of Jmax limitation.
#' @param do_ftemp_kphio (Optional) A logical specifying whether 
#'  temperature-dependence of quantum yield efficiency is used. See \link{ftemp_kphio}
#'  for details. Defaults to \code{TRUE}. Only scalar numbers are accepted. If 
#'  a vector is provided, only the first element will be used.
#' @param do_soilmstress (Optional) A logical specifying whether an empirical 
#' soil moisture stress factor is to be applied to down-scale light use 
#' efficiency (and only light use efficiency). Defaults to \code{FALSE}.
#' @param returnvar (Optional) A character string of vector of character strings
#'  specifying which variables are to be returned (see return below).
#' @param verbose Logical, defines whether verbose messages are printed. 
#'  Defaults to \code{FALSE}.
#'
#' @return A named list of numeric values (including temperature and pressure 
#' dependent parameters of the photosynthesis model, P-model predictions, 
#' including all its corollary). This includes :
#' 
#' \itemize{
#'  \item \code{ca}: Ambient CO2 expressed as partial pressure (Pa)
#'  
#'  \item \code{gammastar}: Photorespiratory compensation point \eqn{\Gamma*}, 
#'   (Pa), see \link{calc_gammastar}.
#'  
#'  \item \code{kmm}: Michaelis-Menten coefficient \eqn{K} for photosynthesis 
#'  (Pa), see \link{calc_kmm}.
#'  
#'  \item \code{ns_star}: Change in the viscosity of water, relative to its 
#'   value at 25 deg C (unitless).
#'   \deqn{\eta* = \eta(T) / \eta(25 deg C)}
#'   This is used to scale the unit cost of transpiration. 
#'   Calculated following Huber et al. (2009).
#'  
#'  \item \code{chi}: Optimal ratio of leaf internal to ambient CO2 (unitless). 
#'   Derived following Prentice et al.(2014) as:
#'  \deqn{
#'   \chi = \Gamma* / ca + (1- \Gamma* / ca) \xi / (\xi + \sqrt D )
#'   }
#'   with
#'   \deqn{
#'    \xi = \sqrt (\beta (K+ \Gamma*) / (1.6 \eta*))
#'   }
#'   \eqn{\beta} is given by argument \code{beta}, \eqn{K} is 
#'   \code{kmm} (see \link{calc_kmm}), \eqn{\Gamma*} is 
#'   \code{gammastar} (see \link{calc_gammastar}). \eqn{\eta*} is \code{ns_star}.
#'   \eqn{D} is the vapour pressure deficit (argument \code{vpd}), \eqn{ca} is 
#'   the ambient CO2 partial pressure in Pa (\code{ca}).
#'   
#'   \item \code{ci}: Leaf-internal CO2 partial pressure (Pa), calculated as \eqn{(\chi ca)}.
#'   
#'   \item \code{lue}: Light use efficiency (g C / mol photons), calculated as
#'                         \deqn{
#'                              LUE = \phi(T) \phi0 m' Mc
#'                         }
#'                         where \eqn{\phi(T)} is the temperature-dependent quantum yield efficiency modifier
#'                         (\link{ftemp_kphio}) if \code{do_ftemp_kphio==TRUE}, and 1 otherwise. \eqn{\phi 0}
#'                         is given by argument \code{kphio}.
#'                         \eqn{m'=m} if \code{method_jmaxlim=="none"}, otherwise
#'                         \deqn{
#'                                m' = m \sqrt( 1 - (c/m)^(2/3) )
#'                         }
#'                         with \eqn{c=0.41} (Wang et al., 2017) if \code{method_jmaxlim=="wang17"}. \eqn{Mc} is
#'                         the molecular mass of C (12.0107 g mol-1). \eqn{m} is given returned variable \code{mj}.
#'                         If \code{do_soilmstress==TRUE}, \eqn{LUE} is multiplied with a soil moisture stress factor,
#'                         calculated with \link{calc_soilmstress}.
#'         \item \code{mj}: Factor in the light-limited assimilation rate function, given by
#'                         \deqn{
#'                             m = (ci - \Gamma*) / (ci + 2 \Gamma*)
#'                        }
#'                        where \eqn{\Gamma*} is given by \code{calc_gammastar}.
#'         \item \code{mc}: Factor in the Rubisco-limited assimilation rate function, given by
#'                         \deqn{
#'                             mc = (ci - \Gamma*) / (ci + K)
#'                        }
#'                        where \eqn{K} is given by \code{calc_kmm}.
#'         \item \code{gpp}: Gross primary production (g C m-2), calculated as
#'                        \deqn{
#'                            GPP = Iabs LUE
#'                        }
#'                        where \eqn{Iabs} is given by \code{fapar*ppfd} (arguments), and is
#'                        \code{NA} if \code{fapar==NA} or \code{ppfd==NA}. Note that \code{gpp} scales with
#'                        absorbed light. Thus, its units depend on the units in which \code{ppfd} is given.
#'         \item \code{iwue}: Intrinsic water use efficiency (iWUE, Pa), calculated as
#'                        \deqn{
#'                              iWUE = ca (1-\chi)/(1.6)
#'                        }
#'         \item \code{gs}: Stomatal conductance (gs, in mol C m-2 Pa-1), calculated as
#'                        \deqn{
#'                             gs = A / (ca (1-\chi))
#'                        }
#'                        where \eqn{A} is \code{gpp}\eqn{/Mc}.
#'         \item \code{vcmax}: Maximum carboxylation capacity \eqn{Vcmax} (mol C m-2) at growth temperature (argument
#'                       \code{tc}), calculated as
#'                       \deqn{
#'                            Vcmax = \phi(T) \phi0 Iabs n
#'                       }
#'                       where \eqn{n} is given by \eqn{n=m'/mc}.
#'         \item \code{vcmax25}: Maximum carboxylation capacity \eqn{Vcmax} (mol C m-2) normalised to 25 deg C
#'                      following a modified Arrhenius equation, calculated as \eqn{Vcmax25 = Vcmax / fv},
#'                      where \eqn{fv} is the instantaneous temperature response by Vcmax and is implemented
#'                      by function \link{ftemp_inst_vcmax}.
#'         \item \code{jmax}: The maximum rate of RuBP regeneration () at growth temperature (argument
#'                       \code{tc}), calculated using
#'                       \deqn{
#'                            A_J = A_C
#'                       }
#'         \item \code{rd}: Dark respiration \eqn{Rd} (mol C m-2), calculated as
#'                      \deqn{
#'                          Rd = b0 Vcmax (fr / fv)
#'                      }
#'                      where \eqn{b0} is a constant and set to 0.015 (Atkin et al., 2015), \eqn{fv} is the
#'                      instantaneous temperature response by Vcmax and is implemented by function
#'                      \link{ftemp_inst_vcmax}, and \eqn{fr} is the instantaneous temperature response
#'                      of dark respiration following Heskel et al. (2016) and is implemented by function
#'                      \link{ftemp_inst_rd}.
#' }
#'
#' Additional variables are contained in the returned list if argument \code{method_jmaxlim=="smith19"}
#' \itemize{
#'  \item \code{omega}: Term corresponding to \eqn{\omega}, defined by Eq. 16 in
#'   Smith et al. (2019), and Eq. E19 in Stocker et al. (2019).
#'   
#'  \item \code{omega_star}: Term corresponding to \eqn{\omega^\ast}, defined by
#'   Eq. 18 in Smith et al. (2019), and Eq. E21 in Stocker et al. (2019).
#'  }patm
#'
#' @references  
#'  Bernacchi, C. J., Pimentel, C., and Long, S. P.:  In vivo temperature response func-tions  of  parameters
#'  required  to  model  RuBP-limited  photosynthesis,  Plant  Cell Environ., 26, 1419–1430, 2003
#'
#   Cai, W., and Prentice, I. C.: Recent trends in gross primary production 
#'  and their drivers: analysis and modelling at flux-site and global scales,
#'  Environ. Res. Lett. 15 124050 https://doi.org/10.1088/1748-9326/abc64e, 2020
#
#'  Heskel,  M.,  O’Sullivan,  O.,  Reich,  P.,  Tjoelker,  M.,  Weerasinghe,  L.,  Penillard,  A.,Egerton, J.,
#'  Creek, D., Bloomfield, K., Xiang, J., Sinca, F., Stangl, Z., Martinez-De La Torre, A., Griffin, K.,
#'  Huntingford, C., Hurry, V., Meir, P., Turnbull, M.,and Atkin, O.:  Convergence in the temperature response
#'  of leaf respiration across biomes and plant functional types, Proceedings of the National Academy of Sciences,
#'  113,  3832–3837,  doi:10.1073/pnas.1520282113,2016.
#'
#'  Huber,  M.  L.,  Perkins,  R.  A.,  Laesecke,  A.,  Friend,  D.  G.,  Sengers,  J.  V.,  Assael,M. J.,
#'  Metaxa, I. N., Vogel, E., Mares, R., and Miyagawa, K.:  New international formulation for the viscosity
#'  of H2O, Journal of Physical and Chemical ReferenceData, 38, 101–125, 2009
#'
#'  Prentice,  I. C.,  Dong,  N.,  Gleason,  S. M.,  Maire,  V.,  and Wright,  I. J.:  Balancing the costs
#'  of carbon gain and water transport:  testing a new theoretical frameworkfor  plant  functional  ecology,
#'  Ecology  Letters,  17,  82–91,  10.1111/ele.12211,http://dx.doi.org/10.1111/ele.12211, 2014.
#'
#'  Wang, H., Prentice, I. C., Keenan, T. F., Davis, T. W., Wright, I. J., Cornwell, W. K.,Evans, B. J.,
#'  and Peng, C.:  Towards a universal model for carbon dioxide uptake by plants, Nat Plants, 3, 734–741, 2017.
#'  Atkin, O. K., et al.:  Global variability in leaf respiration in relation to climate, plant func-tional
#'  types and leaf traits, New Phytologist, 206, 614–636, doi:10.1111/nph.13253,
#'  https://nph.onlinelibrary.wiley.com/doi/abs/10.1111/nph.13253.
#'
#'  Smith, N. G., Keenan, T. F., Colin Prentice, I. , Wang, H. , Wright, I. J., Niinemets, U. , Crous, K. Y.,
#'  Domingues, T. F., Guerrieri, R. , Yoko Ishida, F. , Kattge, J. , Kruger, E. L., Maire, V. , Rogers, A. ,
#'  Serbin, S. P., Tarvainen, L. , Togashi, H. F., Townsend, P. A., Wang, M. , Weerasinghe, L. K. and Zhou, S.
#'  (2019), Global photosynthetic capacity is optimized to the environment. Ecol Lett, 22: 506-517.
#'  doi:10.1111/ele.13210
#'
#'  Stocker, B. et al. Geoscientific Model Development Discussions (in prep.)
#'
#' @export
#'
#' @examples \dontrun{
#'  rpmodel(
#'   tc = 20,
#'   vpd = 1000,
#'   co2 = 400,
#'   ppfd = 30,
#'   elv = 0)
#' }
#'

rpmodel <- function(
  tc,
  vpd,
  co2,
  fapar,
  ppfd,
  patm = NA,
  elv = NA,
  kphio = ifelse(c4, 1.0,
                 ifelse(do_ftemp_kphio,
                        ifelse(do_soilmstress,
                               0.087182,
                               0.081785),
                        0.049977)),
  beta = ifelse(c4, 146/9, 146),
  soilm = stopifnot(!do_soilmstress),
  meanalpha = 1.0,
  apar_soilm = 0.0,
  bpar_soilm = 0.73300,
  c4 = FALSE,
  method_jmaxlim = "wang17",
  do_ftemp_kphio = TRUE,
  do_soilmstress = FALSE,
  returnvar = NULL,
  verbose = FALSE 
  ){

  # Check arguments
  if (identical(NA, elv) && identical(NA, patm)){
    stop(
    "Aborted. Provide either elevation (arugment elv) or 
     atmospheric pressure (argument patm)."
    )
  } else if (!identical(NA, elv) && identical(NA, patm)){
    if (verbose) {
      warning(
      "Atmospheric pressure (patm) not provided. Calculating it as a
      function of elevation (elv), assuming standard atmosphere 
      (101325 Pa at sea level)."
      )
    }
    
    patm <- calc_patm(elv)
  }

  #---- Fixed parameters--------------------------------------------------------
  c_molmass <- 12.0107  # molecular mass of carbon (g)
  kPo <- 101325.0       # standard atmosphere, Pa (Allen, 1973)
  kTo <- 25.0           # base temperature, deg C (Prentice, unpublished)
  rd_to_vcmax <- 0.015  # Ratio of Rdark to Vcmax25, number from Atkin et al., 2015 for C3 herbaceous

  #---- Temperature dependence of quantum yield efficiency----------------------
  ## 'do_ftemp_kphio' is not actually a stress function, but is the temperature-dependency of
  ## the quantum yield efficiency after Bernacchi et al., 2003 PCE
  if (length(do_ftemp_kphio) > 1){
    warning("Argument 'do_ftemp_kphio' has length > 1. Only the first element is used.")
    do_ftemp_kphio <- do_ftemp_kphio[1]
  }
  if (do_ftemp_kphio) {
    kphio <- ftemp_kphio( tc, c4 ) * kphio
  } else {
    if (c4){
      kphio <- ftemp_kphio( 15.0, c4 ) * kphio
    }
  }

  #---- soil moisture stress as a function of soil moisture and mean alpha -----
  if (do_soilmstress) {
    if (length(meanalpha) > 1){
      warning("Argument 'meanalpha' has length > 1. Only the first element is used.")
      meanalpha <- meanalpha[1]
    }
    soilmstress <- calc_soilmstress( soilm, meanalpha, apar_soilm, bpar_soilm )
  }
  else {
    soilmstress <- 1.0
  }

  #---- Photosynthesis parameters depending on temperature, pressure, and CO2. -
  ## ambient CO2 partial pression (Pa)
  ca <- co2_to_ca( co2, patm )

  ## photorespiratory compensation point - Gamma-star (Pa)
  gammastar <- calc_gammastar( tc, patm )

  ## Michaelis-Menten coef. (Pa)
  kmm <- calc_kmm( tc, patm )

  ## viscosity correction factor = viscosity( temp, press )/viscosity( 25 degC, 1013.25 Pa)
  ns      <- viscosity_h2o( tc, patm )  # Pa sc4, 1.0,
  ns25    <- viscosity_h2o( kTo, kPo )  # Pa s
  ns_star <- ns / ns25  # (unitless)

  ##----Optimal ci -------------------------------------------------------------
  ## The heart of the P-model: calculate ci:ca ratio (chi) and additional terms
  out_optchi <- optimal_chi( kmm, gammastar, ns_star, ca, vpd, beta, c4 )
  
  ## leaf-internal CO2 partial pressure (Pa)
  ci <- out_optchi$chi * ca

  #---- Corrolary preditions ---------------------------------------------------
  ## intrinsic water use efficiency (in Pa)
  iwue = ( ca - ci ) / 1.6

  #---- Vcmax and light use efficiency -----------------------------------------
  # Jmax limitation comes in only at this step
  if (c4){

    out_lue_vcmax <- lue_vcmax_c4(
      kphio,
      c_molmass,
      soilmstress
      )

  } else if (method_jmaxlim=="wang17"){

    ## apply correction by Jmax limitation
    out_lue_vcmax <- lue_vcmax_wang17(
      out_optchi,
      kphio,
      c_molmass,
      soilmstress
      )

  } else if (method_jmaxlim=="smith19"){

    out_lue_vcmax <- lue_vcmax_smith19(
      out_optchi,
      kphio,
      c_molmass,
      soilmstress
      )

  } else if (method_jmaxlim=="none"){

    out_lue_vcmax <- lue_vcmax_none(
      out_optchi,
      kphio,
      c_molmass,
      soilmstress
      )

  } else {

    stop("rpmodel(): argument method_jmaxlim not idetified.")

  }

  #---- Corrolary preditions ---------------------------------------------------
  # Vcmax25 (vcmax normalized to 25 deg C)
  ftemp25_inst_vcmax  <- ftemp_inst_vcmax( tc, tc, tcref = 25.0 )
  vcmax25_unitiabs  <- out_lue_vcmax$vcmax_unitiabs / ftemp25_inst_vcmax

  ## Dark respiration at growth temperature
  ftemp_inst_rd <- ftemp_inst_rd( tc )
  rd_unitiabs  <- rd_to_vcmax * (ftemp_inst_rd / ftemp25_inst_vcmax) * out_lue_vcmax$vcmax_unitiabs

  #---- Quantities that scale linearly with absorbed light ---------------------
  iabs <- fapar * ppfd

  # Gross Primary Productivity
  gpp <- iabs * out_lue_vcmax$lue   # in g C m-2 s-1

  # Vcmax per unit ground area is the product of the intrinsic quantum
  # efficiency, the absorbed PAR, and 'n'
  vcmax <- iabs * out_lue_vcmax$vcmax_unitiabs

  ## (vcmax normalized to 25 deg C)
  vcmax25 <- iabs * vcmax25_unitiabs

  ## Dark respiration
  rd <- iabs * rd_unitiabs

  # Jmax using again A_J = A_C, derive the "Jmax limitation factor" 
  fact_jmaxlim <- vcmax * (ci + 2.0 * gammastar) / (kphio * iabs * (ci + kmm))
  
  # use definition of Jmax limitation factor (L in Eq. 13) and solve for Jmax.
  jmax <- 4.0 * kphio * iabs / sqrt( (1.0/fact_jmaxlim)^2 - 1.0 )
  
  # ## Alternatively, Jmax can be calculated from Eq. F10 in Stocker et al., 2020
  # kc <- 0.41
  # jmax_alt <- 4.0 * kphio * iabs * sqrt((out_optchi$mj / kc)^(2/3) - 1.0)
  # fact_jmaxlim_alt <- 1.0 / sqrt(1 + (4.0 * kphio * iabs / jmax_alt)^2)
  
  ftemp25_inst_jmax <- ftemp_inst_jmax( tc, tc, tcref = 25.0 )
  jmax25 <- jmax / ftemp25_inst_jmax

  ## Test: at this stage, verify if A_J = A_C
  if (c4){
    a_j = kphio * iabs * out_optchi$mj * fact_jmaxlim
    a_c = vcmax * out_optchi$mc
  } else {
    a_j <- kphio * iabs * (ci - gammastar)/(ci + 2.0 * gammastar) * fact_jmaxlim
    a_c <- vcmax * (ci - gammastar) / (ci + kmm)
  }
  
  a_j_eq_a_c <- all.equal(a_j, a_c, tol = 0.001)
  if (! isTRUE(a_j_eq_a_c)) {
    warning("rpmodel(): light and Rubisco-limited assimilation rates ",
            "are not identical.\n", a_j_eq_a_c)
  }

  # Assimilation is not returned because it should not be confused with what 
  # is usually measured should use instantaneous assimilation for comparison to
  # measurements. This is returned by inst_rpmodel().
  assim <- ifelse(a_j < a_c , a_j, a_c)
  assim_eq_check <- all.equal(assim, gpp / c_molmass, tol = 0.001)
  if (! isTRUE(assim_eq_check)) {
      warning("rpmodel(): Assimilation and GPP are not identical.\n",
              assim_eq_check)
  }

  ## average stomatal conductance
  gs <- assim / (ca - ci)

  ## construct list for output
  out <- list(
              gpp             = gpp,   # remove this again later
              ca              = ca,
              gammastar       = gammastar,
              kmm             = kmm,
              ns_star         = ns_star,
              chi             = out_optchi$chi,
              xi              = out_optchi$xi,
              mj              = out_optchi$mj,
              mc              = out_optchi$mc,
              ci              = ci,
              iwue            = iwue,
              gs              = gs,
              vcmax           = vcmax,
              vcmax25         = vcmax25,
              jmax            = jmax,
              jmax25          = jmax25,
              rd              = rd
              )

  # if (!is.null(returnvar)) out <- out[returnvar]
  return( out )
}