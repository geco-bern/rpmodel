#' Invokes a P-model function call
#'
#' R implementation of the P-model and its corrolary predictions (Prentice et al., 2014; Han et al., 2017).
#'
#' @param tc Temperature, relevant for photosynthesis (deg C)
#' @param vpd Vapour pressure deficit (Pa)
#' @param co2 Atmospheric CO2 concentration (ppm)
#' @param fapar (Optional) Fraction of absorbed photosynthetically active radiation (unitless, defaults to
#' \code{NA})
#' @param ppfd (Optional) Photosynthetic photon flux density (mol m-2 d-1, defaults to \code{NA}). Note that
#' the units of \code{ppfd} (per area and per time) determine the units of outputs \code{lue}, \code{gpp}, 
#' \code{vcmax}, and \code{rd}. For example, if \code{ppfd} is provided in units of mol m-2 month-1, then 
#' respective output variables are returned as per unit months.
#' @param patm Atmospheric pressure (Pa). When provided, overrides \code{elv}, otherwise \code{patm} 
#' is calculated using standard atmosphere (101325 Pa), corrected for elevation (argument \code{elv}), 
#' using the function \link{calc_patm}.
#' @param elv Elevation above sea-level (m.a.s.l.). Is used only for calculating atmospheric pressure (using 
#' standard atmosphere (101325 Pa), corrected for elevation (argument \code{elv}), using the function 
#' \link{calc_patm}), if argument \code{patm} is not provided. If argument \code{patm} is provided, 
#' \code{elv} is overridden. 
#' @param kphio Apparent quantum yield efficiency (unitless). Defaults to 0.0817 for
#' \code{method_jmaxlim="wang17", do_ftemp_kphio=TRUE, do_soilmstress=FALSE}, 0.0870 for 
#' \code{method_jmaxlim="wang17", do_ftemp_kphio=TRUE, do_soilmstress=TRUE}, and 0.0492 for 
#' \code{method_jmaxlim="wang17", do_ftemp_kphio=FALSE, do_soilmstress=FALSE}, corresponding to the empirically 
#' fitted value as presented in Stocker et al. (2019) Geosci. Model Dev. for model setup 'BRC', 'FULL', and 'ORG'
#' respectively. 
#' @param beta Unit cost ratio. Defaults to 146.0 (see Stocker et al., 2019).
#' @param soilm (Optional, used only if \code{do_soilmstress==TRUE}) Relative soil moisture as a fraction 
#' of field capacity (unitless). Defaults to 1.0 (no soil moisture stress). This information is used to calculate 
#' an empirical soil moisture stress factor (\link{calc_soilmstress}) whereby the sensitivity is determined 
#' by average aridity, defined by the local annual mean ratio of actual over potential evapotranspiration, 
#' supplied by argument \code{meanalpha}.
#' @param meanalpha (Optional, used only if \code{do_soilmstress==TRUE}) Local annual mean ratio of 
#' actual over potential evapotranspiration, measure for average aridity. Defaults to 1.0.
#' @param apar_soilm (Optional, used only if \code{do_soilmstress==TRUE}) Parameter determining the 
#' sensitivity of the empirical soil moisture stress function. Defaults to 0.0, the empirically fitted value
#' as presented in Stocker et al. (2019) Geosci. Model Dev. for model setup 'FULL' (corresponding to a setup
#' with \code{method_jmaxlim="wang17", do_ftemp_kphio=TRUE, do_soilmstress=TRUE}).
#' @param bpar_soilm (Optional, used only if \code{do_soilmstress==TRUE}) Parameter determining the 
#' sensitivity of the empirical soil moisture stress function. Defaults to 0.685, the empirically fitted value
#' as presented in Stocker et al. (2019) Geosci. Model Dev. for model setup 'FULL' (corresponding to a setup
#' with \code{method_jmaxlim="wang17", do_ftemp_kphio=TRUE, do_soilmstress=TRUE}).
#' @param c4 (Optional) A logical value specifying whether the C3 or C4 photosynthetic pathway is followed.
#' Defaults to \code{FALSE}. If \code{TRUE}, the leaf-internal CO2 concentration is assumed to be very large
#' and \eqn{m} (returned variable \code{mj}) tends to 1, and \eqn{m'} tends to 0.669 (with \code{c = 0.41}).
#' @param method_optci (Optional) A character string specifying which method is to be used for calculating
#' optimal ci:ca. Defaults to \code{"prentice14"}.
# Available also \code{"prentice14_num"} for a numerical solution to the same optimization criterium as
# used for \code{"prentice14"}.
#' @param method_jmaxlim (Optional) A character string specifying which method is to be used for factoring in
#' Jmax limitation. Defaults to \code{"wang17"},
#' based on Wang Han et al. 2017 Nature Plants and (Smith 1937). Available is also \code{"smith19"}, following
#' the method by Smith et al., 2019 Ecology Letters,
#' and \code{"none"} for ignoring effects of Jmax limitation.
#' @param do_ftemp_kphio (Optional) A logical specifying whether temperature-dependence of quantum yield
#' efficiency after Bernacchi et al., 2003 is to be accounted for. Defaults to \code{TRUE}.
#' @param do_soilmstress (Optional) A logical specifying whether an empirical soil moisture stress factor
#' is to be applied to down-scale light use efficiency (and only light use efficiency). Defaults to \code{FALSE}.
#' @param returnvar (Optional) A character string of vector of character strings specifying which variables
#' are to be returned (see return below).
#' @param verbose Logical, defines whether verbose messages are printed. Defaults to \code{FALSE}.
#'
#' @return A named list of numeric values (including temperature and pressure dependent parameters of the
#' photosynthesis model, P-model predictions, including all its corollary). This includes :
#' \itemize{
#'         \item \code{ca}: Ambient CO2 expressed as partial pressure (Pa)
#'         \item \code{gammastar}: Photorespiratory compensation point \eqn{\Gamma*}, (Pa), see \link{calc_gammastar}.
#'         \item \code{kmm}: Michaelis-Menten coefficient \eqn{K} for photosynthesis (Pa), see \link{calc_kmm}.
#'         \item \code{ns_star}: Change in the viscosity of water, relative to its value at 25 deg C (unitless).
#'                               \deqn{\eta* = \eta(T) / \eta(25 deg C)}
#'                               This is used to scale the unit cost of transpiration. Calculated following Huber
#'                               et al. (2009).
#'         \item \code{chi}: Optimal ratio of leaf internal to ambient CO2 (unitless). Derived following Prentice et al.
#'                          (2014) as:
#'                          \deqn{
#'                                \chi = \Gamma* / ca + (1- \Gamma* / ca) \xi / (\xi + \sqrt D )
#'                          }
#'                          with
#'                          \deqn{
#'                               \xi = \sqrt (\beta (K+ \Gamma*) / (1.6 \eta*))
#'                          }
#'                          \eqn{\beta} is given by argument \code{beta}, \eqn{K} is \code{kmm} (see \link{calc_kmm}),
#'                          \eqn{\Gamma*} is \code{gammastar} (see \link{calc_gammastar}). \eqn{\eta*} is \code{ns_star}.
#'                          \eqn{D} is the vapour pressure deficit (argument \code{vpd}), \eqn{ca} is the
#'                          ambient CO2 partial pressure in Pa (\code{ca}).
#'         \item \code{ci}: Leaf-internal CO2 partial pressure (Pa), calculated as \eqn{(\chi ca)}.
#'         \item \code{lue}: Light use efficiency (g C / mol photons), calculated as
#'                         \deqn{
#'                              LUE = \phi(T) \phi0 m' Mc
#'                         }
#'                         where \eqn{\phi(T)} is the temperature-dependent quantum yield efficiency modifier
#'                         (\link{calc_ftemp_kphio}) if \code{do_ftemp_kphio==TRUE}, and 1 otherwise. \eqn{\phi 0}
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
#'                        where \eqn{\Gamma*} is given by \code{gammastar}.
#'         \item \code{mc}: Factor in the Rubisco-limited assimilation rate function, given by
#'                         \deqn{
#'                             mc = (ci - \Gamma*) / (ci + K)
#'                        }
#'                        where \eqn{K} is given by \code{kmm}.
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
#'                       where \eqn{n} is given by \eqn{n=m/mc}, or
#'                       \deqn{
#'                           n = (ci + K) / (ci + 2 \Gamma*)
#'                       }
#'         \item \code{vcmax25}: Maximum carboxylation capacity \eqn{Vcmax} (mol C m-2) normalised to 25 deg C
#'                      following a modified Arrhenius equation, calculated as \eqn{Vcmax25 = Vcmax / fv},
#'                      where \eqn{fv} is the instantaneous temperature response by Vcmax and is implemented
#'                      by function \link{calc_ftemp_inst_vcmax}.
#'         \item \code{rd}: Dark respiration \eqn{Rd} (mol C m-2), calculated as
#'                      \deqn{
#'                          Rd = b0 Vcmax (fr / fv)
#'                      }
#'                      where \eqn{b0} is a constant and set to 0.015 (Atkin et al., 2015), \eqn{fv} is the
#'                      instantaneous temperature response by Vcmax and is implemented by function
#'                      \link{calc_ftemp_inst_vcmax}, and \eqn{fr} is the instantaneous temperature response
#'                      of dark respiration following Heskel et al. (2016) and is implemented by function
#'                      \link{calc_ftemp_inst_rd}.
#' }
#' 
#' Additional variables are contained in the returned list if argument \code{method_jmaxlim=="smith19"}
#' \itemize{
#'         \item \code{omega}: Term corresponding to \eqn{\omega}, defined by Eq. 16 in Smith et al. (2019), 
#'         and Eq. E19 in Stocker et al. (2019). 
#'         \item \code{omega_star}: Term corresponding to \eqn{\omega^\ast}, defined by Eq. 18 in Smith et al. 
#'         (2019), and Eq. E21 in Stocker et al. (2019). 
#'         }
#'
#' @references  Bernacchi, C. J., Pimentel, C., and Long, S. P.:  In vivo temperature response func-tions  of  parameters  
#'              required  to  model  RuBP-limited  photosynthesis,  Plant  Cell Environ., 26, 1419–1430, 2003
#' 
#'              Heskel,  M.,  O’Sullivan,  O.,  Reich,  P.,  Tjoelker,  M.,  Weerasinghe,  L.,  Penillard,  A.,Egerton, J.,
#'              Creek, D., Bloomfield, K., Xiang, J., Sinca, F., Stangl, Z., Martinez-De La Torre, A., Griffin, K.,
#'              Huntingford, C., Hurry, V., Meir, P., Turnbull, M.,and Atkin, O.:  Convergence in the temperature response
#'              of leaf respiration across biomes and plant functional types, Proceedings of the National Academy of Sciences,
#'              113,  3832–3837,  doi:10.1073/pnas.1520282113,2016.
#' 
#'              Huber,  M.  L.,  Perkins,  R.  A.,  Laesecke,  A.,  Friend,  D.  G.,  Sengers,  J.  V.,  Assael,M. J.,
#'              Metaxa, I. N., Vogel, E., Mares, R., and Miyagawa, K.:  New international formulation for the viscosity
#'              of H2O, Journal of Physical and Chemical ReferenceData, 38, 101–125, 2009
#'
#'              Prentice,  I. C.,  Dong,  N.,  Gleason,  S. M.,  Maire,  V.,  and Wright,  I. J.:  Balancing the costs
#'              of carbon gain and water transport:  testing a new theoretical frameworkfor  plant  functional  ecology,
#'              Ecology  Letters,  17,  82–91,  10.1111/ele.12211,http://dx.doi.org/10.1111/ele.12211, 2014.
#'
#'              Wang, H., Prentice, I. C., Keenan, T. F., Davis, T. W., Wright, I. J., Cornwell, W. K.,Evans, B. J.,
#'              and Peng, C.:  Towards a universal model for carbon dioxide uptake by plants, Nat Plants, 3, 734–741, 2017.
#'              Atkin, O. K., et al.:  Global variability in leaf respiration in relation to climate, plant func-tional
#'              types and leaf traits, New Phytologist, 206, 614–636, doi:10.1111/nph.13253,
#'              https://nph.onlinelibrary.wiley.com/doi/abs/10.1111/nph.13253.
#'              
#'              Smith, N. G., Keenan, T. F., Colin Prentice, I. , Wang, H. , Wright, I. J., Niinemets, U. , Crous, K. Y., 
#'              Domingues, T. F., Guerrieri, R. , Yoko Ishida, F. , Kattge, J. , Kruger, E. L., Maire, V. , Rogers, A. , 
#'              Serbin, S. P., Tarvainen, L. , Togashi, H. F., Townsend, P. A., Wang, M. , Weerasinghe, L. K. and Zhou, S. 
#'              (2019), Global photosynthetic capacity is optimized to the environment. Ecol Lett, 22: 506-517. 
#'              doi:10.1111/ele.13210
#'              
#'              Stocker, B. et al. Geoscientific Model Development Discussions (in prep.) 
#'
#' @export
#'
#' @examples rpmodel( tc = 20, vpd = 1000, co2 = 400, fapar = 1, ppfd = 300, elv = 0)
#'
rpmodel <- function( tc, vpd, co2, fapar, ppfd, patm = NA, elv = NA, 
                     kphio = ifelse(do_ftemp_kphio, ifelse(do_soilmstress, 0.087182, 0.081785), 0.049977), 
                     beta = 146.0, soilm = 1.0, meanalpha = 1.0, apar_soilm = 0.0, bpar_soilm = 0.73300, 
                     c4 = FALSE, method_optci = "prentice14", method_jmaxlim = "wang17", 
                     do_ftemp_kphio = TRUE, do_soilmstress = FALSE, returnvar = NULL, verbose = FALSE ){
  
  # Check arguments
  if (identical(NA, elv) && identical(NA, patm)){
    rlang::abort("Aborted. Provide either elevation (arugment elv) or atmospheric pressure (argument patm).")
  } else if (!identical(NA, elv) && identical(NA, patm)){
    if (verbose) rlang::warn("Atmospheric pressure (patm) not provided. Calculating it as a function of elevation (elv), assuming standard atmosphere (101325 Pa at sea level).")
    patm <- calc_patm(elv)
  }
  
  #-----------------------------------------------------------------------
  # Fixed parameters
  #-----------------------------------------------------------------------
  c_molmass <- 12.0107  # molecular mass of carbon (g)
  kPo <- 101325.0       # standard atmosphere, Pa (Allen, 1973)
  kTo <- 25.0           # base temperature, deg C (Prentice, unpublished)
  rd_to_vcmax <- 0.015  # Ratio of Rdark to Vcmax25, number from Atkin et al., 2015 for C3 herbaceous

  #-----------------------------------------------------------------------
  # Temperature dependence of quantum yield efficiency
  #-----------------------------------------------------------------------
  ## 'do_ftemp_kphio' is not actually a stress function, but is the temperature-dependency of
  ## the quantum yield efficiency after Bernacchi et al., 2003 PCE
  if (do_ftemp_kphio){
    ftemp_kphio <- calc_ftemp_kphio( tc )
  } else {
    ftemp_kphio <- 1.0
  }

  #-----------------------------------------------------------------------
  # Calculate soil moisture stress as a function of soil moisture and mean alpha
  #-----------------------------------------------------------------------
  if (do_soilmstress) {
    soilmstress <- calc_soilmstress( soilm, meanalpha, apar_soilm, bpar_soilm )
  }
  else {
    soilmstress <- 1.0
  }

  #-----------------------------------------------------------------------
  # Photosynthesis model parameters depending on temperature, pressure, and CO2.
  #-----------------------------------------------------------------------
  ## ambient CO2 partial pression (Pa)
  ca <- co2_to_ca( co2, patm )

  ## photorespiratory compensation point - Gamma-star (Pa)
  gammastar <- calc_gammastar( tc, patm )

  ## Michaelis-Menten coef. (Pa)
  kmm <- calc_kmm( tc, patm )   ## XXX Todo: replace 'NA' here with 'patm'

  ## viscosity correction factor = viscosity( temp, press )/viscosity( 25 degC, 1013.25 Pa)
  ns      <- calc_viscosity_h2o( tc, patm )  # Pa s
  ns25    <- calc_viscosity_h2o( kTo, kPo )  # Pa s
  ns_star <- ns / ns25  # (unitless)

  ##-----------------------------------------------------------------------
  ## Optimal ci
  ## The heart of the P-model: calculate ci:ca ratio (chi) and additional terms
  ##-----------------------------------------------------------------------
  if (c4){

    # "dummy" ci:ca for C4 plants
    out_optchi <- calc_chi_c4()

  } else if (method_optci=="prentice14"){

    ## Full formualation (Gamma-star not zero), analytical solution
    ##-----------------------------------------------------------------------
    out_optchi <- calc_optimal_chi( kmm, gammastar, ns_star, ca, vpd, beta )


  } else {

    rlang::abort("rpmodel(): argument method_optci not idetified.")

  }


  ## leaf-internal CO2 partial pressure (Pa)
  ci <- out_optchi$chi * ca

  ##-----------------------------------------------------------------------
  ## Corrolary preditions
  ##-----------------------------------------------------------------------
  ## intrinsic water use efficiency (in Pa)
  iwue = ( ca - ci ) / 1.6

  ##-----------------------------------------------------------------------
  ## Vcmax and light use efficiency
  ##-----------------------------------------------------------------------
  if (c4){

    out_lue_vcmax <- calc_lue_vcmax_c4(kphio, ftemp_kphio, c_molmass, soilmstress)


  } else if (method_jmaxlim=="wang17"){


    out_lue_vcmax <- calc_lue_vcmax_wang17(out_optchi, kphio, ftemp_kphio, c_molmass, soilmstress)


  } else if (method_jmaxlim=="smith19"){


    out_lue_vcmax <- calc_lue_vcmax_smith19(out_optchi, kphio, ftemp_kphio, c_molmass, soilmstress)


  } else if (method_jmaxlim=="none"){


    out_lue_vcmax <- calc_lue_vcmax_none(out_optchi, kphio, ftemp_kphio, c_molmass, soilmstress)


  } else {

    rlang::abort("rpmodel(): argument method_jmaxlim not idetified.")

  }


  ##-----------------------------------------------------------------------
  ## Corrolary preditions
  ##-----------------------------------------------------------------------
  ## Vcmax25 (vcmax normalized to 25 deg C)
  ftemp25_inst_vcmax  <- calc_ftemp_inst_vcmax( tc, tc, tcref = 25.0 )
  vcmax25_unitiabs  <- out_lue_vcmax$vcmax_unitiabs  / ftemp25_inst_vcmax

  ## Dark respiration at growth temperature
  ftemp_inst_rd <- calc_ftemp_inst_rd( tc )
  rd_unitiabs  <- rd_to_vcmax * (ftemp_inst_rd / ftemp25_inst_vcmax) * out_lue_vcmax$vcmax_unitiabs


  ##-----------------------------------------------------------------------
  ## Quantities that scale linearly with absorbed light
  ##-----------------------------------------------------------------------
  len <- length(out_lue_vcmax[[1]])
  iabs <- rep(fapar * ppfd, len)

  ## Gross primary productivity
  gpp <- ifelse(!is.na(iabs), iabs * out_lue_vcmax$lue, rep(NA, len))   # in g C m-2 s-1

  ## Vcmax per unit ground area is the product of the intrinsic quantum
  ## efficiency, the absorbed PAR, and 'n'
  vcmax <- ifelse(!is.na(iabs), iabs * out_lue_vcmax$vcmax_unitiabs, rep(NA, len))

  ## (vcmax normalized to 25 deg C)
  vcmax25 <- ifelse(!is.na(iabs), iabs * vcmax25_unitiabs, rep(NA, len))

  ## Dark respiration
  rd <- ifelse(!is.na(iabs), iabs * rd_unitiabs, rep(NA, len))

  ## construct list for output
  out <- list(
              ca              = rep(ca, len),
              gammastar       = rep(gammastar, len),
              kmm             = rep(kmm, len),
              ns_star         = rep(ns_star, len),
              chi             = out_optchi$chi,
              mj              = out_optchi$mj,
              mc              = out_optchi$mc,
              ci              = ci,
              lue             = out_lue_vcmax$lue,
              gpp             = gpp,
              iwue            = iwue,
              gs              = (gpp / c_molmass) / (ca - ci),
              vcmax           = vcmax,
              vcmax25         = vcmax25,
              rd              = rd
              )

  if (!is.null(returnvar)) out <- out[returnvar]

  return( out )

}


calc_optimal_chi <- function( kmm, gammastar, ns_star, ca, vpd, beta ){
  #-----------------------------------------------------------------------
  # Input:    - float, 'kmm' : Pa, Michaelis-Menten coeff.
  #           - float, 'ns_star'  : (unitless) viscosity correction factor for water
  #           - float, 'vpd' : Pa, vapor pressure deficit
  # Output:   float, ratio of ci/ca (chi)
  # Features: Returns an estimate of leaf internal to ambient CO2
  #           partial pressure following the "simple formulation".
  # Depends:  - kc
  #           - ns
  #           - vpd
  #-----------------------------------------------------------------------

  ## leaf-internal-to-ambient CO2 partial pressure (ci/ca) ratio
  xi  <- sqrt( (beta * ( kmm + gammastar ) ) / ( 1.6 * ns_star ) )
  chi <- gammastar / ca + ( 1.0 - gammastar / ca ) * xi / ( xi + sqrt(vpd) )

  # Define variable substitutes:
  vdcg <- ca - gammastar
  vacg <- ca + 2.0 * gammastar
  vbkg <- beta * (kmm + gammastar)

  ## wrap if condition in a function to allow vectorization
  calc_mj <- function(ns_star, vpd, vbkg){
    vsr <- sqrt( 1.6 * ns_star * vpd / vbkg ) 

    # Based on the mc' formulation (see Regressing_LUE.pdf)
    mj <- vdcg / ( vacg + 3.0 * gammastar * vsr )
      
    return(mj)
  }
  
  # Check for negatives, vectorized
  mj <- ifelse(ns_star>0 & vpd>0 & vbkg>0, calc_mj(ns_star, vpd, vbkg), rep(NA, length(vpd)))
  
  ## alternative variables
  gamma <- gammastar / ca
  kappa <- kmm / ca

  ## mc
  mc <- (chi - gamma) / (chi + kappa)

  ## mj:mv
  mjoc <- (chi + kappa) / (chi + 2 * gamma)

  out <- list( chi=chi, mc=mc, mj=mj, mjoc=mjoc )
  return(out)
}


calc_lue_vcmax_wang17 <- function(out_optchi, kphio, ftemp_kphio, c_molmass, soilmstress){
  
  ## Include effect of Jmax limitation
  len <- length(out_optchi[[1]])
  mprime <- calc_mprime( out_optchi$mj )

  out <- list(

    ## Light use efficiency (gpp per unit absorbed light)
    lue = kphio * ftemp_kphio * mprime * c_molmass * soilmstress,

    ## Vcmax normalised per unit absorbed PPFD (assuming iabs=1), with Jmax limitation
    vcmax_unitiabs = kphio * ftemp_kphio * out_optchi$mjoc * mprime / out_optchi$mj * soilmstress,

    ## complement for non-smith19 
    omega               = rep(NA, len),
    omega_star          = rep(NA, len)
    )

  return(out)
}


calc_lue_vcmax_smith19 <- function(out_optchi, kphio, ftemp_kphio, c_molmass, soilmstress){

  len <- length(out_optchi[[1]])
  
  # Adopted from Nick Smith's code:
  # Calculate omega, see Smith et al., 2019 Ecology Letters
  calc_omega <- function( theta, c_cost, m ){

    cm <- 4 * c_cost / m                        # simplification term for omega calculation
    v  <- 1/(cm * (1 - theta * cm)) - 4 * theta # simplification term for omega calculation

    # account for non-linearities at low m values
    capP <- (((1/1.4) - 0.7)^2 / (1-theta)) + 3.4
    aquad <- -1
    bquad <- capP
    cquad <- -(capP * theta)
    m_star <- (4 * c_cost) / polyroot(c(aquad, bquad, cquad))

    omega <- ifelse(  m < Re(m_star[1]),
                      -( 1 - (2 * theta) ) - sqrt( (1 - theta) * v),
                      -( 1 - (2 * theta))  + sqrt( (1 - theta) * v)
                      )
    return(omega)
  }

  ## constants
  theta <- 0.85    # should be calibratable?
  c_cost <- 0.05336251


  ## factors derived as in Smith et al., 2019
  omega <- calc_omega( theta = theta, c_cost = c_cost, m = out_optchi$mj )          # Eq. S4
  omega_star <- 1.0 + omega - sqrt( (1.0 + omega)^2 - (4.0 * theta * omega) )       # Eq. 18

  ## Effect of Jmax limitation
  mprime <- out_optchi$mj * omega_star / (8.0 * theta)

  ## Light use efficiency (gpp per unit absorbed light)
  lue <- kphio * ftemp_kphio * mprime * c_molmass * soilmstress

  # calculate Vcmax per unit aborbed light
  vcmax_unitiabs  <- kphio * ftemp_kphio * out_optchi$mjoc * omega_star / (8.0 * theta) * soilmstress   # Eq. 19

  out <- list(
    lue            = lue,
    vcmax_unitiabs = vcmax_unitiabs,
    omega          = omega,
    omega_star     = omega_star
    )

  return(out)
}


calc_lue_vcmax_none <- function(out_optchi, kphio, ftemp_kphio, c_molmass, soilmstress){
  ## Do not include effect of Jmax limitation
  len <- length(out_optchi[[1]])
  
  out <- list(

    ## Light use efficiency (gpp per unit absorbed light)
    lue = kphio * ftemp_kphio * out_optchi$mj * c_molmass * soilmstress,

    ## Vcmax normalised per unit absorbed PPFD (assuming iabs=1), with Jmax limitation
    vcmax_unitiabs = kphio * ftemp_kphio * out_optchi$mjoc * soilmstress,

    ## complement for non-smith19 
    omega               = rep(NA, len),
    omega_star          = rep(NA, len)
    )

  return(out)
}


calc_lue_vcmax_c4 <- function( kphio, ftemp_kphio, c_molmass, soilmstress ){

  len <- length(kphio)
  out <- list(
    ## Light use efficiency (gpp per unit absorbed light)
    lue = kphio * ftemp_kphio * c_molmass * soilmstress,

    ## Vcmax normalised per unit absorbed PPFD (assuming iabs=1), with Jmax limitation
    vcmax_unitiabs = kphio * ftemp_kphio * soilmstress,

    ## complement for non-smith19 
    omega               = rep(NA, len),
    omega_star          = rep(NA, len)
  )

  return(out)
}


calc_chi_c4 <- function(){
  #//////////////////////////////////////////////////////////////////
  # (Dummy-) ci:ca for C4 photosynthesis
  #-----------------------------------------------------------------------
  out <- list( chi=9999, mc=1, mj=1, mjoc=1 )
  return(out)
}


calc_mprime <- function( mc ){
  #-----------------------------------------------------------------------
  # Input:  mc   (unitless): factor determining LUE
  # Output: mpi (unitless): modiefied m accounting for the co-limitation
  #                         hypothesis after Prentice et al. (2014)
  #-----------------------------------------------------------------------
  kc <- 0.41          # Jmax cost coefficient

  mpi <- mc^2 - kc^(2.0/3.0) * (mc^(4.0/3.0))

  # Check for negatives:
  mpi <- ifelse(mpi>0, sqrt(mpi), NA)

  return(mpi)
}


co2_to_ca <- function( co2, patm ){
  #-----------------------------------------------------------------------
  # Input:    - float, annual atm. CO2, ppm (co2)
  #           - float, monthly atm. pressure, Pa (patm)
  # Output:   - ca in units of Pa
  # Features: Converts ca (ambient CO2) from ppm to Pa.
  #-----------------------------------------------------------------------
  ca   <- ( 1.0e-6 ) * co2 * patm         # Pa, atms. CO2
  return( ca )
}


density_h2o <- function( tc, p ){
  #-----------------------------------------------------------------------
  # Input:    - float, air temperature (tc), degrees C
  #           - float, atmospheric pressure (p), Pa
  # Output:   float, density of water, kg/m^3
  # Features: Calculates density of water at a given temperature and
  #           pressure using the Tumlirz Equation
  # Ref:      F.H. Fisher and O.E Dial, Jr. (1975) Equation of state of
  #           pure water and sea water, Tech. Rept., Marine Physical
  #           Laboratory, San Diego, CA.
  #-----------------------------------------------------------------------

  # Calculate lambda, (bar cm^3)/g:
  my_lambda <- 1788.316 +
          21.55053*tc +
        -0.4695911*tc*tc +
     (3.096363e-3)*tc*tc*tc +
    -(7.341182e-6)*tc*tc*tc*tc

  # Calculate po, bar
  po <- 5918.499 +
           58.05267*tc +
         -1.1253317*tc*tc +
     (6.6123869e-3)*tc*tc*tc +
    -(1.4661625e-5)*tc*tc*tc*tc

  # Calculate vinf, cm^3/g
  vinf <- 0.6980547 +
    -(7.435626e-4)*tc +
     (3.704258e-5)*tc*tc +
    -(6.315724e-7)*tc*tc*tc +
     (9.829576e-9)*tc*tc*tc*tc +
   -(1.197269e-10)*tc*tc*tc*tc*tc +
    (1.005461e-12)*tc*tc*tc*tc*tc*tc +
   -(5.437898e-15)*tc*tc*tc*tc*tc*tc*tc +
     (1.69946e-17)*tc*tc*tc*tc*tc*tc*tc*tc +
   -(2.295063e-20)*tc*tc*tc*tc*tc*tc*tc*tc*tc

  # Convert pressure to bars (1 bar <- 100000 Pa)
  pbar <- (1e-5)*p

  # Calculate the specific volume (cm^3 g^-1):
  v <- vinf + my_lambda/(po + pbar)

  # Convert to density (g cm^-3) -> 1000 g/kg; 1000000 cm^3/m^3 -> kg/m^3:
  rho <- (1e3/v)

  return(rho)
}


calc_viscosity_h2o <- function( tc, p ) {
  #-----------------------------------------------------------------------
  # Input:    - float, ambient temperature (tc), degrees C
  #           - float, ambient pressure (p), Pa
  # Return:   float, viscosity of water (mu), Pa s
  # Features: Calculates viscosity of water at a given temperature and
  #           pressure.
  # Depends:  density_h2o
  # Ref:      Huber, M. L., R. A. Perkins, A. Laesecke, D. G. Friend, J. V.
  #           Sengers, M. J. Assael, ..., K. Miyagawa (2009) New
  #           international formulation for the viscosity of H2O, J. Phys.
  #           Chem. Ref. Data, Vol. 38(2), pp. 101-125.
  #-----------------------------------------------------------------------

  # Define reference temperature, density, and pressure values:
  tk_ast  <- 647.096    # Kelvin
  rho_ast <- 322.0      # kg/m^3
  mu_ast  <- 1e-6       # Pa s

  # Get the density of water, kg/m^3
  rho <- density_h2o(tc, p)

  # Calculate dimensionless parameters:
  tbar  <- (tc + 273.15)/tk_ast
  tbarx <- tbar^(0.5)
  tbar2 <- tbar^2
  tbar3 <- tbar^3
  rbar  <- rho/rho_ast

  # Calculate mu0 (Eq. 11 & Table 2, Huber et al., 2009):
  mu0 <- 1.67752 + 2.20462/tbar + 0.6366564/tbar2 - 0.241605/tbar3
  mu0 <- 1e2*tbarx/mu0

  # Create Table 3, Huber et al. (2009):
  h_array <- array(0.0, dim=c(7,6))
  h_array[1,] <- c(0.520094, 0.0850895, -1.08374, -0.289555, 0.0, 0.0)  # hj0
  h_array[2,] <- c(0.222531, 0.999115, 1.88797, 1.26613, 0.0, 0.120573) # hj1
  h_array[3,] <- c(-0.281378, -0.906851, -0.772479, -0.489837, -0.257040, 0.0) # hj2
  h_array[4,] <- c(0.161913,  0.257399, 0.0, 0.0, 0.0, 0.0) # hj3
  h_array[5,] <- c(-0.0325372, 0.0, 0.0, 0.0698452, 0.0, 0.0) # hj4
  h_array[6,] <- c(0.0, 0.0, 0.0, 0.0, 0.00872102, 0.0) # hj5
  h_array[7,] <- c(0.0, 0.0, 0.0, -0.00435673, 0.0, -0.000593264) # hj6

  # Calculate mu1 (Eq. 12 & Table 3, Huber et al., 2009):
  mu1 <- 0.0
  ctbar <- (1.0/tbar) - 1.0
  # print(paste("ctbar",ctbar))
  # for i in xrange(6):
  for (i in 1:6){
    coef1 <- ctbar^(i-1)
    # print(paste("i, coef1", i, coef1))
    coef2 <- 0.0
    for (j in 1:7){
      coef2 <- coef2 + h_array[j,i] * (rbar - 1.0)^(j-1)
    }
    mu1 <- mu1 + coef1 * coef2
  }
  mu1 <- exp( rbar * mu1 )
  # print(paste("mu1",mu1))

  # Calculate mu_bar (Eq. 2, Huber et al., 2009)
  #   assumes mu2 = 1
  mu_bar <- mu0 * mu1

  # Calculate mu (Eq. 1, Huber et al., 2009)
  mu <- mu_bar * mu_ast    # Pa s

  return( mu )
}
