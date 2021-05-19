#' Instantaneous scaling of rates
#'
#' Calculates instantaneous scaling of rates, given acclimated 
#' capacities returned by function \code{rpmodel}.
#'
#' @param x An object (list) returned by function \code{rpmodel}.
#' @param tc Temperature, relevant for photosynthesis (deg C)
#' @param vpd Vapour pressure deficit (Pa)
#' @param co2 Atmospheric CO2 concentration (ppm)
#' @param fapar (Optional) Fraction of absorbed photosynthetically active 
#'  radiation (unitless, defaults to \code{NA})
#' @param ppfd (Optional) Photosynthetic photon flux density (mol m-2 d-1, 
#'  defaults to \code{NA}). Note that the units of 
#'  \code{ppfd} (per area and per time) determine the units of outputs 
#'  \code{lue}, \code{gpp}, \code{vcmax}, and \code{rd}. For example, if 
#'  \code{ppfd} is provided in units of mol m-2 month-1, then
#'  respective output variables are returned as per unit months.
#' @param patm Atmospheric pressure (Pa). When provided, overrides \code{elv},
#'  otherwise \code{patm} is calculated using standard atmosphere (101325 Pa),
#'  corrected for elevation (argument \code{elv}), using the function
#'  \link{patm}.
#' @param elv Elevation above sea-level (m.a.s.l.). Is used only for calculating
#'  atmospheric pressure (using standard atmosphere (101325 Pa), corrected for
#'  elevation (argument \code{elv}), using the function \link{patm}),
#'  if argument \code{patm} is not provided. If argument \code{patm} is 
#'  provided, \code{elv} is overridden.
#' @param kphio Apparent quantum yield efficiency (unitless). Defaults to 0.0817
#'  for \code{method_jmaxlim="wang17", do_ftemp_kphio=TRUE, do_soilmstress=FALSE},
#'  0.0870 for \code{method_jmaxlim="wang17", do_ftemp_kphio=TRUE, do_soilmstress=TRUE},
#'  and 0.0492 for \code{method_jmaxlim="wang17", do_ftemp_kphio=FALSE, do_soilmstress=FALSE},
#'  corresponding to the empirically fitted value as presented in Stocker et al.
#'  (2019) Geosci. Model Dev. for model setup 'BRC', 'FULL', and 'ORG'
#'  respectively.
#' @param verbose provide verbose feedback (default = FALSE)
#' @return A list containing instantaneous rates calculated from acclimated 
#'  photosyntetic capacities (Vcmax25,Jmax25, "average" gs). This includes :
#' \itemize{
#' \item \code{vcmax}
#' \item \code{jmax}
#' \item \code{rd}
#' }
#'
#' @export
#'
#' @examples
#' \dontrun{
#'  dampen_vec(
#'   vec = 20 * (sin(doy*pi/(365)))^2 + rnorm(365, mean = 0, sd = 5),
#'   tau = 40 )
#' }
#'

inst_rpmodel <- function(
  x,
  tc,
  vpd,
  co2,
  fapar,
  ppfd,
  patm = NA,
  elv = NA,
  kphio,
  verbose = FALSE 
  ){

  ## output
  out_inst <- list()

  #---- Check arguments ----
  if (identical(NA, elv) && identical(NA, patm)){
    stop("Aborted. Provide either elevation (arugment elv) or
                 atmospheric pressure (argument patm).")
  } else if (!identical(NA, elv) && identical(NA, patm)){
    if (verbose) warning("Atmospheric pressure (patm) not provided. 
                             Calculating it as a function of elevation (elv),
                             assuming standard atmosphere
                             (101325 Pa at sea level).")
    patm <- patm(elv)
  }

  #---- Parameters ----
  # this parameter value must be the same as the one used in rpmodel. 
  # Should think of a better way to specify parameters that are "public" 
  # across the package.
  rd_to_vcmax <- 0.015  # Ratio of Rdark to Vcmax25, number from Atkin et al., 2015 for C3 herbaceous
  c_molmass <- 12.0107  # molecular mass of carbon (g)

  #---- Required quantities ----
  iabs <- ppfd   # leaf-level quantity - no fapar

  # ambient CO2 partial pression (Pa)
  ca <- co2_to_ca( co2, patm )

  # photorespiratory compensation point - Gamma-star (Pa)
  gammastar <- gammastar( tc, patm )

  # Michaelis-Menten coef. (Pa)
  kmm <- kmm( tc, patm )   ## XXX Todo: replace 'NA' here with 'patm'

  #---- Vcmax, at current temperature ----
  ftemp25_inst_vcmax  <- ftemp_inst_vcmax( tc, tc, tcref = 25.0 )
  vcmax <- x$vcmax25 * ftemp25_inst_vcmax

  #---- Rd: Dark respiration ----
  ftemp_inst_rd <- ftemp_inst_rd( tc )
  rd <- rd_to_vcmax * x$vcmax25 * ftemp_inst_rd

  #---- Jmax, at current temperature ----
  ftemp25_inst_jmax  <- ftemp_inst_jmax( tc, tc, tcref = 25.0 )
  jmax <- x$jmax25 * ftemp25_inst_jmax

  #---- Aj, gs free ----
  L <- 1.0 / sqrt(1.0 + ((4.0 * kphio * iabs)/jmax)^2)
  kv <- (ca - gammastar) / (1 + x$xi / sqrt(vpd))
  ci_j <- ca - kv
  a_j <- L * kphio * iabs * (ci_j - gammastar)/(ci_j + 2 * gammastar)
  gs_j <- a_j / kv

  #---- Ac, gs free ----
  ci_c <- ci_j
  a_c <- vcmax * (ci_c - gammastar)/(ci_c + kmm)
  gs_c <- a_j / kv

  ## Ac
  # A <- -1.0 * x$gs
  # B <- x$gs * ca - x$gs * kmm - vcmax
  # C <- x$gs * ca * kmm + vcmax * gammastar
  #
  # ci_c <- QUADM(A, B, C)
  # a_c <- vcmax * (ci_c - gammastar) / (ci_c + kmm)
  #
  ## AJ 
  # L <- 1.0 / sqrt(1.0 + ((4.0 * kphio * iabs)/jmax)^2)
  # A <- -x$gs
  # B <- x$gs * ca - 2 * gammastar * x$gs - L * kphio * iabs
  # C <- 2 * gammastar * x$gs * ca + L * kphio * iabs * gammastar
  #
  # ci_j <- QUADM(A, B, C)
  # a_j  <- kphio * iabs * (ci_j - gammastar)/(ci_j + 2 * gammastar) * L

  #---- A ----
  # Take minimum of the two assimilation rates and maximum of the two ci
  assim <- ifelse(a_j < a_c, a_j, a_c)
  ci <- ifelse(ci_c > ci_j, ci_c, ci_j)

  #---- GPP ----
  gpp <- assim * c_molmass * fapar

  # construct list for output
  out <- list(
              gpp = gpp,
              assim = assim,
              vcmax = vcmax,
              jmax = jmax,
              rd = rd,
              ci = ci,
              a_j = a_j,
              a_c = a_c,
              ci_j = ci_j,
              ci_c = ci_c
              )

  return( out )
}