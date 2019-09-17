#' Calculates the CO2 compensation point
#'
#' Calculates the photorespiratory CO2 compensation point in absence of dark
#' respiration, \eqn{\Gamma*} (Farquhar, 1980).
#'
#' @param tc Temperature, relevant for photosynthesis (degrees Celsius)
#' @param patm Atmospheric pressure (Pa)
#'
#' @details The temperature and pressure-dependent photorespiratory
#' compensation point in absence of dark respiration \eqn{\Gamma* (T,p)}
#' is calculated from its value at standard temperature (\eqn{T0 = 25 }deg C)
#' and atmospheric pressure (\eqn{p0 = 101325} Pa), referred to as \eqn{\Gamma*0},
#' quantified by Bernacchi et al. (2001) to 4.332 Pa (their value in molar
#' concentration units is multiplied here with 101325 Pa to yield 4.332 Pa).
#' \eqn{\Gamma*0} is modified by temperature following an Arrhenius-type temperature
#' response function \eqn{f(T, \Delta Ha)} (implemented by \link{calc_ftemp_arrh})
#' with activation energy \eqn{\Delta Ha = 37830} J mol-1  and is corrected for
#' atmospheric pressure \eqn{p(z)} (see \link{calc_patm}) at elevation \eqn{z}.
#' \deqn{
#'       \Gamma* = \Gamma*0 f(T, \Delta Ha) p(z) / p_0
#' }
#' \eqn{p(z)} is given by argument \code{patm}.
#'
#' @references Farquhar,  G.  D.,  von  Caemmerer,  S.,  and  Berry,  J.  A.:
#'             A  biochemical  model  of photosynthetic CO2 assimilation in leaves of
#'             C 3 species, Planta, 149, 78–90, 1980.
#'
#'             Bernacchi,  C.  J.,  Singsaas,  E.  L.,  Pimentel,  C.,  Portis,  A.
#'             R.  J.,  and  Long,  S.  P.:Improved temperature response functions
#'             for models of Rubisco-limited photosyn-thesis, Plant, Cell and
#'             Environment, 24, 253–259, 2001
#'
#' @return A numeric value for \eqn{\Gamma*} (in Pa)
#'
#' @examples print("CO2 compensation point at 20 degrees Celsius and standard atmosphere (in Pa):")
#' print(calc_gammastar(20, 101325))
#'
#' @export
#'
calc_gammastar <- function( tc, patm ) {
  #-----------------------------------------------------------------------
  # Input:    float, air temperature, degrees C (tc)
  # Output:   float, gamma-star, Pa (gammastar)
  # Features: Returns the temperature-dependent photorespiratory
  #           compensation point, Gamma star (Pascals), based on constants
  #           derived from Bernacchi et al. (2001) study.
  # Ref:      Bernacchi et al. (2001), Improved temperature response
  #           functions for models of Rubisco-limited photosynthesis,
  #           Plant, Cell and Environment, 24, 253--259.
  #-----------------------------------------------------------------------
  dha    <- 37830       # (J/mol) Activation energy, Bernacchi et al. (2001)
  gs25_0 <- 4.332       # Pa, value based on Bernacchi et al. (2001), converted to Pa by T. Davis assuming elevation of 227.076 m.a.s.l.

  gammastar <- gs25_0 * patm / calc_patm(0.0) * calc_ftemp_arrh( (tc + 273.15), dha=dha )

  return( gammastar )
}
