#' Michaelis Menten coefficient
#'
#' Calculates the Michaelis Menten coefficient of Rubisco-limited assimilation
#' as a function of temperature and atmospheric pressure.
#'
#' @param tc Temperature, relevant for photosynthesis (deg C)
#' @param patm Atmospheric pressure (Pa)
#'
#' @details The Michaelis-Menten coefficient \eqn{K} of Rubisco-limited
#' photosynthesis is determined by the Michalis-Menten constants for
#' O2 and CO2 (Farquhar, 1980) according to:
#' \deqn{
#'   K = Kc ( 1 + pO2 / Ko)
#' }
#' where \eqn{Kc} is the Michaelis-Menten constant for CO2 (Pa), \eqn{Ko} is
#' the Michaelis-Menten constant for O2 (Pa), and \eqn{pO2} is the partial
#' pressure of oxygen (Pa), calculated as \eqn{0.209476 p}, where \eqn{p} is
#' given by argument \code{patm}.  \eqn{Kc} and \eqn{Ko} follow a temperature
#' dependence, given by the Arrhenius Equation \eqn{f}:
#' \deqn{
#'    Kc = Kc25 f(T, \Delta Hkc)
#' }
#' \deqn{
#'    Ko = Ko25 f(T, \Delta Hko)
#' }
#' Values \eqn{\Delta Hkc} (79430 J mol-1), \eqn{\Delta Hko} (36380 J mol-1),
#' \eqn{Kc25} (39.97 Pa), and \eqn{Ko25} (27480 Pa) are taken from Bernacchi
#' et al. (2001) and have been converted from values given therein to units of Pa
#' by multiplication with the standard atmosphere (101325 Pa). \eqn{T} is given
#' by the argument \code{tc}.
#'
#' @return A numeric value for \eqn{K}
#'
#' @examples asdf
#'
#' @export
#'
calc_kmm <- function( tc, patm ) {
  #-----------------------------------------------------------------------
  # Input:    - float, air temperature, deg C (temp)
  #           - float, atmospheric pressure, Pa (patm)
  # Output:   float, Pa (mmk)
  # Features: Returns the temperature & pressure dependent Michaelis-Menten
  #           coefficient, K (Pa).
  # Ref:      Bernacchi et al. (2001), Improved temperature response
  #           functions for models of Rubisco-limited photosynthesis,
  #           Plant, Cell and Environment, 24, 253--259.
  #-----------------------------------------------------------------------
  dhac   <- 79430      # (J/mol) Activation energy, Bernacchi et al. (2001)
  dhao   <- 36380      # (J/mol) Activation energy, Bernacchi et al. (2001)
  kco    <- 2.09476e5  # (ppm) O2 partial pressure, Standard Atmosphere

  ## k25 parameters are not dependent on atmospheric pressure
  kc25 <- 39.97   # Pa, value based on Bernacchi et al. (2001), converted to Pa by T. Davis assuming elevation of 227.076 m.a.s.l.
  ko25 <- 27480   # Pa, value based on Bernacchi et al. (2001), converted to Pa by T. Davis assuming elevation of 227.076 m.a.s.l.

  tk <- tc + 273.15

  kc <- kc25 * calc_ftemp_arrh( tk, dha=dhac )
  ko <- ko25 * calc_ftemp_arrh( tk, dha=dhao )

  po  <- kco * (1e-6) * patm         # O2 partial pressure
  kmm <- kc * (1.0 + po/ko)

  return(kmm)
}
